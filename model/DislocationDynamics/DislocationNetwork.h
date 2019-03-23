/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez <ramirezbrf@gmail.com>
 * Copyright (C) 2011 by Mamdouh Mohamed  <msm07d@fsu.edu>
 * Copyright (C) 2011 by Tamer Crsoby     <tamercrosby@gmail.com>
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>
 * Copyright (C) 2011 by Yinan Cui        <cuiyinan@ucla.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// valgrind --leak-check=full --show-leak-kinds=all ./DDomp

#ifndef model_DISLOCATIONNETWORK_H_
#define model_DISLOCATIONNETWORK_H_

#ifdef _MODEL_MPI_
#define _MODEL_DD_MPI_
#endif


#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif

#include <chrono>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/StdDeque>


#include <model/LoopNetwork/LoopNetwork.h>
#include <model/Utilities/TerminalColors.h>
//#include <model/IO/EigenDataReader.h>
#include <model/DislocationDynamics/DislocationNetworkTraits.h>
#include <model/DislocationDynamics/DislocationNetworkComponent.h>
#include <model/DislocationDynamics/DislocationNode.h>
#include <model/DislocationDynamics/DislocationSegment.h>
#include <model/DislocationDynamics/DislocationLoop.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/DislocationNetworkRemesh.h>
#include <model/DislocationDynamics/Junctions/DislocationJunctionFormation.h>
#include <model/DislocationDynamics/CrossSlip/DislocationCrossSlip.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/DislocationDynamics/IO/DislocationNetworkIO.h>
#include <model/DislocationDynamics/ElasticFields/DislocationParticle.h>
#include <model/DislocationDynamics/ElasticFields/DislocationStress.h>
#include <model/ParticleInteraction/ParticleSystem.h>
#include <model/MPI/MPIcout.h> // defines mode::cout
#include <model/ParticleInteraction/SingleFieldPoint.h>
#include <model/DislocationDynamics/DDtimeIntegrator.h>
#include <model/Threads/EqualIteratorRange.h>
#include <model/DislocationDynamics/BoundingLineSegments.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundaryTransmission.h>
//#include <model/DislocationDynamics/Polycrystals/GrainBoundaryDissociation.h>
#include <model/DislocationDynamics/BVP/BVPsolver.h>
#include <model/DislocationDynamics/Polycrystals/Polycrystal.h>
#include <model/DislocationDynamics/DislocationNodeContraction.h>
#include <model/DislocationDynamics/ElasticFields/EshelbyInclusion.h>
#include <model/IO/TextFileParser.h>


#ifndef ExternalLoadControllerFile
#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/DummyExternalLoadController.h>
#endif
#include ExternalLoadControllerFile

//#include <model/DislocationDynamics/ExternalLoadController.h>
#include <model/DislocationDynamics/DislocationInjector.h>

namespace model
{

template <int dim>
struct DislocationNetworkBase : public GlidePlaneObserver<dim>
{// Class storing objects that need to be destroyed last

	SimplicialMesh<dim>& mesh;
	Polycrystal<dim>& poly;
	BVPsolver<dim,2>& bvpSolver;

	//        DislocationNetworkBase() :
	//        /* init */ mesh(TextFileParser("inputFiles/DD.txt").readScalar<int>("meshID",true)),
	//        /* init */ poly("./inputFiles/polycrystal.txt",mesh,*this),
	//        /* init */ bvpSolver(mesh)
	//        {
	//                            assert(mesh.simplices().size() && "MESH IS EMPTY.");
	//        }

	DislocationNetworkBase(SimplicialMesh<dim>& _mesh,
			Polycrystal<dim>& _poly,
			BVPsolver<dim,2>& _bvpSolver) :
				/* init */ mesh(_mesh),
				/* init */ poly(_poly),
				/* init */ bvpSolver(_bvpSolver)
	{
		//            assert(mesh.simplices().size() && "MESH IS EMPTY.");
	}

};

template <int _dim, short unsigned int _corder, typename InterpolationType>
class DislocationNetwork : public DislocationNetworkBase<_dim>, // must be first in inheritance tree
/* base                 */ public LoopNetwork<DislocationNetwork<_dim,_corder,InterpolationType> >,
/* base                 */ public ParticleSystem<DislocationParticle<_dim> >,
/* base                 */ public std::map<size_t,EshelbyInclusion<_dim>>
{




public:

	enum SimulationType{FINITE_NO_FEM=0,FINITE_FEM=1,PERIODIC=2};

	static constexpr int dim=_dim; // make dim available outside class
	static constexpr int corder=_corder; // make dim available outside class
	typedef DislocationNetwork<dim,corder,InterpolationType> DislocationNetworkType;
	typedef LoopNetwork<DislocationNetworkType> LoopNetworkType;
	typedef DislocationNode<dim,corder,InterpolationType> NodeType;
	typedef DislocationSegment<dim,corder,InterpolationType> LinkType;
	typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
	typedef NetworkComponent<NodeType,LinkType> NetworkComponentType;
	typedef Eigen::Matrix<double,dim,dim>	MatrixDimD;
	typedef Eigen::Matrix<double,dim,1>		VectorDim;
	typedef typename TypeTraits<DislocationNetworkType>::LoopType LoopType;
	typedef GlidePlaneObserver<dim> GlidePlaneObserverType;
	typedef DislocationParticle<_dim> DislocationParticleType;
	typedef typename DislocationParticleType::StressField StressField;
	typedef typename DislocationParticleType::DisplacementField DisplacementField;
	typedef ParticleSystem<DislocationParticleType> ParticleSystemType;
	typedef typename ParticleSystemType::SpatialCellType SpatialCellType;
	typedef SpatialCellObserver<DislocationParticleType,_dim> SpatialCellObserverType;
	typedef BVPsolver<dim,2> BvpSolverType;
	typedef typename BvpSolverType::FiniteElementType FiniteElementType;
	typedef typename FiniteElementType::ElementType ElementType;
	typedef typename LoopNetworkType::IsNodeType IsNodeType;
	typedef DislocationNetworkIO<DislocationNetworkType> DislocationNetworkIOType;
	typedef Polycrystal<dim> PolycrystalType;
	typedef ExternalLoadController<dim> ExternalLoadControllerType;
	typedef std::map<size_t,EshelbyInclusion<_dim>> EshelbyInclusionContainerType;
	typedef NetworkLinkObserver<LinkType> NetworkLinkObserverType;
	typedef typename NetworkLinkObserverType::LinkContainerType NetworkLinkContainerType;

#ifdef DislocationNucleationFile
#include DislocationNucleationFile
#endif

private:

	/**********************************************************************/
	void updateLoadControllers(const long int& runID)
	{/*! Updates bvpSolver using the stress and displacement fields of the
	 *  current DD configuration.
	 */
		const int quadraturePerTriangle=37;
		if(use_bvp)
		{
			if (!(runID%use_bvp))
			{// enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
				model::cout<<"		Updating elastic bvp... "<<std::endl;
				this->bvpSolver.template assembleAndSolve<DislocationNetworkType,quadraturePerTriangle>(*this);
			}
		}
		if (use_externalStress)
		{
			extStressController.update(*this,runID);
		}
	}

	/**********************************************************************/
	void computeNodaVelocities(const long int& runID)
	{

		switch (timeIntegrationMethod)
		{
		case 0:
			DDtimeIntegrator<0>::computeNodaVelocities(*this,runID);
			break;

			//                case 1:
				//                    dt=DDtimeIntegrator<1>::integrate(*this);
			//                    break;

		default:
			assert(0 && "time integration method not implemented");
			break;
		}

		//            if(NodeType::use_velocityFilter)
		//            {
		//                assert(0 && "velocityFilter not implemented yet.");
		//                //            for(auto& node : this->nodes())
		//                //            {
		//                //                node.second->applyVelocityFilter(vMax);
		//                //            }
		//            }


		model::cout<<std::setprecision(3)<<std::scientific<<"		dt="<<dt<<std::endl;
	}


	/**********************************************************************/
	void updateVirtualBoundaryLoops()
	{

		switch (simulationType)
		{
		case FINITE_FEM:
		{
			if(useVirtualExternalLoops)
			{
				const auto t0= std::chrono::system_clock::now();
				model::cout<<"		Updating virtual boundary loops "<<std::flush;

				// First clean up outdated boundary loops
				std::set<size_t> removeLoops;
				for(const auto& loop : this->loops())
				{
					if(loop.second->isVirtualBoundaryLoop() && loop.second->links().size()!=4)
					{// clean up left over loops from topological operations
						removeLoops.insert(loop.second->sID);
					}
				}



				for(const size_t& loopID : removeLoops)
				{// Remove the virtual loops with ID in removeLoops
					this->deleteLoop(loopID);
				}


				// Now reconstruct virtual boundary loops
				std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,size_t>> virtualLoopVector;
				for(const auto& link : this->links())
				{
					if(link.second->isBoundarySegment() && !link.second->hasZeroBurgers())
					{

						virtualLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{link.second->sink,link.second->source,link.second->source->virtualBoundaryNode(),link.second->sink->virtualBoundaryNode()},
								link.second->burgers(),
								(*link.second->grains().begin())->grainID);

					}
				}

				for(const auto& tup : virtualLoopVector)
				{// Insert the new virtual loops
					this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup));
				}

				model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
			}
			break;
		}

		case PERIODIC:
		{
			const auto t0= std::chrono::system_clock::now();

			std::vector<std::tuple<std::vector<std::shared_ptr<NodeType>>,VectorDim,VectorDim,VectorDim,size_t>> DislocationLoopVector;
			for(const auto& link : this->links())
			{

				if(link.second->isBoundarySegment() && !link.second->hasZeroBurgers() && !link.second->sink->isImageBoundaryNode && !link.second->source->isImageBoundaryNode)
				{// First make the check whether the nodes connected to the sinks are not imageNodes and they lie on the boundary
					//Replace all the imageBoundaryNode with the containers
					if (link.second->sink->imageBoundaryNodeContainer().size()==1 && link.second->source->imageBoundaryNodeContainer().size()==1)
					{
						const auto sink_image=link.second->sink->imageBoundaryNodeContainer().begin()->second;
						const auto source_image=link.second->source->imageBoundaryNodeContainer().begin()->second;
//
//					const auto imageLink(this->link(std::min(link.second->sink->imageBoundaryNode()->sID,link.second->source->imageBoundaryNode()->sID),std::max(link.second->sink->imageBoundaryNode()->sID
//							,link.second->source->imageBoundaryNode()->sID)));
					const auto imageLink(this->link(std::min(sink_image->sID, source_image->sID),std::max(sink_image->sID, source_image->sID)));
					//const auto imageLink(this->link(link.second->source->imageBoundaryNode()->sID,link.second->sink->imageBoundaryNode()->sID));

					if(!imageLink.first)
					{// Checking whether the virtual nodes of the original link are already connected or not
						//						std::cout<<"Size of virtual boundary loop from source  ID"<<link.second->source->virtualBoundaryNode()->sID<<"is"<<link.second->source->virtualBoundaryNode()->loops().size()<<std::endl;
						//						std::cout<<"Size of loopLinks from source  ID"<<link.second->source->virtualBoundaryNode()->sID<<"is"<<link.second->source->virtualBoundaryNode()->loopLinks().size()<<std::endl;
						//						std::cout<<"Size of neighbors from source  ID"<<link.second->source->virtualBoundaryNode()->sID<<"is"<<link.second->source->virtualBoundaryNode()->neighbors().size()<<std::endl;
						//						std::cout<<"Size of virtual boundary loop from sink ID"<<link.second->sink->virtualBoundaryNode()->sID<<"is"<<link.second->sink->virtualBoundaryNode()->loops().size()<<std::endl;
						//						std::cout<<"Size of loopLinks from sink  ID"<<link.second->sink->virtualBoundaryNode()->sID<<"is"<<link.second->sink->virtualBoundaryNode()->loopLinks().size()<<std::endl;
						//						std::cout<<"Size of neighbors from sink  ID"<<link.second->sink->virtualBoundaryNode()->sID<<"is"<<link.second->sink->virtualBoundaryNode()->neighbors().size()<<std::endl;
//
//						if(   link.second->  sink->imageBoundaryNode()->loops().size()==0
//								&& link.second->source->imageBoundaryNode()->loops().size()==0)


						if(   (sink_image->loops().size()==0
									&& source_image->loops().size()==0))
//								||  (sink_image->loops().size()>0
//											&& source_image->loops().size()>0))	//the second condition is added later after visit to UCLA
						{// If the virtual nodes are not connected to any loops, create a midnode based on the two virtual nodes and form a loop between them
//							bool temp=true;
//							bool temp1=false;
//							if((sink_image->loops().size()>0
//									&& source_image->loops().size()>0))
//							{	temp1=true;
//								for (const auto& iter:sink_image->neighbors())
//								{
//									model::cout<<"sink image "<<sink_image->sID<<" neighbors are "<<std::get<0>(iter.second)->sID<<std::endl;
//									if (std::get<0>(iter.second)->sID==source_image->sID)
//										{
//										temp=false;
//
//										}
//								}
//							}
//							model::cout<<"Source SID is "<<source_image->sID<<std::endl;
//							model::cout<<"temp for case 1 "<<source_image->sID<<"-> "<<sink_image->sID<<"is ::::"<<temp<<std::endl;
//							model::cout<<"temp1 for case 1 "<<source_image->sID<<"-> "<<sink_image->sID<<"is ::::"<<temp1<<std::endl;
//
//							if (temp)
//							{
								model::cout<<redColor<<"Entered the case 1 "<<source_image->sID<<" ->"<<sink_image->sID<<std::endl;
								VectorDim N = (*link.second->loopLinks().begin())->loop()->glidePlane->unitNormal;
								// Direction along which the new node should be placed
								VectorDim X=((-source_image->get_P()+sink_image->get_P()).cross(N)).normalized();
								//							VectorDim X=link.second->burgers();
								// Make sure the node is placed inside the simulation domain
								if (!(X.dot(link.second->source->bndNormal().normalized()) > 0)) X = -X;
//								if (temp1)			this was for multiple loop case
//								X=X*10;
//								else
								X=X*(10);
								// Create a new node based on two virtual nodes
								std::shared_ptr<NodeType> midNode(new NodeType(this,X+0.5*(source_image->get_P()+sink_image->get_P())
										,0.5*(source_image->get_V()+sink_image->get_V()),1));
								std::cout<<"Creating the mid node"<<midNode->sID<<std::endl;
								// Create a container to store the nodes and the direction of the loop to be inserted
								DislocationLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{source_image,midNode,sink_image},
										link.second->burgers(),N,source_image->get_P(),
										(*link.second->grains().begin())->grainID);
								//							DislocationLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{link.second->sink->imageBoundaryNode(),midNode,link.second->source->imageBoundaryNode()},
								//																link.second->burgers(),N,link.second->source->imageBoundaryNode()->get_P(),
								//																(*link.second->grains().begin())->grainID);
								//                                    std::cout<<"Creating the dislocationloop"<<std::endl;
							}

//						}
						else if((sink_image->loops().size()>0
								&& source_image->loops().size()>0))
						{
							model::cout<<"Case for the merging of the loops needed to be taken into consideration::"<<std::endl;
//							bool temp=true;
//							bool temp1=true;
//							for (const auto& iter:sink_image->neighbors())
//							{
//								model::cout<<"sink image "<<sink_image->sID<<" neighbors are "<<std::get<0>(iter.second)->sID<<std::endl;
//								if (std::get<0>(iter.second)->sID==source_image->sID)
//								{
//									temp=false;
//
//								}
//							}
//
////							model::cout<<"Source SID is "<<source_image->sID<<std::endl;
////							model::cout<<"temp for case 1 "<<source_image->sID<<"-> "<<sink_image->sID<<"is ::::"<<temp<<std::endl;
////							model::cout<<"temp1 for case 1 "<<source_image->sID<<"-> "<<sink_image->sID<<"is ::::"<<temp1<<std::endl;
//
//							if (temp)
//							{
//							model::cout<<redColor<<"Entered the case 1 "<<source_image->sID<<" ->"<<sink_image->sID<<std::endl;
//							VectorDim N = (*link.second->loopLinks().begin())->loop()->glidePlane->unitNormal;
//							// Direction along which the new node should be placed
//							VectorDim X=((-source_image->get_P()+sink_image->get_P()).cross(N)).normalized();
//							//							VectorDim X=link.second->burgers();
//							// Make sure the node is placed inside the simulation domain
//							if (!(X.dot(link.second->source->bndNormal().normalized()) > 0)) X = -X;
//							if (temp1)			//this was for multiple loop case
//							X=X*10;
//							else
//							X=X*(10);
//							// Create a new node based on two virtual nodes
//							std::shared_ptr<NodeType> midNode(new NodeType(this,X+0.5*(source_image->get_P()+sink_image->get_P())
//									,0.5*(source_image->get_V()+sink_image->get_V()),1));
//							std::cout<<"Creating the mid node"<<midNode->sID<<std::endl;
//							// Create a container to store the nodes and the direction of the loop to be inserted
//							DislocationLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{source_image,midNode,sink_image},
//									link.second->burgers(),N,source_image->get_P(),
//									(*link.second->grains().begin())->grainID);
//							}
						}
						else if(  sink_image->loops().size()!=0
								&& source_image->loops().size()==0)
						{
							// Loop over the links associated with the virtual node that has a loop
//							std::cout<<"The node to which expansion is taking place is"<<link.second->source->imageBoundaryNode()->sID<<" at position "<<link.second->source->imageBoundaryNode()->get_P()<<std::endl;
							for(const auto& pair: sink_image->linksByLoopID())
							{
								VectorDim N1 = (*pair.second.begin())->loop()->rightHandedNormal();
								VectorDim B1 = (*pair.second.begin())->loop()->burgers();
								VectorDim N = (*link.second->loopLinks().begin())->loop()->rightHandedNormal();
								VectorDim B = (*link.second->loopLinks().begin())->loop()->burgers();
								model::cout<<magentaColor<<"The value of glide normal condition for source expansion is"<<((N-N1).norm())<<std::endl;
								model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
								model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
								model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
								model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
//								if ((N-N1).norm()<=FLT_EPSILON && (B-B1).norm()<=FLT_EPSILON )
//								{
//									std::cout<<"New loop needs to be merged on source"<<std::endl; //TODO not checked for the slip system for node expansion for the periodic loops
//								}
								// Find the link which is a boundary segment
								if (((N-N1).norm()<=FLT_EPSILON && (B-B1).norm()<=FLT_EPSILON) || ((N+N1).norm()<=FLT_EPSILON && (B+B1).norm()<=FLT_EPSILON))
								{std::cout<<"The source image node to which expansion is taking place is"<<source_image->sID<<" at position "<<source_image->get_P()<<std::endl;
//								adding the condition using loop based approach
								for (const auto& pair2: pair.second)
								{
									if (!(pair2)->pLink->isBoundarySegment() && !(pair2)->pLink->hasZeroBurgers())
									{
										model::cout<<redColor<<"Expanding for the loop with source "<<(pair2)->pLink->source->sID<<" and "
												"sink "<<(pair2)->pLink->sink->sID<<std::endl;
										model::cout<<redColor<<"Checking if boundary segment "<<(pair2)->pLink->isBoundarySegment()<<std::endl;
										model::cout<<redColor<<"Checking if zero burgers segment"<<(pair2)->pLink->hasZeroBurgers()<<std::endl;


										this->expand((pair2)->pLink->source->sID,(pair2)->pLink->sink->sID,source_image);
										//													model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
										//													model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
										//													model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
										//													model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
										break;
									}
								}

//								if ((*pair.second.begin())->pLink->isBoundarySegment())
//								{
//									// Change the link of the other segment by using expand
//									//TODO check if the node has the same slip plane normal (nodes corresponding to rbegin should have real nodes which belong to the same loop as the real node corresponding to new node)
//
//									this->expand((*pair.second.rbegin())->pLink->source->sID,(*pair.second.rbegin())->pLink->sink->sID,source_image);
//									model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
//									model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
//									model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
//									model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
//									break;
//								}
//								else if ((*pair.second.rbegin())->pLink->isBoundarySegment())
//								{
//									this->expand((*pair.second.begin())->pLink->source->sID,(*pair.second.begin())->pLink->sink->sID,source_image);
//									model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
//									model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
//									model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
//									model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
//									break;
//								}
								}
							}

							// Second case - Expand the existing loop on the virtual node instead of creating a new loop
							//                                    // Check is made to see whether one of the node has any loop associated with it
							//                                    if(link.second->source->virtualBoundaryNode()->loops().size()==0)
							//                                    {
							//
							//                                    }
							//
							//                                    else
							//                                    {
							//
							//                                    }

						}
						else if(  sink_image->loops().size()==0
								&& source_image->loops().size()!=0)
						{
							// Same procedure done for the other case
//							std::cout<<"The node to which expansion is taking place is"<<link.second->sink->imageBoundaryNode()->sID<<" at position "<<link.second->source->imageBoundaryNode()->get_P()<<std::endl;
							for(const auto& pair: source_image->linksByLoopID())
							{
								VectorDim N1 = (*pair.second.begin())->loop()->rightHandedNormal();
								VectorDim B1 = (*pair.second.begin())->loop()->burgers();
								VectorDim N = (*link.second->loopLinks().begin())->loop()->rightHandedNormal();
								VectorDim B = (*link.second->loopLinks().begin())->loop()->burgers();
								model::cout<<magentaColor<<"The value for sink expansion of glide normal condition at sink "<<((N-N1).norm())<<std::endl;
								model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
								model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
								model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
								model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
								//									if ((N-N1).norm()==FLT_EPSILON && (B-B1).norm()==FLT_EPSILON )
									//									{
									//										std::cout<<"New loop needs to be merged on sink"<<std::endl;
								//									}
								if (((N-N1).norm()<=FLT_EPSILON && (B-B1).norm()<=FLT_EPSILON) || ((N+N1).norm()<=FLT_EPSILON && (B+B1).norm()<=FLT_EPSILON))
								{std::cout<<"The sink_image node to which expansion is taking place is"<<sink_image->sID<<" at position "<<sink_image->get_P()<<std::endl;
								//								adding the condition using loop based approach
								for (const auto& pair2: pair.second)
								{
									if (!(pair2)->pLink->isBoundarySegment() && !(pair2)->pLink->hasZeroBurgers())
									{
										model::cout<<redColor<<"Expanding for the loop with source "<<(pair2)->pLink->source->sID<<" and "
												"sink "<<(pair2)->pLink->sink->sID<<std::endl;
										model::cout<<redColor<<"Checking if boundary segment "<<(pair2)->pLink->isBoundarySegment()<<std::endl;
										model::cout<<redColor<<"Checking if zero burgers segment"<<(pair2)->pLink->hasZeroBurgers()<<std::endl;


										this->expand((pair2)->pLink->source->sID,(pair2)->pLink->sink->sID,sink_image);
										//													model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
										//													model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
										//													model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
										//													model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
										break;
									}
								}
//								if ((*pair.second.begin())->pLink->isBoundarySegment())
//								{
//									this->expand((*pair.second.rbegin())->pLink->source->sID,(*pair.second.rbegin())->pLink->sink->sID,sink_image);
//									model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
//									model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
//									model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
//									model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
//									break;
//								}
//								else if ((*pair.second.rbegin())->pLink->isBoundarySegment())
//								{
//									this->expand((*pair.second.begin())->pLink->source->sID,(*pair.second.begin())->pLink->sink->sID,sink_image);
//									model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
//									model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
//									model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
//									model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
//									break;
//								}
								}
							}
						}
//						else if(  sink_image->loops().size()!=0
//														&& source_image->loops().size()!=0)
//						{
//							//this case is not working properly
//							model::cout<<greenBoldColor<<"Entering the condition where both source and sinks images have loops "<<std::endl<<" source node"<<source_image->sID<<std::endl;
//							VectorDim N = (*link.second->loopLinks().begin())->loop()->glidePlane->unitNormal;
//							// Direction along which the new node should be placed
//							VectorDim X=((-source_image->get_P()+sink_image->get_P()).cross(N)).normalized();
//							//							VectorDim X=link.second->burgers();
//							// Make sure the node is placed inside the simulation domain
//							if (!(X.dot(link.second->source->bndNormal().normalized()) > 0)) X = -X;
//							X=X*(10);
//							// Create a new node based on two virtual nodes
//							std::shared_ptr<NodeType> midNode(new NodeType(this,X+0.5*(source_image->get_P()+sink_image->get_P())
//									,0.5*(source_image->get_V()+sink_image->get_V()),1));
//							//                                    std::cout<<"Creating the mid node"<<midNode->sID<<std::endl;
//							// Create a container to store the nodes and the direction of the loop to be inserted
//							DislocationLoopVector.emplace_back(std::vector<std::shared_ptr<NodeType>>{source_image,midNode,sink_image},
//									link.second->burgers(),N,source_image->get_P(),
//									(*link.second->grains().begin())->grainID);
//
//						}

						else
						{

							model::cout<<greenBoldColor<<"The loop corresponding to Not considered case is "<<std::endl<<" source node"<<source_image->sID<<std::endl;
							model::cout<<greenBoldColor<<" sink node"<<sink_image->sID<<std::endl;
								assert(0 && "CASE NOT CONSIDERED SO FAR");


						}
					}
				}
					//Here the condition for inserting multiple loops is needed to be inserted
					else if (link.second->sink->imageBoundaryNodeContainer().size()==3 || link.second->source->imageBoundaryNodeContainer().size()==3)
					{

							model::cout<<"Inserting multiple loops into the system for real node "<<std::endl;
//							model::cout<<"sink image node container size "<<link.second->sink->imageBoundaryNodeContainer().size()<<std::endl<<"source "
//									"image node container size "<<link.second->source->imageBoundaryNodeContainer().size();
							for (const auto& imnodeiter:((link.second->sink->imageBoundaryNodeContainer().size()==3)?link.second->sink->imageBoundaryNodeContainer():link.second->source->imageBoundaryNodeContainer()))
							{
								model::cout<<magentaColor<<"entered the loop for the case key "<<imnodeiter.first<<std::endl;
								model::cout<<magentaColor<<"Corresponding image node is "<<imnodeiter.second->sID<<std::endl;
								model::cout<<magentaColor<<"printing the linksbyloopid size "<<imnodeiter.second->linksByLoopID().size()<<std::endl;
//								model::cout<<magentaColor<<"printing the imnode iter.first "<<imnodeiter.first<<std::endl;

//								if (imnodeiter.second->linksByLoopID().size()==0 && imnodeiter.first<=3 )
								if (imnodeiter.second->loops().size()==0 && imnodeiter.first<=3 )
								{

									//start inserting the loops
									for (const auto& rneighnodeiter:(link.second->sink->imageBoundaryNodeContainer().size()==3?link.second->sink->neighbors():link.second->source->neighbors()))
									{	model::cout<<redColor<<"entering the loop before the case for checking the neighbors "<<std::endl;
										model::cout<<redColor<<"Neighbor id is "<<std::get<0>(rneighnodeiter.second)->sID<<std::endl;
										model::cout<<redColor<<"Corresponding image container size is"<<std::get<0>(rneighnodeiter.second)->imageBoundaryNodeContainer().size()<<std::endl;
										if(std::get<0>(rneighnodeiter.second)->isImageBoundaryNode )//case if neighbor node is image node
										{
											model::cout<<redColor<<"Corresponding master node is"<<std::get<0>(rneighnodeiter.second)->masterNode->sID<<std::endl;
											if (std::get<0>(rneighnodeiter.second)->masterNode->imageBoundaryNodeContainer().begin()->first==imnodeiter.first)
											{
												for(const auto& pair1: ((std::get<0>(rneighnodeiter.second)->masterNode)->linksByLoopID()))
												{
													VectorDim N1 = (*pair1.second.begin())->loop()->rightHandedNormal();
													VectorDim B1 = (*pair1.second.begin())->loop()->burgers();
													VectorDim N = (*link.second->loopLinks().begin())->loop()->rightHandedNormal();
													VectorDim B = (*link.second->loopLinks().begin())->loop()->burgers();
													if (((N-N1).norm()<=FLT_EPSILON && (B-B1).norm()<=FLT_EPSILON) || ((N+N1).norm()<=FLT_EPSILON && (B+B1).norm()<=FLT_EPSILON))
													{
														//std::cout<<"The node to which expansion is taking place is"<<sink_image->sID<<" at position "<<source_image->get_P()<<std::endl;
														for (const auto& pair2: pair1.second)
														{
															if (!(pair2)->pLink->isBoundarySegment() && !(pair2)->pLink->hasZeroBurgers())
															{
																model::cout<<redColor<<"Expanding for the loop with source "<<(pair2)->pLink->source->sID<<" and "
																		"sink "<<(pair2)->pLink->sink->sID<<std::endl;
																model::cout<<redColor<<"Checking if boundary segment "<<(pair2)->pLink->isBoundarySegment()<<std::endl;
																model::cout<<redColor<<"Checking if zero burgers segment"<<(pair2)->pLink->hasZeroBurgers()<<std::endl;


																this->expand((pair2)->pLink->source->sID,(pair2)->pLink->sink->sID,imnodeiter.second);
																//													model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
																//													model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
																//													model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
																//													model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
																break;
															}
														}
													}

												}
											}

										}


//										else if (!std::get<0>(rneighnodeiter.second)->isImageBoundaryNode)
										else if (!std::get<0>(rneighnodeiter.second)->isImageBoundaryNode)
										{	model::cout<<blueColor<<"entered the loop for the case for checking the neighbors with ID "<<std::get<0>(rneighnodeiter.second)->sID<<std::endl;
											model::cout<<redColor<<"Inserting loop for the multiple loop case"<<std::endl;
//											model::cout<<redColor<<" Value of the the value corresponding to the key "<<std::get<0>(rneighnodeiter.second)->imageBoundaryNodeContainer().begin()->first <<" rneighnodeiter"
//													".second is "<<std::get<0>(rneighnodeiter.second)->imageBoundaryNodeContainer().begin()->second->linksByLoopID().size()<<std::endl;
											//model::cout<<redColor<<"Data type for loop id is "<<typeid(std::get<0>(rneighnodeiter.second)->imageBoundaryNodeContainer().begin()->second->linksByLoopID()).name()<<std::endl;
											for (const auto& pairtemp: std::get<0>(rneighnodeiter.second)->imageBoundaryNodeContainer())
											{
												if (imnodeiter.first==pairtemp.first)
												{
											for(const auto& pair1: pairtemp.second->linksByLoopID())
											{
												VectorDim N1 = (*pair1.second.begin())->loop()->rightHandedNormal();
												VectorDim B1 = (*pair1.second.begin())->loop()->burgers();
												VectorDim N = (*link.second->loopLinks().begin())->loop()->rightHandedNormal();
												VectorDim B = (*link.second->loopLinks().begin())->loop()->burgers();
												if (((N-N1).norm()<=FLT_EPSILON && (B-B1).norm()<=FLT_EPSILON) || ((N+N1).norm()<=FLT_EPSILON && (B+B1).norm()<=FLT_EPSILON))
												{
													//std::cout<<"The node to which expansion is taking place is"<<sink_image->sID<<" at position "<<source_image->get_P()<<std::endl;
													for (const auto& pair2: pair1.second)
													{
														if (!(pair2)->pLink->isBoundarySegment() && !(pair2)->pLink->hasZeroBurgers())
														{
															model::cout<<redColor<<"Expanding for the loop with source "<<(pair2)->pLink->source->sID<<" and "
																	"sink "<<(pair2)->pLink->sink->sID<<std::endl;
															model::cout<<redColor<<"Checking if boundary segment "<<(pair2)->pLink->isBoundarySegment()<<std::endl;
															model::cout<<redColor<<"Checking if zero burgers segment"<<(pair2)->pLink->hasZeroBurgers()<<std::endl;


															this->expand((pair2)->pLink->source->sID,(pair2)->pLink->sink->sID,imnodeiter.second);
															//													model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
															//													model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
															//													model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
															//													model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
															break;
														}
													}
												}

											}
												}
										}
											//assert(0 && "Simulation stopped by the programmer after insertion of multiple loops");

										}
										else
										{
											model::cout<<"Neighbor node id is "<< std::get<0>(rneighnodeiter.second)->sID<<" Corresponding master node is "<<std::get<0>(rneighnodeiter.second)->masterNode->sID<<std::endl;
											assert (0 && "neighbor node is an image boundary node");
										}
									}
									model::cout<<magentaColor<<"Loops should be inserted for the key "<<imnodeiter.first<<" for the real node ";//<<(link.second->sink->imageBoundaryNodeContainer().size()==4)?link.second->sink->sID:link.second->source->sID;

								}
								model::cout<<blueColor<<"printing the linksbyloopid size "<<imnodeiter.second->linksByLoopID().size()<<std::endl;
								model::cout<<blueColor<<"printing the imnode iter.first "<<imnodeiter.first<<std::endl;
								model::cout<<blueColor<<"image node neighbor size "<<imnodeiter.second->neighbors().size()<<std::endl;
								if (imnodeiter.second->loops().size()==0 && imnodeiter.first>3 )
								{
									//here insert for the diagonally opposite ends

									for(const auto& pair0: ((link.second->sink->imageBoundaryNodeContainer().size()==3)?link.second->sink->imageBoundaryNodeContainer():link.second->source->imageBoundaryNodeContainer()))
									{
										if(pair0.first<=3)
										{
											for(const auto& pair: pair0.second->neighbors())
											{	model::cout<<redColor<<"pair.second "<<std::get<0>(pair.second)->sID<<" is imageboundary node "<< std::get<0>(pair.second)->isImageBoundaryNode<<std::endl;
											model::cout<<redColor<<"pair.second "<<std::get<0>(pair.second)->sID<<" image node container size is "<<std::get<0>(pair.second)->imageBoundaryNodeContainer().size()<<std::endl;
											if(!(std::get<0>(pair.second)->isImageBoundaryNode) && std::get<0>(pair.second)->imageBoundaryNodeContainer().size()==1)
											{
//												if (imnodeiter.second->loops().size()==0)
												if(pair0.first!=std::get<0>(pair.second)->imageBoundaryNodeContainer().begin()->first)
												{for (const auto& pair1:std::get<0>(pair.second)->imageBoundaryNodeContainer().begin()->second->linksByLoopID())
												{
													model::cout<<redColor<<"entered the diagonally opposite loop condition"<<std::endl;
													VectorDim N1 = (*pair1.second.begin())->loop()->rightHandedNormal();
													VectorDim B1 = (*pair1.second.begin())->loop()->burgers();
													VectorDim N = (*link.second->loopLinks().begin())->loop()->rightHandedNormal();
													VectorDim B = (*link.second->loopLinks().begin())->loop()->burgers();
													if (((N-N1).norm()<=FLT_EPSILON && (B-B1).norm()<=FLT_EPSILON) || ((N+N1).norm()<=FLT_EPSILON && (B+B1).norm()<=FLT_EPSILON))
													{
														//std::cout<<"The node to which expansion is taking place is"<<sink_image->sID<<" at position "<<source_image->get_P()<<std::endl;
														//Expand using the loops taking into consideration zero burgers vector segment

														for (const auto& pair2: pair1.second)
														if (!(pair2)->pLink->isBoundarySegment() && !(pair2)->pLink->hasZeroBurgers())
														{
															model::cout<<redColor<<"Expanding for the loop with source "<<(pair2)->pLink->source->sID<<" and "
																	"sink "<<(pair2)->pLink->sink->sID<<std::endl;
															model::cout<<redColor<<"Checking if boundary segment "<<(pair2)->pLink->isBoundarySegment()<<std::endl;
															model::cout<<redColor<<"Checking if zero burgers segment"<<(pair2)->pLink->hasZeroBurgers()<<std::endl;


//															if ((pair2)->pLink->source->isBoundaryNode() && (pair2)->pLink->sink->isBoundaryNode())  //this is done so as to remove anomalous expansion at the diagonal
															{
																this->expand((pair2)->pLink->source->sID,(pair2)->pLink->sink->sID,imnodeiter.second);
																//													model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
																//													model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
																//													model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
																//													model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
																break;
															} //Check here for the condition
//															this->expand((pair2)->pLink->source->sID,(pair2)->pLink->sink->sID,imnodeiter.second);
////															model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
////															model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
////															model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
////															model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
//															break;


														}

//														if ((*pair1.second.begin())->pLink->isBoundarySegment() || (*pair1.second.begin())->pLink->hasZeroBurgers())
//														{
//															this->expand((*pair1.second.rbegin())->pLink->source->sID,(*pair1.second.rbegin())->pLink->sink->sID,imnodeiter.second);
//															//													model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
//															//													model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
//															//													model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
//															//													model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
//															break;
//														}
//														else if ((*pair1.second.rbegin())->pLink->isBoundarySegment() || (*pair1.second.begin())->pLink->hasZeroBurgers())
//														{
//															this->expand((*pair1.second.begin())->pLink->source->sID,(*pair1.second.begin())->pLink->sink->sID,imnodeiter.second);
//															//													model::cout<<redBoldColor<<"N1 is equal to "<<N1<<std::endl;
//															//													model::cout<<greenBoldColor<<"N is equal to "<<N<<std::endl;
//															//													model::cout<<redBoldColor<<"B1 is equal to "<<B1<<std::endl;
//															//													model::cout<<greenBoldColor<<"B is equal to "<<B<<std::endl;
//															break;
//														}
													}
												}

												}
											}
											else
											{

												model::cout<<"Node in consideration is  "<<imnodeiter.second->sID<<std::endl;
												model::cout<<"Case for the neighbor node being on the edge not implemented::"<<std::endl;
//													assert(0 && "Case for the neighbor node also on the edge not considered");
											}

										}
									}
									}
									//}

								}
							}

//

					}
					else
					{
						model::cout<<"Loop in consideration is "<<link.second->source->sID<<" -> "<<link.second->sink->sID<<std::endl;
						model::cout<<"Source->sink image container size "<<link.second->source->imageBoundaryNodeContainer().size()<<" -> "
								<<link.second->sink->imageBoundaryNodeContainer().size()<<std::endl;

						if (link.second->source->imageBoundaryNodeContainer().size()>0 && link.second->sink->imageBoundaryNodeContainer().size() ) //Check this condition for once.
						assert(0 &&"Case for including loops through multiple images not considered so far:::");
					}

				}
			}
			for(const auto& tup : DislocationLoopVector)
			{// Insert the new loops
				this->insertLoop(std::get<0>(tup),std::get<1>(tup),std::get<2>(tup),std::get<3>(tup),std::get<4>(tup));
				//if (runID>390)
//				model::cout<<redColor<<(std::get<0>(tup)).at(0)->sID<<"->"<<(std::get<0>(tup)).at(1)->sID<<"->"<<(std::get<0>(tup)).at(2)->sID<<std::endl;
			}

			model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;

			break;
		}

		default:
		{
			break;
		}
		}



	}



public:

	const int simulationType; // 0 finite sample, no FEM, // 1 finite sample FEM, // 2 periodic prism
	int timeIntegrationMethod;
	int maxJunctionIterations;
	//        long int runID;
	double totalTime;
	double dt;
	double vMax;
	//        size_t Nsteps;
	MatrixDimD _plasticDistortionFromVelocities;
	MatrixDimD _plasticDistortionFromAreas;
	MatrixDimD _plasticDistortionRateFromVelocities;
	MatrixDimD _plasticDistortionRateFromAreas;
	int ddSolverType;
	bool computeDDinteractions;
	int crossSlipModel;
	bool use_boundary;
	unsigned int use_bvp;
	bool useVirtualExternalLoops;
	bool use_externalStress;
	bool use_extraStraightSegments;
	ExternalLoadControllerType extStressController;
	std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> ssdeq;
	std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> straightSegmentsDeq;
	int  outputFrequency;
	bool outputBinary;
	bool outputGlidePlanes;
	bool outputSpatialCells;
	//        bool outputPKforce;
	bool outputElasticEnergy;
	bool outputMeshDisplacement;
	bool outputFEMsolution;
	bool outputDislocationLength;
	bool outputPlasticDistortion;
	bool outputPlasticDistortionRate;
	bool outputQuadraturePoints;
	bool outputLinkingNumbers;
	bool outputLoopLength;
	bool outputSegmentPairDistances;
	unsigned int _userOutputColumn;
	bool use_stochasticForce;
	int dislocationImages_x;
	int dislocationImages_y;
	int dislocationImages_z;
	double surfaceAttractionDistance;
	std::vector<VectorDim,Eigen::aligned_allocator<VectorDim>> periodicShifts;
	bool computePlasticDistortionIncrementally;
	std::string folderSuffix;

	/**********************************************************************/
	DislocationNetwork(int& argc, char* argv[],
			SimplicialMesh<dim>& _mesh,
			Polycrystal<dim>& _poly,
			BVPsolver<dim,2>& _bvpSolver,
			long int& runID) :
				/* init */ DislocationNetworkBase<dim>(_mesh,_poly,_bvpSolver)
				/* init */,simulationType(TextFileParser("inputFiles/DD.txt").readScalar<int>("simulationType",true))
				/* init */,timeIntegrationMethod(TextFileParser("inputFiles/DD.txt").readScalar<int>("timeIntegrationMethod",true))
				/* init */,maxJunctionIterations(TextFileParser("inputFiles/DD.txt").readScalar<int>("maxJunctionIterations",true))
				//        /* init */,runID(TextFileParser("inputFiles/DD.txt").readScalar<int>("startAtTimeStep",true)),
				/* init */,totalTime(0.0),
				/* init */ dt(0.0),
				/* init */ vMax(0.0),
				//        /* init */ Nsteps(TextFileParser("inputFiles/DD.txt").readScalar<size_t>("Nsteps",true)),
				/* init */ _plasticDistortionFromVelocities(MatrixDimD::Zero()),
				/* init */ _plasticDistortionFromAreas(MatrixDimD::Zero()),
				/* init */ _plasticDistortionRateFromVelocities(MatrixDimD::Zero()),
				/* init */ _plasticDistortionRateFromAreas(MatrixDimD::Zero()),
				/* init */ ddSolverType(TextFileParser("inputFiles/DD.txt").readScalar<int>("ddSolverType",true)),
				/* init */ computeDDinteractions(TextFileParser("inputFiles/DD.txt").readScalar<int>("computeDDinteractions",true)),
				/* init */ crossSlipModel(TextFileParser("inputFiles/DD.txt").readScalar<int>("crossSlipModel",true)),
				/* init */ use_boundary(true),
				/* init */ use_bvp(TextFileParser("inputFiles/DD.txt").readScalar<int>("use_bvp",true)),
				/* init */ useVirtualExternalLoops(TextFileParser("inputFiles/DD.txt").readScalar<int>("useVirtualExternalLoops",true)),
				/* init */ use_externalStress(false),
				/* init */ use_extraStraightSegments(false),
				/* init */ outputFrequency(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFrequency",true)),
				/* init */ outputBinary(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputBinary",true)),
				/* init */ outputGlidePlanes(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputGlidePlanes",true)),
				/* init */ outputSpatialCells(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputSpatialCells",true)),
				//        /* init */ outputPKforce(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPKforce",true)),
				/* init */ outputElasticEnergy(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputElasticEnergy",true)),
				/* init */ outputMeshDisplacement(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputMeshDisplacement",true)),
				/* init */ outputFEMsolution(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFEMsolution",true)),
				/* init */ outputDislocationLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputDislocationLength",true)),
				/* init */ outputPlasticDistortion(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortion",true)),
				/* init */ outputPlasticDistortionRate(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortionRate",true)),
				/* init */ outputQuadraturePoints(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputQuadraturePoints",true)),
				/* init */ outputLinkingNumbers(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLinkingNumbers",true)),
				/* init */ outputLoopLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLoopLength",true)),
				/* init */ outputSegmentPairDistances(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputSegmentPairDistances",true)),
				/* init */ _userOutputColumn(3),
				/* init */ use_stochasticForce(TextFileParser("inputFiles/DD.txt").readScalar<int>("use_stochasticForce",true)),
				/* init */ dislocationImages_x(use_boundary? TextFileParser("inputFiles/DD.txt").readScalar<int>("dislocationImages_x",true) : 0),
				/* init */ dislocationImages_y(use_boundary? TextFileParser("inputFiles/DD.txt").readScalar<int>("dislocationImages_y",true) : 0),
				/* init */ dislocationImages_z(use_boundary? TextFileParser("inputFiles/DD.txt").readScalar<int>("dislocationImages_z",true) : 0),
				/* init */ surfaceAttractionDistance(TextFileParser("inputFiles/DD.txt").readScalar<double>("surfaceAttractionDistance",true)),
				/* init */ computePlasticDistortionIncrementally(TextFileParser("inputFiles/DD.txt").readScalar<int>("computePlasticDistortionIncrementally",true)),
				/* init */ folderSuffix("")
				{

		// Some sanity checks
		//            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");

		// Initialize static variables
		LinkType::initFromFile("inputFiles/DD.txt");
		NodeType::initFromFile("inputFiles/DD.txt");
		LoopType::initFromFile("inputFiles/DD.txt");
		DislocationNetworkComponentType::initFromFile("inputFiles/DD.txt");
		DislocationStressBase<dim>::initFromFile("inputFiles/DD.txt");
		DDtimeIntegrator<0>::initFromFile("inputFiles/DD.txt");
		DislocationCrossSlip<DislocationNetworkType>::initFromFile("inputFiles/DD.txt");
		SpatialCellObserverType::setCellSize(TextFileParser("inputFiles/DD.txt").readScalar<double>("dislocationCellSize",true));
		DislocationDisplacement<dim>::use_multipole=TextFileParser("inputFiles/DD.txt").readScalar<double>("use_DisplacementMultipole",true);
		DislocationStress<dim>::use_multipole=TextFileParser("inputFiles/DD.txt").readScalar<double>("use_StressMultipole",true);
		DislocationEnergy<dim>::use_multipole=TextFileParser("inputFiles/DD.txt").readScalar<double>("use_EnergyMultipole",true);


		//            EDR.readScalarInFile(fullName.str(),"use_DisplacementMultipole",DislocationDisplacement<dim>::use_multipole);
		//            EDR.readScalarInFile(fullName.str(),"use_StressMultipole",DislocationStress<dim>::use_multipole);
		//            EDR.readScalarInFile(fullName.str(),"use_EnergyMultipole",DislocationEnergy<dim>::use_multipole);
		//            EDR.readScalarInFile(fullName.str(),"quadPerLength",LinkType::quadPerLength); // quadPerLength
		//            LinkType::quadPerLength=TextFileParser("inputFiles/DD.txt").readScalar<double>("quadPerLength",true);
		//            assert((LinkType::quadPerLength)>=0.0 && "quadPerLength MUST BE >= 0.0");

		//                        EDR.readScalarInFile(fullName.str(),"coreSize",StressField::a); // core-width
		//            StressField::a=TextFileParser("inputFiles/DD.txt").readScalar<double>("coreSize",true)
		//            assert((StressField::a)>0.0 && "coreSize MUST BE > 0.");
		//            StressField::a2=StressField::a*StressField::a;

		//            EDR.readScalarInFile(fullName.str(),"stochasticForceSeed",stochasticForceSeed);
		int stochasticForceSeed=TextFileParser("inputFiles/DD.txt").readScalar<int>("stochasticForceSeed",true);
		if(stochasticForceSeed<0)
		{
			StochasticForceGenerator::init(std::chrono::system_clock::now().time_since_epoch().count());
		}
		else
		{
			StochasticForceGenerator::init(stochasticForceSeed);
		}

		if(argc>1)
		{
			folderSuffix=argv[1];
			std::cout<<"folderSuffix="<<folderSuffix<<std::endl;
		}

		ParticleSystemType::initMPI(argc,argv);





		if(outputPlasticDistortion)
		{
			_userOutputColumn+=9;
		}
		if(outputPlasticDistortionRate)
		{
			_userOutputColumn+=9;
		}

		if(outputDislocationLength)
		{
			_userOutputColumn+=4;
		}

		// IO
		io().read("./","DDinput.txt",runID);

		// Set up periodic shifts
		const VectorDim meshDimensions(use_boundary? (this->mesh.xMax()-this->mesh.xMin()).eval() : VectorDim::Zero());
		model::cout<<"meshDimensions="<<meshDimensions.transpose()<<std::endl;
		for(int i=-dislocationImages_x;i<=dislocationImages_x;++i)
		{
			for(int j=-dislocationImages_y;j<=dislocationImages_y;++j)
			{
				for(int k=-dislocationImages_z;k<=dislocationImages_z;++k)
				{
					const Eigen::Array<int,dim,1> cellID((Eigen::Array<int,dim,1>()<<i,j,k).finished());
					periodicShifts.push_back((meshDimensions.array()*cellID.template cast<double>()).matrix());

				}
			}
		}
		model::cout<<"periodic shift vectors:"<<std::endl;
		for(const auto& shift : periodicShifts)
		{
			model::cout<<shift.transpose()<<std::endl;

		}

		if(simulationType==PERIODIC)
		{// in case of periodic simulation, mesh must be commensurate to lattice periodicity
			assert(this->poly.grains().size()==1 && "ONLY SINGLE-CRYSTAL PERIODIC SIMULATIONS SUPPORTED.");
			model::cout<<"Checking that mesh and lattice are commensurate"<<std::endl;
			const VectorDim meshSize(this->mesh.xMax()-this->mesh.xMin());
			for(int d=0;d<dim;++d)
			{
				VectorDim v(VectorDim::Zero());
				v(d)=1*meshSize(d);
				LatticeVector<dim> lv(v,this->poly.grains().begin()->second);
			}
		}

		// Initializing configuration
		move(0.0);	// initial configuration
				}


	/**********************************************************************/
	void singleStep(const long int& runID)
	{
		//! A simulation step consists of the following:
		//            model::cout<<blueBoldColor<< "runID="<<runID<<" (of "<<Nsteps<<")"
		//            /*                    */<< ", time="<<totalTime
		//            /*                    */<< ": nodes="<<this->nodes().size()
		//            /*                    */<< ", segments="<<this->links().size()
		//            /*                    */<< ", loopSegments="<<this->loopLinks().size()
		//            /*                    */<< ", loops="<<this->loops().size()
		//            /*                    */<< ", components="<<this->components().size()
		//            /*                    */<< defaultColor<<std::endl;

		//! 1- Check that all nodes are balanced
		//            checkBalance();

		//! 2 - Update quadrature points
		updateQuadraturePoints();
		updateStressStraightSegments();

		for(auto& loop : this->loops()) // TODO: PARALLELIZE THIS LOOP
		{// copmute slipped areas and right-handed normal
			loop.second->update();
		}
		updatePlasticDistortionFromAreas();

		//! 3- Calculate BVP correction
		updateLoadControllers(runID);

		//#ifdef DislocationNucleationFile
		//            if(use_bvp && !(runID%use_bvp))
		//            {
		//                nucleateDislocations(); // needs to be called before updateQuadraturePoints()
		//                updateQuadraturePoints();
		//            }
		//#endif

		computeNodaVelocities(runID);


		//! 4- Solve the equation of motion


		//! 5- Compute time step dt (based on max nodal velocity) and increment totalTime
		// make_dt();

		if(outputElasticEnergy)
		{
			typedef typename DislocationParticleType::ElasticEnergy ElasticEnergy;
			this->template computeNeighborField<ElasticEnergy>();
		}

		//! 6- Output the current configuration before changing it
		//            output(runID);
		io().output(runID);


		//! 7- Moves DislocationNodes(s) to their new configuration using stored velocity and dt
		move(dt);

		//! 8- Update accumulated quantities (totalTime and plasticDistortion)
		totalTime+=dt;
		updatePlasticDistortionRateFromVelocities();


		//! 9- Contract segments of zero-length
		//            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();

		//! 10- Cross Slip (needs upated PK force)
		DislocationCrossSlip<DislocationNetworkType>(*this);


		//                        GrainBoundaryTransmission<DislocationNetworkType>(*this).transmit();
		GrainBoundaryTransmission<DislocationNetworkType>(*this).directTransmit();

		//
		//            GrainBoundaryDissociation<DislocationNetworkType>(*this).dissociate();

		//            poly.reConnectGrainBoundarySegments(*this); // this makes stressGauss invalid, so must follw other GB operations


		//! 11- detect loops that shrink to zero and expand as inverted loops
		//            DislocationNetworkRemesh<DislocationNetworkType>(*this).loopInversion(dt);

		//! 12- Form Junctions
		DislocationJunctionFormation<DislocationNetworkType>(*this).formJunctions(maxJunctionIterations,DDtimeIntegrator<0>::dxMax);

		//            // Remesh may contract juncitons to zero lenght. Remove those juncitons:
		//            DislocationJunctionFormation<DislocationNetworkType>(*this).breakZeroLengthJunctions();

		updateVirtualBoundaryLoops();

		DislocationNetworkRemesh<DislocationNetworkType>(*this).remesh(runID);

		//! 13- Node redistribution



		//            mergeLoopsAtNodes();

		//            DislocationInjector<DislocationNetworkType>(*this).insertRandomStraightDislocation();

		//! 9- Contract segments of zero-length
		//            DislocationNetworkRemesh<DislocationNetworkType>(*this).contract0chordSegments();

		//! 14- If BVP solver is not used, remove DislocationSegment(s) that exited the boundary
		//            removeBoundarySegments();

		//            removeSmallComponents(3.0*dx,4);

		//            make_bndNormals();

		//! 16 - Increment runID counter
		//            ++runID;     // increment the runID counter
	}

	//        const VectorDim periodicVector(const Eigen::Array<int,dim,1>& cellID) const
	//        {
	//            return (meshDimensions.array()*cellID.template cast<double>()).matrix();
	//        }

	/**********************************************************************/
	const EshelbyInclusionContainerType& eshelbyInclusions() const
	{
		return *this;
	}

	EshelbyInclusionContainerType& eshelbyInclusions()
	{
		return *this;
	}


	/**********************************************************************/
	bool contract(std::shared_ptr<NodeType> nA,
			std::shared_ptr<NodeType> nB)
	{
		return DislocationNodeContraction<DislocationNetworkType>(*this).contract(nA,nB);
	}

	/**********************************************************************/
	//incorporating imageNodes as well
//	bool contract(std::shared_ptr<NodeType> nA,
//				std::shared_ptr<NodeType> nB)
//		{
//		if (nA->imageBoundaryNodeContainer().size()==1 && nB->imageBoundaryNodeContainer().size()==1)
//		{
////
//			const bool imageNodeCont=(*this).contractSecond(nA->imageBoundaryNodeContainer().at(0),nB->imageBoundaryNodeContainer().at(0));
//			return (*this).contractSecond(nA,nB);
//		}
//		else
//			return DislocationNodeContraction<DislocationNetworkType>(*this).contract(nA,nB);
//
//		}



	/**********************************************************************/
	const double& get_dt() const
	{/*!\returns the current time step increment dt
	 */
		return dt;
	}

	/**********************************************************************/
	void set_dt(const double& dt_in,const double& vMax_in)
	{/*!\param[in] dt_in
	 * Sets the time increment dt to dt_in
	 */
		dt=dt_in;
		vMax=vMax_in;
	}

	/**********************************************************************/
	const double& get_totalTime() const
	{/*! The elapsed simulation time step in dimensionless units
	 */
		return totalTime;
	}

	/**********************************************************************/
	void updateStressStraightSegments()
	{

		const auto t0= std::chrono::system_clock::now();
		model::cout<<"		Collecting StressStraight objects: "<<std::flush;

		straightSegmentsDeq.clear();
		//                std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> straightSegmentsDeq;
		size_t currentSize=0;
		//            if(computeDDinteractions)
		//            {
		for(const auto& link : this->networkLinks())
		{
			link.second->addToStressStraight(straightSegmentsDeq);
		}

		currentSize=straightSegmentsDeq.size();

		const VectorDim meshSize(this->mesh.xMax()-this->mesh.xMin());

		for(int i=-dislocationImages_x;i<=dislocationImages_x;++i)
		{
			for(int j=-dislocationImages_y;j<=dislocationImages_y;++j)
			{
				for(int k=-dislocationImages_z;k<=dislocationImages_z;++k)
				{

					const Eigen::Matrix<int,3,1> cellID((Eigen::Matrix<int,3,1>()<<i,j,k).finished());

					if( cellID.squaredNorm()!=0) //skip current cell
							{
						for (size_t c=0;c<currentSize;++c)
						{
							const VectorDim P0=straightSegmentsDeq[c].P0+(meshSize.array()*cellID.cast<double>().array()).matrix();
							const VectorDim P1=straightSegmentsDeq[c].P1+(meshSize.array()*cellID.cast<double>().array()).matrix();

							straightSegmentsDeq.emplace_back(P0,P1,straightSegmentsDeq[c].b);
						}
							}
				}
			}
		}


		//            }
		model::cout<< straightSegmentsDeq.size()<<" straight segments ("<<currentSize<<"+"<<straightSegmentsDeq.size()-currentSize<<" images)"<<std::flush;
		model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;

	}

	/**********************************************************************/
	void assembleAndSolve(const long int& runID)
	{/*! Performs the following operatons:
	 */
#ifdef _OPENMP
		const size_t nThreads = omp_get_max_threads();
#else
		const size_t nThreads = 1;
#endif

		//! -1 Compute the interaction StressField between dislocation particles


		if(corder==0)
		{// For straight segments use analytical expression of stress field
			//                const auto t0= std::chrono::system_clock::now();
			//                model::cout<<"		Collecting StressStraight objects: "<<std::flush;
			//
			////                std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> straightSegmentsDeq;
			//                size_t currentSize=0;
			//                if(computeDDinteractions)
			//                {
			//                    for(const auto& link : this->networkLinks())
			//                    {
			//                        link.second->addToStressStraight(straightSegmentsDeq);
			//                    }
			//
			//                    currentSize=straightSegmentsDeq.size();
			//
			//                    const VectorDim meshSize(this->mesh.xMax()-this->mesh.xMin());
			//
			//                    for(int i=-dislocationImages_x;i<=dislocationImages_x;++i)
			//                    {
			//                        for(int j=-dislocationImages_y;j<=dislocationImages_y;++j)
			//                        {
			//                            for(int k=-dislocationImages_z;k<=dislocationImages_z;++k)
			//                            {
			//
			//                                const Eigen::Matrix<int,3,1> cellID((Eigen::Matrix<int,3,1>()<<i,j,k).finished());
			//
			//                                if( cellID.squaredNorm()!=0) //skip current cell
			//                                {
			//                                    for (size_t c=0;c<currentSize;++c)
			//                                    {
			//                                        const VectorDim P0=straightSegmentsDeq[c].P0+(meshSize.array()*cellID.cast<double>().array()).matrix();
			//                                        const VectorDim P1=straightSegmentsDeq[c].P1+(meshSize.array()*cellID.cast<double>().array()).matrix();
			//
			//                                        straightSegmentsDeq.emplace_back(P0,P1,straightSegmentsDeq[c].b);
			//                                    }
			//                                }
			//                            }
			//                        }
			//                    }
			//
			//
			//                }
			//                model::cout<< straightSegmentsDeq.size()<<" straight segments ("<<currentSize<<"+"<<straightSegmentsDeq.size()-currentSize<<" images)"<<std::flush;
			//                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
			//
			const auto t1= std::chrono::system_clock::now();
			model::cout<<"		Computing analytical stress field at quadrature points ("<<nThreads<<" threads) "<<std::flush;
#ifdef _OPENMP
			//                const size_t nThreads = omp_get_max_threads();
			EqualIteratorRange<typename NetworkLinkContainerType::iterator> eir(this->links().begin(),this->links().end(),nThreads);
			//                SOMETHING WRONG HERE? CPU USE SAYS 100% ONLY?

#pragma omp parallel for
			for(size_t thread=0;thread<eir.size();thread++)
			{
				for (auto& linkIter=eir[thread].first;linkIter!=eir[thread].second;linkIter++)
				{
					linkIter->second->assemble(straightSegmentsDeq);
				}
			}
#else
	for (auto& linkIter : this->links())
	{
		linkIter.second->assemble(straightSegmentsDeq);
	}
#endif
	model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;






	//#ifdef _OPENMP
	//#pragma omp parallel for
	//#endif
	//                for (unsigned int k=0; k<this->particleSystem().size();++k)
	//                {
	//                    if(this->particleSystem()[k].template fieldPointBase<StressField>().enabled)
	//                    {
	//                        this->particleSystem()[k].template fieldPointBase<StressField>() += StressField::addSourceContribution(this->particleSystem()[k],straightSegmentsDeq);
	//                    }
	//                }
	//                model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
		}
		else
		{// For curved segments use quandrature integration of stress field
			//                assert(0 && "RE-ENABLE THIS. New Material class is incompatible with old implmentation using static objects");

			const auto t0= std::chrono::system_clock::now();

			if(computeDDinteractions)
			{

				if(dislocationImages_x!=0 || dislocationImages_y!=0 || dislocationImages_z!=0)
				{
					assert(0 && "FINISH HERE");
				}

				model::cout<<"		Computing numerical stress field at quadrature points ("<<nThreads<<" threads)..."<<std::flush;
				if (use_extraStraightSegments)
				{
					this->template computeNeighborField<StressField>(ssdeq);
				}
				else
				{
					this->template computeNeighborField<StressField>();
				}
				model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
			}

			//! -2 Loop over DislocationSegments and assemble stiffness matrix and force vector
			const auto t1= std::chrono::system_clock::now();
			model::cout<<"		Computing segment stiffness matrices and force vectors ("<<nThreads<<" threads)..."<<std::flush;
			typedef void (LinkType::*LinkMemberFunctionPointerType)(void); // define type of Link member function
			assert(0 && "THIS CASE MUST BE REWORKED, since LinkType::assemble DOES NOT EXIST ANYMORE");
			//LinkMemberFunctionPointerType Lmfp(&LinkType::assemble); // Lmfp is a member function pointer to Link::assemble
			//this->parallelExecute(Lmfp);
			model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]."<<defaultColor<<std::endl;


		}

		//! -3 Loop over DislocationSubNetworks, assemble subnetwork stiffness matrix and force vector, and solve
		const auto t2= std::chrono::system_clock::now();
		model::cout<<"		Assembling NetworkComponents and solving "<<std::flush;



		switch (ddSolverType)
		{
		case 1: // iterative solver
		{
			model::cout<<"(MINRES "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP // MINRES is not multi-threaded. So parallelize loop over NetworkComponents.
#pragma omp parallel for
			for (unsigned int k=0;k<this->components().size();++k)
			{
				auto snIter(this->components().begin());
				std::advance(snIter,k);
				DislocationNetworkComponentType(*snIter->second).iterativeSolve();
			}
#else
	for (const auto& networkComponent : this->components())
	{
		DislocationNetworkComponentType(*networkComponent.second).iterativeSolve();
	}
#endif
break;
		}

		case 2: // direct solver
		{
#ifdef _MODEL_PARDISO_SOLVER_
			model::cout<<"(PardisoLDLT "<<nThreads<<" threads)..."<<std::flush;
			for (const auto& networkComponent : this->components())
			{
				DislocationNetworkComponentType(*networkComponent.second).directSolve();
			}
#else
	model::cout<<"(SimplicialLDLT "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP // SimplicialLDLT is not multi-threaded. So parallelize loop over NetworkComponents.
#pragma omp parallel for
	for (unsigned int k=0;k<this->components().size();++k)
	{
		auto snIter(this->components().begin());
		std::advance(snIter,k);
		DislocationNetworkComponentType(*snIter->second).directSolve();
	}
#else
	for (const auto& networkComponent : this->components())
	{
		DislocationNetworkComponentType(*networkComponent.second).directSolve();
	}
#endif
#endif
break;
		}

		default: // lumped solver
		{
			model::cout<<"(lumpedSolver "<<nThreads<<" threads)..."<<std::flush;
#ifdef _OPENMP
#pragma omp parallel for
			for (unsigned int k=0;k<this->components().size();++k)
			{
				auto snIter(this->components().begin());
				std::advance(snIter,k);
				DislocationNetworkComponentType(*snIter->second).lumpedSolve(runID);
			}
#else
	for (const auto& networkComponent : this->components())
	{
		DislocationNetworkComponentType(*networkComponent.second).lumpedSolve(runID);
	}
#endif
break;
		}
		}

		//            if(DislocationNetworkComponentType::use_directSolver)
		//            {
		//
		//
		//            }
		//            else // iterative solver
		//            {
		//
		//            }

		model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t2)).count()<<" sec]."<<defaultColor<<std::endl;
	}

	/**********************************************************************/
	void updateQuadraturePoints()
	{
		if(corder>0)
		{// quadrature points only used for curved segments
			model::cout<<"		Updating quadrature points... "<<std::flush;
			const auto t0=std::chrono::system_clock::now();

			// Clear DislocationParticles
			this->clearParticles(); // this also destroys all cells

			// Populate DislocationParticles
			for(auto& linkIter : this->links())
			{
				linkIter.second->updateQuadraturePoints(*this);
			}
			model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
		}
	}


	/**********************************************************************/
	void move(const double & dt_in)
	{/*! Moves all nodes in the DislocationNetwork using the stored velocity and current dt
	 */
		model::cout<<"		Moving DislocationNodes (dt="<<dt_in<< ")... "<<std::flush;
		const auto t0= std::chrono::system_clock::now();
		for (auto& nodeIter : this->nodes())
		{
			nodeIter.second->move(dt_in,DDtimeIntegrator<0>::dxMax);
		}
		model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
	}

	//        /**********************************************************************/
	//        void runSteps() __attribute__ ((deprecated))
	//        {/*! Runs Nsteps simulation steps
	//          */
	//            const auto t0= std::chrono::system_clock::now();
	//            while (runID<Nsteps)
	//            {
	//                model::cout<<std::endl; // leave a blank line
	//                singleStep();
	//            }
	//            updateQuadraturePoints(); // necessary if quadrature data are necessary in main
	//            model::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<Nsteps<< " simulation steps completed in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" [sec]"<<defaultColor<<std::endl;
	//        }

	/**********************************************************************/
	void updatePlasticDistortionRateFromVelocities()
	{
		if(computePlasticDistortionIncrementally)
		{
			_plasticDistortionRateFromVelocities.setZero();
			for (const auto& linkIter : this->links())
			{
				_plasticDistortionRateFromVelocities+= linkIter.second->plasticDistortionRate();
			}
			_plasticDistortionFromVelocities += _plasticDistortionRateFromVelocities*dt;
		}
	}

	/**********************************************************************/
	void updatePlasticDistortionFromAreas()
	{
		if(!computePlasticDistortionIncrementally)
		{
			const MatrixDimD old(_plasticDistortionFromAreas);
			_plasticDistortionFromAreas.setZero();
			for(const auto& loop : this->loops())
			{
				_plasticDistortionFromAreas+= loop.second->plasticDistortion();
			}
			_plasticDistortionRateFromAreas=(_plasticDistortionFromAreas-old)/dt;
		}
	}

	/**********************************************************************/
	const MatrixDimD& plasticDistortionRate() const
	{
		return computePlasticDistortionIncrementally? _plasticDistortionRateFromVelocities : _plasticDistortionRateFromAreas;
	}

	/**********************************************************************/
	const MatrixDimD& plasticDistortion() const
	{
		return computePlasticDistortionIncrementally? _plasticDistortionFromVelocities : _plasticDistortionFromAreas;
	}

	/**********************************************************************/
	MatrixDimD plasticStrainRate() const
	{/*!\returns the plastic distortion rate tensor generated by this
	 * DislocaitonNetwork.
	 */
		//const MatrixDimD temp(plasticDistortionRate());
		return (plasticDistortionRate()+plasticDistortionRate().transpose())*0.5;
	}

	/**********************************************************************/
	std::tuple<double,double,double,double> networkLength() const
        		{/*!\returns the total line length of the DislocationNetwork. The return
        		 * value is a tuple, where the first value is the length of bulk glissile
        		 * dislocations, the second value is the length of bulk sessile
        		 * dislocations, and the third value is the length accumulated on
        		 * the mesh boundary.
        		 */
		double bulkGlissileLength(0.0);
		double bulkSessileLength(0.0);
		double boundaryLength(0.0);
		double grainBoundaryLength(0.0);

		for(auto& loop : this->loops())
		{
			for(const auto& loopLink : loop.second->links())
			{
				if(!loopLink.second->pLink->hasZeroBurgers())
				{
					if(loopLink.second->pLink->isBoundarySegment())
					{
						boundaryLength+=loopLink.second->pLink->chord().norm();
					}
					else if(loopLink.second->pLink->isGrainBoundarySegment())
					{
						grainBoundaryLength+=loopLink.second->pLink->chord().norm();
					}
					else
					{
						if(loopLink.second->pLink->isSessile())
						{
							bulkSessileLength+=loopLink.second->pLink->chord().norm()/loopLink.second->pLink->loopLinks().size();
						}
						else
						{
							bulkGlissileLength+=loopLink.second->pLink->chord().norm()/loopLink.second->pLink->loopLinks().size();
						}
					}
				}
			}
		}

		//
		//
		//            for (const auto& linkIter : this->links())
			//            {
			//                const double temp(linkIter.second->arcLength());
			//                if(linkIter.second->isBoundarySegment())
				//                {
				//                    boundaryLength+=temp;
				//                }
			//                else
				//                {
				//                    if(linkIter.second->isSessile())
					//                    {
					//                        bulkSessileLength+=temp;
					//                    }
		//                    else
		//                    {
		//                        bulkGlissileLength+=temp;
		//
		//                    }
		//                }
		//
		//            }
		return std::make_tuple(bulkGlissileLength,bulkSessileLength,boundaryLength,grainBoundaryLength);
        		}

	//        /**********************************************************************/
	//        const long int& runningID() const
	//        {/*! The current simulation step ID.
	//          */
	//            return runID;
	//        }


	/**********************************************************************/
	const unsigned int& userOutputColumn()  const
	{
		return _userOutputColumn;
	}

	/**********************************************************************/
	MatrixDimD stress(const VectorDim& P) const
	{/*!\param[in] P position vector
	 * \returns The stress field generated by the DislocationNetwork at P
	 *
	 * Note:
	 */

		assert(false && "REWORK THIS FOR STRAIGHT SEGMENTS");

		SingleFieldPoint<StressField> fieldPoint(P,true);
		this->template computeField<SingleFieldPoint<StressField>,StressField>(fieldPoint);
		return fieldPoint.field();
	}

	/**********************************************************************/
	std::pair<bool,const Simplex<dim,dim>*> pointIsInsideMesh(const VectorDim& P0, const Simplex<dim,dim>* const guess) const
        		{/*!\param[in] P0 position vector
        		 * \param[in] guess pointer of the Simplex where the search starts
        		 * \returns true if P0 is inside the mesh
        		 */
		std::pair<bool,const Simplex<dim,dim>*> temp(true,NULL);
		if (use_boundary)
		{
			temp=this->mesh.searchWithGuess(P0,guess);
		}
		return temp;
        		}

	/**********************************************************************/
	const ParticleSystemType& particleSystem() const
	{
		return *this;
	}

	/**********************************************************************/
	ParticleSystemType& particleSystem()
	{
		return *this;
	}

	/**********************************************************************/
	DislocationNetworkIOType io()
	{
		return DislocationNetworkIOType(*this,folderSuffix);
	}

	/**********************************************************************/
	DislocationNetworkIOType io() const
	{
		return DislocationNetworkIOType(*this,folderSuffix);
	}

};

}
#endif


