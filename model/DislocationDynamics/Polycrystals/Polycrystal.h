/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Polycrystal_H_
#define model_Polycrystal_H_

#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>
#include <utility>
#include <tuple>
#include <map>
#include <vector>
#include <deque>
#include <tuple>
#include <chrono>
#include <random>
#include <Eigen/Core>
#include <model/MPI/MPIcout.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundary.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/LatticeMath/LatticeVector.h>
#include <model/DislocationDynamics/Materials/PeriodicElement.h>
#include <model/IO/EigenDataReader.h>
#include <model/IO/SequentialOutputFile.h>
#include <model/DislocationDynamics/ElasticFields/StressStraight.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundaryType.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/DislocationDynamics/Materials/MaterialsLibrary.h>

namespace model
{
    
    
    
    template <int dim>
    class Polycrystal :
    /* base */ private std::map<size_t,Grain<dim>>,
    /* base */ private std::map<std::pair<size_t,size_t>,GrainBoundary<dim>>,
    /* base */ public std::vector<StressStraight<dim> >
    {
        
//        static constexpr int dim=NetworkType::dim;
        typedef SimplicialMesh<dim> SimplicialMeshType;
        typedef MeshRegion<Simplex<dim,dim> > MeshRegionType;
        typedef MeshRegionObserver<MeshRegionType> MeshRegionObserverType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Grain<dim> GrainType;
        typedef GrainBoundary<dim> GrainBoundaryType;
        
        static constexpr PeriodicElement<13,Isotropic> Al=PeriodicElement<13,Isotropic>();
        static constexpr PeriodicElement<28,Isotropic> Ni=PeriodicElement<28,Isotropic>();
        static constexpr PeriodicElement<29,Isotropic> Cu=PeriodicElement<29,Isotropic>();
        static constexpr PeriodicElement<74,Isotropic>  W=PeriodicElement<74,Isotropic>();
        
        

        
        unsigned int materialZ;
        

        
    public:
        
        MaterialsLibrary materials;
        
        const SimplicialMeshType& mesh;
        
        /**********************************************************************/
        Polycrystal(const SimplicialMeshType& mesh_in) :
        /* init */ materialZ(0),
        /* init */ mesh(mesh_in)
        {
            model::cout<<"Creating Polycrystal"<<std::endl;
            
            
        }
        
        /**********************************************************************/
//        template<typename NetworkType>
        void init(GlidePlaneObserver<dim>& dn, const std::string& fullName)
        {
            
            
            EigenDataReader EDR;
            
            EDR.readScalarInFile(fullName,"material",materialZ); // material by atomic number Z
            
            // Construct Grains
            model::cout<<"Creating Grains"<<std::endl;
            for(const auto& rIter : mesh.regions())
            {
                
                // Read grain orientation
                
                int C2GfromSlipSystem=-1;
                EDR.readScalarInFile(fullName,"C2G"+std::to_string(rIter.second->regionID)+"fromSlipSystem",C2GfromSlipSystem);
                
                Eigen::Matrix<double,dim,dim> C2Gtemp(Eigen::Matrix<double,dim,dim>::Identity());
                
                if(C2GfromSlipSystem>=0)
                {
                
                    grains().emplace(std::piecewise_construct,
                                     std::forward_as_tuple(rIter.second->regionID),
                                     std::forward_as_tuple(*(rIter.second),
                                                           materialZ,
                                                           C2Gtemp));
                    
                    if(C2GfromSlipSystem<grains().at(rIter.second->regionID).slipSystems().size())
                    {
                        const SlipSystem& slipSystem(grains().at(rIter.second->regionID).slipSystems()[C2GfromSlipSystem]);
                        C2Gtemp.row(2)=slipSystem.unitNormal;
                        C2Gtemp.row(0)=slipSystem.s.cartesian().normalized();
                        C2Gtemp.row(1)=C2Gtemp.row(2).cross(C2Gtemp.row(0));
                        grains().at(rIter.second->regionID).rotate(C2Gtemp);
                    }
                    else
                    {
                        model::cout<<"C2G"+std::to_string(rIter.second->regionID)+"fromSlipSystem"<<" exceedes the number of SlipSystems in Grain "<<rIter.second->regionID<<" ("<<grains().at(rIter.second->regionID).slipSystems().size()<<")."<<std::endl;
                        assert(0 && "WRONG SLIPSYSTEM ID");
                    }
                    
                }
                else
                {
                    EDR.readMatrixInFile(fullName,"C2G"+std::to_string(rIter.second->regionID),C2Gtemp); // crystal-to-global orientation
                    
                    grains().emplace(std::piecewise_construct,
                                     std::forward_as_tuple(rIter.second->regionID),
                                     std::forward_as_tuple(*(rIter.second),
                                                           materialZ,
                                                           C2Gtemp));
                }
                
                
//                if()
//                {
//                
//                }
                
//                for(int i=0;i<3;i++)
//                {
//                    double c2gNorm(C2Gtemp.row(i).norm());
//                    for(int j=0;j<3;j++)
//                    {
//                        assert(C2Gtemp(i,j)-std::round(C2Gtemp(i,j))==0 && "User must pass C2G matrix in integer form to limit rounding errors");
//                        C2Gtemp(i,j)=C2Gtemp(i,j)/c2gNorm;
//                    }
//                }
                

                
            }
            
            // Construct GrainsBoundaries
            model::cout<<"Creating GrainsBoundaries"<<std::endl;
            grainBoundaryDislocations().clear();
//            int fileID=1;
            for(const auto& rgnBnd : mesh.regionBoundaries())
            {
                
//                const auto pr=rgnBnd.second.unsortedBoundary();
//                std::ofstream myfile1;
//                myfile1.open ("rb"+std::to_string(fileID)+".txt");
//                //myfile << "Writing this to a file.\n";
//                for(const auto& edge : pr)
//                {
//                    myfile1<<edge->child(0).P0.transpose()<<" "<<edge->child(1).P0.transpose()<<"\n";
//                }
//                myfile1.close();
                
                
                grainBoundaries().emplace(std::piecewise_construct,
                                          std::forward_as_tuple(rgnBnd.first),
                                          std::forward_as_tuple(rgnBnd.second,
                                                                grain(rgnBnd.first.first),
                                                                grain(rgnBnd.first.second),
                                                                dn,
                                                                mesh));
                
//                std::ofstream myfile;
//                myfile.open ("gb"+std::to_string(fileID)+".txt");
//                //myfile << "Writing this to a file.\n";
//                for(const auto& pair : grainBoundary(rgnBnd.first.first,rgnBnd.first.second).meshIntersections)
//                {
//                    myfile<<" "<<pair.second.transpose()<<"\n";
//                }
//                myfile.close();
//                fileID++;
            }
            
            
            if(grainBoundaryDislocations().size())
            {
                model::SequentialOutputFile<'B',1>::set_count(0);
                model::SequentialOutputFile<'B',1>::set_increment(1);
                model::SequentialOutputFile<'B',true> stressStraightFile;
                size_t n=0;
                for(const auto& sStraight : grainBoundaryDislocations())
                {
                    stressStraightFile<<n<<"\t"<<sStraight.P0.transpose()<<"\t"<<sStraight.P1.transpose()<<"\t"<<sStraight.b.transpose()<<std::endl;
                    n++;
                }
            }
        }
        
        /**********************************************************************/
        const std::deque<GrainBoundaryType>& grainBoundaryTypes() const
        {
            switch (materialZ)
            {
                case Al.Z:
                    return PeriodicElement<Al.Z,Isotropic>::grainBoundaryTypes;
                    break;
                case Ni.Z:
                    return PeriodicElement<Ni.Z,Isotropic>::grainBoundaryTypes;
                    break;
                case Cu.Z:
                    return PeriodicElement<Cu.Z,Isotropic>::grainBoundaryTypes;
                    break;
                default:
                    assert(0 && "grainBoundaryTypes not implemented.");
                    break;
            }
        }
        
        /**********************************************************************/
        GrainType& grain(const size_t& k)
        {
            return grains().at(k);
        }
        
        /**********************************************************************/
        const GrainType& grain(const size_t& k) const
        {
            return grains().at(k);
        }
        
        /**********************************************************************/
        const std::map<size_t,GrainType>& grains() const
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<size_t,GrainType>& grains()
        {
            return *this;
        }
        
        /**********************************************************************/
        std::map<std::pair<size_t,size_t>,GrainBoundaryType>& grainBoundaries()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::map<std::pair<size_t,size_t>,GrainBoundaryType>& grainBoundaries() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const GrainBoundaryType& grainBoundary(const size_t& i,
                                                        const size_t& j) const
        {
            assert(i!=j && "GrainBoundary IDs cannot be the same.");
            return (i<j)? grainBoundaries().at(std::make_pair(i,j)) : grainBoundaries().at(std::make_pair(j,i));
        }
        
        /**********************************************************************/
        GrainBoundaryType& grainBoundary(const size_t& i,
                                          const size_t& j)
        {
            assert(i!=j && "GrainBoundary IDs cannot be the same.");
            return (i<j)? grainBoundaries().at(std::make_pair(i,j)) : grainBoundaries().at(std::make_pair(j,i));
        }
        
        /**********************************************************************/
        LatticeVectorType latticeVectorFromPosition(const VectorDim& p,
                                                    const Simplex<dim,dim>* const guess) const
        {
            const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
            assert(temp.first && "Position not found in mesh");
            return grain(temp.second->region->regionID).latticeVector(p);
        }
        
        /**********************************************************************/
        LatticeVectorType latticeVectorFromPosition(const VectorDim& p) const
        {
            return latticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDim& p,
                                                                        const Simplex<dim,dim>* const guess) const
        {
            const std::pair<bool,const Simplex<dim,dim>*> temp(mesh.searchWithGuess(p,guess));
            assert(temp.first && "Position not found in mesh");
            return grain(temp.second->region->regionID).reciprocalLatticeVector(p);
        }
        
        /**********************************************************************/
        ReciprocalLatticeVectorType reciprocalLatticeVectorFromPosition(const VectorDim& p) const
        {
            return reciprocalLatticeVectorFromPosition(p,&(mesh.simplices().begin()->second));
        }
        
        /**********************************************************************/
        std::vector<StressStraight<dim>>& grainBoundaryDislocations()
        {
            return *this;
        }
        
        /**********************************************************************/
        const std::vector<StressStraight<dim>>& grainBoundaryDislocations() const
        {
            return *this;
        }
        
        /**********************************************************************/
        VectorDim randomPoint() const
        {
            std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
            std::uniform_real_distribution<double> distribution(0.0,1.0);

            VectorDim P0(VectorDim::Zero());
            
            P0 << mesh.xMin(0)+distribution(generator)*(mesh.xMax(0)-mesh.xMin(0)),
            /* */ mesh.xMin(1)+distribution(generator)*(mesh.xMax(1)-mesh.xMin(1)),
            /* */ mesh.xMin(2)+distribution(generator)*(mesh.xMax(2)-mesh.xMin(2));
            
            return P0;
            
        }
        
        /**********************************************************************/
        std::pair<LatticeVector<dim>,int> randomLatticePointInMesh() const
        {
            const VectorDim P0=randomPoint();
            auto searchResult=mesh.search(P0);
            if(searchResult.first)
            {// point inside
                const LatticeVector<dim> L0 = grain(searchResult.second->region->regionID).snapToLattice(P0);
                searchResult=mesh.searchRegionWithGuess(L0.cartesian(),searchResult.second);
                if(searchResult.first)
                {// point inside
                    return std::make_pair(L0,searchResult.second->region->regionID);
                }
                else
                {
                    return randomLatticePointInMesh();
                }
            }
            else
            {
                return randomLatticePointInMesh();
            }
            
        }
        
        //        /**********************************************************************/
        //        template <typename DislocationNetworkType>
        //        void reConnectGrainBoundarySegments(DislocationNetworkType& DN) const
        //        {
        //            const auto t0= std::chrono::system_clock::now();
        //            model::cout<<"		re-connecting GrainBoundarySegments ("<<std::flush;
        //            std::map<std::pair<size_t,size_t>,LatticeVectorType> disconnectMap;
        //            for (const auto& link : DN.links())
        //            {
        ////                if(link.second->source->isGrainBoundaryNode() && link.second->sink->isGrainBoundaryNode())
        ////                {
        //                    if(link.second->source->rID2()>=0 &&
        //                       link.second->sink->rID2()>=0 &&
        //                       link.second->source->rID2()==link.second->sink->rID2())
        //                    {
        //                        disconnectMap.emplace(std::piecewise_construct,
        //                                              std::forward_as_tuple(link.second->source->sID,link.second->sink->sID),
        //                                              std::forward_as_tuple(link.second->flow));
        //                    }
        ////                }
        //            }
        //
        //            model::cout<<disconnectMap.size()<<std::flush;
        //
        //            for(const auto& segment : disconnectMap)
        //            {
        //                DN.template disconnect<false>(segment.first.first,segment.first.second);
        //                DN.template connect(segment.first.first,segment.first.second,segment.second);
        ////                assert(DN.link(std::get<0>(tuple),std::get<1>(tuple)).second->isGrainBoundarySegment() && "SEGMENT IS NOT A GB SEGMENT");
        //            }
        //            
        //            model::cout<<std::setprecision(3)<<magentaColor<<" reconnected) ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
        //        }
        
    };
    
    
} // end namespace
#endif

