/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po             <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez       <ramirezbrf@gmail.com>
 * Copyright (C) 2011 by Mamdouh Mohamed        <msm07d@fsu.edu>
 * Copyright (C) 2011 by Tamer Crsoby           <tamercrosby@gmail.com>
 * Copyright (C) 2011 by Can Erel               <canerel55@gmail.com>
 * Copyright (C) 2011 by Yinan Cui              <cuiyinan@ucla.edu>
 * Copyright (C) 2017 by Sabyasachi Chatterjee  <sabyasac@andrew.cmu.edu>
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

// template header cpp
// https://www.codeproject.com/Articles/48575/How-to-Define-a-Template-Class-in-a-h-File-and-Imp

#ifdef _OPENMP
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE // disable Eigen Internal openmp Parallelization
#endif

#include <vector>
#include <chrono>
#include <map>
#include <memory>
#include <numbers>

#include <Eigen/Dense>

#include <TerminalColors.h>
#include <DislocationNetworkTraits.h>
#include <DislocationDynamicsModule.h>
#include <LoopNetwork.h>
//#include <DislocationNetworkComponent.h>
#include <DislocationLoop.h>
#include <DislocationLoopNode.h>
#include <DislocationLoopLink.h>
#include <DislocationNode.h>
#include <DislocationSegment.h>
#include <Hermite.h>

#include <DislocationNetworkRemesh.h>
#include <DislocationJunctionFormation.h>
// #include <DislocationCrossSlip.h>
////#include <Material.h>
//#include <UniqueOutputFile.h>
#include <DislocationNetworkIO.h>
//#include <DislocationParticle.h>
#include <DislocationFieldBase.h>
////#include <ParticleSystem.h>
//
////#include <SingleFieldPoint.h>
#include <DDtimeIntegrator.h>
//#include <EqualIteratorRange.h>
////#include <BoundingLineSegments.h>
//#include <GrainBoundaryTransmission.h>
////#include <GrainBoundaryDissociation.h>
#include <DefectiveCrystalParameters.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <BVPsolver.h>
#include <ExternalLoadControllerBase.h>
#include <GlidePlaneModule.h>

#include <DislocationNodeContraction.h>
#include <EshelbyInclusionBase.h>
#include <SphericalInclusion.h>
#include <PolyhedronInclusion.h>

#include <DDconfigIO.h>
#include <DislocationGlideSolver.h>
#include <CrossSlipModels.h>

//#include <TextFileParser.h>
////#include <DisplacementPoint.h>


////#include <ExternalLoadController.h>
//#include <DislocationInjector.h>
//#include <PeriodicDislocationLoop.h>
////#include <PeriodicLoopObserver.h>
//
////#include <PeriodicDislocationSuperLoop.h>
////#include <PlanarDislocationSuperLoop.h>
//
//#ifdef _MODEL_GREATWHITE_
//#include <MooseSolution.h>
//#endif

#ifndef NDEBUG
#define VerboseDislocationNetwork(N,x) if(verboseDislocationNetwork>=N){std::cout<<x;}
#else
#define VerboseDislocationNetwork(N,x)
#endif

namespace model
{
    template <int dim, short unsigned int corder>
    class DislocationNetwork :public LoopNetwork<DislocationNetwork<dim,corder> >
    /*                     */,public std::map<size_t,std::shared_ptr<EshelbyInclusionBase<dim>>>
    /*                     */,public std::map<size_t,PolyhedronInclusionNodeIO<dim>>

    {
        
    public:
        
        typedef TypeTraits<DislocationNetwork<dim,corder>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef typename TraitsType::VectorLowerDim VectorLowerDim;
        typedef typename TraitsType::MatrixDim MatrixDim;
        typedef DislocationNetworkIO<LoopNetworkType> DislocationNetworkIOType;
        typedef std::map<size_t,std::shared_ptr<EshelbyInclusionBase<dim>>> EshelbyInclusionContainerType;
        typedef std::map<size_t,PolyhedronInclusionNodeIO<dim>> PolyhedronInclusionNodeContainerType;
        typedef BVPsolver<dim,2> BvpSolverType;
        typedef typename BvpSolverType::FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;

        
        static int verboseDislocationNetwork;


    public:

        const DefectiveCrystalParameters& simulationParameters;
        const SimplicialMesh<dim>& mesh;
        const Polycrystal<dim>& poly;
        GlidePlaneFactory<dim> glidePlaneFactory;
        std::shared_ptr<PeriodicGlidePlaneFactory<dim>> periodicGlidePlaneFactory;
        const std::unique_ptr<BVPsolver<dim,2>>& bvpSolver;
        const std::unique_ptr<ExternalLoadControllerBase<dim>>& externalLoadController;
        const std::vector<VectorDim>& periodicShifts;
        DislocationNetworkRemesh<LoopNetworkType> networkRemesher;
        DislocationJunctionFormation<DislocationNetwork<dim,corder>> junctionsMaker;
        const std::shared_ptr<BaseCrossSlipModel<DislocationNetwork<dim,corder>>> crossSlipModel;
        DislocationCrossSlip<DislocationNetwork<dim,corder>> crossSlipMaker;
        DislocationNodeContraction<LoopNetworkType> nodeContractor;
        DDtimeIntegrator timeIntegrator;
        std::shared_ptr<StochasticForceGenerator> stochasticForceGenerator;
        DislocationNetworkIO<LoopNetworkType> networkIO;
        int ddSolverType;
        bool computeDDinteractions;
        int  outputFrequency;
        bool outputBinary;
        bool outputMeshDisplacement;
        bool outputFEMsolution;
        bool outputQuadraturePoints;
        bool outputLinkingNumbers;
        bool outputLoopLength;
        bool outputSegmentPairDistances;
        const bool computeElasticEnergyPerLength;
        double alphaLineTension;
        std::set<const LoopNodeType*> danglingBoundaryLoopNodes;
        const bool use_velocityFilter;
        const double velocityReductionFactor;
        const int verboseDislocationNode;
        bool capMaxVelocity;

        DislocationNetwork(const DefectiveCrystalParameters& _simulationParameters,
                           const SimplicialMesh<dim>& _mesh,
                           const Polycrystal<dim>& _poly,
                           const std::unique_ptr<BVPsolver<dim,2>>& _bvpSolver,
                           const std::unique_ptr<ExternalLoadControllerBase<dim>>& _externalLoadController,
                           const std::vector<VectorDim>& _periodicShifts,
                           long int& runID);
        
        void setConfiguration(const DDconfigIO<dim>&);
        MatrixDim plasticDistortionRate() const;
        MatrixDim plasticDistortion() const;
        MatrixDim plasticStrain() const;
        std::map<std::pair<int,int>,double> slipSystemPlasticDistortion() const;
        MatrixDim plasticStrainRate() const;
        void updateGeometry();//
        void updateRates();//
        DislocationNetworkIOType& io();
        const DislocationNetworkIOType& io() const;
        std::tuple<double,double,double,double> networkLength() const;
        const EshelbyInclusionContainerType& eshelbyInclusions() const;
        EshelbyInclusionContainerType& eshelbyInclusions();
        const PolyhedronInclusionNodeContainerType& polyhedronInclusionNodes() const;
        PolyhedronInclusionNodeContainerType& polyhedronInclusionNodes();
        VectorDim displacement(const VectorDim& x) const;
        void displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>& fieldPoints) const;
        MatrixDim stress(const VectorDim& x) const;
        void stress(std::deque<FEMfaceEvaluation<ElementType,dim,dim>>& fieldPoints) const;
        void assembleGlide(const long int& runID, const double& maxVelocity);
        void solveGlide(const long int& runID);
        void moveGlide(const double & dt_in);
        void storeSingleGlideStepDiscreteEvents(const long int& runID);
        void executeSingleGlideStepDiscreteEvents(const long int& runID);
        void updateBoundaryNodes();
        bool contract(std::shared_ptr<NetworkNodeType> nA,std::shared_ptr<NetworkNodeType> nB);
        
    };
    
}
#endif
