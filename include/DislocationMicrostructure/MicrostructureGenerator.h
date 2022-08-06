/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGenerator_H_
#define model_MicrostructureGenerator_H_

#include <filesystem>
#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <assert.h>
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <limits>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
//#include <PlaneMeshIntersection.h>
//#include <DislocationNodeIO.h>
//#include <DislocationLoopIO.h>
//#include <DislocationLoopLinkIO.h>
//#include <DislocationLoopNodeIO.h>
#include <DDtraitsIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationLinkingNumber.h>
#include <TextFileParser.h>
#include <DislocationInjector.h>
#include <MeshBoundarySegment.h>
//#include <ConfinedDislocationObject.h>
#include <GlidePlaneModule.h>
#include <MeshModule.h>
#include <Plane.h>
#include <MicrostructureGeneratorBase.h>


namespace model
{

    

    class MicrostructureGenerator : public std::map<std::string,MicrostructureGeneratorBase* const>
    {
        constexpr static int dim=MicrostructureGeneratorBase::dim;
        typedef typename  MicrostructureGeneratorBase::VectorDimD VectorDimD;
//        typedef typename  MicrostructureGeneratorBase::VectorDimI VectorDimI;
        typedef typename  MicrostructureGeneratorBase::DislocationLoopType DislocationLoopType;


        DDtraitsIO traitsIO;
        DDconfigIO<3> configIO;
        DDauxIO<dim> auxIO;

    public:
        
//        const TextFileParser parser;
        const bool outputBinary;
//        const int meshID;
        const std::set<int> periodicFaceIDs;
//        const std::string meshFilename;
        const SimplicialMesh<dim> mesh;
        const double minSize;
        const double maxSize;
        const Polycrystal<dim> poly;
        GlidePlaneFactory<dim> glidePlaneFactory;
        PeriodicGlidePlaneFactory<dim> periodicGlidePlaneFactory;
        std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>> uniqueNetworkNodeMap;
        
        MicrostructureGenerator(const std::string& folderName);
        
        const DDconfigIO<3>& config() const;
        const DDauxIO<3>& aux() const;
        size_t insertLoop(const VectorDimD& b,const VectorDimD& unitNormal,const VectorDimD& P0,const size_t& grainID,const DislocationLoopType& loopType);
        size_t insertLoopNode(const size_t& loopID,const VectorDimD& loopNodePos,const size_t& networkNodeID,const VectorDimD& loopNodeShift,const std::pair<short int,short int>& periodicEdgeIDs);
        std::vector<size_t> insertLoopLinks(const size_t& loopID,const std::vector<size_t>& loopNodeIDs);
        size_t insertNetworkNode(const VectorDimD& networkNodePos);
        size_t insertInclusion(const VectorDimD& pos,const double& R, const Eigen::Matrix<double,dim,dim>& eT, const double& vrc,const int&type);
        void writeConfigFiles(const size_t& fileID);
        
        void insertJunctionLoop(const std::vector<VectorDimD>& loopNodePos,
                                const std::shared_ptr<PeriodicGlidePlane<3>>& periodicPlane,
                                const VectorDimD& b,
                                const VectorDimD& unitNormal,
                                const VectorDimD& P0,
                                const size_t& grainID,
                                const DislocationLoopIO<dim>::DislocationLoopType& loopType);

        bool allPointsInGrain(const std::vector<VectorDimD>& points,const int& grainID);


        
    };
//
//    struct PolyPoint
//    {
//
//        PeriodicPlanePatch<3>* periodicPlanePatch() const;
//
//
//
//    };
//
//    class MicrostructureGenerator : public MicrostructureGeneratorBase
//    {
//        constexpr static int dim=3;
//        typedef Eigen::Matrix<double,dim,1> VectorDimD;
//        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
//        typedef LatticeDirection<dim> LatticeDirectionType;
//
//        typedef Eigen::Matrix<double,dim,dim>    MatrixDimD;
//        typedef Eigen::Matrix<long int,dim,dim>    MatrixDimI;
//        typedef PolycrystallineMaterial<dim,Isotropic> MaterialType;
//
//        typedef BoundingMeshSegments<dim> MeshBoundaryContainerType;
//
//
////        static double min(const double& a,const double& b);
////
////
////        static double max(const double& a,const double& b);
//
//
////        static VectorDimD randomOrthogonalUnitVector(VectorDimD v);
////
////
//////        size_t insertLoopNode(const size_t& loopNodeID,const size_t& loopID,const VectorDimD& loopNodePos,const size_t& networkNodeID,const VectorDimD& loopNodeShift,const std::pair<short int,short int>& periodicEdgeIDs);
////
////        bool addSingleLoop(const bool randomizeBurgersSense,
////                           const std::vector<size_t>& nodeIDs,
////                           const std::vector<size_t>& loopNodeIDs,
////                           VectorDimD b,
////                           const VectorDimD& unitNormal,
////                           const VectorDimD& P0,
////                           const int& grainID,
////                           const int& loopType,
////                           const std::vector<VectorDimD>& loopNodeShifts,
////                           const std::vector<std::pair<short int,short int>>& periodicEdgeIDs // default -1 for non-periodic nodes
////                           );
////
////        bool addSingleLoop(const bool randomizeBurgersSense,
////                           const std::vector<VectorDimD>& nodePos,
////                           const std::vector<VectorDimD>& loopNodePos,
////                           VectorDimD b,
////                           const VectorDimD& unitNormal,
////                           const VectorDimD& P0,
////                           const int& grainID,
////                           const int& loopType,
////                           const std::vector<VectorDimD>& loopNodeShifts,
////                           const std::vector<std::pair<short int,short int>>& periodicEdgeIDs // default -1 for non-periodic nodes
////                           );
////
////        bool addSingleLoopwithJunction(const bool randomizeBurgersSense,
////                           const std::vector<VectorDimD>& nodePos,
////                           const std::vector<VectorDimD>& loopNodePos,
////                           VectorDimD b,
////                           const VectorDimD& unitNormal,
////                           const VectorDimD& P0,
////                           const int& grainID,
////                           const int& loopType,
////                           const std::vector<VectorDimD>& loopNodeShifts,
////                           const std::vector<std::pair<short int,short int>>& periodicEdgeIDs // default -1 for non-periodic nodes
////                           );
////
////
////        void addStraightDislocations();
////
////
////        void addFrankReadSources();
////
////
////        void addSingleArmDislocations();
////
////
////        void addPrismaticLoops();
////
////        void addIndividualStraightDislocations();
////
////
////        void addNonPlanarLoops();
//
//
////        void addStatisticallyHomegeneousPeriodicLoops();
//        void addIndividualPeriodicStraighDislocation();
//
////        void addFrankLoops();
////
////
////        void addEshelbyInclusions();
////
////
////        std::mt19937 generator;
////        size_t nodeID;
//////        size_t snID;
////        size_t loopNodeID;
////        size_t loopID;
////        std::list<PeriodicGlidePlane<dim>> periodicGlidePlaneContainer;
////
////        std::deque<std::vector<VectorDimD>> loopPoints;
////        std::deque<VectorDimD> loopBurgers;
////        const bool enforceMonotonicHelicity;
////        double helicity;
//
//    public:
//
////        DDconfigIO<3> configIO;
////        DDauxIO<dim> auxIO;
//
//
//        const bool outputBinary;
////        const int meshID;
//        const std::set<int> periodicFaceIDs;
//        const std::string meshFilename;
//        const SimplicialMesh<dim> mesh;
//        const double minSize;
//        const double maxSize;
//        Polycrystal<dim> poly;
//         GlidePlaneFactory<dim> glidePlaneFactory;
//        PeriodicGlidePlaneFactory<dim> pgpf;
//
////        // Straight Dislocations
////        const double targetStraightDislocationDensity;
////        const double fractionSessileStraightDislocationDensity;
////
////        // Frank-Read sources
////        const double targetFrankReadDislocationDensity;
////        const double FrankReadSizeMean;
////        const double FrankReadSizeStd;
////        const double FrankReadAspectRatioMean;
////        const double FrankReadAspectRatioStd;
////
////        // Single-arm sources
////        const double targetSingleArmDislocationDensity;
////
////        // Prismatic loops
////        const double targetPrismaticLoopDensity;
////
////        // Individual dislocations
////        const std::vector<int> straightDislocationsSlipSystemIDs;
////        const std::vector<double> straightDislocationsAngleFromScrewOrientation;
////        const Eigen::Matrix<double,Eigen::Dynamic,dim> pointsAlongStraightDislocations;
////
////        // Frank  Loops
////        const double targetFrankLoopDensity;
////        const double frankLoopRadiusMean;
////        const double frankLoopRadiusStd;
////        const double frankLoopSides;
////
////        // NonPlanar  Loops
////        const double targetNonPlanarLoopDensity;
////        const double nonPlanarLoopRadiusMean;
////        const double nonPlanarLoopRadiusStd;
////        const double nonPlanarLoopSides;
////
////        // Periodic  Loops
////        const double targetPeriodicLoopDensity;
////        const double periodicLoopRadiusMean;
////        const double periodicLoopRadiusStd;
////        const double periodicLoopSides;
//
//
////        // Junction  Loops
////        const size_t targetJunctionLoops;
////        const double targetJunctionLoopsSize;
//
////        // PlanarDipolar  Loops
////        const double targetPlanarDipolarLoopDensity;
////        const double planarDipolarLoopMean;
////        const double planarDipolarLoopStd;
////        const double planarDipolarLoopAspectRatioMean;
////        const double planarDipolarLoopAspectRatioStd;
//
////        // Irradiation Loops
////        const double targetIrradiationLoopDensity;
//////        const double averageLoopSize;
////        const double irradiationLoopsDiameterLognormalDistribution_M;
////        const double irradiationLoopsDiameterLognormalDistribution_S;
////        const double irradiationLoopsDiameterLognormalDistribution_A;
////        // const double fraction111Loops;  // fraction of [111] glissile loop;
////        // const bool mobile111Loops;
////
////        // SFTs
////        const double targetSFTdensity;
////
////        // Inclusions
////        const std::vector<double> targetInclusionDensities;
////        const std::vector<double> inclusionsDiameterLognormalDistribution_M;
////        const std::vector<double> inclusionsDiameterLognormalDistribution_S;
////        const std::vector<double> inclusionsDiameterLognormalDistribution_A;
////        const Eigen::Matrix<double,Eigen::Dynamic,dim*dim> inclusionsTransformationStrains;
////        const Eigen::Matrix<double,Eigen::Dynamic,dim> inclusionsPatterns;
////        //        const std::vector<double> inclusionsMobilityReduction;
////
////        bool isInclusionsUsed() const;
////
////
////        MicrostructureGenerator(int argc, char* argv[]) ;
////
////
////        void writeConfigFiles(const size_t& fileID);
////
////        void addStackingFaultTetrahedra();
////
////
////        bool allPointsInGrain(const std::vector<VectorDimD>& points,const int& grainID);
////
////
////        void addIrradiationLoopsFCC();
////
////
////        void addIrradiationLoopsHCP();
////
////
////
////        void addIrradiationLoopsBCC();
////
////
////        void addIrradiationLoops();
////
////
////        double deltaHelicity(const std::vector<VectorDimD>& newPoints,
////                             const VectorDimD& newBurgers) const;
////
////
////        std::pair<LatticeVector<dim>,int> randomPointInMesh() const;
////
////        double randomSize();
////
////
////        int randomSign();
//
//
//        static std::map<double,VectorDimD> boundaryProjection(const VectorDimD& P0,
//                                                              const VectorDimD& P1,
//                                                              const VectorDimD& D,
////                                                              const PlaneMeshIntersectionContainerType& pp
//                                                              const MeshBoundaryContainerType& pp);
//
//
//        static std::pair<int,VectorDimD> boundaryProjection(const VectorDimD& P,
//                                                            const VectorDimD& D,
////                                                            const PlaneMeshIntersectionContainerType& pp,
//                                                            const MeshBoundaryContainerType& pp);
//
//
//
//    };

}
#endif
