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
#include <memory>

//#include <Simplex.h>
#include <SimplicialMesh.h>
#include <Polycrystal.h>
#include <PolycrystallineMaterialBase.h>
#include <LatticeModule.h>
#include <DDtraitsIO.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <DislocationDynamicsBase.h>
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


//        DDtraitsIO traitsIO;
        DDconfigIO<3> configIO;
        DDauxIO<3> auxIO;

    public:
        
//        const TextFileParser parser;
        DislocationDynamicsBase<3>& ddBase;
        const bool outputBinary;
        const double minSize;
        const double maxSize;
        std::map<VectorDimD,size_t,CompareVectorsByComponent<double,dim,float>> uniqueNetworkNodeMap;
        
        MicrostructureGenerator(DislocationDynamicsBase<3>& ddBase_in);
        
        void readMicrostructureFile();
        const DDtraitsIO& traits() const;
        const DDconfigIO<3>& config() const;
        const DDauxIO<3>& aux() const;
        DDconfigIO<3>& config();
        DDauxIO<3>& aux();

        size_t insertLoop(const VectorDimD& b,const VectorDimD& unitNormal,const VectorDimD& P0,const size_t& grainID,const DislocationLoopType& loopType);
        size_t insertLoopNode(const size_t& loopID,const VectorDimD& loopNodePos,const size_t& networkNodeID,const VectorDimD& loopNodeShift,const std::pair<short int,short int>& periodicEdgeIDs);
        std::vector<size_t> insertLoopLinks(const size_t& loopID,const std::vector<size_t>& loopNodeIDs);
        size_t insertNetworkNode(const VectorDimD& networkNodePos);
        size_t insertInclusion(const VectorDimD& pos,const double& R, const Eigen::Matrix<double,dim,dim>& eT, const double& vrc,const int&type);

        size_t insertInclusion(const std::map<size_t,Eigen::Vector3d>& nodes,const std::map<size_t,std::vector<size_t>>& faceMap, const Eigen::Matrix<double,dim,dim>& eT, const double& vrc,const int&type);


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


}
#endif
