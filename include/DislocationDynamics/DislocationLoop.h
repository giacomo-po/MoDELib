/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2019 by Yash Pachaury      <ypachaur@purdue.edu>
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoop_h_
#define model_DislocationLoop_h_


#include <DislocationDynamicsModule.h>


#ifndef NDEBUG
#define VerboseDislocationLoop(N,x) if(verboseDislocationLoop>=N){std::cout<<x;}
#else
#define VerboseDislocationLoop(N,x)
#endif

namespace model
{
    template <int _dim, short unsigned int corder>
    class DislocationLoop : public Loop<DislocationLoop<_dim,corder>>
    {

        
        
    public:
        
        typedef TypeTraits<DislocationLoop<_dim,corder>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        typedef typename TraitsType::VectorDim VectorDim;
        typedef typename TraitsType::MatrixDim MatrixDim;
        typedef typename TraitsType::GrainType GrainType;
        typedef typename TraitsType::GlidePlaneType GlidePlaneType;
        typedef typename TraitsType::PeriodicGlidePlaneType PeriodicGlidePlaneType;
        typedef typename TraitsType::ReciprocalLatticeDirectionType ReciprocalLatticeDirectionType;
        static int verboseDislocationLoop;

        const std::shared_ptr<GlidePlaneType> glidePlane;
        const std::shared_ptr<PeriodicGlidePlaneType> periodicGlidePlane;
        const GrainType& grain;
        const int loopType;

        
    private:
        
        VectorDim nA;
        VectorDim nAR;
        double _slippedArea;
        double _slippedAreaRate;

        VectorDim _rightHandedUnitNormal;
        ReciprocalLatticeDirectionType _rightHandedNormal;
        std::shared_ptr<SlipSystem> _slipSystem;

        
    public:
        
        DislocationLoop(LoopNetworkType* const,
                        const VectorDim&,
                        const std::shared_ptr<GlidePlaneType>& glidePlane_in);
//        DislocationLoop(LoopNetworkType *const ,
//                        const VectorDim &,
//                        const int &,
//                        const int &);
        ~DislocationLoop();
        std::shared_ptr<LoopType> clone() const;
        const double& slippedArea() const;
        const double& slippedAreaRate() const;
        const VectorDim& rightHandedUnitNormal() const;
        const ReciprocalLatticeDirectionType& rightHandedNormal() const;
        std::tuple<double,double,double> loopLength() const;
        std::shared_ptr<SlipSystem> searchSlipSystem() const;
        void updateSlipSystem();
        void updateGeometry();
        MatrixDim plasticDistortion() const;
        MatrixDim plasticDistortionRate() const;
        VectorDim burgers() const;
        const std::shared_ptr<SlipSystem>&  slipSystem() const;
        bool isVirtualBoundaryLoop() const;
        double solidAngle(const VectorDim& x) const;
        
        static void initFromFile(const std::string&);
        static double planarSolidAngle(const VectorDim& x,const VectorDim& planePoint,const VectorDim& rhN,const std::vector<std::pair<VectorDim,VectorDim>>& polygonSegments);
        template <typename T> static int sgn(const T& val);
        
        //        typedef DislocationLoop<_dim,corder> LoopType;

//
//        constexpr static int dim=_dim;
//
////        std::shared_ptr<SlipSystem> slipSystem;
////
////        /**********************************************************************/
////        void updateSlipSystem()
////        {
////            slipSystem=nullptr;
//////            const ReciprocalLatticeDirection<dim> rightHandedDir(this->grain.reciprocalLatticeDirection(this->rightHandedNormal()));
////            for(const auto& ss : this->grain.slipSystems())
////            {
////                if(  ((this->flow()-ss->s).squaredNorm()==0 && (this->rightHandedNormal()-ss->n).squaredNorm()==0)
////                   ||((this->flow()+ss->s).squaredNorm()==0 && (this->rightHandedNormal()+ss->n).squaredNorm()==0))
////                {
////                    slipSystem=ss;
////                }
////            }
////        }
//
//    public:
//
//        typedef DislocationLoop<dim,corder> DislocationLoopType;
//        typedef PlanarDislocationLoop<DislocationLoopType> BaseLoopType;
//        typedef typename TypeTraits<DislocationLoopType>::LoopNetworkType LoopNetworkType;
//        typedef PeriodicDislocationLoop<LoopNetworkType> PeriodicDislocationLoopType;
//        typedef Eigen::Matrix<double,dim,1> VectorDim;
//        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
//
//
//
//        /**********************************************************************/
//        DislocationLoop(LoopNetworkType* const dn,
//                        const VectorDim& B,
//                        const std::shared_ptr<GlidePlane<dim>>& glidePlane) :
////        /* init */ BaseLoopType(dn,glidePlane->grain.latticeVector(B),glidePlane)
//        /* init */ BaseLoopType(dn,glidePlane->grain.rationalLatticeDirection(B),glidePlane)
//        //        /* init */,slipSystem(nullptr)
//        //        /* init */,isGlissile(this->flow().dot(this->glidePlane->n)==0)
//        {
//        }
//
//        /**********************************************************************/
//        DislocationLoop(LoopNetworkType* const dn,
//                        const VectorDim& B,
//                        const std::shared_ptr<GlidePlane<dim>>& glidePlane,
////                        const std::shared_ptr<PeriodicDislocationLoopType>& pLoop,
//                        const VectorDim& shift) :
////        /* init */ BaseLoopType(dn,glidePlane->grain.latticeVector(B),glidePlane,pLoop,shift)
//        /* init */ BaseLoopType(dn,glidePlane->grain.rationalLatticeDirection(B),glidePlane,shift)
//        //        /* init */,slipSystem(nullptr)
//        //        /* init */,isGlissile(this->flow().dot(this->glidePlane->n)==0)
//        {
//        }
//
//        /**********************************************************************/
//        DislocationLoop(LoopNetworkType* const dn,
//                        const VectorDim& B,
//                        const int& grainID,
//                        const int& _loopType) :
//        /* init */ BaseLoopType(dn,dn->poly.grain(grainID).rationalLatticeDirection(B),grainID,_loopType)
////        /* init */,slipSystem(nullptr)
////        /* init */,isGlissile(false)
//        {// Virtual dislocation loop
//        }
//
////        /**********************************************************************/
////        DislocationLoop(const DislocationLoop& other) :
////        /* base init */ BaseLoopType(other)
////        /* init */,slipSystem(other.slipSystem)
//////        /* init */,isGlissile(other.isGlissile)
////        {
////        }
//

        /******************************************************************/


        /**********************************************************************/
//
        /**********************************************************************/
//
//
//        /**********************************************************************/
//        template <class T>
//        friend T& operator << (T& os, const DislocationLoopType& dL)
//        {
//            os<< DislocationLoopIO<dim>(dL);
//            return os;
//        }
//
   };
    
}
#endif
