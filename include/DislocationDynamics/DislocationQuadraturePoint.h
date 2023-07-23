/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

//Implements non-uniform open quadrature rule

#ifndef model_DislocationQuadraturePoint_H_
#define model_DislocationQuadraturePoint_H_

#include <Eigen/Dense>
#include <SplineBase.h>
#include <QuadratureDynamic.h>
#include <QuadPowDynamic.h>
#include <StressStraight.h>
#include <SegmentSegmentDistance.h>
#include <StraightDislocationSegment.h>
#include <EshelbyInclusionBase.h>
#include <DefectiveCrystalParameters.h>
#include <Polygon2D.h>
#include <CatmullRomSplineSegment.h>
#include <DislocationSegment.h>

namespace model
{
    template<int dim,int corder>
    struct DislocationQuadraturePoint
    {
        
        static constexpr int Ncoeff= SplineBase<dim,corder>::Ncoeff;
        static constexpr int pOrder= SplineBase<dim,corder>::pOrder;
        static constexpr int Ndof= SplineBase<dim,corder>::Ndof;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<float,dim,dim> MatrixDimF;
        typedef Eigen::Matrix<double,dim,Ndof> MatrixDimNdof;
        typedef Eigen::Matrix<double,Ncoeff,Ncoeff> MatrixNcoeff;
        typedef Eigen::Matrix<double,Ncoeff,dim> MatrixNcoeffDim;
        typedef DislocationSegment<dim,corder> LinkType;
        typedef   QuadratureDynamic<1,UniformOpen,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadratureDynamicType;
        typedef QuadPowDynamic<pOrder,UniformOpen,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadPowDynamicType;

//        typedef   QuadratureDynamic<1,NonUniformOpen,3,5,7,9,11,13,15,17,19,21,23,25,27,29> QuadratureDynamicType; // this causes some line wiggles
//        typedef QuadPowDynamic<pOrder,NonUniformOpen,3,5,7,9,11,13,15,17,19,21,23,25,27,29> QuadPowDynamicType; // this causes some line wiggles

        const size_t sourceID;
        const size_t sinkID;
        const int qID;
        const Eigen::Matrix<double,1,Ncoeff> SF;    // Spline shape-functions at this quadrature point
        const VectorDim r;                          // position
        const VectorDim ru;                         // parametric tangent dr/du with u in [0:1]
        const double j;                             // jacobian dl/du
        const VectorDim rl;                         // unit tangent dr/dl
        const double dL;                            // line length corresponding to this quadrature point
        
        MatrixDim stress;
        VectorDim pkForce;
        VectorDim stackingFaultForce;
        VectorDim lineTensionForce;
        VectorDim glideVelocity;
        double elasticEnergyPerLength;
        double coreEnergyPerLength;
        int inclusionID;
        
        DislocationQuadraturePoint(const LinkType& parentSegment,
                                   const int& q,const int& qOrder,
                                   const MatrixNcoeff& SFCH,
                                   const MatrixNcoeffDim& qH);
        DislocationQuadraturePoint(const size_t sourceID_in,
                                   const size_t sinkID_in,
                                   const int& q,const int& qOrder,
                                   const MatrixNcoeff& SFCH,
                                   const MatrixNcoeffDim& qH);
        DislocationQuadraturePoint();
        MatrixDimNdof SFgaussEx() const;
        static VectorDim  getGlideVelocity(const LinkType& parentSegment,
                                           const VectorDim& r,
                                           const VectorDim& fPK,
                                           const MatrixDim& S,
                                           const VectorDim& rl,
                                           const double& dL,
                                           const int& inclusionID
                                           );
        
 
        void updateForcesAndVelocities(const LinkType& parentSegment);
        MatrixDim forceToStress(const VectorDim& force,const VectorDim& unitTangent,const LinkType& parentSegment) const;
        
    };
    
    template<int dim,int corder>
    class DislocationQuadraturePointContainer : public std::deque<DislocationQuadraturePoint<dim,corder>>
    
    {
        static constexpr int Ncoeff= SplineBase<dim,corder>::Ncoeff;
        static constexpr int Ndof= SplineBase<dim,corder>::Ndof;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,Ndof,1> VectorNdof;
        typedef Eigen::Matrix<double,Ndof,Ndof> MatrixNdof;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,Ndof> MatrixDimNdof;
        typedef DislocationQuadraturePointContainer<dim,corder> DislocationQuadraturePointContainerType;
        typedef DislocationQuadraturePoint<dim,corder> DislocationQuadraturePointType;
        typedef std::deque<DislocationQuadraturePointType> BaseContainerType;
        typedef Eigen::Matrix<double,Ncoeff,Ncoeff> MatrixNcoeff;
        typedef Eigen::Matrix<double,Ncoeff,dim> MatrixNcoeffDim;
        typedef typename DislocationQuadraturePointType::QuadratureDynamicType QuadratureDynamicType;
        typedef typename DislocationQuadraturePointType::QuadPowDynamicType QuadPowDynamicType;
        typedef DislocationSegment<dim,corder> LinkType;

        VectorNdof nodalVelocityLinearKernel(const int& k) const;
        MatrixNdof nodalVelocityBilinearKernel(const int& k) const;
        VectorDim pkKernel(const int& k) const;
        MatrixDim stressKernel(const int& k) const;
        VectorDim glideVelocityKernel(const int& k) const;
        void updateForcesAndVelocitiesKernel(const LinkType& parentSegment,const StraightDislocationSegment<dim>& ss,const double& L0,const VectorDim& c);

    public:
        
        void updateForcesAndVelocities(const LinkType& parentSegment,
                                       const double&,
                                       const double& );
        void updateQuadraturePoints(const LinkType& parentSegment,
                                    const double& quadPerLength,
                                    const bool& isClimbStep);
        const DislocationQuadraturePointContainerType& quadraturePoints() const;
        DislocationQuadraturePointContainerType& quadraturePoints();
        const DislocationQuadraturePointType& quadraturePoint(const int& k) const;
        VectorNdof nodalVelocityVector(const LinkType& parentSegment) const;
        MatrixNdof nodalVelocityMatrix(const LinkType& parentSegment) const;
        VectorDim glideVelocityIntegral() const;
        VectorDim pkIntegral() const;
        MatrixDim stressIntegral() const;
        
    };
    
}
#endif
