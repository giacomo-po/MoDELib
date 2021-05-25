/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>, 
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DislocationNETWORKTRAITS_H_
#define model_DislocationNETWORKTRAITS_H_

#include <Eigen/Dense>
#include <TypeTraits.h>
#include <MPIcout.h>
#include <Quadrature.h>
#include <QuadratureDynamic.h>
#include <QuadPowDynamic.h>
#include <LatticeMath.h>
#include <SplineBase.h>
#include <GlidePlaneModule.h>

namespace model
{
	

    
	/************************************************************/	
	/*	Class Predeclarations ***********************************/
	template <int dim, short unsigned int corder, typename InterpolationType>
	class DislocationNetwork;
		
    template <int dim, short unsigned int corder, typename InterpolationType>
    class DislocationLoop;

    template <int dim, short unsigned int corder, typename InterpolationType>
    class DislocationLoopNode;

    template <int dim, short unsigned int corder, typename InterpolationType>
    class DislocationLoopLink;

	template <int dim, short unsigned int corder, typename InterpolationType>
	class DislocationNode;
		
	template <int dim, short unsigned int corder, typename InterpolationType>
	struct DislocationSegment;
    
	
	/********************************************************************/	
	/*	DislocationNetworkTraitsBase: a base class for Dislocation Network Traits */
	template <int _dim, short unsigned int _corder, typename InterpolationType>
	struct DislocationNetworkTraitsBase
    {
        static constexpr int dim=_dim;
        static constexpr int corder=_corder;
        typedef DislocationNetwork   <dim,corder,InterpolationType>	LoopNetworkType;
        typedef DislocationLoop      <dim,corder,InterpolationType> LoopType;
        typedef DislocationLoopNode  <dim,corder,InterpolationType>	LoopNodeType;
        typedef DislocationLoopLink  <dim,corder,InterpolationType> LoopLinkType;
        typedef DislocationNode      <dim,corder,InterpolationType> NetworkNodeType;
        typedef DislocationSegment   <dim,corder,InterpolationType>	NetworkLinkType;
        
        
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        typedef RationalLatticeDirection<dim> RationalLatticeDirectionType;
        
        typedef RationalLatticeDirectionType                         FlowType;
//        typedef double                         FlowType;
        typedef Eigen::Matrix<double,dim,1>                         VectorDim;
        typedef Eigen::Matrix<double,dim-1,1>                         VectorLowerDim;
        typedef Eigen::Matrix<double,dim,dim>                       MatrixDim;
        typedef Grain<dim>                       GrainType;
        typedef GlidePlane<dim>                       GlidePlaneType;
        typedef PeriodicGlidePlane<dim>                       PeriodicGlidePlaneType;

//        static constexpr FlowType zeroFlow=FlowType::Zero();
//        typedef QuadratureDynamic<1,UniformOpen,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadratureDynamicType;
//        typedef QuadPowDynamic<SplineBase<dim,corder>::pOrder,UniformOpen,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadPowDynamicType;
        
        enum MeshLocation{outsideMesh=-1, insideMesh=0, onMeshBoundary=1, onRegionBoundary=2};

    
    };

	template <int dim, short unsigned int corder, typename InterpolationType>
	struct TypeTraits<DislocationNetwork<dim,corder,InterpolationType> > :
	public DislocationNetworkTraitsBase <dim,corder,InterpolationType>{};

    template <int dim, short unsigned int corder, typename InterpolationType>
    struct TypeTraits<DislocationLoop<dim,corder,InterpolationType> > :
    public DislocationNetworkTraitsBase <dim,corder,InterpolationType>{};

    template <int dim, short unsigned int corder, typename InterpolationType>
    struct TypeTraits<DislocationLoopNode<dim,corder,InterpolationType> > :
    public DislocationNetworkTraitsBase <dim,corder,InterpolationType>{};

    template <int dim, short unsigned int corder, typename InterpolationType>
    struct TypeTraits<DislocationLoopLink<dim,corder,InterpolationType> > :
    public DislocationNetworkTraitsBase <dim,corder,InterpolationType>{};

    template <int dim, short unsigned int corder, typename InterpolationType>
	struct TypeTraits<DislocationNode<dim,corder,InterpolationType> > :
	public DislocationNetworkTraitsBase <dim,corder,InterpolationType>{};

    template <int dim, short unsigned int corder, typename InterpolationType>
	struct TypeTraits<DislocationSegment<dim,corder,InterpolationType> > :
	public DislocationNetworkTraitsBase <dim,corder,InterpolationType>{};
    
    
//    template<>
//    struct NullFlow<RationalLatticeDirection<3>>
//    {
//        static bool isZero(const RationalLatticeDirection<3>& flow)
//        {
//            return flow.squaredNorm()<FLT_EPSILON;
//        }
//    };

} // namespace model
#endif
