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
#include <model/Utilities/TypeTraits.h>
#include <model/Geometry/Splines/SplineConsts.h>
#include <model/Quadrature/Quadrature.h>
#include <model/Quadrature/QuadratureDynamic.h>
#include <model/Quadrature/QuadPowDynamic.h>
#include <model/LatticeMath/LatticeMath.h>

namespace model
{
	
	/************************************************************/	
	/*	Class Predeclarations ***********************************/
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ template <short unsigned int, size_t> class QuadratureRule>
	class DislocationNetwork;
		
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ template <short unsigned int, size_t> class QuadratureRule>
	class DislocationNode;
		
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ template <short unsigned int, size_t> class QuadratureRule>
	class DislocationSegment;
	
	
	/********************************************************************/	
	/*	DislocationNetworkTraitsBase: a base class for Dislocation Network Traits */
	template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ template <short unsigned int, size_t> class QuadratureRule>
	struct DislocationNetworkTraitsBase
    {
        enum{dim=_dim};
		typedef DislocationNetwork   <dim,corder,InterpolationType,QuadratureRule>	NetworkType;
		typedef DislocationNode      <dim,corder,InterpolationType,QuadratureRule>	NodeType;
		typedef DislocationSegment   <dim,corder,InterpolationType,QuadratureRule>	LinkType;
//		typedef Eigen::Matrix<double,dim,1>													FlowType;
        typedef LatticeVector<3>                                                            FlowType;
        typedef QuadratureDynamic<1,QuadratureRule,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadratureDynamicType;
        typedef    QuadPowDynamic<3,QuadratureRule,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadPowDynamicType;
    };

	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ template <short unsigned int, size_t> class QuadratureRule>
	struct TypeTraits<DislocationNetwork<dim,corder,InterpolationType,QuadratureRule> > :
	public DislocationNetworkTraitsBase <dim,corder,InterpolationType,QuadratureRule>{};

    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ template <short unsigned int, size_t> class QuadratureRule>
	struct TypeTraits<DislocationNode<dim,corder,InterpolationType,QuadratureRule> > :
	public DislocationNetworkTraitsBase <dim,corder,InterpolationType,QuadratureRule>{};

    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ template <short unsigned int, size_t> class QuadratureRule>
	struct TypeTraits<DislocationSegment<dim,corder,InterpolationType,QuadratureRule> > :
	public DislocationNetworkTraitsBase <dim,corder,InterpolationType,QuadratureRule>{};

} // namespace model
#endif
