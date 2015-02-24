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

namespace model
{
	
	/************************************************************/	
	/*	Class Predeclarations ***********************************/
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationNetwork;
		
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationNode;
		
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationSegment;
	
	
	/********************************************************************/	
	/*	DislocationNetworkTraitsBase: a base class for Dislocation Network Traits */
	template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct DislocationNetworkTraitsBase
    {
        enum{dim=_dim};
		typedef DislocationNetwork   <dim,corder,InterpolationType,qOrder,QuadratureRule>	NetworkType;
		typedef DislocationNode      <dim,corder,InterpolationType,qOrder,QuadratureRule>	NodeType;
		typedef DislocationSegment   <dim,corder,InterpolationType,qOrder,QuadratureRule>	LinkType;
		typedef Eigen::Matrix<double,dim,1>													FlowType;
	};

	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct TypeTraits<DislocationNetwork<dim,corder,InterpolationType,qOrder,QuadratureRule> > :
	public DislocationNetworkTraitsBase <dim,corder,InterpolationType,qOrder,QuadratureRule>{};

    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct TypeTraits<DislocationNode<dim,corder,InterpolationType,qOrder,QuadratureRule> > :
	public DislocationNetworkTraitsBase <dim,corder,InterpolationType,qOrder,QuadratureRule>{};

    template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct TypeTraits<DislocationSegment<dim,corder,InterpolationType,qOrder,QuadratureRule> > :
	public DislocationNetworkTraitsBase <dim,corder,InterpolationType,qOrder,QuadratureRule>{};

} // namespace model
#endif
