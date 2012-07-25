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
#include <model/Dislocations/Materials/Copper.h>
#include <model/Quadrature/Quadrature.h>
#include <model/Quadrature/GaussLegendre.h>

namespace model {
	
	/************************************************************/	
	/*	Class Predeclarations ***********************************/
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
	/*	   */ typename MaterialType>
	class DislocationNetwork;
	
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
	/*	   */ typename MaterialType>
	class DislocationSubNetwork;
	
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
	/*	   */ typename MaterialType>
	class DislocationNode;
		
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
	/*	   */ typename MaterialType>
	class DislocationSegment;
	
	
	/********************************************************************/	
	/*	DislocationNetworkTraitsBase: a base class for Dislocation Network Traits */
	template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule, 
	/*	   */ typename _MaterialType>
	struct DislocationNetworkTraitsBase{
        enum{dim=_dim};
        typedef _MaterialType MaterialType;
		typedef DislocationNetwork   <dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType>	NetworkType;
		typedef DislocationSubNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType>	SubNetworkType;
		typedef DislocationNode      <dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType>	NodeType;
		typedef DislocationSegment   <dim,corder,InterpolationType,alpha,qOrder,QuadratureRule,MaterialType>	LinkType;
		typedef Eigen::Matrix<double,dim,1>																		FlowType;
	};
	
	/************************************************************************************/	
	/*	TypeTraits for: dim=3, corder=1, alpha=centripetal, InterpolationType=CatmullRom 
	 qOrder=16, QuadratureRule=GaussLegendre
	 */
	template<>
	struct TypeTraits<DislocationNetwork	<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper>{};
	
	template<>
	struct TypeTraits<DislocationSubNetwork	<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper>{};
	
	template<>
	struct TypeTraits<DislocationNode		<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper>{};
	
	template<>
	struct TypeTraits<DislocationSegment	<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,GaussLegendre,Copper>{};
	Eigen::Matrix<double,1,16> AvoidMacOsXBug1=Quadrature<1,16,GaussLegendre>::abscissas; // is this static initialization fiasco?
	Eigen::Matrix<double,1,16> AvoidMacOsXBug2=Quadrature<1,16,GaussLegendre>::weights;   // is this static initialization fiasco?
	/************************************************************************************/	
	
	/************************************************************************************/		
	/************************************************************************************/	
	/*	TypeTraits for: dim=3, corder=1, alpha=centripetal, InterpolationType=CatmullRom 
	 qOrder=8, QuadratureRule=GaussLegendre
	 */
	template<>
	struct TypeTraits<DislocationNetwork	<3,1,CatmullRom,centripetal,8,GaussLegendre,Copper> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,8,GaussLegendre,Copper>{};
	
	template<>
	struct TypeTraits<DislocationSubNetwork	<3,1,CatmullRom,centripetal,8,GaussLegendre,Copper> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,8,GaussLegendre,Copper>{};
	
	template<>
	struct TypeTraits<DislocationNode		<3,1,CatmullRom,centripetal,8,GaussLegendre,Copper> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,8,GaussLegendre,Copper>{};
	
	template<>
	struct TypeTraits<DislocationSegment	<3,1,CatmullRom,centripetal,8,GaussLegendre,Copper> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,8,GaussLegendre,Copper>{};
	Eigen::Matrix<double,1,8> AvoidMacOsXBug3=Quadrature<1,8,GaussLegendre>::abscissas; // is this static initialization fiasco?
	Eigen::Matrix<double,1,8> AvoidMacOsXBug4=Quadrature<1,8,GaussLegendre>::weights;   // is this static initialization fiasco?
	/************************************************************************************/	
	
	
	
	
	
	
} // namespace model
#endif
