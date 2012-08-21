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
//#include <model/Quadrature/GaussLegendre.h>

namespace model {
	
	/************************************************************/	
	/*	Class Predeclarations ***********************************/
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationNetwork;
	
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationSubNetwork;
	
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationNode;
		
	template <short unsigned int dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	class DislocationSegment;
	
	
	/********************************************************************/	
	/*	DislocationNetworkTraitsBase: a base class for Dislocation Network Traits */
	template <short unsigned int _dim, short unsigned int corder, typename InterpolationType,
	/*	   */ double & alpha, short unsigned int qOrder, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct DislocationNetworkTraitsBase{
        enum{dim=_dim};
        //typedef _MaterialType MaterialType;
		typedef DislocationNetwork   <dim,corder,InterpolationType,alpha,qOrder,QuadratureRule>	NetworkType;
		typedef DislocationSubNetwork<dim,corder,InterpolationType,alpha,qOrder,QuadratureRule>	SubNetworkType;
		typedef DislocationNode      <dim,corder,InterpolationType,alpha,qOrder,QuadratureRule>	NodeType;
		typedef DislocationSegment   <dim,corder,InterpolationType,alpha,qOrder,QuadratureRule>	LinkType;
		typedef Eigen::Matrix<double,dim,1>																		FlowType;
	};
	
	/************************************************************************************/	
	/*	TypeTraits for: dim=3, corder=1, alpha=centripetal, InterpolationType=CatmullRom 
	 qOrder=16, QuadratureRule=GaussLegendre
	 */
	template<>
	struct TypeTraits<DislocationNetwork	<3,1,CatmullRom,centripetal,16,GaussLegendre> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,GaussLegendre>{};
	
	template<>
	struct TypeTraits<DislocationSubNetwork	<3,1,CatmullRom,centripetal,16,GaussLegendre> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,GaussLegendre>{};
	
	template<>
	struct TypeTraits<DislocationNode		<3,1,CatmullRom,centripetal,16,GaussLegendre> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,GaussLegendre>{};
	
	template<>
	struct TypeTraits<DislocationSegment	<3,1,CatmullRom,centripetal,16,GaussLegendre> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,GaussLegendre>{};
	Eigen::Matrix<double,1,16> AvoidMacOsXBug1=Quadrature<1,16,GaussLegendre>::abscissas; // is this static initialization fiasco?
	Eigen::Matrix<double,1,16> AvoidMacOsXBug2=Quadrature<1,16,GaussLegendre>::weights;   // is this static initialization fiasco?
	/************************************************************************************/	
	
	/************************************************************************************/		
	/************************************************************************************/	
	/*	TypeTraits for: dim=3, corder=1, alpha=centripetal, InterpolationType=CatmullRom 
	 qOrder=8, QuadratureRule=GaussLegendre
	 */
	template<>
	struct TypeTraits<DislocationNetwork	<3,1,CatmullRom,centripetal,8,GaussLegendre> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,8,GaussLegendre>{};
	
	template<>
	struct TypeTraits<DislocationSubNetwork	<3,1,CatmullRom,centripetal,8,GaussLegendre> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,8,GaussLegendre>{};
	
	template<>
	struct TypeTraits<DislocationNode		<3,1,CatmullRom,centripetal,8,GaussLegendre> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,8,GaussLegendre>{};
	
	template<>
	struct TypeTraits<DislocationSegment	<3,1,CatmullRom,centripetal,8,GaussLegendre> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,8,GaussLegendre>{};
	Eigen::Matrix<double,1,8> AvoidMacOsXBug3=Quadrature<1,8,GaussLegendre>::abscissas; // is this static initialization fiasco?
	Eigen::Matrix<double,1,8> AvoidMacOsXBug4=Quadrature<1,8,GaussLegendre>::weights;   // is this static initialization fiasco?
	/************************************************************************************/	
    
    /************************************************************************************/	
	/*	TypeTraits for: dim=3, corder=1, alpha=centripetal, InterpolationType=CatmullRom 
	 qOrder=16, QuadratureRule=UniformOpen
	 */
	template<>
	struct TypeTraits<DislocationNetwork	<3,1,CatmullRom,centripetal,16,UniformOpen> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationSubNetwork	<3,1,CatmullRom,centripetal,16,UniformOpen> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationNode		<3,1,CatmullRom,centripetal,16,UniformOpen> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationSegment	<3,1,CatmullRom,centripetal,16,UniformOpen> > : 
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,16,UniformOpen>{};
	Eigen::Matrix<double,1,16> AvoidMacOsXBug5=Quadrature<1,16,UniformOpen>::abscissas; // is this static initialization fiasco?
	Eigen::Matrix<double,1,16> AvoidMacOsXBug6=Quadrature<1,16,UniformOpen>::weights;   // is this static initialization fiasco?
	/************************************************************************************/	

	
    /************************************************************************************/
	/*	TypeTraits for: dim=3, corder=1, alpha=centripetal, InterpolationType=CatmullRom
	 qOrder=32, QuadratureRule=UniformOpen
	 */
	template<>
	struct TypeTraits<DislocationNetwork	<3,1,CatmullRom,centripetal,32,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,32,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationSubNetwork	<3,1,CatmullRom,centripetal,32,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,32,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationNode		<3,1,CatmullRom,centripetal,32,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,32,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationSegment	<3,1,CatmullRom,centripetal,32,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,32,UniformOpen>{};
	Eigen::Matrix<double,1,32> AvoidMacOsXBug7=Quadrature<1,32,UniformOpen>::abscissas; // is this static initialization fiasco?
	Eigen::Matrix<double,1,32> AvoidMacOsXBug8=Quadrature<1,32,UniformOpen>::weights;   // is this static initialization fiasco?
	/************************************************************************************/

    /************************************************************************************/
	/*	TypeTraits for: dim=3, corder=1, alpha=centripetal, InterpolationType=CatmullRom
	 qOrder=64, QuadratureRule=UniformOpen
	 */
	template<>
	struct TypeTraits<DislocationNetwork	<3,1,CatmullRom,centripetal,64,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,64,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationSubNetwork	<3,1,CatmullRom,centripetal,64,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,64,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationNode		<3,1,CatmullRom,centripetal,64,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,64,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationSegment	<3,1,CatmullRom,centripetal,64,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,centripetal,64,UniformOpen>{};
	Eigen::Matrix<double,1,64> AvoidMacOsXBug9=Quadrature<1,64,UniformOpen>::abscissas; // is this static initialization fiasco?
	Eigen::Matrix<double,1,64> AvoidMacOsXBug10=Quadrature<1,64,UniformOpen>::weights;   // is this static initialization fiasco?
	/************************************************************************************/

    /************************************************************************************/
	/*	TypeTraits for: dim=3, corder=1, alpha=centripetal, InterpolationType=CatmullRom
	 qOrder=64, QuadratureRule=UniformOpen
	 */
	template<>
	struct TypeTraits<DislocationNetwork	<3,1,CatmullRom,chordal,64,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,chordal,64,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationSubNetwork	<3,1,CatmullRom,chordal,64,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,chordal,64,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationNode		<3,1,CatmullRom,chordal,64,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,chordal,64,UniformOpen>{};
	
	template<>
	struct TypeTraits<DislocationSegment	<3,1,CatmullRom,chordal,64,UniformOpen> > :
	public DislocationNetworkTraitsBase		<3,1,CatmullRom,chordal,64,UniformOpen>{};
	Eigen::Matrix<double,1,64> AvoidMacOsXBug11=Quadrature<1,64,UniformOpen>::abscissas; // is this static initialization fiasco?
	Eigen::Matrix<double,1,64> AvoidMacOsXBug12=Quadrature<1,64,UniformOpen>::weights;   // is this static initialization fiasco?
	/************************************************************************************/

	
	
	
	
	
} // namespace model
#endif
