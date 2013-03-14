/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_QUADRATUREDYNAMIC_H_
#define model_QUADRATUREDYNAMIC_H_

#include <assert.h>
#include <Eigen/Dense>
//#include <model/Quadrature/GaussLegendre.h>
#include <model/Quadrature/Quadrature.h>



namespace model {
	
	
	/*******************************************************************************************************/
	/*******************************************************************************************************/
	/*! \brief A static class template for generating quadrature abcsissas and weights for arbitrary-order
	 *  and QuadratureRule and dimension.
	 */
    
    template<short unsigned int dim, short unsigned int _qOrderMax, template <short unsigned int, short unsigned int> class QuadratureRule = GaussLegendre>
	struct QuadratureDynamic;
    
    template<short unsigned int dim, short unsigned int _qOrderMax, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct QuadratureDynamic {
        
        enum{qOrderMax=_qOrderMax};   // make qOrderMax available externally		
        typedef typename VectorDimTypeSelector<dim,_qOrderMax>::VectorDim VectorDim;
		
		/* weight ********************************************/
		static const double& weight(const short unsigned int & N, const short unsigned int & k) {
            assert(N <= qOrderMax && "N EXCEEDS qOrderMax");
            assert(k < N          && "k must be less than N");
			return (N==qOrderMax)? Quadrature<dim,qOrderMax,QuadratureRule>::weights(k) : QuadratureDynamic<dim,qOrderMax-1,QuadratureRule>::weights(N,k);
		}
		
		/* abscissa ******************************************/
		static const VectorDim abscissa(const short unsigned int & N, const short unsigned int & k) {
            assert(N <= qOrderMax && "N EXCEEDS qOrderMax");
            assert(k < N          && "k must be less than N");
			return (N==qOrderMax)? Quadrature<dim,qOrderMax,QuadratureRule>::abscissa(k) : QuadratureDynamic<dim,qOrderMax-1,QuadratureRule>::abscissa(N,k);
		}
        
		/* integrate_function *******************************/
		template<typename IntegrandType>
		static void integrate(const short unsigned int & N,IntegrandType (*fp)(const double&), IntegrandType &intgrl) {
            assert(N <= qOrderMax && "N EXCEEDS qOrderMax");
            if (N==qOrderMax){
                Quadrature<dim,qOrderMax,QuadratureRule>::integrate(fp,intgrl);
            }
            else{
                QuadratureDynamic<dim,qOrderMax-1,QuadratureRule>::integrate(N,fp,intgrl);
            }
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename IntegrandType, typename ...Args>
		static void integrate(const short unsigned int & N, const AnyClass* const C, IntegrandType &intgrl, IntegrandType (AnyClass::*mfp)(const VectorDim&, const Args&...) const, const Args&...args) {
            assert(N <= qOrderMax && "N EXCEEDS qOrderMax");
            if (N==qOrderMax){
                Quadrature<dim,qOrderMax,QuadratureRule>::integrate(C,intgrl,mfp,args...);
            }
            else{
                QuadratureDynamic<dim,qOrderMax-1,QuadratureRule>::integrate(N,C,intgrl,mfp,args...);
            }
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename IntegrandType, typename ...Args>
		static void integrate(const short unsigned int & N,const AnyClass* const C, IntegrandType &intgrl, IntegrandType (AnyClass::*mfp)(const int&, const Args&... ) const, const Args&...args) {
            assert(N <= qOrderMax && "N EXCEEDS qOrderMax");
            if (N==qOrderMax){
                Quadrature<dim,qOrderMax,QuadratureRule>::integrate(C,intgrl,mfp,args...);
            }
            else{
                QuadratureDynamic<dim,qOrderMax-1,QuadratureRule>::integrate(N,C,intgrl,mfp,args...);
            }
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename ...Args>
		static void execute(const short unsigned int & N,AnyClass* const C, void (AnyClass::*mfp)(const VectorDim&, const Args&...), const Args&...args) {
            assert(N <= qOrderMax && "N EXCEEDS qOrderMax");
            if (N==qOrderMax){
                Quadrature<dim,qOrderMax,QuadratureRule>::execute(C,mfp,args...);
            }
            else{
                QuadratureDynamic<dim,qOrderMax-1,QuadratureRule>::execute(N,C,mfp,args...);
            }
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename ...Args>
		static void execute(const short unsigned int & N,AnyClass* const C, void (AnyClass::*mfp)(const int&, const Args&...), const Args&...args) {
            assert(N <= qOrderMax && "N EXCEEDS qOrderMax");
            if (N==qOrderMax){
                Quadrature<dim,qOrderMax,QuadratureRule>::execute(C,mfp,args...);
            }
            else{
                QuadratureDynamic<dim,qOrderMax-1,QuadratureRule>::execute(N,C,mfp,args...);
            }
		}
		
		
	};
    
    
    template<short unsigned int dim, template <short unsigned int, short unsigned int> class QuadratureRule>
	struct QuadratureDynamic<dim,1,QuadratureRule> : public Quadrature<dim,1,QuadratureRule>{
        
        enum{qOrderMax=1};
        typedef typename VectorDimTypeSelector<dim,qOrderMax>::VectorDim VectorDim;

        
        /* weight ********************************************/
		static const double& weight(const short unsigned int & N, const short unsigned int & k) {
            assert(N ==1 && "N EXCEEDS qOrderMax");
            assert(k < N          && "k must be less than N");
			return Quadrature<dim,qOrderMax,QuadratureRule>::weights(k);
		}
		
		/* abscissa ******************************************/
		static const VectorDim abscissa(const short unsigned int & N, const short unsigned int & k) {
            assert(N ==1 && "N EXCEEDS qOrderMax");
            assert(k < N          && "k must be less than N");
			return Quadrature<dim,qOrderMax,QuadratureRule>::abscissa(k);
		}
		
        
		/* integrate_function *******************************/
		template<typename IntegrandType>
		static void integrate(const short unsigned int & N,IntegrandType (*fp)(const double&), IntegrandType &intgrl) {
            assert(N ==1 && "N EXCEEDS qOrderMax");
            Quadrature<dim,qOrderMax,QuadratureRule>::integrate(fp,intgrl);
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename IntegrandType, typename ...Args>
		static void integrate(const short unsigned int & N, const AnyClass* const C, IntegrandType &intgrl, IntegrandType (AnyClass::*mfp)(const VectorDim&, const Args&...) const, const Args&...args) {
            assert(N ==1 && "N EXCEEDS qOrderMax");
            Quadrature<dim,qOrderMax,QuadratureRule>::integrate(C,intgrl,mfp,args...);
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename IntegrandType, typename ...Args>
		static void integrate(const short unsigned int & N,const AnyClass* const C, IntegrandType &intgrl, IntegrandType (AnyClass::*mfp)(const int&, const Args&... ) const, const Args&...args) {
            assert(N ==1 && "N EXCEEDS qOrderMax");
            Quadrature<dim,qOrderMax,QuadratureRule>::integrate(C,intgrl,mfp,args...);
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename ...Args>
		static void execute(const short unsigned int & N,AnyClass* const C, void (AnyClass::*mfp)(const VectorDim&, const Args&...), const Args&...args) {
            assert(N ==1 && "N EXCEEDS qOrderMax");
            Quadrature<dim,qOrderMax,QuadratureRule>::execute(C,mfp,args...);
		}
		
		/* integrate_function *******************************/
		template <typename AnyClass, typename ...Args>
		static void execute(const short unsigned int & N, AnyClass* const C, void (AnyClass::*mfp)(const int&, const Args&...), const Args&...args) {
            assert(N ==1 && "N EXCEEDS qOrderMax");
            Quadrature<dim,qOrderMax,QuadratureRule>::execute(C,mfp,args...);
		}
        
        
    };
    
	
    
	//////////////////////////////////////////////////////////////
} // namespace model
#endif

