/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_QUADRATURE_H_
#define model_QUADRATURE_H_

#include <assert.h>
#include <Eigen/Dense>
#include <GaussLegendre.h>
#include <UniformOpen.h>
#include <NonUniformOpen.h>

namespace model
{
    
    template<short unsigned int dim, size_t qOrder>
    struct VectorDimTypeSelector
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        static const VectorDim abscissa(const short unsigned int & k, const Eigen::Matrix<double,dim,qOrder>& abscissas)
        {
            return abscissas.col(k);
        }
    };
    
    template<size_t qOrder>
    struct VectorDimTypeSelector<1,qOrder>
    {
        typedef double VectorDim;
        
        static const VectorDim& abscissa(const short unsigned int & k, const Eigen::Matrix<double,1,qOrder>& abscissas)
        {
            return abscissas(k);
        }
    };
    
    
    /**************************************************************************/
    /**************************************************************************/
    /*! \brief A static class template for generating quadrature abcsissas and 
     * weights for arbitrary-order and QuadratureRule and dimension.
     */
    template<short unsigned int dim, size_t qOrder, template <short unsigned int, size_t> class QuadratureRule = GaussLegendre>
    struct Quadrature
    {
        
        constexpr static size_t quadratureOrder=qOrder;
        
        
        typedef typename VectorDimTypeSelector<dim,qOrder>::VectorDim VectorDim;
        
        //! the static matrix storing abcsissa values
//        static const Eigen::Matrix<double,dim,qOrder>	abscissas;
        static const Eigen::MatrixXd	abscissas;
        
        //! the static vector storing weight values
//        static const Eigen::Matrix<double,1,  qOrder>	weights;
        static const Eigen::MatrixXd	weights;
        
        /**********************************************************************/
        static const double& weight(const short unsigned int & k)
        {/*\returns the k-th weight
          */
            return weights(k);
        }
        
        /**********************************************************************/
        static const VectorDim abscissa(const short unsigned int & k)
        {/*\returns the k-th abscissa
          */
            return VectorDimTypeSelector<dim,qOrder>::abscissa(k,abscissas);
        }
        
        /**********************************************************************/
        template<typename IntegrandType>
        static void integrate(IntegrandType (*fp)(const double&), IntegrandType &intgrl)
        {/*! Integrates function in range [0,1]
          * @param[in]  fp     the function pointer.
          * @param[out] intgrl the returned integral (values are accumulated)
          *
          * Sample code: \code
          * #include <iostream>
          * #include <Eigen/Dense>
          * #include <Quadrature.h>
          *
          * Eigen::Matrix<double,3,3> integrand(const double& u) 
          * {
          *     return Eigen::Matrix<double,3,3>::Ones()*u*u;
          * }
          *
          * int main()
          * {
          *     Eigen::Matrix<double,3,3> I=Eigen::Matrix<double,3,3>::Zero();
          *     Q.integrate(&integrand,I);
          * 	std::cout<< I <<std::endl;
          * 	return 0;
          * }
          * \endcode
          */
            for(size_t k=0;k<qOrder;++k)
            {
                intgrl+=(*fp)(k) * weights(k);
            }
        }
        
        /**********************************************************************/
        template <typename AnyClass, typename IntegrandType, typename ...Args>
        static void integrate(const AnyClass* const C, IntegrandType &intgrl, IntegrandType (AnyClass::*mfp)(const VectorDim&, const Args&...) const, const Args&...args)
        //		static void integrate(AnyClass* const C, IntegrandType &intgrl, IntegrandType (AnyClass::*mfp)(const VectorDim&, const Args&...), const Args&...args)
        {/*! Integrates non-static member function of arbitrary class in range [0,1]
          * @param[in]  C    The pointer to the class instance
          * @param[out] intgrl the returned integral (values are accumulated)
          * @param[in]  mfp    the member function pointer.
          *
          * Sample code: \code
          * #include <iostream>
          * #include <Eigen/Dense>
          * #include <Quadrature.h>
          *
          * class SomeClass
          * {
          *
          *	Eigen::Matrix<double,3,3> integrand(const double& u) const 
          * {
          *		return Eigen::Matrix<double,3,3>::Ones()*u*u;
          *	}
          *
          *	public:
          *
          *	template <size_t qOrder>
          *	void compute_integral()
          * {
          *		Eigen::Matrix<double,3,3> I=Eigen::Matrix<double,3,3>::Zero();
          *		model::Quadrature<1,qOrder,model::GaussLegendre>::integrate(this,&SomeClass::integrand,I);
          *		std::cout<< I <<std::endl;
          *	}
          *
          * };
          *
          * int main()
          * {
          * 	SomeClass sc;
          *     sc.compute_integral<16>();
          * 	return 0;
          * }
          * \endcode
          */
            for(size_t k=0;k<qOrder;++k)
            {
                intgrl+=(C->*mfp)(abscissa(k),args...) * weights(k);
            }
        }
        
        /**********************************************************************/
        template <typename AnyClass, typename IntegrandType, typename ...Args>
        static void integrate(const AnyClass* const C, IntegrandType &intgrl, IntegrandType (AnyClass::*mfp)(const int&, const Args&... ) const, const Args&...args)
        {/*! Integrates non-static member function of arbitrary class in range [0,1]
          * @param[out] intgrl The memory area to copy to.
          * @param[in]  mfp    The member function pointer.
          * @param[in]  C    The pointer to ...
          *
          * Sample code: \code
          * #include <iostream>
          * #include <Eigen/Dense>
          * #include <Quadrature.h>
          *
          * class SomeClass{
          *
          *	Eigen::Matrix<double,3,3> integrand(const double& u) const {
          *		return Eigen::Matrix<double,3,3>::Ones()*u*u;
          *	}
          *
          *	public:
          *
          *	template <size_t qOrder>
          *	void compute_integral(){
          *		Eigen::Matrix<double,3,3> I=Eigen::Matrix<double,3,3>::Zero();
          *		model::Quadrature<1,qOrder,model::GaussLegendre>::integrate(this,&SomeClass::integrand,I);
          *		std::cout<< I <<std::endl;
          *	}
          *
          * };
          *
          * int main(){
          * 	SomeClass sc;
          *  sc.compute_integral<16>();
          * 	return 0;
          * }
          * \endcode
          */
            
            for(size_t k=0;k<qOrder;++k)
            {
                intgrl+=(C->*mfp)(k,args...) * weights(k);
            }
        }
        
        /**********************************************************************/
        template <typename AnyClass, typename ...Args>
        static void execute(AnyClass* const C, void (AnyClass::*mfp)(const VectorDim&, const Args&...), const Args&...args)
        {/*! Execute non-static member function of arbitrary class, at quadrature points in range [0,1]
          
          */
            for(size_t k=0;k<qOrder;++k)
            {
                (C->*mfp)(abscissa(k),args...);
            }
        }
        
        /**********************************************************************/
        template <typename AnyClass, typename ...Args>
        static void execute(AnyClass* const C, void (AnyClass::*mfp)(const int&, const Args&...), const Args&...args)
        {/*! Execute non-static member function of arbitrary class, at quadrature points in range [0,1]
          
          */
            for(size_t k=0;k<qOrder;++k)
            {
                (C->*mfp)(k,args...);
            }
        }
        
        
    };
    
    
    //! The static const member data abscissas is initialized only once when Quadrature<1,qOrder,QuadratureRule> is first encountered.
    template<short unsigned int dim, size_t qOrder,  template <short unsigned int, size_t> class QuadratureRule>
//    const Eigen::Matrix<double,dim,qOrder> Quadrature<dim,qOrder,QuadratureRule>::abscissas=QuadratureRule<dim,qOrder>::abcsissasAndWeights().template block<dim,qOrder>(0,0);
    const Eigen::MatrixXd Quadrature<dim,qOrder,QuadratureRule>::abscissas=QuadratureRule<dim,qOrder>::abcsissasAndWeights().template block<dim,qOrder>(0,0);
    
    //! The static const member data abscissas is initialized only once when Quadrature<1,qOrder,QuadratureRule> is first encountered.
    template<short unsigned int dim, size_t qOrder,  template <short unsigned int, size_t> class QuadratureRule>
//    const Eigen::Matrix<double,1,qOrder> Quadrature<dim,qOrder,QuadratureRule>::weights=QuadratureRule<dim,qOrder>::abcsissasAndWeights().template block<1,qOrder>(dim,0);
    const Eigen::MatrixXd Quadrature<dim,qOrder,QuadratureRule>::weights=QuadratureRule<dim,qOrder>::abcsissasAndWeights().template block<1,qOrder>(dim,0);
    
} // namespace model
#endif
