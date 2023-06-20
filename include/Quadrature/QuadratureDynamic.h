/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_QUADRATUREDYNAMIC_H_
#define model_QUADRATUREDYNAMIC_H_

#include <math.h>       /* round, floor, ceil, trunc */
#include <assert.h>
#include <set>
#include <Eigen/Dense>
#include <Quadrature.h>

namespace model
{
    /**************************************************************************/
    /**************************************************************************/
    template<
    short unsigned int dim,
    template <short unsigned int, size_t> class QuadratureRule,
    int...Orders>
    struct QuadratureDynamic;
    
    /**************************************************************************/
    /**************************************************************************/
    template<
    short unsigned int dim,
    template <short unsigned int, size_t> class QuadratureRule,
    int order,
    int... otherOrders>
    struct QuadratureDynamic<dim,QuadratureRule,order,otherOrders...>
    {
        typedef typename VectorDimTypeSelector<dim,order>::VectorDim VectorDim;
        
        /**********************************************************************/
        static const double& weight(const short unsigned int & N,
                                    const short unsigned int & k)
        {
//            assert(k < N          && "k must be less than N");
            if(k>=N)
            {
                throw std::runtime_error("k must be less than N");
            }
            return (N==order)? Quadrature<dim,order,QuadratureRule>::weight(k) : QuadratureDynamic<dim,QuadratureRule,otherOrders...>::weight(N,k);
        }
        
        /**********************************************************************/
        static const VectorDim abscissa(const short unsigned int & N,
                                        const short unsigned int & k)
        {
//            std::cout<<"order="<<order<<std::endl;
//            assert(k < N          && "k must be less than N");
            if(k>=N)
            {
                throw std::runtime_error("k must be less than N");
            }
            return (N==order)? Quadrature<dim,order,QuadratureRule>::abscissa(k) : QuadratureDynamic<dim,QuadratureRule,otherOrders...>::abscissa(N,k);
        }
        
        /**********************************************************************/
        template<typename IntegrandType>
        static void integrate(const short unsigned int & N,
                              IntegrandType (*fp)(const double&),
                              IntegrandType &intgrl)
        {
            if (N==order)
            {
                Quadrature<dim,order,QuadratureRule>::integrate(fp,intgrl);
            }
            else
            {
                QuadratureDynamic<dim,QuadratureRule,otherOrders...>::integrate(N,fp,intgrl);
            }
        }
        
        /**********************************************************************/
        template <typename AnyClass, typename IntegrandType, typename ...Args>
        static void integrate(const short unsigned int & N,
                              const AnyClass* const C, IntegrandType &intgrl,
                              IntegrandType (AnyClass::*mfp)(const VectorDim&, const Args&...) const,
                              const Args&...args)
        {
            if (N==order)
            {
                Quadrature<dim,order,QuadratureRule>::integrate(C,intgrl,mfp,args...);
            }
            else
            {
                QuadratureDynamic<dim,QuadratureRule,otherOrders...>::integrate(N,C,intgrl,mfp,args...);
            }
        }
        
        /**********************************************************************/
        template <typename AnyClass, typename IntegrandType, typename ...Args>
        static void integrate(const short unsigned int & N,
                              const AnyClass* const C, IntegrandType &intgrl,
                              IntegrandType (AnyClass::*mfp)(const int&, const Args&... ) const,
                              const Args&...args)
        {
            if (N==order)
            {
                Quadrature<dim,order,QuadratureRule>::integrate(C,intgrl,mfp,args...);
            }
            else
            {
                QuadratureDynamic<dim,QuadratureRule,otherOrders...>::integrate(N,C,intgrl,mfp,args...);
            }
        }
        
        template <typename AnyClass, typename ...Args>
        static void execute(const short unsigned int & N,AnyClass* const C,
                            void (AnyClass::*mfp)(const VectorDim&, const Args&...),
                            const Args&...args)
        {
            if (N==order)
            {
                Quadrature<dim,order,QuadratureRule>::execute(C,mfp,args...);
            }
            else
            {
                QuadratureDynamic<dim,QuadratureRule,otherOrders...>::execute(N,C,mfp,args...);
            }
        }
        
        template <typename AnyClass, typename ...Args>
        static void execute(const short unsigned int & N,AnyClass* const C,
                            void (AnyClass::*mfp)(const int&, const Args&...),
                            const Args&...args)
        {
            if (N==order)
            {
                Quadrature<dim,order,QuadratureRule>::execute(C,mfp,args...);
            }
            else
            {
                QuadratureDynamic<dim,QuadratureRule,otherOrders...>::execute(N,C,mfp,args...);
            }
        }
        
        static const std::set<int> orderSet;
        
        static std::set<int> fillSet()
        {
            std::set<int> temp;
            QuadratureDynamic<dim,QuadratureRule,otherOrders...>::fillSet(temp);
            temp.emplace(order);
            return temp;
        }
        
        static std::set<int> fillSet(std::set<int>& temp)
        {
            QuadratureDynamic<dim,QuadratureRule,otherOrders...>::fillSet(temp);
            temp.emplace(order);
            return temp;
        }
        
        static int lowerOrder(const int& k)
        {
            std::set<int>::const_iterator temp=orderSet.lower_bound(k);
            return temp==orderSet.end()? *orderSet.rbegin() : *temp;
        }
        
        static int lowerOrder(const double& d)
        {
            return lowerOrder(int(round(d)));
        }
        
        
    };
    
    template<
    short unsigned int dim,
    template <short unsigned int, size_t> class QuadratureRule,
    int order,
    int... otherOrders>
    const std::set<int> QuadratureDynamic<dim,QuadratureRule,order,otherOrders...>::orderSet=QuadratureDynamic<dim,QuadratureRule,order,otherOrders...>::fillSet();
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<
    short unsigned int dim,
    template <short unsigned int, size_t> class QuadratureRule,
    int order>
    struct QuadratureDynamic<dim,QuadratureRule,order>
    {
        typedef typename VectorDimTypeSelector<dim,order>::VectorDim VectorDim;
        
        /* weight ********************************************/
        static const double& weight(const short unsigned int & N,
                                    const short unsigned int & k)
        {
            //            assert(k < N          && "k must be less than N");
            //            assert(N==order && "quadrature order N not found in QuadratureDynamic");
            if(k>=N)
            {
                throw std::runtime_error("k must be less than N");
            }
            if(N!=order)
            {
                throw std::runtime_error("quadrature order N not found in QuadratureDynamic");
            }
            return Quadrature<dim,order,QuadratureRule>::weight(k);
        }
        
        /* abscissa ******************************************/
        static const VectorDim abscissa(const short unsigned int & N,
                                        const short unsigned int & k)
        {
//            assert(k < N          && "k must be less than N");
//            assert(N==order && "quadrature order N not found in QuadratureDynamic");
            if(k>=N)
            {
                throw std::runtime_error("k must be less than N");
            }
            if(N!=order)
            {
                throw std::runtime_error("quadrature order N not found in QuadratureDynamic");
            }
            return Quadrature<dim,order,QuadratureRule>::abscissa(k);
        }
        
        template<typename IntegrandType>
        static void integrate(const short unsigned int & N,
                              IntegrandType (*fp)(const double&),
                              IntegrandType &intgrl)
        {
//            assert(N==order && "quadrature order N not found in QuadratureDynamic");
            if(N!=order)
            {
                throw std::runtime_error("quadrature order N not found in QuadratureDynamic");
            }
            Quadrature<dim,order,QuadratureRule>::integrate(fp,intgrl);
        }
        
        template <typename AnyClass, typename IntegrandType, typename ...Args>
        static void integrate(const short unsigned int & N,
                              const AnyClass* const C, IntegrandType &intgrl,
                              IntegrandType (AnyClass::*mfp)(const VectorDim&, const Args&...) const,
                              const Args&...args)
        {
//            assert(N==order && "quadrature order N not found in QuadratureDynamic");
            if(N!=order)
            {
                throw std::runtime_error("quadrature order N not found in QuadratureDynamic");
            }
            Quadrature<dim,order,QuadratureRule>::integrate(C,intgrl,mfp,args...);
        }
        
        template <typename AnyClass, typename IntegrandType, typename ...Args>
        static void integrate(const short unsigned int & N,
                              const AnyClass* const C, IntegrandType &intgrl,
                              IntegrandType (AnyClass::*mfp)(const int&, const Args&... ) const,
                              const Args&...args)
        {
//            assert(N==order && "quadrature order N not found in QuadratureDynamic");
            if(N!=order)
            {
                throw std::runtime_error("quadrature order N not found in QuadratureDynamic");
            }
            Quadrature<dim,order,QuadratureRule>::integrate(C,intgrl,mfp,args...);
        }
        
        template <typename AnyClass, typename ...Args>
        static void execute(const short unsigned int & N,AnyClass* const C,
                            void (AnyClass::*mfp)(const VectorDim&, const Args&...),
                            const Args&...args)
        {
//            assert(N==order && "quadrature order N not found in QuadratureDynamic");
            if(N!=order)
            {
                throw std::runtime_error("quadrature order N not found in QuadratureDynamic");
            }
            Quadrature<dim,order,QuadratureRule>::execute(C,mfp,args...);
        }
        
        template <typename AnyClass, typename ...Args>
        static void execute(const short unsigned int & N,AnyClass* const C,
                            void (AnyClass::*mfp)(const int&, const Args&...),
                            const Args&...args)
        {
//            assert(N==order && "quadrature order N not found in QuadratureDynamic");
            if(N!=order)
            {
                throw std::runtime_error("quadrature order N not found in QuadratureDynamic");
            }
            Quadrature<dim,order,QuadratureRule>::execute(C,mfp,args...);
        }
        
        static const std::set<int> orderSet;
        
        static std::set<int> fillSet(std::set<int>& temp)
        {
            temp.emplace(order);
            return temp;
        }
        
        static std::set<int> fillSet()
        {
            std::set<int> temp;
            temp.emplace(order);
            return temp;
        }
        
    };
    
    template<
    short unsigned int dim,
    template <short unsigned int, size_t> class QuadratureRule,
    int order>
    const std::set<int> QuadratureDynamic<dim,QuadratureRule,order>::orderSet=QuadratureDynamic<dim,QuadratureRule,order>::fillSet();
    
    
} // namespace model
#endif

