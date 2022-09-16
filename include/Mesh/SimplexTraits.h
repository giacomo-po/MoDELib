/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_SimplexTraits_H_
#define  model_SimplexTraits_H_

#include <iostream>
#include <memory> // shared_ptr
#include <array> // std::array
#include <set>
#include <Eigen/Dense>

namespace model
{
    
    
    template<short int dim, short int order>
    class Simplex;
    
    template<short int dim,short int order>
    class SimplexChild;
    
    
    template<short int nVertices>
    struct SimplexID : public std::array<size_t,nVertices>
    {
        typedef size_t ScalarIDType;
        typedef typename std::array<size_t,nVertices>::size_type size_type;
        
        SimplexID(const std::set<size_t>& set)
        {
            if(set.size()==nVertices)
            {
                int k(0);
                for (const size_t& m : set)
                {
                    this->operator[](k)=m;
                    ++k;
                }
            }
            else
            {
                std::cout<<"set= ";
                for (const size_t& m : set)
                {
                    std::cout<<m<<" ";
                }
                std::cout<<std::endl;
                throw std::runtime_error("Cannot create SimplexID from this set.");
            }
        }
        
        ScalarIDType& operator()(const size_type& k)
        {
            return this->operator[](k);
        }

        const ScalarIDType& operator()(const size_type& k) const
        {
            return this->operator[](k);
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const SimplexID& xID)
        {
            for(const auto& val : xID)
            {
                os  << val<<" ";
            }
            return os;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim, short int order>
    struct SimplexTraits
    {
        
        static_assert(dim>0,"dim must be > 0");
        static_assert(order>0,"order must be >= 0");
        static_assert(order<=dim,"order must be <= dim");
        
        static constexpr int nVertices=SimplexTraits<dim,order-1>::nVertices+1;
        static constexpr int nEdges=SimplexTraits<dim,order-1>::nEdges+SimplexTraits<dim,order-1>::nVertices;
        static constexpr int nFaces=nVertices;
        
        typedef Simplex<dim,order>   SimplexType;
        typedef Simplex<dim,order+1> ParentSimplexType;

        typedef std::array<std::shared_ptr<Simplex<dim,order-1> >,nFaces> BaseArrayType;

        typedef typename SimplexID<nVertices>::ScalarIDType ScalarIDType;
        typedef SimplexID<nVertices> SimplexIDType;
        
        /**********************************************************************/
        static SimplexIDType sortID(const SimplexIDType& vIN)
        {
            std::set<size_t> set;
            for (int k=0;k<nVertices;++k)
            {
                set.insert(vIN(k));
//
//                const bool couldInsert(set.insert(vIN(k)).second);
//                assert(couldInsert && "VERTEX IDs ARE NOT UNIQUE.");
            }
            return SimplexIDType(set);
        }
        
        /**********************************************************************/
        static typename SimplexTraits<dim,order-1>::SimplexIDType faceID(const SimplexIDType& xID,
                                                                         const size_t& j)
        {/*!@param[in] xID the ID of the parent Simplex
          * @param[in] j the j-th face of the parent Simplex
          * \returns the ID of the child Simplex which is the j-th face of xID
          */
            assert(j<nVertices && "REQUESTING NON-EXISTING FACE");
            std::set<size_t> set;
            for (size_t k=0;k<nVertices;++k)
            {
                if(k!=j)
                {
                    set.insert(xID[k]);
                }
            }
            return typename SimplexTraits<dim,order-1>::SimplexIDType(set);
        }
        
        /**********************************************************************/
        static std::array<std::shared_ptr<typename SimplexTraits<dim,order-1>::SimplexIDType>,nFaces> faceIDs(const SimplexIDType& xID)
        {
            std::array<std::shared_ptr<typename SimplexTraits<dim,order-1>::SimplexIDType>,nFaces> temp;
            for (size_t j=0;j<nFaces;++j)
            {
                temp[j].reset(new typename SimplexTraits<dim,order-1>::SimplexIDType(faceID(xID,j)));
            }
            return temp;
        }
        
        /**************************************************************************/
        static size_t faceOrder(const SimplexIDType& parentID,
                                const typename SimplexTraits<dim,order-1>::SimplexIDType& childID)
        {
            const std::array<std::shared_ptr<typename SimplexTraits<dim,order-1>::SimplexIDType>,nFaces> faceIDS(faceIDs(parentID));
            
            int temp=-1;
            for(size_t k=0;k<nFaces;++k)
            {
                if(*faceIDS[k]==childID)
                {
                    temp=k;
                    break;
                }
            }
            assert(temp>=0 && "FACE NOT FOUND");
            return temp;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    struct SimplexTraits<dim,0>
    {
        enum {order=0};
        enum {nVertices=1};
        enum {nEdges=0};
        enum {nFaces=nVertices};
        
        typedef Simplex<dim,order>   SimplexType;
        typedef Simplex<dim,order+1> ParentSimplexType;
        typedef typename SimplexID<nVertices>::ScalarIDType ScalarIDType;
        typedef SimplexID<nVertices> SimplexIDType;

        /**********************************************************************/
        static SimplexIDType sortID(const SimplexIDType& vIN)
        {
            return vIN;
        }
        
    };
    
    
}
#endif
