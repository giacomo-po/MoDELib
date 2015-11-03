/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_SimplexTraits_H_
#define  model_SimplexTraits_H_

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
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim, short int order>
    struct SimplexTraits
    {
        
        static_assert(dim>0,"dim must be > 0");
        static_assert(order>0,"order must be >= 0");
        static_assert(order<=dim,"order must be <= dim");
        
        
        enum {nVertices=SimplexTraits<dim,order-1>::nVertices+1};
        enum {nEdges=SimplexTraits<dim,order-1>::nEdges+SimplexTraits<dim,order-1>::nVertices};
        enum {nFaces=nVertices};
        
        typedef std::array<std::shared_ptr<Simplex<dim,order-1> >,nFaces> BaseArrayType;
        
        typedef size_t ScalarIDType;
        
        typedef Eigen::Matrix<ScalarIDType,1,nVertices> SimplexIDType;
        
        /**********************************************************************/
        static SimplexIDType sortID(const SimplexIDType& vIN)
        {
            std::set<size_t> set;
            for (int k=0;k<nVertices;++k)
            {
                const bool couldInsert(set.insert(vIN(k)).second);
                assert(couldInsert && "VERTEX IDs ARE NOT UNIQUE.");
            }
            
            SimplexIDType temp;
            int k(0);
            //			for (std::set<size_t>::const_iterator iter=set.begin();iter!=set.end();++iter)
            //            {
            //				temp(k)=(*iter);
            //				++k;
            //			}
            for (const size_t& m : set)
            {
                temp(k)=m;
                ++k;
            }
            return temp;
        }
        
        /**********************************************************************/
        static typename SimplexTraits<dim,order-1>::SimplexIDType faceID(const SimplexIDType& xID,
                                                                         const size_t& j)
        {/*!@param[in] xID the ID of the parent Simplex
          * @param[in] j the j-th face of the parent Simplex
          * \returns the ID of the child Simplex which is the j-th face of xID
          */
            assert(j<nVertices && "REQUESTING NON-EXISTING FACE");
            
            typename SimplexTraits<dim,order-1>::SimplexIDType temp;
            size_t m(0);
            for (size_t k=0;k<nVertices;++k)
            {
                if(k!=j)
                {
                    temp(m)=(xID(k));
                    ++m;
                }
            }
            return temp;
        }
        
        /**********************************************************************/
        static std::array<typename SimplexTraits<dim,order-1>::SimplexIDType,nFaces> faceIDs(const SimplexIDType& xID)
        {
            std::array<typename SimplexTraits<dim,order-1>::SimplexIDType,nFaces> temp;
            for (size_t j=0;j<nFaces;++j)
            {
                temp[j]=faceID(xID,j);
            }
            return temp;
        }
        
        /**************************************************************************/
        static size_t faceOrder(const SimplexIDType& parentID,
                                const typename SimplexTraits<dim,order-1>::SimplexIDType& childID)
        {
            const std::array<typename SimplexTraits<dim,order-1>::SimplexIDType,nFaces> faceIDS(faceIDs(parentID));
            
            int temp=-1;
            for(size_t k=0;k<nFaces;++k)
            {
                if(faceIDS[k]==childID)
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
        
        typedef size_t ScalarIDType;
        
        typedef Eigen::Matrix<ScalarIDType,1,nVertices> SimplexIDType;
        
        /**********************************************************************/
        static SimplexIDType sortID(const SimplexIDType& vIN)
        {
            return vIN;
        }
        
    };
    
    
}
#endif
