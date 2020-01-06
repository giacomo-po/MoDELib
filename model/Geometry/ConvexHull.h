/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ConvexHull_H_
#define model_ConvexHull_H_

#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <vector>
#include <cfloat>

namespace model
{
    
    
    template <typename T>
    class reverse_range
    {
        const T &x;
        
    public:
        reverse_range(T &x) : x(x) {}
        
        auto begin() const -> decltype(this->x.rbegin())
        {
            return x.rbegin();
        }
        
        auto end() const -> decltype(this->x.rend())
        {
            return x.rend();
        }
    };
    
    template <typename T>
    reverse_range<T> reverse_iterate(T &x)
    {
        return reverse_range<T>(x);
    }
    
    template<int dim,typename T>
    struct HullPoint : public std::array<double,dim>
    {
        const T* const t ;
        
        HullPoint() :
        /* init */t(nullptr)
        {
            
        }
        
        HullPoint(const std::array<double,dim>& p_in,
                  const T* const t_in) :
        /* init */ std::array<double,dim>(p_in)
        /* init */,t(t_in)
        {
            
        }
        
//        HullPoint(const HullPoint& other)
//        {}
        
    };
    
    
    template <int dim,typename T>
    struct ConvexHull
    {
        
        static_assert("ConvexHull not implemented for thid dimension");
    };
    
    template <typename T>
    struct ConvexHull<2,T> :  public std::set<HullPoint<2,T>>
    {/*! Constructs the Convex-Hull of a planar set of points
      * algorithm adapted from https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
      */
        typedef HullPoint<2,T> HullPointType;
        typedef std::vector<HullPointType> HullPointContainerType;
        static constexpr int dim=2;
        static constexpr double tol=FLT_EPSILON;
        
    private:
        
        static double cross(const HullPointType& O, const HullPointType& A, const HullPointType& B)
        {
            return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
        }
        
    public:
        
        HullPointContainerType getPoints()
        {
            HullPointContainerType H;
            if(this->size()<=3)
            {
                for(const auto& pt : *this)
                {
                    H.push_back(pt);
                }
            }
            else
            {
                // Build lower hull
                for (const auto& pt : *this)
                {
                    while(H.size()>=2 && cross(H[H.size()-2], H[H.size()-1], pt)<= tol)
                    {
                        H.pop_back();
                    }
                    H.push_back(pt);
                }
                H.pop_back();

                
                const size_t lowerSize(H.size());
                
                for (const auto& pt : reverse_iterate(*this))
                {
                    while(H.size()-lowerSize>=2 && cross(H[H.size()-2], H[H.size()-1], pt)<= tol)
                    {
                        H.pop_back();
                    }
                    H.push_back(pt);
                }
                H.pop_back();

            }
            return H;
        }


        
    };
    
}
#endif

