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
#include <vector>
#include <cfloat>

namespace model
{
    
    template<int dim,typename T>
    struct HullPoint : public std::array<double,dim>
    {
        const T*  t ;
        
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
    struct ConvexHull<2,T> :  public std::vector<HullPoint<2,T>>
    {
        typedef HullPoint<2,T> HullPointType;
        typedef std::vector<HullPointType> HullPointContainerType;
        static constexpr int dim=2;
        static constexpr double tol=FLT_EPSILON;
        
    private:
        
        double cross(const HullPointType& O, const HullPointType& A, const HullPointType& B)
        {
            return (A[0] - O[0]) * (B[1] - O[1]) - (A[1] - O[1]) * (B[0] - O[0]);
        }
        
    public:
        
        HullPointContainerType getPoints()
        {
            if(this->size()<=3)
            {
                return *this;
            }
            else
            {
                std::sort(this->begin(), this->end());
                HullPointContainerType H(*this);
                size_t k(0);

                // Build lower hull
                for (size_t i = 0; i < this->size(); ++i)
                {
                    while (k >= 2 && cross(H[k-2], H[k-1], this->operator[](i)) <= tol) k--;
                    H[k++] = this->operator[](i);
                }

                // Build upper hull
                for (size_t i = this->size()-1, t = k+1; i > 0; --i)
                {
                    while (k >= t && cross(H[k-2], H[k-1], this->operator[](i-1)) <= tol) k--;
                    H[k++] = this->operator[](i-1);
                }

                H.resize(k-1);
                return H;
            }
        }

//        HullPointContainerType getPoints()
//        {
//            if(this->size()<=3)
//            {
//                return *this;
//            }
//            else
//            {
//                std::sort(this->begin(), this->end());
//                HullPointContainerType H;
//                size_t k(0);
//
//                // Build lower hull
//                for (size_t i = 0; i < this->size(); ++i)
//                {
//                    while (k >= 2 && cross(H[k-2], H[k-1], this->operator[](i)) <= tol) k--;
//                    k++;
//                    H.push_back(this->operator[](i));
//                }
//
//                // Build upper hull
//                for (size_t i = this->size()-1, t = k+1; i > 0; --i)
//                {
//                    while (k >= t && cross(H[k-2], H[k-1], this->operator[](i-1)) <= tol) k--;
//                    k++;
//                    H.push_back(this->operator[](i-1));
//                }
//
////                H.resize(k-1);
//                return H;
//            }
//        }
        
    };
    
}
#endif
