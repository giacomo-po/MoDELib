/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexChild_H_
#define model_SimplexChild_H_

#include <set>
#include <model/Mesh/SimplexTraits.h>
#include <model/Mesh/SimplexCompare.h>
#include <model/Mesh/BoundarySimplex.h>

namespace model
{
    
    template<short int dim,short int order>
    class SimplexChild : private std::set<const Simplex<dim,order+1>*,SimplexCompare<dim,order+1> >
    {
        
    public:
        //        typedef SimplexChild<dim,order> SimplexChildType;
        typedef Simplex<dim,order+1> ParentSimplexType;
        typedef typename SimplexTraits<dim,order+1>::SimplexIDType ParentSimplexIDType;
        typedef typename SimplexTraits<dim,order+1>::ScalarIDType ScalarIDType;
        typedef std::set<const ParentSimplexType*,SimplexCompare<dim,order+1> >  ParentContainerType;
        
    public:
        
        /**********************************************************************/
        void addToParents(const ParentSimplexType* const pP)
        {/*!@param[in] pID the ID of the parent Simplex
          * @param[in] pP the pointer to the parent Simplex
          *
          * Adds pP to the parentContainer.
          */
            //            const bool couldInsert(parentContainer.insert(pP).second);
            const bool couldInsert(this->insert(pP).second);
            assert(couldInsert && "COULD NOT INSERT PARENT IN parentContainer");
            
            // HERE WE SHOULD LOOP OVER PARENTS AND ADD pP TO THEIR NEIGHBORS
            
        }
        
        /**********************************************************************/
        void removeFromParents(const ParentSimplexType* const pP)
        {/*!@param[in] pID the ID of the parent Simplex
          *
          * Removes pID from the parentContainer.
          */
            //            const int nRemoved(parentContainer.erase(pP));
            const int nRemoved(this->erase(pP));
            assert(nRemoved==1 && "COULD NOT REMOVE PARENT IN parentContainer");
            
            
            // HERE WE SHOULD LOOP OVER PARENTS AND REMOVE pP FROM THEIR NEIGHBORS
            
        }
        
        /**********************************************************************/
        typename ParentContainerType::const_iterator parentBegin() const
        {
            return this->begin();
        }
        
        /**********************************************************************/
        typename ParentContainerType::const_iterator parentEnd() const
        {
            return this->end();
        }
        
        /**********************************************************************/
        const ParentContainerType& parents() const
        {
            return *this;
        }
        
        /**********************************************************************/
        bool isBoundarySimplex() const
        {
            return BoundarySimplex<dim,dim-order>::isBoundarySimplex(*this);
            
        }
        
        /**********************************************************************/
        bool isRegionBoundarySimplex() const
        {
            return BoundarySimplex<dim,dim-order>::isRegionBoundarySimplex(*this);
        }
        
        /**********************************************************************/
        std::set<int> regionIDs() const
        {
            std::set<int> temp;
            for(const auto& parent : parents())
            {
                const std::set<int> parentIDs(parent->regionIDs());
                for(const int& regionID : parentIDs)
                {
                    temp.insert(regionID);
                }
            }
            return temp;
        }
        
//        /**********************************************************************/
//        Eigen::Matrix<double,dim,1> outNormal() const
//        {
//            return BoundarySimplex<dim,dim-order>::outNormal(*this);
//        }
        
        //        /**********************************************************************/
        //        Eigen::Matrix<double,dim,1> outNormal() const
        //        {
        //            Eigen::Matrix<double,dim,1> temp(Eigen::Matrix<double,dim,1>::Zero());
        //            for(auto ele : *this)
        //            {
        //                const Eigen::Matrix<double,dim+1,1> bary(ele->simplex.pos2bary(P0));
        //                for(int k=0;k<dim+1;++k)
        //                {
        //                    if (std::fabs(bary(k))<FLT_EPSILON && ele->simplex.child(k).isBoundarySimplex())
        //                    {
        //                        temp += ele->simplex.nda.col(k).normalized();
        //                    }
        //                }
        //            }
        //            const double tempNorm(temp.norm());
        //            return tempNorm>FLT_EPSILON? (temp/tempNorm).eval() : Eigen::Matrix<double,dim,1>::Zero();
        //        }
        
        
    };
    
}	// close namespace
#endif
