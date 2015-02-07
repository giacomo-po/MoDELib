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
        typedef Simplex<dim,order+1> ParentSimplexType;
        typedef typename SimplexTraits<dim,order+1>::SimplexIDType ParentSimplexIDType;
        typedef typename SimplexTraits<dim,order+1>::ScalarIDType ScalarIDType;
        typedef std::set<const ParentSimplexType*,SimplexCompare<dim,order+1> >  ParentContainerType;

//        typedef std::set<const ParentSimplexType*,SimplexCompare<dim,order+1> >  ParentContainerType;

        
//    private:
//        ParentContainerType parentContainer;

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
            //            return parentContainer.begin();
        }
        
        /**********************************************************************/
        typename ParentContainerType::const_iterator parentEnd() const
        {
            return this->end();
//            return parentContainer.end();
        }
        
        /**********************************************************************/
        const ParentContainerType& parents() const
        {
            return *this;
//            return parentContainer;
        }
        
        /**********************************************************************/
        bool isBoundarySimplex() const
        {
            return BoundarySimplex<dim,dim-order>::template isBoundarySimplex<SimplexChild>(*this);
        }
        
        /**********************************************************************/
        bool isRegionBoundarySimplex() const
        {
            return BoundarySimplex<dim,dim-order>::template isRegionBoundarySimplex<SimplexChild>(*this);
        }
        
    };
    
}	// close namespace
#endif