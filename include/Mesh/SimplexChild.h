/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexChild_H_
#define model_SimplexChild_H_

#include <memory>
#include <map>
#include <SimplexTraits.h>
//#include <SimplexCompare.h>
#include <BoundarySimplex.h>
// gcc8 incomplete type problem https://www.reddit.com/r/cpp/comments/8joovp/gcc_8_cant_forward_declare_with_associative/
namespace model
{
    
    template<short int dim,short int order>
    class SimplexChild : private std::map<typename SimplexTraits<dim,order+1>::SimplexIDType,const Simplex<dim,order+1>* const>
    /*               */, private std::map<typename SimplexTraits<dim,order  >::SimplexIDType,const Simplex<dim,order  >* const>

    {
        
    public:
        typedef typename SimplexTraits<dim,order>::SimplexType   SimplexType;
        typedef std::map<typename SimplexTraits<dim,order  >::SimplexIDType,const Simplex<dim,order  >* const>  SiblingsContainerType;
        typedef typename SimplexTraits<dim,order+1>::SimplexType   ParentSimplexType;
        typedef std::map<typename SimplexTraits<dim,order+1>::SimplexIDType,const Simplex<dim,order+1>* const>  ParentContainerType;

        typedef typename SimplexTraits<dim,order+1>::SimplexIDType ParentSimplexIDType;
        
    private:
        
        void make_siblings()
        {
            siblings().clear();
            for(const auto& parent : parents())
            {
                for(int c=0; c<ParentSimplexType::nFaces;++c)
                {
                    siblings().emplace(parent.second->child(c).xID,&parent.second->child(c));
                }
            }
        }
        
        
    public:
                
        void addToParents(const ParentSimplexType* const pP)
        {/*!@param[in] pID the ID of the parent Simplex
          * @param[in] pP the pointer to the parent Simplex
          *
          * Adds pP to the parentContainer.
          */
            //            const bool couldInsert(parentContainer.insert(pP).second);
            const bool couldInsert(parents().emplace(pP->xID,pP).second);
            if(!couldInsert)
            {
                throw std::runtime_error("COULD NOT INSERT SIMPLEX IN SIPLEX-MAP");
            }

            // HERE WE SHOULD LOOP OVER PARENTS AND ADD pP TO THEIR NEIGHBORS
            //update();
            make_siblings();
        }
        
        /**********************************************************************/
        void removeFromParents(const ParentSimplexType* const pP)
        {/*!@param[in] pID the ID of the parent Simplex
          *
          * Removes pID from the parentContainer.
          */
            const int nRemoved(parents().erase(pP->xID));
            if(nRemoved!=1)
            {
                throw std::runtime_error("COULD NOT REMOVE SIMPLEX FROM SIPLEX-MAP");
            }

            
            // HERE WE SHOULD LOOP OVER PARENTS AND REMOVE pP FROM THEIR NEIGHBORS
            //update();
            make_siblings();
            
        }
        
        ParentContainerType& parents()
        {
            return *this;
        }

        const ParentContainerType& parents() const
        {
            return *this;
        }
        
        const SiblingsContainerType& siblings() const
        {
            return *this;
        }
        
        SiblingsContainerType& siblings()
        {
            return *this;
        }
        
        bool isBoundarySimplex() const
        {
            return BoundarySimplex<dim,dim-order>::isBoundarySimplex(*this);
            
        }
        
        bool isRegionBoundarySimplex() const
        {
            return BoundarySimplex<dim,dim-order>::isRegionBoundarySimplex(*this);
        }

        bool isInRegion(const size_t& rID) const
        {
            const std::set<size_t> rIDs=regionIDs();
            return rIDs.find(rID)!=rIDs.end();
        }
        
        std::set<size_t> regionIDs() const
        {
            std::set<size_t> temp;
            for(const auto& parent : parents())
            {
                const std::set<size_t> parentIDs(parent.second->regionIDs());
                for(const size_t& regionID : parentIDs)
                {
                    temp.insert(regionID);
                }
            }
            return temp;
        }
        
        
    };
    
}
#endif
