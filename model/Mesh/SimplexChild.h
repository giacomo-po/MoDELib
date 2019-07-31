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
//    class SimplexChild : private std::set<const Simplex<dim,order+1>*,SimplexCompare<dim,order+1> >
//    /*               */, private std::set<const Simplex<dim,order  >*,SimplexCompare<dim,order  > >
//    class SimplexChild : private std::set<const typename SimplexTraits<dim,order+1>::SimplexType*,SimplexCompare<dim,order+1> >
//    /*               */, private std::set<const typename SimplexTraits<dim,order  >::SimplexType*,SimplexCompare<dim,order  > >
    class SimplexChild : private std::map<typename SimplexTraits<dim,order+1>::SimplexIDType,const Simplex<dim,order+1>* const>
    /*               */, private std::map<typename SimplexTraits<dim,order  >::SimplexIDType,const Simplex<dim,order  >* const>

    {
        
    public:
        //        typedef SimplexChild<dim,order> SimplexChildType;
        typedef typename SimplexTraits<dim,order>::SimplexType   SimplexType;
        typedef std::map<typename SimplexTraits<dim,order  >::SimplexIDType,const Simplex<dim,order  >* const>  SiblingsContainerType;

//        typedef Simplex<dim,order+1> ParentSimplexType;
        typedef typename SimplexTraits<dim,order+1>::SimplexType   ParentSimplexType;
        typedef std::map<typename SimplexTraits<dim,order+1>::SimplexIDType,const Simplex<dim,order+1>* const>  ParentContainerType;

        typedef typename SimplexTraits<dim,order+1>::SimplexIDType ParentSimplexIDType;
//        typedef typename SimplexTraits<dim,order+1>::ScalarIDType ScalarIDType;
        
    private:
        
        //        std::set<int> _regionIDs;
//        SiblingsContainerType _siblings;
        
        //        /**********************************************************************/
        //        void make_regionIDs()
        //        {
        ////            std::set<int> temp;
        //            _regionIDs.clear();
        //            for(const auto& parent : parents())
        //            {
        //                const std::set<int> parentIDs(parent->regionIDs());
        //                for(const int& regionID : parentIDs)
        //                {
        //                    _regionIDs.insert(regionID);
        //                }
        //            }
        ////            return temp;
        //        }
        
        /**********************************************************************/
        void make_siblings()
        {
            //            SiblingsContainerType temp;
            siblings().clear();
            for(const auto& parent : parents())
            {
                for(int c=0; c<ParentSimplexType::nFaces;++c)
                {
                    siblings().emplace(parent.second->child(c).xID,&parent.second->child(c));
                }
            }
            //            return temp;
        }
        
        
    public:
        
        
        //        /**********************************************************************/
        //        SimplexChild()
        //        {
        //            make_regionIDs();
        //            make_siblings();
        //        }
        //
        ////        void update()
        ////        {
        ////            make_regionIDs();
        ////            make_siblings();
        ////        }
        
        
        /**********************************************************************/
        void addToParents(const ParentSimplexType* const pP)
        {/*!@param[in] pID the ID of the parent Simplex
          * @param[in] pP the pointer to the parent Simplex
          *
          * Adds pP to the parentContainer.
          */
            //            const bool couldInsert(parentContainer.insert(pP).second);
            const bool couldInsert(parents().emplace(pP->xID,pP).second);
            assert(couldInsert && "COULD NOT INSERT PARENT IN parentContainer");
            
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
            //            const int nRemoved(parentContainer.erase(pP));
            const int nRemoved(parents().erase(pP->xID));
            assert(nRemoved==1 && "COULD NOT REMOVE PARENT IN parentContainer");
            
            
            // HERE WE SHOULD LOOP OVER PARENTS AND REMOVE pP FROM THEIR NEIGHBORS
            //update();
            make_siblings();
            
        }
        
//        /**********************************************************************/
//        typename ParentContainerType::const_iterator parentBegin() const
//        {
//            return parents().begin();
//        }
//
//        /**********************************************************************/
//        typename ParentContainerType::const_iterator parentEnd() const
//        {
//            return parents().end();
//        }

        /**********************************************************************/
        ParentContainerType& parents()
        {
            return *this;
        }

        /**********************************************************************/
        const ParentContainerType& parents() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const SiblingsContainerType& siblings() const
        {
            return *this;
        }
        
        /**********************************************************************/
        SiblingsContainerType& siblings()
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
        bool isInRegion(const int& rID) const
        {
            const std::set<int> rIDs=regionIDs();
            return rIDs.find(rID)!=rIDs.end();
        }
        
        /**********************************************************************/
        std::set<int> regionIDs() const
        {
            std::set<int> temp;
            for(const auto& parent : parents())
            {
                const std::set<int> parentIDs(parent.second->regionIDs());
                for(const int& regionID : parentIDs)
                {
                    temp.insert(regionID);
                }
            }
            return temp;
        }
        
        
    };
    
}
#endif
