/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_NetworkComponent_h_
#define model_NetworkComponent_h_

#include <iomanip>
#include <assert.h>
#include <map>
#include <utility>
//#include <boost/ptr_container/ptr_map.hpp>
//#include <boost/tuple/tuple.hpp>
//#include <boost/utility.hpp>
//#include <model/Network/NetworkLink.h>
#include <model/Network/Operations/includeNetworkOperations.h>
#include <model/Utilities/AddressBook.h>
#include <model/Utilities/TypeTraits.h>
#include <model/Utilities/NonCopyable.h>

namespace model
{
    // Class pre-declaration
    template<typename Derived>
    class NetworkLink;
    
    /************************************************************/
    /************************************************************/
    /*! \brief Class template containing the network component
     */
    template <typename NodeType, typename LinkType>
    class NetworkComponent : public NonCopyable,
// boost::noncopyable,
    /*              */ private	std::map<size_t,NodeType* const>,
    /*              */ private	std::map<std::pair<size_t,size_t>,LinkType* const>,
    /*              */ public	AddressBook<NetworkComponent<NodeType,LinkType> >
    {
        
        
        
        //			#include <model/Network/NetworkTypedefs.h>
        
        friend class NetworkLink<LinkType>; // allow NetworkLink to access std::map<size_t,NodeType* const>
        
        
    public:
        
        typedef std::map<size_t,NodeType* const>						NetworkComponentNodeContainerType;
        typedef std::map<std::pair<size_t,size_t>,LinkType* const>		NetworkComponentLinkContainerType;
        
        
        /************************************************************/
        /* Constructor with pointer to Node *************************/
        NetworkComponent(NodeType* const pN){
            add(pN);
        }
        
//        /************************************************************/
//        /* Constructor with pointer to Link *************************/
//        NetworkComponent(LinkType* const pL){
//            add(pL);
//        }
        
        /**********************************************************************/
        ~NetworkComponent(){
            assert(NetworkComponentNodeContainerType::empty() && "Destroying non-empty NetworkComponent! NetworkComponent contains Vertices.");	// make sure that only empty NetworkComponents are deleted
            assert(NetworkComponentLinkContainerType::empty() && "Destroying non-empty NetworkComponent! NetworkComponent contains Edges.");		// make sure that only empty NetworkComponents are deleted
        }
        
        /**********************************************************************/
        size_t nodeOrder() const
        {
            return NetworkComponentNodeContainerType::size();
        }
        
        /**********************************************************************/
        typename NetworkComponentNodeContainerType::iterator nodeBegin()
        {
            return NetworkComponentNodeContainerType::begin();
        }
        
        /**********************************************************************/
        typename NetworkComponentNodeContainerType::const_iterator nodeBegin() const
        {
            return NetworkComponentNodeContainerType::begin();
        }
        
        /**********************************************************************/
        typename NetworkComponentNodeContainerType::iterator nodeEnd()
        {
            return NetworkComponentNodeContainerType::end();
        }
        
        /**********************************************************************/
        typename NetworkComponentNodeContainerType::const_iterator nodeEnd() const
        {
            return NetworkComponentNodeContainerType::end();
        }
        
        /**********************************************************************/
        size_t snID(const NodeType* const & pN) const
        {
            return std::distance(NetworkComponentNodeContainerType::begin(), NetworkComponentNodeContainerType::find(pN->sID) );
        }
        
        /**********************************************************************/
        void add(NodeType* const pN)
        {
            assert(NetworkComponentNodeContainerType::insert(std::make_pair(pN->sID,pN)).second);
        }
        
        /**********************************************************************/
        void remove(NodeType* const pN)
        {
            assert(NetworkComponentNodeContainerType::erase(pN->sID)==1);
        }
        
        /**********************************************************************/
        size_t linkOrder() const
        {
            return NetworkComponentLinkContainerType::size();
        }
        
        /**********************************************************************/
        typename NetworkComponentLinkContainerType::iterator linkBegin()
        {
            return NetworkComponentLinkContainerType::begin();
        }
        
        typename NetworkComponentLinkContainerType::const_iterator linkBegin() const {
            return NetworkComponentLinkContainerType::begin();
        }
        
        /**********************************************************************/
        typename NetworkComponentLinkContainerType::iterator linkEnd()
        {/*! @return An iterator to the end of the link container
          */
            return NetworkComponentLinkContainerType::end();
        }
        
        typename NetworkComponentLinkContainerType::const_iterator linkEnd() const
        {
            return NetworkComponentLinkContainerType::end();
        }
        
        /**********************************************************************/
        size_t snID(const LinkType* const & pL) const
        {/*! @param[in] pL
          *  @return The ID of pL in this NetworkComponent
          */
            return std::distance(NetworkComponentLinkContainerType::begin(), NetworkComponentLinkContainerType::find(pL->nodeIDPair  ) );
        }
        
        /**********************************************************************/
        void add(LinkType* const pL)
        {/*! @param[in] pL
          */
            assert(NetworkComponentLinkContainerType::insert(std::make_pair(pL->nodeIDPair ,pL)).second);
        }
        
        /**********************************************************************/
        void remove(LinkType* const pL)
        {/*! @param[in] pL
          */
            assert(NetworkComponentLinkContainerType::erase(pL->nodeIDPair )==1);
        }
        
        /**********************************************************************/
        const std::map<std::pair<size_t,size_t>,LinkType* const>& links() const
        {
            return *this;
        }

        /**********************************************************************/
        const std::map<size_t,NodeType* const>& nodes() const
        {
            return *this;
        }
        
    };
} // namespace model
#endif
