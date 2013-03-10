/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_SUBNETWORK_H_
#define model_SUBNETWORK_H_

#include <iomanip>
#include <assert.h>
#include <map>
#include <utility>

#include <boost/ptr_container/ptr_map.hpp>
//#include <boost/tuple/tuple.hpp>
#include <boost/utility.hpp>

#include <model/Network/Operations/includeNetworkOperations.h>


#include <model/Utilities/AddressBook.h>
#include <model/Utilities/TypeTraits.h>

namespace model {

		/************************************************************/
		/************************************************************/
		template <typename Derived>
		class SubNetwork : boost::noncopyable,
		/*              */ private	std::map<size_t,typename TypeTraits<Derived>::NodeType* const>,
		/*              */ private	std::map<std::pair<size_t,size_t>,typename TypeTraits<Derived>::LinkType* const>,
		/*              */ public	AddressBook<Derived>{
			
			/*! \brief CRTP-based class template which 
			 */
			
			#include <model/Network/NetworkTypedefs.h>
			
		public:	
			
			/************************************************************/
			/* Constructor with pointer to Node *************************/
			SubNetwork(NodeType* const pN){
				add(pN);
			}
			
			/************************************************************/
			/* Constructor with pointer to Link *************************/
			SubNetwork(LinkType* const pL){
				add(pL);
			}
			
			/************************************************************/
			// Destructor
			~SubNetwork(){
				assert(SubNetworkNodeContainerType::empty() && "Destroying non-empty SubNetwork! Subnetwork contains Vertices.");	// make sure that only empty subnetworks are deleted
				assert(SubNetworkLinkContainerType::empty() && "Destroying non-empty SubNetwork! Subnetwork contains Edges.");		// make sure that only empty subnetworks are deleted
			}
			
			/************************************************************/
			// NODE OPERATIONS
			/************************************************************/
			// nodeOrder
			size_t nodeOrder() const {
				return SubNetworkNodeContainerType::size();
			}
			
			/************************************************************/
			// nodeBegin
			typename SubNetworkNodeContainerType::iterator nodeBegin() {
				return SubNetworkNodeContainerType::begin();
			}
			
			/************************************************************/
			// nodeBegin
			typename SubNetworkNodeContainerType::const_iterator nodeBegin() const {
				return SubNetworkNodeContainerType::begin();
			}
			
			/************************************************************/
			// nodeEnd
			typename SubNetworkNodeContainerType::iterator nodeEnd() {
				return SubNetworkNodeContainerType::end();
			}
			
			/************************************************************/
			// nodeEnd
			typename SubNetworkNodeContainerType::const_iterator nodeEnd() const {
				return SubNetworkNodeContainerType::end();
			}
			
			/************************************************************/
			// snID
			size_t snID(const NodeType* const & pN) const {
				return std::distance(SubNetworkNodeContainerType::begin(), SubNetworkNodeContainerType::find(pN->sID) );
			}
			
			/************************************************************/
			void add(NodeType* const pN){
				assert(SubNetworkNodeContainerType::insert(std::make_pair(pN->sID,pN)).second);
			}
			
			/************************************************************/
			void remove(NodeType* const pN){
				assert(SubNetworkNodeContainerType::erase(pN->sID)==1);
			}
			
			/************************************************************/
			// LINK OPERATIONS
			/************************************************************/
			// linkOrder
			size_t linkOrder() const {
				return SubNetworkLinkContainerType::size();
			}
			
			/************************************************************/
			// linkBegin
			typename SubNetworkLinkContainerType::iterator linkBegin()  {
				return SubNetworkLinkContainerType::begin();
			}
			
			typename SubNetworkLinkContainerType::const_iterator linkBegin() const {
				return SubNetworkLinkContainerType::begin();
			}
			
			/************************************************************/
			// linkEnd
			typename SubNetworkLinkContainerType::iterator linkEnd()
            {/*! @return An iterator to the end of the link container
              */
				return SubNetworkLinkContainerType::end();
			}
			
			typename SubNetworkLinkContainerType::const_iterator linkEnd() const
            {
				return SubNetworkLinkContainerType::end();
			}
			
			/************************************************************/
			// snID (link)
			size_t snID(const LinkType* const & pL) const
            {/*! @param[in] pL
              *  @return The ID of pL in this SubNetwork
              */
				return std::distance(SubNetworkLinkContainerType::begin(), SubNetworkLinkContainerType::find(pL->nodeIDPair  ) );
			}
			
			/************************************************************/
			void add(LinkType* const pL)
            {/*! @param[in] pL
              */
				assert(SubNetworkLinkContainerType::insert(std::make_pair(pL->nodeIDPair ,pL)).second);
			}
			
			/************************************************************/
			void remove(LinkType* const pL)
            {/*! @param[in] pL
              */
				assert(SubNetworkLinkContainerType::erase(pL->nodeIDPair )==1);
			}
			
		};
		
		/************************************************************/
		/************************************************************/
} // namespace model
#endif
