/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VERTEXCONNECTION_H_
#define model_VERTEXCONNECTION_H_

#include <memory>   // for auto_ptr
#include <utility>  // for std::pair
#include <assert.h>
#include <map>
#include <tuple>


namespace model
{
	
	template <typename VertexType, typename EdgeType>
	class VertexConnection
    {
		
        typedef std::map<size_t,VertexType> NetworkVertexMapType;
        typedef std::map<std::pair<size_t,size_t>,EdgeType> NetworkEdgeMapType;

        //! A reference to the network vertex map
		NetworkVertexMapType& networkVertexMapRef;
		
        //! A reference to the network Edge map
		NetworkEdgeMapType& networkEdgeMapRef;
		
	public:
		
		/**********************************************************************/
		VertexConnection(NetworkVertexMapType& networkVertexMapRef_in,
		/*            */ NetworkEdgeMapType&     networkEdgeMapRef_in) :
        /* init list  */ networkVertexMapRef(networkVertexMapRef_in),
		/* init list  */ networkEdgeMapRef(networkEdgeMapRef_in){}
		
		/**********************************************************************/
        template <typename ...EdgeArgTypes>
		bool connect(const size_t& i, const size_t& j, const EdgeArgTypes&... edgeArgs)
        {/*! Connects vertex i to vertex j creating the edge i->j.
          * Requirements are:
          * - i and j are distinct (asserts otherwise)
          * - vertex i exists (asserts otherwise)
          * - vertex j exists (asserts otherwise)
          * - edge i->j does not exist (asserts otherwise)
          * - edge j->i does not exist (asserts otherwise)
          */
			assert(i!=j && "TRYING TO CONNECT A VERTEX TO ITSELF.");
			typename NetworkVertexMapType::iterator Ni(networkVertexMapRef.find(i));
			assert(Ni!=networkVertexMapRef.end() && "TRYING TO CONNECT (i->j): VERTEX i DOES NOT EXIST.");
			typename NetworkVertexMapType::iterator Nj(networkVertexMapRef.find(j));
			assert(Nj!=networkVertexMapRef.end() && "TRYING TO CONNECT (i->j): VERTEX j DOES NOT EXIST.");			
			assert(networkEdgeMapRef.find(std::make_pair(i,j))==networkEdgeMapRef.end() && "EDGE ALREADY EXISTS. CANNOT CONNECT");
			assert(networkEdgeMapRef.find(std::make_pair(j,i))==networkEdgeMapRef.end() && "OPPOSITE EDGE ALREADY EXISTS. CANNOT CONNECT");
            const bool success=networkEdgeMapRef.emplace(std::piecewise_construct,
                                                         std::make_tuple(i,j),
                                                         std::make_tuple(std::make_pair(&Ni->second,&Nj->second), edgeArgs...) ).second;
            assert(success && "CANNOT INSERT EDGE IN networkEdgeMapRef.");
            return success;
		}
		
		/**********************************************************************/
        template<bool removeIsolatedNodes>
		bool disconnect(const size_t& i, const size_t& j)
        {
			typename NetworkEdgeMapType::iterator iterIJ(networkEdgeMapRef.find(std::make_pair(i,j))); // look for edge (i->j)
			if (iterIJ!=networkEdgeMapRef.end())
            {   // edge (i->j) is in the edgeMap
				networkEdgeMapRef.erase(iterIJ);	// remove (i->j)
			}
			else
            { // edge (i->j) is not in the edgeMap
				typename NetworkEdgeMapType::iterator iterJI(networkEdgeMapRef.find(std::make_pair(j,i))); // look for edge (j->i)
				assert(iterJI!=networkEdgeMapRef.end() && "NEITHER (i->j) nor (j->i) ARE IN networkEdgeMapRef.");
				networkEdgeMapRef.erase(iterJI);	// remove (j->i)
			}
			if (removeIsolatedNodes)
            {
				typename NetworkVertexMapType::iterator Ni(networkVertexMapRef.find(i));
				assert(Ni!=networkVertexMapRef.end() && "NODE i DOES NOT EXIST");
				if (Ni->second.is_isolated())
                {
					networkVertexMapRef.erase(Ni); // WARNING: erasing by key, i.e. networkVertexMapRef.erase(i), gives a bug
					//	WHY IS THIS TRIGGERED????				assert(networkVertexMapRef.find(i)==networkVertexMapRef.end() && "NODE i IS STILL IN networkVertexMapRef AFTER ERASE.");
				}
				typename NetworkVertexMapType::iterator Nj(networkVertexMapRef.find(j));
				assert(Nj!=networkVertexMapRef.end() && "NODE j DOES NOT EXIST");
				if (Nj->second.is_isolated())
                {
					networkVertexMapRef.erase(Nj); // WARNING: erasing by key, i.e. networkVertexMapRef.erase(j), gives a bug
					// 	WHY IS THIS TRIGGERED????					assert(networkVertexMapRef.find(j)==networkVertexMapRef.end() && "NODE j IS STILL IN networkVertexMapRef AFTER ERASE.");
				}
			}
			return true;
		}
		
		/**********************************************************************/
        template<bool removeIsolatedNodes>
		bool remove(const size_t& k)
        {
			
			typename NetworkVertexMapType::iterator Nk(networkVertexMapRef.find(k));
			assert(Nk!=networkVertexMapRef.end() && "REMOVING NON-EXISTING VERTEX.");
						
			// Define BoolLinkMFPsize_t as a member-function-pointer of class EdgeType with size_t input and bool output
			typedef bool (EdgeType::*BoolLinkMFPsize_t)(const size_t &) const;
			BoolLinkMFPsize_t MFP = &EdgeType::isIncident;
			disconnect_if<removeIsolatedNodes>(MFP,k);
			//	disconnect_if<removeIsolatedNodes>(&EdgeType::isIncident,k);	// WHY IS THIS NOT WORKING????????
			
			typename NetworkVertexMapType::iterator Nkk(networkVertexMapRef.find(k));	
			if (Nkk!=networkVertexMapRef.end())
            {
				networkVertexMapRef.erase(Nkk); // WARNING: erasing by key, i.e. networkVertexMapRef.erase(j), gives a bug
			}
			assert(networkVertexMapRef.find(k)==networkVertexMapRef.end() && "NODE k IS STILL IN networkVertexMapRef AFTER ERASE.");
						
			return true;
		}
		
		
		/**********************************************************************/
		template<bool removeIsolatedNodes>
		size_t disconnect_if(bool (EdgeType::*Lfptr)(void) const)
        {
			size_t count = 0;
			for(typename NetworkEdgeMapType::iterator edgeIter=networkEdgeMapRef.begin(); edgeIter!=networkEdgeMapRef.end();)
            {
				if ( (edgeIter->second.*Lfptr)() )
                {
					typename NetworkEdgeMapType::iterator toBeDisconnected(edgeIter);
					++edgeIter; // increment iterator before erasing
					count+=disconnect<removeIsolatedNodes>(toBeDisconnected->second.source->sID,toBeDisconnected->second.sink->sID);
				}
				else
                {
					++edgeIter;
				}
			}
			return count;
		}
		
		/**********************************************************************/
		template<bool removeIsolatedNodes, typename T>
		size_t disconnect_if(bool (EdgeType::*Lfptr)(const T &) const, const T & input)
        {
			size_t count = 0;
			for(typename NetworkEdgeMapType::iterator edgeIter=networkEdgeMapRef.begin(); edgeIter!=networkEdgeMapRef.end();)
            {
				if ( (edgeIter->second.*Lfptr)(input) )
                {
					typename NetworkEdgeMapType::iterator toBeDisconnected(edgeIter);
					++edgeIter; // increment iterator before erasing
					count+=disconnect<removeIsolatedNodes>(toBeDisconnected->second.source->sID,toBeDisconnected->second.sink->sID);
				}
				else
                {
					++edgeIter;
				}
			}
			return count;
		}
		
	};
	
} // namespace model
#endif
