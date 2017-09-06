/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexObserver_H_
#define model_SimplexObserver_H_

#include <memory> // shared_ptr
#include <map>
#include <array>
#include <model/Mesh/SimplexTraits.h>
#include <model/Utilities/CompareVectorsByComponent.h>

namespace model
{
	
    /**************************************************************************/
    /**************************************************************************/
	template<short int dim,short int order>
	struct SimplexObserver
    {
		        
        typedef Simplex<dim,order> SimplexType;
		typedef std::map<typename SimplexTraits<dim,order>::SimplexIDType, // key
        /*            */ const Simplex<dim,order>* const, // value
        /*            */ CompareVectorsByComponent<typename SimplexTraits<dim,order>::ScalarIDType,
        /*                                      */ SimplexTraits<dim,order>::nVertices> // key compare
        /*            */ >  SimplexMapType;
		
        typedef typename SimplexMapType::const_iterator const_iterator;
        typedef typename SimplexMapType::iterator iterator;
        
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        typedef std::shared_ptr<SimplexType> SharedPtrType;
		
        
		
        /**********************************************************************/
        static const SimplexMapType& simplices()
        {
            return simplexMap;
        }
        
		/**********************************************************************/
		static size_t size()
        {
			return simplexMap.size();
		}
		
        /**********************************************************************/
		static typename SimplexMapType::const_iterator simplexBegin()
        {
			return simplexMap.begin();
		}
		
		/**********************************************************************/
		static typename SimplexMapType::const_iterator simplexEnd()
        {
			return simplexMap.end();
		}
        
		/**********************************************************************/
		static SharedPtrType pSimplex(const SimplexIDType& xID)
        {
			typename SimplexMapType::const_iterator iter(simplexMap.find(xID));
			return (iter!=simplexMap.end())? (*(iter->second->parentBegin()))->child(xID) : SharedPtrType(new SimplexType(xID));
		}
        
        /**********************************************************************/
        static const SimplexType& simplex(const SimplexIDType& xID)
        {
            typename SimplexMapType::const_iterator iter(simplexMap.find(xID));
            assert(iter!=simplexMap.end() && "Simplex not in simplexMap");
            return *(iter->second);
//            return (iter!=simplexMap.end())? (*(iter->second->parentBegin()))->child(xID) : SharedPtrType(new SimplexType(xID));
        }
		
        
        /**********************************************************************/
        static std::array<std::shared_ptr<Simplex<dim,order-1> >,SimplexTraits<dim,order>::nFaces> faces(const SimplexIDType& xID)
        {
            std::array<std::shared_ptr<Simplex<dim,order-1> >,SimplexTraits<dim,order>::nFaces> temp;
            for (size_t j=0;j<SimplexTraits<dim,order>::nFaces;++j)
            {
                temp[j]=SimplexObserver<dim,order-1>::pSimplex(SimplexTraits<dim,order>::faceID(xID,j));
            }
            return temp;
        }
        
        /**********************************************************************/
        static void insertSimplex(const Simplex<dim,order>& s)
        {            
            const bool couldInsert(simplexMap.insert(std::make_pair(s.xID,&s)).second);
            assert(couldInsert && "COULD NOT INSERT SIMPLEX IN SIPLEX-MAP");
        }
        
        /**********************************************************************/
        static void removeSimplex(const Simplex<dim,order>& s)
        {            
            const int nRemoved(simplexMap.erase(s.xID));
            assert(nRemoved==1 && "COULD NOT REMOVE SIMPLEX FROM SIPLEX-MAP");
        }
        
        /**********************************************************************/
        static typename SimplexMapType::const_iterator find(const SimplexIDType& key)
        {
            return simplexMap.find(key);
        }
        
        
	private:
		static  SimplexMapType simplexMap;
        
	};
    
	// Declare static data members
	template<short int dim,short int order>
    std::map<typename SimplexTraits<dim,order>::SimplexIDType, // key
    /*            */ const Simplex<dim,order>* const, // value
    /*            */ CompareVectorsByComponent<typename SimplexTraits<dim,order>::ScalarIDType,
    /*                                      */ SimplexTraits<dim,order>::nVertices> // key compare
    /*            */ >  SimplexObserver<dim,order>::simplexMap;
    
}	// close namespace
#endif
