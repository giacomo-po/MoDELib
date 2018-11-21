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
#include <SimplexTraits.h>
#include <CompareVectorsByComponent.h>
#include <MPIcout.h>

namespace model
{

    template<int dim>
    struct SimplicialMesh;

    template<short int dim,short int order>
    struct SimplexObserverBase;
    
    template<short int dim,short int order>
    struct SimplexObserverRecursion;
    
    /**************************************************************************/
    /**************************************************************************/
    template<short int dim>
    struct SimplexObserver : public SimplexObserverRecursion<dim,dim>
    {
        
        template<short int o>
        SimplexObserverBase<dim,o>& observer()
        {
            static_assert(o>=0 && o<=dim,"order must be >=0 and <= dim");
            return *this;
        }
        
        /**********************************************************************/
        template<short int o>
        const SimplexObserverBase<dim,o>& observer() const
        {
            static_assert(o>=0 && o<=dim,"order must be >=0 and <= dim");
            return *this;
        }
        
        /**********************************************************************/
        void info() const
        {
            SimplexObserverRecursion<dim,dim>::info();
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template< short int dim,short int order>
    struct SimplexObserverRecursion : public SimplexObserverBase<dim,order>,
    /*                             */ public SimplexObserverRecursion<dim,order-1>
    {
        
        /**********************************************************************/
        void info() const
        {
            SimplexObserverRecursion<dim,order-1>::info();
            
            size_t nT(0);
            size_t nB(0);
            double volT=0.0;
            double volB=0.0;
            
            for (auto& pSimplex : SimplexObserverBase<dim,order>::simplices())
            {
                nT++;
                volT+=pSimplex.second->vol0;
                
                if(pSimplex.second->isBoundarySimplex())
                {
                    nB++;
                    volB+=pSimplex.second->vol0;
                }
            }
            
            model::cout<<"    Simplex<"<<dim<<","<<order  <<"> #="<<nT<<", vol="<<volT;
            model::cout<<"     (boundary #="<<nB<<", vol="<<volB<<")\n";
            
        }
        
        
    };
    
    template< short int dim>
    struct SimplexObserverRecursion<dim,0> : public SimplexObserverBase<dim,0>
    {
        
        /**********************************************************************/
        void info() const
        {
            
            size_t nT(0);
            size_t nB(0);
            
            for (auto& pSimplex : SimplexObserverBase<dim,0>::simplices())
            {
                nT++;
                
                if(pSimplex.second->isBoundarySimplex())
                {
                    nB++;
                    assert(fabs(pSimplex.second->outNormal().norm()-1.0)<FLT_EPSILON);
                }
                else
                {
                    assert(pSimplex.second->outNormal().norm()<FLT_EPSILON);
                }
            }
            
            model::cout<<"    Simplex<"<<dim<<","<<0  <<"> #="<<nT;
            model::cout<<"     (boundary #="<<nB<<")\n";
        }
        
    };

    
    
    /**************************************************************************/
    /**************************************************************************/
	template<short int dim,short int order>
    struct SimplexObserverBase : public std::map<typename SimplexTraits<dim,order>::SimplexIDType, // key
    /*                                        */ const Simplex<dim,order>* const, // value
    /*                                        */ CompareVectorsByComponent<typename SimplexTraits<dim,order>::ScalarIDType,
    /*                                        */ SimplexTraits<dim,order>::nVertices>> // key compare
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
        const SimplexMapType& simplices() const
        {
            return *this;
        }
        
//		/**********************************************************************/
//		static size_t size()
//        {
//			return simplexMap.size();
//		}
		
//        /**********************************************************************/
//		static typename SimplexMapType::const_iterator simplexBegin()
//        {
//			return simplexMap.begin();
//		}
//		
//		/**********************************************************************/
//		static typename SimplexMapType::const_iterator simplexEnd()
//        {
//			return simplexMap.end();
//		}
        
		/**********************************************************************/
        SharedPtrType sharedSimplex(SimplicialMesh<dim>* const mesh,
                                           const SimplexIDType& xID)
        {
			typename SimplexMapType::const_iterator iter(this->find(xID));
			return (iter!=this->end())? (*(iter->second->parentBegin()))->child(xID) : SharedPtrType(new SimplexType(mesh,xID));
		}
        
        /**********************************************************************/
        const SimplexType& simplex(const SimplexIDType& xID) const
        {
            typename SimplexMapType::const_iterator iter(this->find(xID));
            assert(iter!=this->end() && "Simplex not in simplexMap");
            return *(iter->second);
//            return (iter!=simplexMap.end())? (*(iter->second->parentBegin()))->child(xID) : SharedPtrType(new SimplexType(xID));
        }
		
        
        /**********************************************************************/
        std::array<std::shared_ptr<Simplex<dim,order-1> >,SimplexTraits<dim,order>::nFaces> faces(SimplicialMesh<dim>* const mesh,
                                                                                                         const SimplexIDType& xID)
        {
            std::array<std::shared_ptr<Simplex<dim,order-1> >,SimplexTraits<dim,order>::nFaces> temp;
            for (size_t j=0;j<SimplexTraits<dim,order>::nFaces;++j)
            {
                temp[j]=mesh->template observer<order-1>().sharedSimplex(mesh,SimplexTraits<dim,order>::faceID(xID,j));
            }
            return temp;
        }
        
        /**********************************************************************/
        void insertSimplex(const Simplex<dim,order>& s)
        {            
            const bool couldInsert(this->emplace(s.xID,&s).second);
            assert(couldInsert && "COULD NOT INSERT SIMPLEX IN SIPLEX-MAP");
        }
        
        /**********************************************************************/
        void removeSimplex(const Simplex<dim,order>& s)
        {            
            const int nRemoved(this->erase(s.xID));
            assert(nRemoved==1 && "COULD NOT REMOVE SIMPLEX FROM SIPLEX-MAP");
        }
        
//        /**********************************************************************/
//        static typename SimplexMapType::const_iterator find(const SimplexIDType& key)
//        {
//            return simplexMap.find(key);
//        }
//        
//        
//	private:
//		static  SimplexMapType simplexMap;
        
	};
    
//	// Declare static data members
//	template<short int dim,short int order>
//    std::map<typename SimplexTraits<dim,order>::SimplexIDType, // key
//    /*            */ const Simplex<dim,order>* const, // value
//    /*            */ CompareVectorsByComponent<typename SimplexTraits<dim,order>::ScalarIDType,
//    /*                                      */ SimplexTraits<dim,order>::nVertices> // key compare
//    /*            */ >  SimplexObserver<dim,order>::simplexMap;
    
}	// close namespace
#endif
