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
#include <model/MPI/MPIcout.h>

namespace model
{

//    template<int dim>
//    struct SimplicialMesh;

    
    template<typename Derived,short int dim,short int order>
    struct SimplexObserverBase;
    
    template<typename Derived,short int dim,short int order>
    struct SimplexObserverRecursion;
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename Derived,short int dim>
    struct SimplexObserver : public SimplexObserverRecursion<Derived,dim,dim>
    {
        
        template<short int o>
        SimplexObserverBase<Derived,dim,o>& observer()
        {
            static_assert(o>=0 && o<=dim,"order must be >=0 and <= dim");
            return *this;
        }
        
        template<short int o>
        const SimplexObserverBase<Derived,dim,o>& observer() const
        {
            static_assert(o>=0 && o<=dim,"order must be >=0 and <= dim");
            return *this;
        }
        
        /**********************************************************************/
        void stats() const
        {
            SimplexObserverRecursion<Derived,dim,dim>::stats();
            
            
            
        }
    
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template<typename Derived, short int dim,short int order>
    struct SimplexObserverRecursion : public SimplexObserverBase<Derived,dim,order>,
    /*                             */ public SimplexObserverRecursion<Derived,dim,order-1>
    {
    
        /**********************************************************************/
        void stats() const
        {
            SimplexObserverRecursion<Derived,dim,order-1>::stats();
            
            size_t nT(0);
            size_t nB(0);
            double volT=0.0;
            double volB=0.0;
            
            for (auto& pSimplex : SimplexObserverBase<Derived,dim,order>::simplices())
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
        
//        /**********************************************************************/
//        std::array<std::shared_ptr<Simplex<dim,order-1> >,SimplexTraits<dim,order>::nFaces> faces(const SimplexIDType& xID) const
//        {
//            std::array<std::shared_ptr<Simplex<dim,order-1> >,SimplexTraits<dim,order>::nFaces> temp;
//            for (size_t j=0;j<SimplexTraits<dim,order>::nFaces;++j)
//            {
//                temp[j]=SimplexObserver<dim,order-1>::pSimplex(SimplexTraits<dim,order>::faceID(xID,j));
//            }
//            return temp;
//        }
        
    };

    template<typename Derived, short int dim>
    struct SimplexObserverRecursion<Derived,dim,0> : public SimplexObserverBase<Derived,dim,0>
    {
        
        /**********************************************************************/
        void stats() const
        {
            
            size_t nT(0);
            size_t nB(0);
            
            for (auto& pSimplex : SimplexObserverBase<Derived,dim,0>::simplices())
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
	template<typename Derived, short int dim,short int order>
    class SimplexObserverBase : public std::map<typename SimplexTraits<dim,order>::SimplexIDType, // key
    /*                                        */ const Simplex<dim,order>* const, // value
    /*                                        */ CompareVectorsByComponent<typename SimplexTraits<dim,order>::ScalarIDType,
    /*                                        */ SimplexTraits<dim,order>::nVertices> // key compare
    /*                                        */ >
    {
        
//        const Derived* derived() const
//        {
//            return static_cast<const Derived* const >(this);
//        }

        Derived* derived()
        {
            return static_cast< Derived* const >(this);
        }


        
    public:
        
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
		
        
        ~SimplexObserverBase()
        {
            std::cout<<"Destroying SimplexObserverBase<"<<dim<<","<<order<<">"<<std::endl;
        }
		
        /**********************************************************************/
        const SimplexMapType& simplices() const
        {
            return *this;
        }
        
//		/**********************************************************************/
//        size_t size() const
//        {
//			return this->size();
//		}
		
//        /**********************************************************************/
//        typename SimplexMapType::const_iterator simplexBegin() const
//        {
//			return this->begin();
//		}
//		
//		/**********************************************************************/
//        typename SimplexMapType::const_iterator simplexEnd() const
//        {
//			return this->end();
//		}
        
//		/**********************************************************************/
//        SharedPtrType pSimplex(const SimplexIDType& xID) const
//        {
//			typename SimplexMapType::const_iterator iter(this->find(xID));
//			return (iter!=this->end())? (*(iter->second->parentBegin()))->child(xID) : SharedPtrType(new SimplexType(xID));
//		}

        /**********************************************************************/
        SharedPtrType pSimplex(const SimplexIDType& xID)
        {
            
//            std::cout<<"pSimplex function"<<std::endl;
//            SharedPtrType b(new SimplexType(derived(),xID));
//                            std::cout<<"b created"<<std::endl;
//            std::cout<<xID<<std::endl;
//            std::cout<<this<<std::endl;

            typename SimplexMapType::const_iterator iter(this->find(xID));

//            std::cout<<"here"<<std::endl;
//
//            std::cout<<"iter found"<<(iter!=this->end())<<std::endl;

            
//            return (iter!=this->end())? (*(iter->second->parentBegin()))->child(xID) : SharedPtrType(new SimplexType(const_cast<Derived* const>(derived()),xID));
            return (iter!=this->end())? (*(iter->second->parentBegin()))->child(xID) : SharedPtrType(new SimplexType(derived(),xID));

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
        std::array<std::shared_ptr<Simplex<dim,order-1> >,SimplexTraits<dim,order>::nFaces> faces(const SimplexIDType& xID)
        {
//            std::cout<<"faces start"<<std::endl;
//            std::cout<<"derived="<<derived()<<std::endl;
            std::array<std::shared_ptr<Simplex<dim,order-1> >,SimplexTraits<dim,order>::nFaces> temp;
            for (size_t j=0;j<SimplexTraits<dim,order>::nFaces;++j)
            {
//                std::cout<<j<< "of" << SimplexTraits<dim,order>::nFaces<<std::endl;
//                temp[j]=SimplexObserver<dim,order-1>::pSimplex(SimplexTraits<dim,order>::faceID(xID,j));
                temp[j]=derived()->template observer<order-1>().pSimplex(SimplexTraits<dim,order>::faceID(xID,j));
//                std::cout<<"done"<<std::endl;

            }
//            std::cout<<"faces end"<<std::endl;

            return temp;
        }
        
        /**********************************************************************/
        void insertSimplex(const Simplex<dim,order>& s)
        {            
            const bool couldInsert(this->insert(std::make_pair(s.xID,&s)).second);
            assert(couldInsert && "COULD NOT INSERT SIMPLEX IN SIPLEX-MAP");
        }
        
        /**********************************************************************/
        void removeSimplex(const Simplex<dim,order>& s)
        {            
            const int nRemoved(this->erase(s.xID));
            assert(nRemoved==1 && "COULD NOT REMOVE SIMPLEX FROM SIPLEX-MAP");
        }
        
//        /**********************************************************************/
//        typename SimplexMapType::const_iterator find(const SimplexIDType& key) const
//        {
//            return this->find(key);
//        }
        
        
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
