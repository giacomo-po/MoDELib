/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexObserver_H_
#define model_SimplexObserver_H_

#include <iostream>
#include <memory> // shared_ptr
#include <map>
#include <array>
#include <SimplexTraits.h>
#include <CompareVectorsByComponent.h>


namespace model
{

    template<int dim>
    class SimplicialMesh;

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
        
        template<short int o>
        const SimplexObserverBase<dim,o>& observer() const
        {
            static_assert(o>=0 && o<=dim,"order must be >=0 and <= dim");
            return *this;
        }
        
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
            
            std::cout<<"    Simplex<"<<dim<<","<<order  <<"> #="<<nT<<", vol="<<volT;
            std::cout<<"     (boundary #="<<nB<<", vol="<<volB<<")\n";
            
        }
        
        
    };
    
    template< short int dim>
    struct SimplexObserverRecursion<dim,0> : public SimplexObserverBase<dim,0>
    {
        
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
            
            std::cout<<"    Simplex<"<<dim<<","<<0  <<"> #="<<nT;
            std::cout<<"     (boundary #="<<nB<<")\n";
        }
        
    };

    
    
    /**************************************************************************/
    /**************************************************************************/
	template<short int dim,short int order>
    struct SimplexObserverBase : public std::map<typename SimplexTraits<dim,order>::SimplexIDType, // key
    /*                                        */ const Simplex<dim,order>* const> // value
    {
        typedef typename SimplexTraits<dim,order>::SimplexIDType SimplexIDType;
        typedef Simplex<dim,order> SimplexType;
		typedef std::map<typename SimplexTraits<dim,order>::SimplexIDType, // key
        /*            */ const Simplex<dim,order>* const> SimplexMapType;// value
        typedef typename SimplexMapType::const_iterator const_iterator;
        typedef typename SimplexMapType::iterator iterator;
        
        typedef std::shared_ptr<SimplexType> SharedPtrType;
		
		
        /**********************************************************************/
        const SimplexMapType& simplices() const
        {
            return *this;
        }
        

        
		/**********************************************************************/
        SharedPtrType sharedSimplex(SimplicialMesh<dim>* const mesh,
                                           const SimplexIDType& xID)
        {
			typename SimplexMapType::const_iterator iter(this->find(xID));
//			return (iter!=this->end())? (*(iter->second->parentBegin()))->child(xID) : SharedPtrType(new SimplexType(mesh,xID));
            return (iter!=this->end())? (iter->second->parents().begin())->second->child(xID) : SharedPtrType(new SimplexType(mesh,xID));
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
            if(!couldInsert)
            {
                throw std::runtime_error("COULD NOT INSERT SIMPLEX IN SIPLEX-MAP");
            }
        }
        
        /**********************************************************************/
        void removeSimplex(const Simplex<dim,order>& s)
        {            
            const int nRemoved(this->erase(s.xID));
            if(nRemoved!=1)
            {
                throw std::runtime_error("COULD NOT REMOVE SIMPLEX FROM SIPLEX-MAP");
            }
        }
                
	};
        
}	// close namespace
#endif

