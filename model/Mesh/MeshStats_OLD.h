/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshStats_H_
#define model_MeshStats_H_

#include <cmath>
#include <cfloat>
#include <iostream>
#include <model/Mesh/SimplexObserver.h>
#include <model/MPI/MPIcout.h> // defines mode::cout

namespace model
{
    
    
    /**************************************************************************/
	/**************************************************************************/
	template<int dim,int k>
	struct MeshStats
    {
        
        typedef SimplexObserver<dim,k> SimplexObserverType;
        typedef typename SimplexObserverType::SimplexMapType SimplexMapType;
        
        /**********************************************************************/
        static void stats()
        {
            MeshStats<dim,k-1>::stats();
            
            size_t nT(0);
            size_t nB(0);
            double volT=0.0;
            double volB=0.0;
            
            for (auto& pSimplex : SimplexObserverType::simplices())
            {
                nT++;
                volT+=pSimplex.second->vol0;
                
                if(pSimplex.second->isBoundarySimplex())
                {
                    nB++;
                    volB+=pSimplex.second->vol0;
                    assert(fabs(pSimplex.second->outNormal().norm()-1.0)<FLT_EPSILON);
                }
                else
                {
                    assert(pSimplex.second->outNormal().norm()<FLT_EPSILON);
                }
            }
            
            model::cout<<"    Simplex<"<<dim<<","<<k  <<"> #="<<nT<<", vol="<<volT;
            model::cout<<"     (boundary #="<<nB<<", vol="<<volB<<")\n";

        }

	};
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    struct MeshStats<dim,dim>
    {
        enum{k=dim};
        typedef SimplexObserver<dim,k> SimplexObserverType;
        typedef typename SimplexObserverType::SimplexMapType SimplexMapType;
        
        /**********************************************************************/
        static void stats()
        {
            MeshStats<dim,k-1>::stats();
            
            size_t nT(0);
            size_t nB(0);
            double volT=0.0;
            double volB=0.0;
            
            for (auto& pSimplex : SimplexObserverType::simplices())
            {
                nT++;
                volT+=pSimplex.second->vol0;
                
                if(pSimplex.second->isBoundarySimplex())
                {
                    nB++;
                    volB+=pSimplex.second->vol0;
                }
            }
            
            model::cout<<"    Simplex<"<<dim<<","<<k  <<"> #="<<nT<<", vol="<<volT;
            model::cout<<"     (boundary #="<<nB<<", vol="<<volB<<")\n";
            
        }
        
    };
    
    /**************************************************************************/
	/**************************************************************************/
	template<int dim>
	struct MeshStats<dim,0>
    {
        enum{k=0};
        typedef SimplexObserver<dim,k> SimplexObserverType;
        typedef typename SimplexObserverType::SimplexMapType SimplexMapType;
        
        /**********************************************************************/
        static void stats()
        {

            size_t nT(0);
            size_t nB(0);
            
            for (auto& pSimplex : SimplexObserverType::simplices())
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
            
            model::cout<<"    Simplex<"<<dim<<","<<k  <<"> #="<<nT;
            model::cout<<"     (boundary #="<<nB<<")\n";
        }
	};
    
    
}	// close namespace
#endif
