/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MeshStats_H_
#define model_MeshStats_H_

#include <iostream>
#include <model/Mesh/SimplexObserver.h>
#include <model/MPI/MPIcout.h> // defines mode::cout

namespace model {
    
    
    /**************************************************************************/
	/**************************************************************************/
	template<int dim,int k>
	struct MeshStats
    {
        
        typedef SimplexObserver<dim,k> SimplexObserverType;
        typedef typename SimplexObserverType::SimplexMapType SimplexMapType;
        
        /**********************************************************************/
        static void stats(const bool& countBoundarySimplices)
        {
            MeshStats<dim,k-1>::stats(countBoundarySimplices);
            
            model::cout<<"    Simplex<"<<dim<<","<<k  <<">: "<<SimplexObserverType::size();
            if (countBoundarySimplices)
            {
                size_t nB(0);
                for (typename SimplexMapType::const_iterator sIter=SimplexObserverType::simplexBegin();sIter!=SimplexObserverType::simplexEnd();++sIter)
                {
                    nB+=sIter->second->isBoundarySimplex();
                }
                model::cout<<" ("<<nB<<" boundary)";
            }
            model::cout<<"\n";
            
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
        static void stats(const bool& countBoundarySimplices)
        {
            model::cout<<"    Simplex<"<<dim<<","<<k  <<">: "<<SimplexObserverType::size();
            if (countBoundarySimplices)
            {
                size_t nB(0);
                for (typename SimplexMapType::const_iterator sIter=SimplexObserverType::simplexBegin();sIter!=SimplexObserverType::simplexEnd();++sIter)
                {
                    nB+=sIter->second->isBoundarySimplex();
                }
                                model::cout<<" ("<<nB<<" boundary)";
            }
            model::cout<<"\n";
        }
        
	};
    
    
}	// close namespace
#endif