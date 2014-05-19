/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexBase_H_
#define model_SimplexBase_H_

#include <model/Mesh/SimplexTraits.h>
#include <model/Mesh/SimplexObserver.h>

namespace model
{
	
	
	/**************************************************************************/
	/**************************************************************************/	
	template<short int dim, short int order>
	struct SimplexBase 
    {
        
        typedef SimplexTraits<dim,order> SimplexTraitsType;
        typedef typename SimplexTraitsType::SimplexIDType SimplexIDType;
        
        const SimplexIDType xID;
        
		/**********************************************************************/
        SimplexBase(const SimplexIDType& vIN) : xID(SimplexTraitsType::sortID(vIN))
        {/*! 
          */            
        }
        
		/**********************************************************************/
        std::array<const Simplex<dim,0>*, SimplexTraits<dim,order>::nVertices> vertices() const
        {
            std::array<const Simplex<dim,0>*, SimplexTraits<dim,order>::nVertices> temp;
            for (int v=0;v<SimplexTraits<dim,order>::nVertices;++v)
            {
                typename SimplexTraits<dim,0>::SimplexIDType vID;
                vID<<xID(v);
                temp[v]=SimplexObserver<dim,0>::pSimplex(vID).get();
            }
            return temp;
        }
        
	};
    
}	// close namespace
#endif