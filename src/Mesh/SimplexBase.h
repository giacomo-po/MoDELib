/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexBase_H_
#define model_SimplexBase_H_

#include <StaticID.h>
#include <NonCopyable.h>
#include <SimplexTraits.h>
#include <SimplexObserver.h>

namespace model
{
    
	/**************************************************************************/
	/**************************************************************************/	
	template<short int _dim, short int order>
	struct SimplexBase :
    /* inheritance    */ public NonCopyable,
    /* inheritance    */ public StaticID<SimplexBase<_dim,order> >
    {
     
        static constexpr short int dim=_dim;
        typedef SimplexTraits<dim,order> SimplexTraitsType;
        typedef typename SimplexTraitsType::SimplexIDType SimplexIDType;
        
        SimplicialMesh<dim>* const mesh;
        const SimplexIDType xID;
        
		/**********************************************************************/
        SimplexBase(SimplicialMesh<dim>* const m,
                    const SimplexIDType& vIN) :
        /* init */ mesh(m)
        /* init */,xID(SimplexTraitsType::sortID(vIN))
        {/*!
          */

        }
        
        /**********************************************************************/
        SimplexObserverBase<dim,order>& observer()
        {
            return mesh->template observer<order>();
        }
        
		/**********************************************************************/
        std::array<const Simplex<dim,0>*, SimplexTraits<dim,order>::nVertices> vertices() const
        {
            std::array<const Simplex<dim,0>*, SimplexTraits<dim,order>::nVertices> temp;
            for (int v=0;v<SimplexTraits<dim,order>::nVertices;++v)
            {
                const typename SimplexTraits<dim,0>::SimplexIDType vID(std::set<size_t>{xID[v]});
                temp[v]=&mesh->template observer<0>().simplex(vID);
            }
            return temp;
        }
        
	};
    
}	// close namespace
#endif
