/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexCompare_H_
#define model_SimplexCompare_H_

#include <SimplexTraits.h>
#include <CompareVectorsByComponent.h>


namespace model {
    
    
    /**************************************************************************/
	/**************************************************************************/
	template<short int dim,short int order>
	struct SimplexCompare
    {
        typedef typename SimplexTraits<dim,order>::ScalarIDType Scalar;
        enum {N=SimplexTraits<dim,order>::nVertices};
        typedef CompareVectorsByComponent<Scalar,N> CompareType;
        
        /**********************************************************************/

        bool operator() (const Simplex<dim,order>* const lhs,
                         const Simplex<dim,order>* const rhs) const
        {
            return CompareType().operator()(lhs->xID,rhs->xID);
        }
        
        /**********************************************************************/
        bool operator() (const Simplex<dim,order>& lhs,
                         const Simplex<dim,order>& rhs) const
        {
            return CompareType().operator()(lhs.xID,rhs.xID);
        }
        
	};
    
    
}	// close namespace
#endif
