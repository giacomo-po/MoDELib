/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Hermitian_H_
#define model_Hermitian_H_


namespace model {

//    /**************************************************************************/
//	/**************************************************************************/
//	template<int dim,int order>
//	struct HermitianElement
//    {
//    };
	
    /**************************************************************************/
	/**************************************************************************/
	template<int dim,int order>
	struct Hermitian
    {
        
        enum{nNodes=dim+1};
        enum{nDofPerComponent=1+dim}; // scalar + derivatives
        
    public:

        

    
    };
    
    
}	// close namespace
#endif