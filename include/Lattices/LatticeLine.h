/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeLine_h_
#define model_LatticeLine_h_

#include <LatticeVector.h>
#include <LatticeDirection.h>
#include <math.h>       /* round, floor, ceil, trunc */

namespace model
{
    struct LatticeLine
    {
        typedef LatticeVector<3>    LatticeVectorType;
        typedef LatticeDirection<3> LatticeDirectionType;
        typedef typename LatticeVectorType::VectorDimD VectorDimD;
        
        const LatticeVectorType P;
        const LatticeDirectionType d;
        
        LatticeLine(const LatticeVectorType& P_in,const LatticeVectorType& d_in);
        
        /**********************************************************************/
        LatticeVectorType snapToLattice(const LatticeVectorType& P0) const;
        /**********************************************************************/
        LatticeVectorType snapToLattice(const VectorDimD& P0) const;
        
        /**********************************************************************/
        VectorDimD snapToLine(const VectorDimD& P0) const;
        
        /**********************************************************************/
        bool contains(const LatticeVectorType& P0) const;
        
    };
    
} // end namespace
#endif