/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlane_h_
#define model_LatticePlane_h_

#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticePlaneBase.h>

namespace model
{
//    template <int dim>
    class LatticePlane
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        typedef LatticeVector<dim>    LatticeVectorType;
        
       
    public:
        
        const LatticeVectorType P;
        const LatticePlaneBase n;
        
        /**********************************************************************/
        LatticePlane(const LatticeVectorType& P_in,const LatticePlaneBase& n_in) :
        /* init */ P(P_in),
        /* init */ n(n_in)
        {
            assert(&P.lattice==&n.lattice && "LatticeVectors have different bases.");
        }
        
        /**********************************************************************/
        LatticeVectorType snapToLattice(const VectorDimD& P0) const
        {
            return P+n.snapToLattice(P0-P.cartesian());
        }
        
        /**********************************************************************/
        VectorDimD snapToPlane(const VectorDimD& P0) const
        {
            return P.cartesian()+n.snapToPlane(P0-P.cartesian());
        }
        
        /**********************************************************************/
        bool contains(const LatticeVectorType& L) const
        {
            assert(&P.lattice==&L.lattice && "LatticeVectors have different bases.");
            return (L-P).dot(n)==0;
        }
        
        /**********************************************************************/
        bool contains(const Eigen::Matrix<double,3,1>& P0) const
        {
            return fabs((P0-P.cartesian()).dot(n.cartesian()))<FLT_EPSILON;
        }
        
    };
    
} // end namespace
#endif
