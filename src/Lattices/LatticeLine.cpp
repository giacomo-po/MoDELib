/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticeLine_cpp_
#define model_LatticeLine_cpp_

#include <LatticeModule.h>
#include <LatticeLine.h>
#include <math.h> /* round, floor, ceil, trunc */

namespace model
{
 
    LatticeLine::LatticeLine(const LatticeVectorType& P_in,const LatticeVectorType& d_in) :
    /* init */ P(P_in),
    /* init */ d(d_in)
    {
        assert(&P.lattice==&d.lattice && "LatticeVectors have different bases.");

    }

    /**********************************************************************/
    typename LatticeLine::LatticeVectorType LatticeLine::snapToLattice(const LatticeVectorType &P0) const
    {
        return snapToLattice(P0.cartesian());
    }

    /**********************************************************************/
    typename LatticeLine::LatticeVectorType LatticeLine::snapToLattice(const VectorDimD &P0) const
    {
        const VectorDimD dc(d.cartesian());
        const VectorDimD Pc(P.cartesian());
        const double n = (P0 - Pc).dot(dc) / dc.squaredNorm();
        return P + d * lround(n);
    }

    /**********************************************************************/
    typename LatticeLine::VectorDimD LatticeLine::snapToLine(const VectorDimD &P0) const
    {
        VectorDimD dc(d.cartesian());
        const VectorDimD Pc(P.cartesian());
        const double dNorm(dc.norm());
        assert(dNorm > FLT_EPSILON);
        dc /= dNorm;
        return Pc + (P0 - Pc).dot(dc) * dc;
    }

    /**********************************************************************/
    bool LatticeLine::contains(const LatticeVectorType &P0) const
    {
        assert(&P.lattice == &P0.lattice && "LatticeVectors have different bases.");

        //            assert(&P.covBasis==&P0.covBasis && "LatticeVectors have different bases.");
        //            assert(&P.contraBasis==&P0.contraBasis && "LatticeVectors have different bases.");
        return LatticeDirectionType(LatticeVectorType(P0 - P)).cross(d).squaredNorm() == 0;
    }

} // end namespace
#endif