/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlaneBase_cpp_
#define model_LatticePlaneBase_cpp_

#include <map>
#include <LatticeModule.h>
#include <LatticePlaneBase.h>

namespace model
{
    
    LatticePlaneBase::LatticePlaneBase(const LatticeVectorType& v1_in, const LatticeVectorType& v2_in) :
    /* init base */ ReciprocalLatticeDirectionType(v1_in,v2_in),
    /* init */ primitiveVectors(std::make_pair(v1_in,v2_in))
    {
        assert(primitiveVectors.first.squaredNorm()>0);
        assert(primitiveVectors.second.squaredNorm()>0);
        assert(primitiveVectors.first.cross(primitiveVectors.second).squaredNorm()>0);
    }
        
    /**********************************************************************/
    typename LatticePlaneBase::LatticeVectorType LatticePlaneBase::snapToLattice(const Eigen::Matrix<double, 3, 1> &P) const
    {
        const Eigen::Matrix<double, 3, 2> B((Eigen::Matrix<double, 2, 3>() << primitiveVectors.first.cartesian().transpose(), primitiveVectors.second.cartesian().transpose()).finished().transpose());
        const Eigen::Matrix<double, 2, 1> nd = (B.transpose() * B).inverse() * B.transpose() * P;
        //            const Eigen::Matrix<double,3,1>  p=B*RoundEigen<double,2>::round(nd);
        const Eigen::Matrix<double, 3, 1> p = B * nd.array().round().matrix();
        //            return LatticeVectorType(p,primitiveVectors.first.covBasis,primitiveVectors.second.contraBasis);
        return LatticeVectorType(p, primitiveVectors.first.lattice);
    }

    /**********************************************************************/
    typename LatticePlaneBase::VectorDimD LatticePlaneBase::snapToPlane(const Eigen::Matrix<double, 3, 1> &P) const
    {
        //            return LatticeVectorType(p,primitiveVectors.first.covBasis,primitiveVectors.second.contraBasis);
        VectorDimD cn(this->cartesian());
        const double cnNorm(cn.norm());
        assert(cnNorm > FLT_EPSILON);
        cn /= cnNorm;
        return P - P.dot(cn) * cn;
    }

} // end namespace
#endif
