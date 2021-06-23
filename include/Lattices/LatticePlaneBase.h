/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LatticePlaneBase_h_
#define model_LatticePlaneBase_h_

#include <map>
#include <LatticeVector.h>
#include <LatticeDirection.h>
#include <ReciprocalLatticeDirection.h>

namespace model
{
    struct LatticePlaneBase : public ReciprocalLatticeDirection<3>
    {
        typedef Eigen::Matrix<double,3,1> VectorDimD;
        typedef Eigen::Matrix<long int,3,1> VectorDimI;
        typedef LatticeVector<3>    LatticeVectorType;
        typedef LatticeDirection<3>    LatticeDirectionType;
        typedef ReciprocalLatticeDirection<3> ReciprocalLatticeDirectionType;
        
        const std::pair<LatticeDirectionType,LatticeDirectionType> primitiveVectors;
  
        
        /**********************************************************************/
        LatticePlaneBase(const LatticeVectorType& v1_in, const LatticeVectorType& v2_in);
        

        
        /**********************************************************************/
        LatticeVectorType snapToLattice(const Eigen::Matrix<double,3,1>& P) const;
        
        /**********************************************************************/
        VectorDimD snapToPlane(const Eigen::Matrix<double,3,1>& P) const;
        
    };
    
} // end namespace
#endif