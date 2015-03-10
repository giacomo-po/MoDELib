/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SNAPTOLATTICE_H_
#define model_SNAPTOLATTICE_H_

#include <Eigen/Dense>
#include <model/Math/RoundEigen.h>
#include <model/DislocationDynamics/Materials/CrystalOrientation.h>

namespace model
{
//    template <int...N>
//    int gcd(const int& a,const N& n...)
//    {
//        return gcd(a,gcd(b,n...));
//    }
//    
//    int gcd(const int& a,const int& b)
//    {
//        return b>0? gcd(b, a % b) : a;
//    }
    
    /**********************************************************************/
    /**********************************************************************/
    template <short unsigned int dim>
    struct SnapToLattice
    {

        /**********************************************************************/
        int gcd(const size_t& a,const size_t& b)
        {
            return b>0? gcd(b, a % b) : a;
        }
        
        /**********************************************************************/
        int gcd(const size_t& a,const size_t& b,const size_t& c)
        {
            return gcd(a,gcd(b,c));
        }
        
//        /**********************************************************************/
//        VectorDim round(const VectorDim& P)
//        {
//            return CrystalOrientation<dim>::lattMat()*RoundEigen<double,dim>::round(CrystalOrientation<dim>::invLattMat()*P);
//        }
        
    };
    
} // end namespace
#endif

