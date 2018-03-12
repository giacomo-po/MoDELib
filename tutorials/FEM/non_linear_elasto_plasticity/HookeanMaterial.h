/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

// The following line is needed for shape functions or order 5 or more
//#define EIGEN_STACK_ALLOCATION_LIMIT 1000000

#ifndef model_HookeanMaterial_H_
#define model_HookeanMaterial_H_

#include <Eigen/Dense>

namespace model
{
    template <int dim>
    struct HookeanMaterial
    {
        
        static constexpr int StrainVoigtSize=dim*(dim+1)/2;
        typedef Eigen::Matrix<double,StrainVoigtSize,1> StrainVoigtVector;
        typedef Eigen::Matrix<double,StrainVoigtSize,StrainVoigtSize> StrainVoigtMatrix;
        
        const StrainVoigtMatrix C0;
        
        /**************************************************************************/
        HookeanMaterial(const StrainVoigtMatrix& C_in) :
        /* init */ C0(C_in)
        {
            
        }
        
        /**************************************************************************/
        double w(const StrainVoigtVector& E) const
        {
            return 0.5*E.transpose()*C0*E;
        }
        
        /**************************************************************************/
        StrainVoigtVector S(const StrainVoigtVector& E) const
        {/*!\param[in] E the column vector of strain in Voigt form
          *\returns the column vector of the derivative dw/dE
          */
            return C0*E;
        }
        
        /**************************************************************************/
        StrainVoigtMatrix& C(const StrainVoigtVector& E) const
        {/*!\param[in] E the column vector of strain in Voigt form
          *\returns the column vector of the derivative dw/dE
          */
            return C0;
        }
        
        
    };
}
#endif




