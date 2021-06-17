/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CrackSystem_H_
#define model_CrackSystem_H_

#include <Eigen/Dense>

namespace model
{
    template <int dim>
    class CrackSystem
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
    public:
        
        /**********************************************************************/
        VectorDim displacement(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The displacement field in the DefectiveCrystal at P
          */

            return VectorDim::Zero();
        }
        
        /**********************************************************************/
        template<typename ElementType>
        void displacement(std::vector<FEMnodeEvaluation<ElementType,dim,1>>& fieldPoints) const
        {

        }
        
        /**********************************************************************/
        MatrixDim stress(const VectorDim& x) const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */

            return MatrixDim::Zero();
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortion() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */

            return MatrixDim::Zero();
        }
        
        /**********************************************************************/
        MatrixDim plasticDistortionRate() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            return MatrixDim::Zero();
        }
        
    };
}
#endif
