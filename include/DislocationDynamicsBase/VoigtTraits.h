/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VoigtTraits_H_
#define model_VoigtTraits_H_


#include <Eigen/Dense>

namespace model
{
    
    
    template <int dim>
    struct SymmetricVoigtTraits
    {
        static constexpr int voigtSize=dim*(dim+1)/2;
        typedef Eigen::Matrix<size_t,voigtSize,2> VoigtSizeMatrixType;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;

        
        const VoigtSizeMatrixType tensorIndex; // i=tensorIndex(k,0), j=tensorIndex(k,1)
        const Eigen::Matrix<size_t,dim,dim> voigtIndex; // k=voigtIndex(i,j)

        SymmetricVoigtTraits(const VoigtSizeMatrixType& voigtOrder_in);
        
        MatrixDim v2m(const Eigen::Matrix<double,voigtSize,1>& voigtvector, const bool& is_strain) const ;
        Eigen::Matrix<double,voigtSize,1> m2v(const MatrixDim& input_matrix, const bool& is_strain) const;
        
        
        static Eigen::Matrix<size_t,dim,dim> getVoigtIndex(const VoigtSizeMatrixType& ti);
        
    };
        
}
#endif
