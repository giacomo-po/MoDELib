/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VoigtTraits_cpp_
#define model_VoigtTraits_cpp_


#include <VoigtTraits.h>

namespace model
{
    
template <int dim>
SymmetricVoigtTraits<dim>::SymmetricVoigtTraits(const VoigtSizeMatrixType& tensorIndex_in):
/* init */ tensorIndex(tensorIndex_in)
/* init */,voigtIndex(getVoigtIndex(tensorIndex))
{
    
}

template <int dim>
Eigen::Matrix<size_t,dim,dim> SymmetricVoigtTraits<dim>::getVoigtIndex(const VoigtSizeMatrixType& ti)
{
    
    Eigen::Matrix<size_t,dim,dim> temp(Eigen::Matrix<size_t,dim,dim>::Zero());
    for(int k=0;k<ti.rows();++k)
    {
        temp(ti(k,0),ti(k,1))=k;
        temp(ti(k,1),ti(k,0))=k;
    }
    return temp;
}


template <int dim>
typename SymmetricVoigtTraits<dim>::MatrixDim SymmetricVoigtTraits<dim>::v2m(const Eigen::Matrix<double,voigtSize,1>& voigtvector, const bool& is_strain) const
{
    //from voigt to matrix format
    MatrixDim temp(MatrixDim::Zero());
    for (size_t i=0;i<voigtSize;i++)
    {
        if(is_strain && tensorIndex.row(i)[0] != tensorIndex.row(i)[1])
        {
            temp.row(tensorIndex.row(i)[0])[tensorIndex.row(i)[1]]=0.5*voigtvector[i];
            temp.row(tensorIndex.row(i)[1])[tensorIndex.row(i)[0]]=0.5*voigtvector[i];
        }
        else
        {
            temp.row(tensorIndex.row(i)[0])[tensorIndex.row(i)[1]]=voigtvector[i];
            temp.row(tensorIndex.row(i)[1])[tensorIndex.row(i)[0]]=voigtvector[i];
        }
    }
    return temp;
}

template <int dim>
typename SymmetricVoigtTraits<dim>::VectorVoigt SymmetricVoigtTraits<dim>::m2v(const MatrixDim& input_matrix, const bool& is_strain) const
{
    //from matrix to voigt format
    VectorVoigt temp(VectorVoigt::Zero());
    for (size_t i=0;i<voigtSize;i++)
    {
        if(is_strain && tensorIndex.row(i)[0] != tensorIndex.row(i)[1])
        {
            temp.row(i)[0]=2.0*input_matrix.row(tensorIndex.row(i)[0])[tensorIndex.row(i)[1]];
        }
        else
        {
            temp.row(i)[0]=input_matrix.row(tensorIndex.row(i)[0])[tensorIndex.row(i)[1]];
        }
    }
    return temp;
}

//template class SymmetricVoigtTraits<2>;
template struct SymmetricVoigtTraits<3>;

}
#endif
