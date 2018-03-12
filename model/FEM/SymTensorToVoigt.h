/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SymTensorToVoigt_H_
#define model_SymTensorToVoigt_H_

#include <Eigen/Dense>

namespace model
{
    template <int dim>
    struct SymTensorToVoigt
    {
        
        static int t2v(const int& i,const int& j)
        {
            if(i>=j)
            {/*!Lower diagonals have equation i=j+k. The k-th diagaonl has dim-k elements.
              * So, letting k=i-j, the return is k*dim-k*(k-1)/2+j
              */
                assert(j>=0 && i<dim);
                const int k(i-j);
                return k*dim-k*(k-1)/2+j;
            }
            else
            {
                return t2v(j,i);
            }
        }
        
    };
}
#endif




