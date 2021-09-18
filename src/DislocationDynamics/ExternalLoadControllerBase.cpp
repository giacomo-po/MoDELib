/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *                       Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExternalLoadControllerBase_cpp_
#define model_ExternalLoadControllerBase_cpp_


#include <ExternalLoadControllerBase.h>


namespace model
{
    
 
        /**************************************************************************/
        template <int dim>
        ExternalLoadControllerBase<dim>::ExternalLoadControllerBase(const std::string& _inputFileName) :
        /* init list */ inputFileName(_inputFileName)
        {
            std::cout<<greenBoldColor<<"Reading ExternalLoadController file: "<<inputFileName<<defaultColor<<std::endl;
        }
        
        /**************************************************************************/
        template <int dim>
        ExternalLoadControllerBase<dim>::~ExternalLoadControllerBase()
        {
        }

        /**************************************************************************/
        template <int dim>
        typename ExternalLoadControllerBase<dim>::MatrixDim ExternalLoadControllerBase<dim>::stress(const VectorDim&) const = 0;
        
        /*************************************************************************/
        template <int dim>
        void ExternalLoadControllerBase<dim>::update(const long int& runID) = 0;
        
        /**************************************************************************/
        template <int dim>
        void ExternalLoadControllerBase<dim>::output(const long int& runID,
                            UniqueOutputFile<'F'>& f_file,
                            std::ofstream& F_labels) const = 0;
        
        template class ExternalLoadControllerBase<3>;
}
#endif
