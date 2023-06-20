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
        /* init list */,ExternalStress(MatrixDim::Zero())
        /* init list */,ExternalStressRate(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStressRate",dim,dim,true))
        /* init list */,ExternalStrain(MatrixDim::Zero())
        /* init list */,ExternalStrainRate(TextFileParser(this->inputFileName).readMatrix<double>("ExternalStrainRate",dim,dim,true))
        /* init list */,plasticStrain(MatrixDim::Zero())
        /* init list */,MachineStiffnessRatio(TextFileParser(this->inputFileName).readMatrix<double>("MachineStiffnessRatio",1,voigtSize,true))
        /* init list */,lambda(1.0)
        /* init list */,nu_use(0.12)
        /* init list */,stressmultimachinestiffness(Eigen::Matrix<double,voigtSize,voigtSize>::Zero())
        /* init list */,strainmultimachinestiffness(Eigen::Matrix<double,voigtSize,voigtSize>::Zero())
        {
            std::cout<<greenBoldColor<<"Reading ExternalLoadController file: "<<inputFileName<<defaultColor<<std::endl;
            if((ExternalStressRate-ExternalStressRate.transpose()).norm()>DBL_EPSILON)
            {
                throw std::runtime_error("ExternalStressRate is not symmetric.");
            }
            if((ExternalStrainRate-ExternalStrainRate.transpose()).norm()>DBL_EPSILON)
            {
                throw std::runtime_error("ExternalStressRate is not symmetric.");
            }
//            assert((ExternalStressRate-ExternalStressRate.transpose()).norm()<DBL_EPSILON && "ExternalStressRate is not symmetric.");
//            assert((ExternalStrainRate-ExternalStrainRate.transpose()).norm()<DBL_EPSILON && "ExternalStrainRate is not symmetric.");
        }
        
        /**************************************************************************/
        template <int dim>
        ExternalLoadControllerBase<dim>::~ExternalLoadControllerBase()
        {
        }

//        /**************************************************************************/
//        template <int dim>
//        typename ExternalLoadControllerBase<dim>::MatrixDim ExternalLoadControllerBase<dim>::stress(const VectorDim&) const = 0;
//        
//        /*************************************************************************/
//        template <int dim>
//        void ExternalLoadControllerBase<dim>::update(const long int& runID) = 0;
//        
//        /**************************************************************************/
//        template <int dim>
//        void ExternalLoadControllerBase<dim>::output(const long int& runID,
//                            UniqueOutputFile<'F'>& f_file,
//                            std::ofstream& F_labels) const = 0;
        
        template class ExternalLoadControllerBase<3>;
}
#endif
