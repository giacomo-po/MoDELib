/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *                       Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ExternalLoadControllerBase_H_
#define model_ExternalLoadControllerBase_H_

#include <iostream>
#include <fstream>      // std::stringstream
#include <Eigen/Dense>
#include <TerminalColors.h>
//#include <UniqueOutputFile.h>
//#include <cmath>
//#include <cfloat>
//#include <Material.h>
//#include <EigenDataReader.h>
//#include <IDreader.h>


namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
    class ExternalLoadControllerBase
    {
        
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        
    public:
        
        const std::string inputFileName;

        
        /**************************************************************************/
        ExternalLoadControllerBase(const std::string& _inputFileName) ;
        
        /**************************************************************************/
        virtual ~ExternalLoadControllerBase();
        /**************************************************************************/
        virtual MatrixDim stress(const VectorDim&) const = 0;
        virtual MatrixDim strain(const VectorDim&) const = 0;

        
        /*************************************************************************/
        virtual void update(const long int& runID) = 0;
        
        /**************************************************************************/
        virtual void output(const long int& runID,
                            std::ofstream& f_file,
                            std::ofstream& F_labels) const = 0;
        
    };
}
#endif
