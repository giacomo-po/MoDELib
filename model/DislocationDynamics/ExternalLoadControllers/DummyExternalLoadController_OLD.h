/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2017 by Giacomo Po <gpo@ucla.edu>.
 *                       Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DummyExternalLoadController_H_
#define model_DummyExternalLoadController_H_

#include <iostream>
#include <sstream>      // std::stringstream
#include <cmath>
#include <cfloat>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/IO/EigenDataReader.h>
#include <model/IO/IDreader.h>


namespace model
{
    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
    class DummyExternalLoadController;
    
    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
    using ExternalLoadController=DummyExternalLoadController<dim>;

    /**************************************************************************/
    /**************************************************************************/
    template <int dim>
    class DummyExternalLoadController
    {
        
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;


    public:
        /**************************************************************************/
        DummyExternalLoadController()
        {
        }
        
        /**************************************************************************/
        template <typename DislocationNetworkType>
        void init(const DislocationNetworkType& )
        {
        }
        
        /**************************************************************************/
        static MatrixDim externalStress(const VectorDim&)
        {
            return MatrixDim::Zero();
        }
        
        /*************************************************************************/
        template <typename DislocationNetworkType>
        void update(const DislocationNetworkType& )
        {
        }
        
        /**************************************************************************/
        std::string output() const
        {
            return std::string();
        }
        
    };
}
#endif
