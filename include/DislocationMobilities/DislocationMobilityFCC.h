/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by David Rivera <drivera2@ucla.edu>.
 * Copyright (C) 2013 by Giacomo Po   <gpo@ucla.edu>.
 * Copyright (C) 2013 by Tamer Crosby <tcrosby@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationMobilityFCC_h_
#define _model_DislocationMobilityFCC_h_

#include <DislocationMobilityBase.h>

namespace model
{
    
    struct DislocationMobilityFCC : public DislocationMobilityBase
    {
        
        typedef Eigen::Matrix<double,3,3> MatrixDim;
        typedef Eigen::Matrix<double,3,1> VectorDim;
        
        static constexpr double kB_SI=1.38064852e-23; // [J/K]
        
        const double B0e;
        const double B1e;
        const double B0s;
        const double B1s;
        const double kB;
        
        /**********************************************************************/
        DislocationMobilityFCC(const PolycrystallineMaterialBase& material) ;
        /**********************************************************************/
        double velocity(const MatrixDim& S,
                        const VectorDim& b,
                        const VectorDim& xi,
                        const VectorDim& n,
                        const double& T,
                        const double& dL,
                        const double& dt,
                        const std::shared_ptr<StochasticForceGenerator>& sfg) override;
        
    };
    
}
#endif
