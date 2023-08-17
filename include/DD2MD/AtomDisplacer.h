/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_AtomDisplacerH_
#define model_AtomDisplacer_H_

#include <string>


#include <Eigen/Dense>

#include <DislocationDynamicsBase.h>
#include <MicrostructureGenerator.h>
#include <DDconfigFields.h>

namespace model
{

    

    struct AtomDisplacer
    {
        
        typedef Eigen::Matrix<double,3,1> VectorDim;
    
        AtomDisplacer(const std::string& folderName);
        
        DislocationDynamicsBase<3> ddBase;
//        DDconfigIO<3> configIO;
        MicrostructureGenerator microstructureGenerator;
        DDconfigFields<3> configFields;
        
        const DDconfigIO<3>& config() const;
        DDconfigIO<3>& config();
        void readMicrostructure();
        void readConfiguration(const size_t& runID);
        void writeConfiguration(const size_t& runID);
        double solidAngle(const VectorDim& x) const;
        double solidAngle(const double& x,const double& y,const double& z) const;
        VectorDim dislocationPlasticDisplacement(const VectorDim& x) const;
//        pybind11::array_t<double, pybind11::array::c_style> dislocationPlasticDisplacement(const double& x,const double& y,const double& z) const;
        VectorDim dislocationPlasticDisplacement(const double& x,const double& y,const double& z) const;


    };


}
#endif
