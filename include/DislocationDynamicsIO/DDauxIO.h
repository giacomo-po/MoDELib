/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDauxIO_H_
#define model_DDauxIO_H_

#include <vector>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <cfloat>      // std::ifstream
#include <DDbaseIO.h>
#include <PeriodicGlidePlane.h>
//#include <GlidePlaneIO.h>
#include <GlidePlaneModule.h>
#include <DislocationQuadraturePointIO.h>
#include <MeshNodeIO.h>
#include <DefectiveCrystalParameters.h>

namespace model
{
        
    template <int dim>
    struct DDauxIO : public DDbaseIO
    /*            */,private std::vector<MeshNodeIO<dim>>
//    /*            */,private std::vector<GlidePlaneIO<dim>>
    /*            */,private std::vector<PeriodicPlanePatchIO<dim>>
    /*            */,private std::vector<DislocationQuadraturePointIO<dim>>
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        DDauxIO(const std::string& folderName,const std::string& suffix="");
        void clear();
        void addPeriodicGlidePlane(const PeriodicGlidePlane<dim>& pgp);
        const std::vector<MeshNodeIO<dim>>& meshNodes() const;
        std::vector<MeshNodeIO<dim>>& meshNodes();
//        const std::vector<GlidePlaneIO<dim>>& glidePlanes() const;
//        std::vector<GlidePlaneIO<dim>>& glidePlanes();
        const std::vector<PeriodicPlanePatchIO<dim>>& periodicGlidePlanePatches() const;
        std::vector<PeriodicPlanePatchIO<dim>>& periodicGlidePlanePatches();
        const std::vector<DislocationQuadraturePointIO<dim>>& quadraturePoints() const;
        std::vector<DislocationQuadraturePointIO<dim>>& quadraturePoints();
        void write(const size_t& runID,const bool& outputBinary);
        void writeTxt(const long int& runID);
        void writeBin(const size_t& runID);
        void readBin(const size_t& runID);
        void readTxt(const size_t& runID);
        void printLog(const std::chrono::time_point<std::chrono::system_clock>& t0) const;
        void read(const size_t& runID);
        
    };
    
}
#endif

