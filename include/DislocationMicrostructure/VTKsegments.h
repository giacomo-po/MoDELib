/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_VTKsegments_H_
#define model_VTKsegments_H_

#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <cmath>
#include <list>
#include <vector>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif


#include <TextFileParser.h>
#include <GlidePlaneModule.h>
#include <DDtraitsIO.h>
#include <StressStraight.h>
#include <DislocationQuadraturePoint.h>
#include <QuadratureDynamic.h>
#include <SplineSegmentBase.h>
#include <TerminalColors.h>



namespace model
{

    class VTKsegments : public std::vector<Eigen::Matrix<double,3,1>>
    /*               */,public std::vector<StressStraight<3>>
    /*               */,public std::vector<DislocationQuadraturePoint<3,0>>
    {

        static constexpr int dim = 3; // 3D
        static constexpr int corder = 0; // linear segments
        typedef SplineSegmentBase<dim,corder> SplineSegmentBaseType;
        static constexpr int Ncoeff= SplineSegmentBaseType::Ncoeff; // # of spline coefficients in 1D
        typedef typename SplineSegmentBaseType::MatrixNcoeffDim MatrixNcoeffDim;
        typedef typename SplineSegmentBaseType::VectorDim VectorDim;
        typedef typename SplineSegmentBaseType::MatrixDim MatrixDim;
        typedef Eigen::Matrix<int,dim,1> VectorDimI;

//        typedef Eigen::Matrix<double,3,1> VectorDim;
        typedef   QuadratureDynamic<1,UniformOpen,1,2,3,4,5,6,7,8,16,32,64,128,256,512,1024> QuadratureDynamicType;

        
        void updateQuadraturePoints(const std::string& vtkFilePrefix);

//        VectorDimI getPbcFlags(const std::string& filename) const;
        
        MatrixDim cellMatrix;
        VectorDimI pbcFlags;
        std::vector<VectorDim> periodicShifts;
        
    public:
        
        VTKsegments(const std::string& folderName);
        
//        const DDtraitsIO traitsIO;
        TextFileParser perser;
        const PolycrystallineMaterialBase material;
        const double quadPerLength;
//        const VectorDimI pbcFlags;
        
        
        void writeVTK(const std::string& vtkFilePrefix) const;
        void readVTK(const std::string& vtkFilePrefix);
        
        const std::vector<VectorDim>& nodes() const;
        std::vector<VectorDim>& nodes();
        const std::vector<StressStraight<3>>& segments() const;
        std::vector<StressStraight<3>>& segments();
        const std::vector<DislocationQuadraturePoint<3,0>>& quadraturePoints() const;
        std::vector<DislocationQuadraturePoint<3,0>>& quadraturePoints();


    };

}
#endif
