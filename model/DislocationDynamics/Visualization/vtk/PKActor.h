/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PKActor_H_
#define model_PKActor_H_

#include <Eigen/Dense>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
//#include <vtkMath.h>
//#include <vtkProperty.h>
//#include <vtkTubeFilter.h>
//#include <vtkTransformPolyDataFilter.h>
//#include <vtkTransform.h>

//#include <model/DislocationDynamics/Visualization/vtk/DislocationSegmentActor.h>

// VTK documentation
// http://vtk.1045678.n5.nabble.com/VTK-slow-to-display-300-vtkArrowSource-in-real-time-td5740730.html
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/Visualization/ScaleGlyphs
// http://hiraku0n.blogspot.com/2012/10/writing-scalar-and-vector-field-in-vtk.html
// http://vtk.1045678.n5.nabble.com/Glyphing-vtkImageData-scalars-3D-as-arrows-td3199837.html
namespace model
{
    
    /************************************************************************/
    /************************************************************************/
    struct PKActor
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<float,dim,1>  VectorDim;
        typedef Eigen::Matrix<float,dim,dim,1>  MatrixDim;
        typedef IDreader<'P',3,6,double> PKContainerType;
        
        const PKContainerType& pkContainer;
        
        
        vtkSmartPointer<vtkPoints> points;
        vtkSmartPointer<vtkDoubleArray> vectors;
        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkArrowSource> arrowSource;
        vtkSmartPointer<vtkGlyph3D> glyphs;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;
        
        
        
        /************************************************************************/
        PKActor(const PKContainerType& pkContainer_in) :
        /* init */ pkContainer(pkContainer_in),
        /* init */ points(vtkSmartPointer<vtkPoints>::New()),
        /* init */ vectors(vtkSmartPointer<vtkDoubleArray>::New()),
        /* init */ polyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ arrowSource(vtkSmartPointer<vtkArrowSource>::New()),
        /* init */ glyphs(vtkSmartPointer<vtkGlyph3D>::New()),
        //        /* init */ transform(vtkSmartPointer<vtkTransform>::New()),
        //        /* init */ transformPD(vtkSmartPointer<vtkTransformPolyDataFilter>::New()),
        //        /* init */ arrowSource(vtkSmartPointer<vtkArrowSource>::New()),
        /* init */ mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ actor(vtkSmartPointer<vtkActor>::New())
        {
            
            vectors->SetNumberOfComponents(3);
            vectors->SetName("pkForce");
            
            
            for(const auto& pk : pkContainer)
            {
                Eigen::Map<const Eigen::Matrix<double,1,6>> val(pk.second.data());
                const Eigen::Matrix<float,3,1> P(val.segment<3>(0).template cast<float>());
                const Eigen::Matrix<float,3,1> F(val.segment<3>(3).template cast<float>());
                points->InsertNextPoint(P.data());
                vectors->InsertNextTuple(F.data());
            }
            
            polyData->SetPoints(points);
            polyData->GetPointData()->SetVectors(vectors);
            polyData->Modified();
            
            glyphs->SetSourceConnection(arrowSource->GetOutputPort());
            glyphs->SetInputData(polyData);
            glyphs->ScalingOn();
            glyphs->SetScaleModeToScaleByVector();
            glyphs->SetScaleFactor(1000);
            glyphs->OrientOn();
            glyphs->ClampingOff();
            glyphs->SetVectorModeToUseVector();
            glyphs->SetIndexModeToOff();
            
            mapper->SetInputConnection(glyphs->GetOutputPort());
            mapper->ScalarVisibilityOff();
            
            actor->SetMapper(mapper);
            
            
        }
        
        
        
        /************************************************************************/
        void modify()
        {
            //arrowSource->SetShaftRadius(1.0);
            //arrowSource->SetTipLength(1.0);
            //arrowSource->Update();
            
        }
        
        
    };
    
} // namespace model
#endif







