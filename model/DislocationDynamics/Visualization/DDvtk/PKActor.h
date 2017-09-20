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

#include <model/IO/IDreader.h>


// VTK documentation
// http://vtk.1045678.n5.nabble.com/VTK-slow-to-display-300-vtkArrowSource-in-real-time-td5740730.html
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/Visualization/ScaleGlyphs
// http://hiraku0n.blogspot.com/2012/10/writing-scalar-and-vector-field-in-vtk.html
// http://vtk.1045678.n5.nabble.com/Glyphing-vtkImageData-scalars-3D-as-arrows-td3199837.html

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    struct PKActor : public IDreader<'P',3,6,double>
    {
        static bool showPK;
        static float pkFactor;
        static constexpr int dim=3;
        typedef Eigen::Matrix<float,dim,1>  VectorDim;
        typedef Eigen::Matrix<float,dim,dim,1>  MatrixDim;
        typedef IDreader<'P',3,6,double> ReaderType;
        
//        const PKContainerType& pkContainer;
        
        vtkRenderer* const renderer;

        
        vtkSmartPointer<vtkPoints> points;
        vtkSmartPointer<vtkDoubleArray> vectors;
        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkArrowSource> arrowSource;
        vtkSmartPointer<vtkGlyph3D> glyphs;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;
        
        /**********************************************************************/
        ReaderType& reader()
        {
            return *this;
        }
        
        /**********************************************************************/
        void read(const size_t& frameID)
        {
            // Read data
            if(reader().isGood(frameID,false))
            {
                reader().read(frameID,false);
            }
            else
            {
                reader().read(frameID,true);
            }
            
            
            // insert each datapoint
            for(const auto& pk : reader())
            {
                Eigen::Map<const Eigen::Matrix<double,1,6>> val(pk.second.data());
                const Eigen::Matrix<float,3,1> P(val.segment<3>(0).template cast<float>());
                const Eigen::Matrix<float,3,1> F(val.segment<3>(3).template cast<float>());
                points->InsertNextPoint(P.data());  // origin of arrow
                vectors->InsertNextTuple(F.data()); // arrow vactor
            }
            
            // Fill polyData
            polyData->SetPoints(points);
            polyData->GetPointData()->SetVectors(vectors);
            polyData->Modified();

        }
        
        /**********************************************************************/
//        PKActor(const PKContainerType& pkContainer_in) :
        PKActor(const size_t& frameID,vtkRenderer* const ren) :
//        /* init */ pkContainer(pkContainer_in),
        /* init */ renderer(ren),
        /* init */ points(vtkSmartPointer<vtkPoints>::New()),
        /* init */ vectors(vtkSmartPointer<vtkDoubleArray>::New()),
        /* init */ polyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ arrowSource(vtkSmartPointer<vtkArrowSource>::New()),
        /* init */ glyphs(vtkSmartPointer<vtkGlyph3D>::New()),
        /* init */ mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ actor(vtkSmartPointer<vtkActor>::New())
        {
            
            // Prepare vectors
            vectors->SetNumberOfComponents(3);
            vectors->SetName("pkForce");

            // Read data
            read(frameID);
            
            // Set up glyph
            glyphs->SetSourceConnection(arrowSource->GetOutputPort());
            glyphs->SetInputData(polyData);
            glyphs->ScalingOn();
            glyphs->SetScaleModeToScaleByVector();
//            glyphs->SetScaleFactor(pkFactor);
            glyphs->OrientOn();
            glyphs->ClampingOff();
            glyphs->SetVectorModeToUseVector();
            glyphs->SetIndexModeToOff();
            
            // Set up mapper
            mapper->SetInputConnection(glyphs->GetOutputPort());
            mapper->ScalarVisibilityOff();
            
            // Set up actor
            actor->SetMapper(mapper);
            
            // Add actor to renderer
            renderer->AddActor(actor);


            modify();
        }
        
        /**********************************************************************/
        ~PKActor()
        {
            renderer->RemoveActor(actor);

        }
        
        
        /**********************************************************************/
        void modify()
        {
            //arrowSource->SetShaftRadius(1.0);
            //arrowSource->SetTipLength(1.0);
            //arrowSource->Update();
            glyphs->SetScaleFactor(pkFactor);
            
            if(showPK)
            {
                actor->VisibilityOn();
                
            }
            else
            {
                actor->VisibilityOff();
                
            }

        }
        
    };
    
    bool  PKActor::showPK=false;
    float PKActor::pkFactor=1000.0;
    
} // namespace model
#endif







