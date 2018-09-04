/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionActor_H_
#define model_InclusionActor_H_

#include <Eigen/Dense>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkLookupTable.h>


#include <model/IO/IDreader.h>


// VTK documentation
// http://vtk.1045678.n5.nabble.com/VTK-slow-to-display-300-vtkSphereSource-in-real-time-td5740730.html
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/Visualization/Scaleglyph3D
// http://hiraku0n.blogspot.com/2012/10/writing-scalar-and-vector-field-in-vtk.html
// http://vtk.1045678.n5.nabble.com/Glyphing-vtkImageData-scalars-3D-as-arrows-td3199837.html

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    struct InclusionActor : public IDreader<'E',1,14,float>
    {
        static bool showInclusions;
        //        static float pkFactor;
        static constexpr int dim=3;
        typedef Eigen::Matrix<float,dim,1>  VectorDim;
        typedef Eigen::Matrix<float,dim,dim,1>  MatrixDim;
        typedef IDreader<'E',1,14,float> ReaderType;
        
        
        vtkRenderer* const renderer;
        
        
        vtkSmartPointer<vtkPoints> points;
        vtkSmartPointer<vtkFloatArray> diameters;
        vtkSmartPointer<vtkFloatArray> colors;
        vtkSmartPointer<vtkUnstructuredGrid> grid;
        vtkSmartPointer<vtkSphereSource> sphereSource;
        vtkSmartPointer<vtkGlyph3D> glyph3D;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;
        vtkSmartPointer<vtkLookupTable> lookUpColors;
        
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
            
            if(false)
            {// Color by typeID
                std::set<int> typeIDs;
                for(const auto& pk : reader())
                {
                    Eigen::Map<const Eigen::Matrix<float,1,14>> row(pk.second.data());
                    typeIDs.insert(row(13));
                }
                
                Eigen::Matrix<float,3,1> color00((Eigen::Matrix<float,3,1>()<<38.0,38.0,38.0).finished()/255.0);
                Eigen::Matrix<float,3,1> color05((Eigen::Matrix<float,3,1>()<<153.0,115.0,0.0).finished()/255.0);
                Eigen::Matrix<float,3,1> color10((Eigen::Matrix<float,3,1>()<<77.0,0.0,0.0).finished()/255.0);
                
                lookUpColors->SetNumberOfTableValues(typeIDs.size());
                for(int f=0;f<typeIDs.size();++f)
                {
                    const float u= double(f)/(typeIDs.size()-1);
                    
                    Eigen::Matrix<float,3,1> color=color00;
                    if(u<=0.5)
                    {
                        color=color00*(1.0-u/0.5)+color05*u/0.5;
                    }
                    else
                    {
                        color=color05*(1.0-(u-0.5)/0.5)+color10*(u-0.5)/0.5;
                    }
                    lookUpColors->SetTableValue(f ,color(0) ,color(1) ,color(2) ,0.5); // red
                    
                }
                
                if(typeIDs.size()>1)
                {
                    mapper->SetScalarRange(0, typeIDs.size()-1); // to scale color label (without, col should be between 0 and 1)
                }
            }
            
            // insert each datapoint
            for(const auto& pk : reader())
            {
                Eigen::Map<const Eigen::Matrix<float,1,14>> row(pk.second.data());
                points->InsertNextPoint(row(0), row(1), row(2));
                diameters->InsertNextValue(2.0*row(3));  // origin of arrow
                colors->InsertNextValue(row(13));
                
            }
            
            grid->SetPoints(points);
            grid->GetPointData()->AddArray(diameters);
            grid->GetPointData()->SetActiveScalars("diameters"); // to set radius first
            grid->GetPointData()->AddArray(colors);
            
        }
        
        /**********************************************************************/
        InclusionActor(const size_t& frameID,vtkRenderer* const ren) :
        /* init */ renderer(ren),
        /* init */ points(vtkSmartPointer<vtkPoints>::New()),
        /* init */ diameters(vtkSmartPointer<vtkFloatArray>::New()),
        /* init */ colors(vtkSmartPointer<vtkFloatArray>::New()),
        /* init */ grid(vtkSmartPointer<vtkUnstructuredGrid>::New()),
        /* init */ sphereSource(vtkSmartPointer<vtkSphereSource>::New()),
        /* init */ glyph3D(vtkSmartPointer<vtkGlyph3D>::New()),
        /* init */ mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ actor(vtkSmartPointer<vtkActor>::New()),
        /* init */ lookUpColors(vtkSmartPointer<vtkLookupTable>::New())
        {
            
            
            diameters->SetName("diameters");
            colors->SetName("colors");
            
            
            // Read data
            read(frameID);
            
            sphereSource->SetPhiResolution(20);
            sphereSource->SetThetaResolution(20);
            
            // Set up glyph
            glyph3D->SetInputData(grid);
            glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
            
            // Set up mapper
            mapper->SetInputConnection(glyph3D->GetOutputPort());
            //            mapper->ScalarVisibilityOff();
            
            if(false)
            {
                mapper->SetScalarModeToUsePointFieldData(); // without, color are displayed regarding radius and not color label
                mapper->SelectColorArray("colors"); // !!!to set color (nevertheless you will have nothing)
                mapper->SetLookupTable(lookUpColors);
                
            }
            else
            {
                mapper->ScalarVisibilityOff();
                actor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
                actor->GetProperty()->SetOpacity(0.4); // Give some color to the mesh. (1,1,1) is white
            }
            
            // Set up actor
            actor->SetMapper(mapper);
            
            
            // Add actor to renderer
            renderer->AddActor(actor);
            
            
            modify();
        }
        
        /**********************************************************************/
        ~InclusionActor()
        {
            renderer->RemoveActor(actor);
            
        }
        
        
        /**********************************************************************/
        void modify()
        {
            //            glyph3D->SetScaleFactor(pkFactor);
            
            if(showInclusions)
            {
                actor->VisibilityOn();
                
            }
            else
            {
                actor->VisibilityOff();
                
            }
            
        }
        
    };
    
    bool  InclusionActor::showInclusions=true;
    //    float InclusionActor::pkFactor=1000.0;
    
}
#endif
