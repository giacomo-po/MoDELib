/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneActor_H_
#define model_GlidePlaneActor_H_

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
#include <vtkTriangle.h>


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
    struct GlidePlaneActor : public IDreader<'G',2,5,double>
    {
        static bool showGlidePlane;
        static float opacity;

//        static float pkFactor;
        static constexpr int dim=3;
        typedef Eigen::Matrix<float,dim,1>  VectorDim;
//        typedef Eigen::Matrix<float,dim,dim,1>  MatrixDim;
        typedef IDreader<'G',2,5,double> ReaderType;
        
//        const PKContainerType& pkContainer;
        
        vtkRenderer* const renderer;

        

        vtkSmartPointer<vtkPoints> points;
        
        vtkSmartPointer<vtkCellArray> triangles;
        vtkSmartPointer<vtkPolyData> trianglePolyData;

        
//        vtkSmartPointer<vtkDoubleArray> vectors;
//        vtkSmartPointer<vtkPolyData> polyData;
//        vtkSmartPointer<vtkArrowSource> arrowSource;
//        vtkSmartPointer<vtkGlyph3D> glyphs;
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
            //VectorDim C(VectorDim::Zero())
            
            std::map<size_t,VectorDim> centers;
//            std::map<size_t,size_t> centersSize;
            std::map<size_t,std::deque<VectorDim>> loopPts;
            
            for(const auto& pt : reader())
            {
                const size_t loopID=pt.first[0];
                
                if(centers.find(loopID)==centers.end())
                {
                    centers[loopID]=VectorDim::Zero();
//                    centersSize[loopID]=0;

                }
//                else
//                {
                Eigen::Matrix<float,dim,1>  P=(Eigen::Map<const Eigen::Matrix<double,5,1>>(pt.second.data())).segment<dim>(2).cast<float>();
                
                centers[loopID]+=P;
                
                loopPts[loopID].push_back(P);
                
//                    centersSize[loopID]++;
//                }
//                
//                
//                Eigen::Map<const Eigen::Matrix<double,1,6>> val(pk.second.data());
//                const Eigen::Matrix<float,3,1> P(val.segment<3>(0).template cast<float>());
//                const Eigen::Matrix<float,3,1> F(val.segment<3>(3).template cast<float>());
//                points->InsertNextPoint(P.data());  // origin of arrow
//                vectors->InsertNextTuple(F.data()); // arrow vactor
            }
            
            size_t ptID=0;
            for(const auto& pair : loopPts)
            {
                centers[pair.first]/=pair.second.size();
                
                for(size_t k=0;k<pair.second.size();++k)
                {
                    size_t next = (k<pair.second.size()-1)? k+1 : 0;
                    
//                    vtkSmartPointer<vtkPoints> points(vtkSmartPointer<vtkPoints>::New());
                    points->InsertNextPoint(pair.second[k](0),pair.second[k](1),pair.second[k](2));
                    points->InsertNextPoint(pair.second[next](0),pair.second[next](1),pair.second[next](2));
                    points->InsertNextPoint(centers[pair.first](0),centers[pair.first](1),centers[pair.first](2));

                    
                    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                    triangle->GetPointIds()->SetId ( 0, ptID+0 );
                    triangle->GetPointIds()->SetId ( 1, ptID+1 );
                    triangle->GetPointIds()->SetId ( 2, ptID+2 );
                    triangles->InsertNextCell ( triangle );
                    
                    ptID+=3;


                }
                
                
            }
            
            trianglePolyData->SetPoints ( points );
            trianglePolyData->SetPolys ( triangles );
                        trianglePolyData->Modified();

            
//            // Fill polyData
//            polyData->SetPoints(points);
//            polyData->GetPointData()->SetVectors(vectors);
//            polyData->Modified();

        }
        
        /**********************************************************************/
//        PKActor(const PKContainerType& pkContainer_in) :
        GlidePlaneActor(const size_t& frameID,vtkRenderer* const ren) :
        //        /* init */ pkContainer(pkContainer_in),
        /* init */ renderer(ren),
        /* init */ points(vtkSmartPointer<vtkPoints>::New()),
        /* init */ triangles(vtkSmartPointer<vtkCellArray>::New()),
        /* init */ trianglePolyData(vtkSmartPointer<vtkPolyData>::New()),

//        /* init */ vectors(vtkSmartPointer<vtkDoubleArray>::New()),
//        /* init */ polyData(vtkSmartPointer<vtkPolyData>::New()),
//        /* init */ arrowSource(vtkSmartPointer<vtkArrowSource>::New()),
//        /* init */ glyphs(vtkSmartPointer<vtkGlyph3D>::New()),
        /* init */ mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ actor(vtkSmartPointer<vtkActor>::New())
        {
            

            // Read data
            read(frameID);
            
            mapper->SetInputData(trianglePolyData);
            actor->SetMapper(mapper);
            modify();
            
            renderer->AddActor(actor);

        }
        
        /**********************************************************************/
        ~GlidePlaneActor()
        {
            renderer->RemoveActor(actor);
        }
        
        
        /**********************************************************************/
        void modify()
        {
            
            actor->GetProperty()->SetColor(0.5,0.5,0.5);
            actor->GetProperty()->SetOpacity(opacity);
            //arrowSource->SetShaftRadius(1.0);
            //arrowSource->SetTipLength(1.0);
            //arrowSource->Update();
//            glyphs->SetScaleFactor(pkFactor);
//            
            if(showGlidePlane)
            {
                actor->VisibilityOn();
                
            }
            else
            {
                actor->VisibilityOff();
                
            }

        }
        
    };
    
    bool  GlidePlaneActor::showGlidePlane=false;
    float  GlidePlaneActor::opacity=0.25;

    //    float PKActor::pkFactor=1000.0;
    
} // namespace model
#endif







