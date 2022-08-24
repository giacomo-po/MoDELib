/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplicialMeshActor_H_
#define model_SimplicialMeshActor_H_

#include <random>
#include <algorithm>

#include <QGroupBox>
#include <QGridLayout>
#include <QCheckBox>
#include <QLabel>
#include <QWidget>
#include <QSlider>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkTriangle.h>
#include <vtkCellData.h>
#include <vtkPlane.h>
#include <vtkClipPolyData.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyData.h>
#include <vtkStripper.h>
#include <vtkFeatureEdges.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkSphereSource.h>
#include <vtkPlane.h>
#include <vtkCommand.h>
#include <vtkImplicitPlaneWidget2.h>
#include <vtkImplicitPlaneRepresentation.h>
#include <vtkNamedColors.h>
#include <vtkLookupTable.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPlane.h>
#include <vtkCommand.h>
#include <vtkImplicitPlaneWidget2.h>
#include <vtkImplicitPlaneRepresentation.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkCubeAxesActor.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkDoubleArray.h>

#include <TextFileParser.h>

#include <SimplicialMesh.h>
//#include <MeshNodeIO.h>

//#include <VertexReader.h>

// See this for plane interactor
// https://www.vtk.org/Wiki/VTK/Examples/Cxx/Widgets/ImplicitPlaneWidget2

// camera
// https://www.vtk.org/pipermail/vtkusers/2014-January/082864.html

namespace model
{


//struct FieldActorCallback : public vtkCommand
//
//{
//    vtkPlane *plane;
////    TriangularMesh planeTriangulation;
//
//
//    FieldActorCallback():
//    /* init */ plane(0)
//    {
//    }
//
//    static FieldActorCallback *New()
//    {
//        return new FieldActorCallback;
//
//    }
//
//    virtual void Execute(vtkObject *caller, unsigned long, void*)
//    {
//        vtkImplicitPlaneWidget2 *planeWidget = reinterpret_cast<vtkImplicitPlaneWidget2*>(caller);
//        vtkImplicitPlaneRepresentation  *rep = reinterpret_cast<vtkImplicitPlaneRepresentation*>(planeWidget->GetRepresentation());
//        rep->GetPlane(this->plane);
//    }
//};

    struct SimplicialMeshActor : public QWidget
    {
        
        Q_OBJECT
        
    private slots:
        
        void modify();

        
    public:
        QGridLayout* mainLayout;
        QCheckBox* showMesh;
        QCheckBox* showFaceBoundaries;
//        QCheckBox* showGrainColors;
        QCheckBox* showRegionBoundaries;
        QSlider* sliderRegionBoundaries;
        QCheckBox* showClipPlane;
        QSlider* sliderClipPlane;
        QCheckBox* showAxes;
 //       QSlider* sliderAxes;

        QCheckBox* showPeriodicityVectors;

        
        vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow;
        vtkRenderer* const renderer;
        const SimplicialMesh<3>& mesh;
        
        double dispCorr;

        std::map<size_t,std::pair<size_t,Eigen::Matrix<double,3,1>>> sIDtoVtkPointsMap;
        
        vtkSmartPointer<vtkPoints> pts;
        vtkSmartPointer<vtkPolyData> polydata;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;

        vtkSmartPointer<vtkPoints> facePts;
        vtkSmartPointer<vtkPolyData> facePolydata;
        vtkSmartPointer<vtkPolyDataMapper> faceMapper;
        vtkSmartPointer<vtkActor> faceActor;

        
        vtkSmartPointer<vtkUnsignedCharArray> gbColors;
        vtkSmartPointer<vtkPoints> gbPoints;
        vtkSmartPointer<vtkCellArray> gbTriangles;
        vtkSmartPointer<vtkPolyData> gbTrianglePolyData;
        vtkSmartPointer<vtkPolyDataMapper> gbMapper;
        vtkSmartPointer<vtkActor> gbActor;
        
        vtkSmartPointer<vtkPlane> clipPlane;
        vtkSmartPointer<vtkClipPolyData> clipper;
        vtkSmartPointer<vtkPolyData> clippedPolyData;
        vtkSmartPointer<vtkDataSetMapper> clipMapper;
        vtkSmartPointer<vtkActor> clipActor;
        
        vtkSmartPointer<vtkCubeAxesActor> cubeAxesActor;
//        vtkSmartPointer<vtkAxesActor> axes;
//        vtkSmartPointer<vtkOrientationMarkerWidget> axesWidget;

        vtkSmartPointer<vtkPolyData> periodicityPolydata;
        vtkSmartPointer<vtkGlyph3D> periodicityGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> periodicityMapper;
        vtkSmartPointer<vtkActor> periodicityActor;


//        vtkSmartPointer<FieldActorCallback> myCallback;
//        vtkSmartPointer<vtkImplicitPlaneRepresentation> rep;
//        vtkSmartPointer<vtkImplicitPlaneWidget2> planeWidget;

                

        SimplicialMeshActor(vtkSmartPointer<vtkGenericOpenGLRenderWindow>,vtkRenderer* const renderer_in, const SimplicialMesh<3>& mesh_in);
        void modifyPts();
        
        
    };
    
} // namespace model
#endif







