/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDFieldWidget_H_
#define model_DDFieldWidget_H_

#include <Eigen/Dense>

#include <memory>
//#include <random>
//#include <algorithm>
//
#include <QGroupBox>
#include <QGridLayout>
//#include <QCheckBox>
//#include <QLabel>
#include <QWidget>
#include <QSpinBox>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
//
//#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
//#include <vtkProperty.h>
//#include <vtkRenderer.h>
//#include <vtkTriangle.h>
//#include <vtkCellData.h>
//#include <vtkPlane.h>
//#include <vtkClipPolyData.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
//#include <vtkActor.h>
//#include <vtkCamera.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkDataSetMapper.h>
//#include <vtkPolyData.h>
//#include <vtkStripper.h>
//#include <vtkFeatureEdges.h>
//#include <vtkActor.h>
//#include <vtkProperty.h>
//#include <vtkPolyDataMapper.h>
//#include <vtkClipPolyData.h>
//#include <vtkPlane.h>
//#include <vtkSphereSource.h>
//#include <vtkPlane.h>
//#include <vtkCommand.h>
//#include <vtkImplicitPlaneWidget2.h>
//#include <vtkImplicitPlaneRepresentation.h>
//#include <vtkNamedColors.h>
//#include <vtkLookupTable.h>
//#include <vtkUnsignedCharArray.h>
//#include <vtkPointData.h>
//#include <vtkGenericOpenGLRenderWindow.h>
//#include <vtkPlane.h>
//#include <vtkCommand.h>
//#include <vtkImplicitPlaneWidget2.h>
//#include <vtkImplicitPlaneRepresentation.h>
//#include <vtkAxesActor.h>
//#include <vtkOrientationMarkerWidget.h>
//#include <vtkCubeAxesActor.h>
//#include <vtkGlyph3D.h>
//#include <vtkArrowSource.h>
//#include <vtkDoubleArray.h>
//
//#include <TextFileParser.h>
//
#include <SimplicialMesh.h>
#include <MeshPlane.h>
#include <TriangularMesh.h>

namespace model
{

    struct FieldDataPnt
    {
        
        const Eigen::Matrix<double,3,1> P;
        double solidAngle;
//        Eigen::Matrix<double,3,1> displacementDD;

        FieldDataPnt(const Eigen::Matrix<double,3,1>& Pin);
        
    };


    struct DDPlaneField : public QWidget
                        , public TriangularMesh
                        , public std::deque<FieldDataPnt>
    {
        
        Q_OBJECT

        typedef Eigen::Matrix<double,3,1> VectorDim;

        QGridLayout* mainLayout;
        QGroupBox* groupBox;
        
        QLineEdit* posEdit;
        QLineEdit* normalEdit;
        QLineEdit* meshSizeEdit;
        
        vtkSmartPointer<vtkPolyData> meshPolydata;
        vtkSmartPointer<vtkPolyDataMapper> meshMapper;
        vtkSmartPointer<vtkActor> meshActor;

        vtkGenericOpenGLRenderWindow* const renWin;
        vtkRenderer* const renderer;
        const SimplicialMesh<3>& mesh;
        std::shared_ptr<MeshPlane<3>> plane;

        public slots:
            
            void resetPlane();
        void modify();

        public:

        DDPlaneField(vtkGenericOpenGLRenderWindow* const renWin_in,vtkRenderer* const renderer_in,const SimplicialMesh<3>& mesh_in);
        ~DDPlaneField();
        const std::deque<FieldDataPnt>& dataPnts() const;
        std::deque<FieldDataPnt>& dataPnts();
        void compute();

        
//        DDField(vtkRenderer* const renderer, const SimplicialMesh<3>& mesh_in);
    };


    struct DDFieldWidget : public QWidget
    {
        
        Q_OBJECT
        
        QGridLayout* mainLayout;
        QLabel* boxLabel;
        QSpinBox* spinBox;
        QPushButton* computeButton;
        QGroupBox* groupBox;

        vtkGenericOpenGLRenderWindow* const renWin;
        vtkRenderer* const renderer;
        const SimplicialMesh<3>& mesh;
//        std::vector<std::shared_ptr<DDPlaneField>> planes;
        
    private slots:
        
        void clearLayout(QLayout *layout);
        void resetPlanes();
        void compute();

    public:
 
        DDFieldWidget(vtkGenericOpenGLRenderWindow* const renWin_in,vtkRenderer* const renderer_in, const SimplicialMesh<3>& mesh_in);
        
        
        
    };
    
} // namespace model
#endif







