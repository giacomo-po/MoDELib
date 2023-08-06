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
#include <QGroupBox>
#include <QGridLayout>
#include <QWidget>
#include <QSpinBox>
#include <QLineEdit>
#include <QLabel>
#include <QPushButton>
#include <QComboBox>
#include <QCheckBox>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkPointData.h>


#include <SimplicialMesh.h>
#include <MeshPlane.h>
#include <TriangularMesh.h>
#include <DDconfigIO.h>
#include <Polycrystal.h>
//#include <DDconfigVtk.h>
#include <NetworkLinkActor.h>
#include <InclusionActor.h>
#include <NetworkLoopActor.h>

namespace model
{

    struct FieldDataPnt
    {
        const Eigen::Matrix<double,3,1> P;
        double solidAngle;
        //        Eigen::Matrix<double,3,1> displacementDD;
        Eigen::Matrix<double,3,3> stressDD;
        Eigen::Matrix<double,3,3> stressIN;

        FieldDataPnt(const Eigen::Matrix<double,3,1>& Pin);
        
        double value(const int& valID,const bool& useDD,const bool& useIN) const;
    };


    struct DDPlaneField : public QWidget
    , public TriangularMesh
    , public std::deque<FieldDataPnt>
    {
        
        Q_OBJECT
        
    public:

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
        const Polycrystal<3>& poly;
        std::shared_ptr<MeshPlane<3>> plane;
        
    public slots:
        
        void resetPlane();
        void modify();
        
    public:
        
        DDPlaneField(vtkGenericOpenGLRenderWindow* const renWin_in,vtkRenderer* const renderer_in,const Polycrystal<3>& poly_in);
        ~DDPlaneField();
        const std::deque<FieldDataPnt>& dataPnts() const;
        std::deque<FieldDataPnt>& dataPnts();
        void compute(const DDconfigIO<3>& configIO,const NetworkLoopActor& loops,const NetworkLinkActor& segments,const InclusionActor& inclusions);
        void plotField(const int& valID,const bool& useDD,const bool& useIN,const vtkSmartPointer<vtkLookupTable>& lut);

    };


    struct DDFieldWidget : public QWidget
    {
        
        Q_OBJECT
        
        QGridLayout* mainLayout;
        QLabel* boxLabel;
        QSpinBox* spinBox;
        QGroupBox* groupBox;
        QPushButton* computeButton;
        QComboBox* fieldComboBox;
        QCheckBox* dislocationsCheck;
        QCheckBox* inclusionsCheck;
        QGroupBox* customScaleBox;
        QLineEdit* minScale;
        QLineEdit* maxScale;
        vtkSmartPointer<vtkLookupTable> lut;
        vtkSmartPointer<vtkScalarBarActor> scalarBar;

        
        vtkGenericOpenGLRenderWindow* const renWin;
        vtkRenderer* const renderer;
        const Polycrystal<3>& poly;
        const DDconfigIO<3>& configIO;
        const NetworkLoopActor& loops;
        const NetworkLinkActor& segments;
        const InclusionActor& inclusions;
        
    private slots:
        void clearLayout(QLayout *layout);
        void resetPlanes();
//        void toggleAutoscale();
        void plotField();

   public slots:
        void compute();

    public:
        
        DDFieldWidget(vtkGenericOpenGLRenderWindow* const renWin_in,
                      vtkRenderer* const renderer_in,
                      const Polycrystal<3>& poly_in,
                      const DDconfigIO<3>& configIO_in,
                      const NetworkLoopActor& loops_in,
                      const NetworkLinkActor& segments_in,
                      const InclusionActor& inclusions_in);
        
    };

} // namespace model
#endif







