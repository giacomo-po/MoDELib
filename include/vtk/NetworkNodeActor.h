/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkNodeActor_h_
#define model_NetworkNodeActor_h_


#include <deque>
#include <string>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>
#include <QLineEdit>
#include <QSlider>


#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPolyLine.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLabeledDataMapper.h>
#include <vtkFloatArray.h>
#include <IDreader.h>
#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
//#include <DDconfigVtkBase.h>

namespace model
{
    struct NetworkNodeActor : public QWidget
//: public DDconfigVtkBase
    {
        
        Q_OBJECT
        private slots:
            void modify();

        private:
        vtkGenericOpenGLRenderWindow* const renderWindow;

        QGridLayout* mainLayout;
        QCheckBox* showNodes;
        QSlider* sliderNodeRadius;
        QCheckBox* showNodeLabels;
        QCheckBox* showSpecificNodeLabel;
        QLineEdit* showSpecificNodeLabelEdit;
        QCheckBox* showVelocities;
        QLineEdit* velocityScaleEdit;

        
        public:
                
        vtkSmartPointer<vtkPolyData> nodePolyData;
        vtkSmartPointer<vtkGlyph3D> nodeGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> nodeMapper;
        vtkSmartPointer<vtkActor> nodeActor;
        
        vtkSmartPointer<vtkPolyData> labelPolyData;
        vtkSmartPointer<vtkLabeledDataMapper> labelMapper;
        vtkSmartPointer<vtkActor2D> labelActor;

        vtkSmartPointer<vtkPolyData> velocityPolyData;
        vtkSmartPointer<vtkGlyph3D> velocityGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> velocityMapper;
        vtkSmartPointer<vtkActor> velocityActor;
        
        vtkSmartPointer<vtkPolyData> specificNodeLabelPolyData;
        vtkSmartPointer<vtkLabeledDataMapper> specificNodeLabelMapper;
        vtkSmartPointer<vtkActor2D> specificNodeLabelActor;
        
//        size_t singleNodeID;
        unsigned char nodeClr[4][3];
        
        NetworkNodeActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const);
        void updateConfiguration(const DDconfigIO<3>& configIO);
        
        
    };
    
} // namespace model
#endif
