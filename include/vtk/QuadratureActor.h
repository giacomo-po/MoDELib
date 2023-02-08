/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_QuadratureActor_h_
#define model_QuadratureActor_h_


#include <deque>
#include <string>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>
#include <QLineEdit>
#include <QComboBox>
#include <QGroupBox>
#include <QVBoxLayout>

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
#include <vtkRenderer.h>
#include <vtkDataSetMapper.h>

#include <IDreader.h>
#include <PlanarPolygon.h>
#include <DDauxIO.h>
#include <MeshPlane.h>
#include <Polycrystal.h>
#include <GlidePlaneFactory.h>


namespace model
{
    struct QuadratureActor : public QWidget
    {
        
        Q_OBJECT
        private slots:
            void modify();

        private:
        vtkGenericOpenGLRenderWindow* const renderWindow;
        vtkRenderer* const renderer;


        QGridLayout* mainLayout;

        QGroupBox* pkBox;
        QLineEdit* pkScaleEdit;

        QGroupBox* sfBox;
        QLineEdit* sfScaleEdit;

        QGroupBox* ltBox;
        QLineEdit* ltScaleEdit;

        
        QGroupBox* velocityBox;
        QLineEdit* velocityScaleEdit;

        
        vtkSmartPointer<vtkPolyData> pkPolydata;
        vtkSmartPointer<vtkGlyph3D> pkGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> pkMapper;
        vtkSmartPointer<vtkActor> pkActor;

        vtkSmartPointer<vtkPolyData> sfPolydata;
        vtkSmartPointer<vtkGlyph3D> sfGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> sfMapper;
        vtkSmartPointer<vtkActor> sfActor;

        vtkSmartPointer<vtkPolyData> ltPolydata;
        vtkSmartPointer<vtkGlyph3D> ltGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> ltMapper;
        vtkSmartPointer<vtkActor> ltActor;

        
        vtkSmartPointer<vtkPolyData> velocityPolydata;
        vtkSmartPointer<vtkGlyph3D> velocityGlyphs;
        vtkSmartPointer<vtkPolyDataMapper> velocityMapper;
        vtkSmartPointer<vtkActor> velocityActor;

        
        public:
                        
        QuadratureActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const,const Polycrystal<3>& poly,const DDtraitsIO& traitsIO);
        void updateConfiguration(const DDauxIO<3>& auxIO);
        
        
    };
    
} // namespace model
#endif
