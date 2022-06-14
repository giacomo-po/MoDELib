/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneActor_h_
#define model_GlidePlaneActor_h_


#include <deque>
#include <string>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>
#include <QLineEdit>
#include <QComboBox>

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
#include <GlidePlaneNoise.h>

//#include <DDconfigVtkBase.h>

namespace model
{
    struct GlidePlaneActor : public QWidget
//: public DDconfigVtkBase
    {
        
        Q_OBJECT
        private slots:
            void modify();

        private:
        vtkGenericOpenGLRenderWindow* const renderWindow;
        vtkRenderer* const renderer;


        QGridLayout* mainLayout;
        QCheckBox* showGlidePlanes;
        QCheckBox* showGlidePlanesNoise;
        QComboBox* glidePlanesNoiseBox;
        
//        vtkSmartPointer<vtkPoints> glidePlanePoints;
        vtkSmartPointer<vtkPolyData> glidePlanePolydata;
        vtkSmartPointer<vtkPolyDataMapper> glidePlaneMapper;
        vtkSmartPointer<vtkActor> glidePlaneActor;
        
        std::vector<vtkSmartPointer<vtkDataSetMapper>> noiseMappers;
        std::vector<vtkSmartPointer<vtkActor>> noiseActors;


        GlidePlaneFactory<3> glidePlaneFactory;
        GlidePlaneNoise planeNoise;
        
        public:
                
//        vtkSmartPointer<vtkPolyData> nodePolyData;
//        vtkSmartPointer<vtkGlyph3D> nodeGlyphs;
//        vtkSmartPointer<vtkPolyDataMapper> nodeMapper;
//        vtkSmartPointer<vtkActor> nodeActor;
//
//        vtkSmartPointer<vtkPolyData> labelPolyData;
//        vtkSmartPointer<vtkLabeledDataMapper> labelMapper;
//        vtkSmartPointer<vtkActor2D> labelActor;
//
//        vtkSmartPointer<vtkPolyData> velocityPolyData;
//        vtkSmartPointer<vtkGlyph3D> velocityGlyphs;
//        vtkSmartPointer<vtkPolyDataMapper> velocityMapper;
//        vtkSmartPointer<vtkActor> velocityActor;
//
//        vtkSmartPointer<vtkPolyData> singleNodeLabelPolyData;
//        vtkSmartPointer<vtkLabeledDataMapper> singleNodeLabelMapper;
//        vtkSmartPointer<vtkActor2D> singleNodeLabelActor;
//
//        size_t singleNodeID;
//        unsigned char nodeClr[4][3];
        
        GlidePlaneActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const,const Polycrystal<3>& poly,const DDtraitsIO& traitsIO);
        void updateConfiguration(const DDauxIO<3>& auxIO);
        
        
    };
    
} // namespace model
#endif