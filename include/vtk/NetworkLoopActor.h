/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLoopActor_h_
#define model_NetworkLoopActor_h_


#include <deque>
#include <string>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>
#include <QSlider>
#include <QColorDialog>
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
#include <vtkTriangleFilter.h>
#include <IDreader.h>
#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <PeriodicGlidePlaneFactory.h>
#include <DislocationLoopPatches.h>

namespace model
{
    struct NetworkLoopActor : public QWidget
    /*                     */,private std::map<size_t,DislocationLoopPatches<3>>
    {
//        static constexpr int dim=3;
//        enum ColorScheme {colorBurgers=0,colorSessile=1,colorNormal=2,colorEdgeScrew=3,colorComponent=4};
//        typedef Eigen::Matrix<double,dim,1>  VectorDim;
        
        Q_OBJECT
        private slots:
            void modify();

        private:
        vtkGenericOpenGLRenderWindow* const renderWindow;

        QGridLayout* mainLayout;
        QCheckBox* showLoops;
        
        QGroupBox* slippedAreaBox;
//        QCheckBox* showSlippedArea;
        QSlider* sliderSlippedArea;
//        QColorDialog* colorSlippedArea;

//        double tubeRadius;
//        ColorScheme clr;


//        vtkSmartPointer<vtkPoints> loopPoints;
        vtkSmartPointer<vtkPolyData> loopPolyData;
        vtkSmartPointer<vtkPolyDataMapper> loopMapper;
        vtkSmartPointer<vtkActor> loopActor;

        
        vtkSmartPointer<vtkPolyData> areaPolyData;
        vtkSmartPointer<vtkTriangleFilter> areaTriangleFilter; // needed to handle non-convex area polygons
        vtkSmartPointer<vtkPolyDataMapper> areaMapper;
        vtkSmartPointer<vtkActor> areaActor;

        
        public:
                
        const Polycrystal<3>& poly;
        PeriodicGlidePlaneFactory<3>& periodicGlidePlaneFactory;
        
        NetworkLoopActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const,
                         const Polycrystal<3>& poly_in,
                         PeriodicGlidePlaneFactory<3>& pgf);
        void updateConfiguration(const DDconfigIO<3>& configIO);
        const std::map<size_t,DislocationLoopPatches<3>>& loopPatches() const;
        std::map<size_t,DislocationLoopPatches<3>>& loopPatches();

//        Eigen::Matrix<int,3,1> computeColor(const VectorDim& burgers, const VectorDim& chord, const VectorDim& planeNormal) const;
        
    };
    
} // namespace model
#endif
