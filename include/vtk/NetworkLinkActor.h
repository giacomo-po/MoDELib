/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLinkActor_h_
#define model_NetworkLinkActor_h_


#include <deque>
#include <string>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>
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
#include <IDreader.h>
#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <DDconfigFields.h>

namespace model
{
    struct NetworkLinkActor : public QWidget
    {
        static constexpr int dim=3;
        enum ColorScheme {colorBurgers=0,colorSessile=1,colorNormal=2,colorEdgeScrew=3,colorComponent=4};
        typedef Eigen::Matrix<double,dim,1>  VectorDim;
        
        Q_OBJECT
        private slots:
        void modify();
        Eigen::Matrix<int,3,1> vector2Clr(VectorDim) const;
        
        private:
        
        
        vtkGenericOpenGLRenderWindow* const renderWindow;

        QGridLayout* mainLayout;
        QCheckBox* showLinks;
        QSlider* sliderLinksRadius;

        QCheckBox* showZeroLinks;
        QComboBox* linksColorBox;

//        double tubeRadius;
        ColorScheme clr;


        vtkSmartPointer<vtkPolyData> polyData;
        vtkSmartPointer<vtkPolyData> polyDataBnd;
        vtkSmartPointer<vtkPolyData> polyData0;
        vtkSmartPointer<vtkTubeFilter> tubeFilter;
        vtkSmartPointer<vtkTubeFilter> tubeFilterBnd;
        vtkSmartPointer<vtkTubeFilter> tubeFilter0;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper;
        vtkSmartPointer<vtkActor> tubeActor;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapperBnd;
        vtkSmartPointer<vtkActor> tubeActorBnd;
        vtkSmartPointer<vtkPolyDataMapper> tubeMapper0;
        vtkSmartPointer<vtkActor> tubeActor0;
        
        public:
        
        const DDconfigFields<3>& configFields;
                
        NetworkLinkActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const,const DDconfigFields<3>& configFields_in);
        void updateConfiguration(vtkPolyData* const nodePolyData);

        Eigen::Matrix<int,3,1> computeColor(const VectorDim& burgers, const VectorDim& chord, const VectorDim& planeNormal) const;
    };
    
} // namespace model
#endif
