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

namespace model
{
    struct NetworkLinkActor : public QWidget
    {
        
        Q_OBJECT
        private slots:
            void modify();

        private:
        vtkGenericOpenGLRenderWindow* const renderWindow;

        QGridLayout* mainLayout;
        QCheckBox* showLinks;

        
        public:
                
        
        NetworkLinkActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const);
        void updateConfiguration(const DDconfigIO<3>& configIO);
        
        
    };
    
} // namespace model
#endif
