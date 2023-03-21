/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionActor_H_
#define model_InclusionActor_H_

#include <Eigen/Dense>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkLookupTable.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkFloatArray.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkDataSetMapper.h>

#include <DDconfigIO.h>


// VTK documentation
// http://vtk.1045678.n5.nabble.com/VTK-slow-to-display-300-vtkSphereSource-in-real-time-td5740730.html
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/Visualization/Scaleglyph3D
// http://hiraku0n.blogspot.com/2012/10/writing-scalar-and-vector-field-in-vtk.html
// http://vtk.1045678.n5.nabble.com/Glyphing-vtkImageData-scalars-3D-as-arrows-td3199837.html

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    struct InclusionActor : public QWidget
    {
        
        static constexpr int dim=3;
        typedef Eigen::Matrix<float,dim,1>  VectorDim;
        typedef Eigen::Matrix<float,dim,dim,1>  MatrixDim;

        
        Q_OBJECT
        private slots:
            void modify();

        private:
        vtkGenericOpenGLRenderWindow* const renderWindow;
        QGridLayout* mainLayout;
        QCheckBox* showInclusions;

        vtkSmartPointer<vtkUnstructuredGrid> grid;
        vtkSmartPointer<vtkSphereSource> sphereSource;
        vtkSmartPointer<vtkGlyph3D> glyph3D;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;
        vtkSmartPointer<vtkLookupTable> lookUpColors;
        
        vtkSmartPointer<vtkDataSetMapper> polyhedronMapper;
        vtkSmartPointer<vtkActor> polyhedronActor;
        
//        vtkSmartPointer<vtkPoints> polyhedronPoints;
        
    public:
        InclusionActor(vtkGenericOpenGLRenderWindow* const,vtkRenderer* const);
        void updateConfiguration(const DDconfigIO<3>& configIO);
        
    };
        
}
#endif
