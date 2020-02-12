/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeActor_H_
#define model_DislocationNodeActor_H_

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>

#include <model/DislocationDynamics/Visualization/vtk/DislocationSegmentActor.h>


namespace model
{
    
    /************************************************************************/
    /************************************************************************/
    struct DislocationNodeActor
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<float,dim,1>  VectorDim;
        
        vtkSmartPointer<vtkSphereSource> sphereSource;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;

        /************************************************************************/
        DislocationNodeActor(const Eigen::Matrix<float,dim,1>& P) :
        /* init */ sphereSource(vtkSmartPointer<vtkSphereSource>::New()),
        /* init */ mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ actor(vtkSmartPointer<vtkActor>::New())
        {
            
            sphereSource->SetCenter(P(0), P(1), P(2));
            sphereSource->SetRadius(DislocationSegmentActor::tubeRadius*1.2);
            mapper->SetInputConnection(sphereSource->GetOutputPort());
            actor->SetMapper(mapper);
            
        }
        
        /************************************************************************/
        void modify()
        {
            sphereSource->SetRadius(DislocationSegmentActor::tubeRadius*1.2);
//            sphereSource->Modified();
        }

    };
    
} // namespace model
#endif







