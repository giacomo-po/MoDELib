/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PKActor_H_
#define model_PKActor_H_

#include <Eigen/Dense>

#include <vtkVersion.h>
#include <vtkArrowSource.h>
#include <vtkTransform.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkTransformPolyDataFilter.h>

#include <model/DislocationDynamics/Visualization/vtk/DislocationSegmentActor.h>


namespace model
{
    
    /************************************************************************/
    /************************************************************************/
    struct PKActor
    {
        static constexpr int dim=3;
        typedef Eigen::Matrix<float,dim,1>  VectorDim;
        typedef Eigen::Matrix<float,dim,dim,1>  MatrixDim;
        
        vtkSmartPointer<vtkTransform> transform;
        vtkSmartPointer<vtkTransformPolyDataFilter> transformPD;
        vtkSmartPointer<vtkArrowSource> arrowSource;
        vtkSmartPointer<vtkPolyDataMapper> mapper;
        vtkSmartPointer<vtkActor> actor;

        /************************************************************************/
        PKActor(const VectorDim& P,
                const VectorDim& v) :
        /* init */ transform(vtkSmartPointer<vtkTransform>::New()),
        /* init */ transformPD(vtkSmartPointer<vtkTransformPolyDataFilter>::New()),
        /* init */ arrowSource(vtkSmartPointer<vtkArrowSource>::New()),
        /* init */ mapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ actor(vtkSmartPointer<vtkActor>::New())
        {
            
            const double length=v.norm();
            
            if(length>FLT_EPSILON)
            {
                MatrixDim R(MatrixDim::Zero());
                R.col(0)=v/length;
                
                while(R.col(1).squaredNorm()<FLT_EPSILON)
                {
                    R.col(1)=VectorDim::Random().cross(R.col(0));
                }
                R.col(1).normalize();
                R.col(2)=R.col(0).cross(R.col(1));

                
                                double startPoint[3];
                
                vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
                for (unsigned int i = 0; i < 3; i++)
                {
                    matrix->SetElement(i, 0, R(i,0));
                    matrix->SetElement(i, 1, R(i,1));
                    matrix->SetElement(i, 2, R(i,2));
                    
                    startPoint[i]=P(i);
                }
                
                std::cout<<startPoint[0]<<","<<startPoint[1]<<","<<startPoint[2]<<std::endl;

                // Apply the transforms
                transform->Translate(startPoint);
                //transform->Concatenate(matrix);
                transform->Scale(length*10000, length*10000, length*10000);

                // Transform the polydata
                transformPD->SetTransform(transform);
                transformPD->SetInputConnection(arrowSource->GetOutputPort());
             
                //Update mapper and actor for the arrow
                mapper->SetInputConnection(arrowSource->GetOutputPort());
//                mapper->SetInputConnection(transformPD->GetOutputPort());
                actor->SetMapper(mapper);

                
            }
            

            

            
            
            

            
//            sphereSource->SetCenter(P(0), P(1), P(2));
//            sphereSource->SetRadius(DislocationSegmentActor::tubeRadius*1.2);
//            mapper->SetInputConnection(sphereSource->GetOutputPort());
//            actor->SetMapper(mapper);
//            actor->setPosition();
        }
        
        /************************************************************************/
        void modify()
        {
            //arrowSource->SetShaftRadius(1.0);
            //arrowSource->SetTipLength(1.0);
            arrowSource->Update();

        }
        

    };
    
} // namespace model
#endif







