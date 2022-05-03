/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_NetworkLinkActor_cpp_
#define model_NetworkLinkActor_cpp_

#include <iostream>
#include <deque>
#include <string>

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
#include <vtkTextProperty.h>
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
//#include <IDreader.h>
//#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <NetworkLinkActor.h>
#include <vtkRenderer.h>

namespace model
{
    
        
        /**********************************************************************/
        NetworkLinkActor::NetworkLinkActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer) :
        /* init */ renderWindow(renWin)
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,showLinks(new QCheckBox(this))
        {
            showLinks->setChecked(true);
            showLinks->setText("show links");

            mainLayout->addWidget(showLinks,0,0,1,1);
            this->setLayout(mainLayout);

            connect(showLinks,SIGNAL(stateChanged(int)), this, SLOT(modify()));

            

        }
        

        
        /**********************************************************************/
        void NetworkLinkActor::updateConfiguration(const DDconfigIO<3>& configIO)
        {// https://stackoverflow.com/questions/6878263/remove-individual-points-from-vtkpoints
            std::cout<<"Updating links..."<<std::flush;
            const auto t0= std::chrono::system_clock::now();

            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void NetworkLinkActor::modify()
        {
            
//            nodeActor->SetVisibility(showLinks->isChecked());
//            labelActor->SetVisibility(showLinks->isChecked());
//            velocityActor->SetVisibility(showLinks->isChecked());
            
//            nodeGlyphs->SetScaleFactor(2.0*this->tubeRadius*1.2);
//            
//            if(this->showVelocities)
//            {
//                velocityActor->VisibilityOn();
//            }
//            else
//            {
//                velocityActor->VisibilityOff();
//            }
//            
//            if(this->showNodeIDs)
//            {
//                labelActor->VisibilityOn();
//                
//            }
//            else
//            {
//                labelActor->VisibilityOff();
//            }
//            
//            velocityGlyphs->SetScaleFactor(this->velocityFactor);
//            
//            
//            if(this->showSingleNode)
//            {
//                // HERE WE SHOULD CHANGE THE NODE POSITION BASED ON NODE ID
//                // OTHERWISE THE SELECTED NODE WILL BE VISIBLE ONLY UPON LOADING A NEW FRAME
//                std::cout<<"RELOAD FRAME TO SHOW SELECTED NODE"<<std::endl;
//                singleNodeLabelActor->VisibilityOn();
//            }
//            else
//            {
//                singleNodeLabelActor->VisibilityOff();
//            }
//            
//            if(this->showLinks)
//            {
//                nodeActor->VisibilityOn();
//            }
//            else
//            {
//                nodeActor->VisibilityOff();
//            }

            renderWindow->Render();
        }
        
    
} // namespace model
#endif
