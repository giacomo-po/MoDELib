/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigVtk_cpp_
#define model_DDconfigVtk_cpp_

#include <algorithm>
#include <QString>

#include <DDconfigVtk.h>
#include <SimplicialMesh.h>

//#include <QSize>
//#include <QPixmap>
#include <QImage>

namespace model
{
        
        /**********************************************************************/
        DDconfigVtk::DDconfigVtk(const DDtraitsIO& traitsIO_in,vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const ren,
                                 QVTKOpenGLStereoWidget* const qvtkGLwidget_in,
                                 const Polycrystal<3>& poly,PeriodicGlidePlaneFactory<3>& pgpf) :
        /* init */ DDconfigIO<3>(traitsIO_in.evlFolder,"")
        /* init */,DDauxIO<3>(traitsIO_in.auxFolder,"")
        /* init */,renderWindow(renWin)
        /* init */,qvtkGLwidget(qvtkGLwidget_in)
        /* init */,traitsIO(traitsIO_in)
        /* init */,nodes(new NetworkNodeActor(renWin,ren))
        /* init */,segments(new NetworkLinkActor(renWin,ren))
        /* init */,loops(new NetworkLoopActor(renWin,ren,poly,pgpf))
        /* init */,inclusions(new InclusionActor(renWin,ren))
        /* init */,glidePlanes(new GlidePlaneActor(renWin,ren,poly,traitsIO))
        /* init */,quadrature(new QuadratureActor(renWin,ren,poly,traitsIO))
        /* init */,chartActor(new ChartActor(traitsIO,renderWindow,ren))
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,frameIDedit(new QLineEdit("0"))
        /* init */,plusFrameButton(new QPushButton(">"))
        /* init */,minusFrameButton(new QPushButton("<"))
        /* init */,frameIncrementEdit(new QLineEdit("1"))
        /* init */,saveImage(new QCheckBox(this))
        /* init */,tabWidget(new QTabWidget())
        {
            
            tabWidget->addTab(nodes, tr(std::string("Nodes").c_str()));
            tabWidget->addTab(segments, tr(std::string("Segments").c_str()));
            tabWidget->addTab(loops, tr(std::string("Loops").c_str()));
            tabWidget->addTab(inclusions, tr(std::string("Inclusions").c_str()));
            tabWidget->addTab(glidePlanes, tr(std::string("GlidePlanes").c_str()));
            tabWidget->addTab(quadrature, tr(std::string("Quadrature").c_str()));
            tabWidget->addTab(chartActor, tr(std::string("Chart").c_str()));
            
            saveImage->setText("save PNG");


            mainLayout->addWidget(frameIDedit,0,0,1,1);
            mainLayout->addWidget(minusFrameButton,0,1,1,1);
            mainLayout->addWidget(plusFrameButton,0,2,1,1);
            mainLayout->addWidget(frameIncrementEdit,0,3,1,1);
            mainLayout->addWidget(saveImage,1,0,1,1);


            mainLayout->addWidget(tabWidget,2,0,1,4);
//            controlsBox->setLayout(mainLayout);
            this->setLayout(mainLayout);

            
            updateConfiguration(0);
            
            
            connect(frameIDedit,SIGNAL(returnPressed()), this, SLOT(updateConfiguration()));
            connect(plusFrameButton,SIGNAL(pressed()), this, SLOT(nextConfiguration()));
            connect(minusFrameButton,SIGNAL(pressed()), this, SLOT(prevConfiguration()));
            

//            connect(frameIDedit,SIGNAL(keyPressEvent(QKeyEvent*)), this, SLOT(updateConfiguration()));
        }

        const DDconfigIO<3>& DDconfigVtk::configIO() const
        {
            return *this;
        }

        DDconfigIO<3>& DDconfigVtk::configIO()
        {
            return *this;
        }

        const DDauxIO<3>& DDconfigVtk::auxIO() const
        {
            return *this;
        }

        DDauxIO<3>& DDconfigVtk::auxIO()
        {
            return *this;
        }

        void DDconfigVtk::nextConfiguration()
        {
            const long int currentFrameID(std::atoi(frameIDedit->text() .toStdString().c_str()));
            const long int currentIncrement(std::atoi(frameIncrementEdit->text() .toStdString().c_str()));
            const long int nextFrameID(std::max((long int)0,currentFrameID+currentIncrement));
            frameIDedit->setText(QString::fromStdString(std::to_string(nextFrameID)));
            if(!updateConfiguration())
            {
                frameIDedit->setText(QString::fromStdString(std::to_string(currentFrameID)));
            }
        }

        void DDconfigVtk::prevConfiguration()
        {
            const long int currentFrameID(std::atoi(frameIDedit->text() .toStdString().c_str()));
            const long int currentIncrement(std::atoi(frameIncrementEdit->text() .toStdString().c_str()));
            const long int nextFrameID(std::max((long int)0,currentFrameID-currentIncrement));
            frameIDedit->setText(QString::fromStdString(std::to_string(nextFrameID)));
            if(!updateConfiguration())
            {
                frameIDedit->setText(QString::fromStdString(std::to_string(currentFrameID)));
            }
        }

        bool DDconfigVtk::updateConfiguration()
        {
            const size_t frameID(std::atoi(frameIDedit->text() .toStdString().c_str()));
            if(updateConfiguration(frameID))
            {
                renderWindow->Render();
                return true;
            }
            else
            {
                return false;
            }
        }

        bool DDconfigVtk::updateConfiguration(const size_t& frameID)
        {
            try
            {
                configIO().read(frameID);
                auxIO().read(frameID);
                nodes->updateConfiguration(*this);
                segments->updateConfiguration(*this,nodes->nodePolyData);
                loops->updateConfiguration(*this,nodes->nodePolyData);
                inclusions->updateConfiguration(*this);
                glidePlanes->updateConfiguration(*this);
                quadrature->updateConfiguration(*this);
                chartActor->updateConfiguration(frameID);
                
                if(saveImage->isChecked())
                {
                    QImage img=qvtkGLwidget->grabFramebuffer();
                    img.save(QString::fromStdString(traitsIO.evlFolder+"/img_"+std::to_string(frameID)+".png"), "PNG", -1);

                }

                
                return true;
            }
            catch(const std::exception& e)
            {
                std::cout<<e.what()<<std::endl;
                return false;
            }
        }
        
        void DDconfigVtk::modify()
        {
//            nodes.modify();
//            segments.modify();
//            loops.modify();
        }
        
}
#endif
