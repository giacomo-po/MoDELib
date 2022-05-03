/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigVtk_h_
#define model_DDconfigVtk_h_

#include <QGroupBox>
#include <QGridLayout>
#include <QCheckBox>
#include <QLabel>
#include <QTabWidget>
#include <QTextEdit>
#include <QLineEdit>
#include <QPushButton>

#include <vtkRenderer.h>
#include <vtkGenericOpenGLRenderWindow.h>

#include <DDconfigIO.h>
//#include <DDconfigVtkBase.h>
#include <SimplicialMesh.h>
#include <NetworkNodeActor.h>
//#include <DislocationSegmentActor.h>
//#include <DislocationLoopActor.h>

namespace model
{
    struct DDconfigVtk : public QWidget,
                         public DDconfigIO<3>
//: public DDconfigVtkBase
    /*                */
    {
        
        Q_OBJECT
        
        vtkGenericOpenGLRenderWindow* const renderWindow;
        NetworkNodeActor* nodes;
//        DislocationSegmentActor segments;
//        DislocationLoopActor loops;
        private slots:
        void updateConfiguration();
        void nextConfiguration();
        void prevConfiguration();

        public:

        
//        QGroupBox* controlsBox;
        QGridLayout* mainLayout;
        QLineEdit* frameIDedit;
        QPushButton* plusFrameButton;
        QPushButton* minusFrameButton;
        QLineEdit* frameIncrementEdit;

        QTabWidget* tabWidget;


        
        /**********************************************************************/
        DDconfigVtk(const std::string& folderName,vtkGenericOpenGLRenderWindow* const, vtkRenderer* const ren,const SimplicialMesh<3>& mesh);
        void updateConfiguration(const size_t& frameID);

        void modify();
        
    };
    
}
#endif
