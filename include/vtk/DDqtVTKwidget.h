/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDqtVTKwidget_h_
#define model_DDqtVTKwidget_h_

#include <QGridLayout>

#include <QVTKOpenGLStereoWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>

#include <TextFileParser.h>
#include <SimplicialMesh.h>
#include <SimplicialMeshActor.h>
#include <DDconfigVtk.h>

namespace model
{
    

    

struct DDqtVTKwidget : public QWidget
//public QVTKOpenGLStereoWidget
{
    Q_OBJECT
    

    
private:
    
    QGridLayout* mainLayout;
    QTabWidget* tabWidget;
    QVTKOpenGLStereoWidget* openglWidget;
    
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderer> renderer;
    
    std::string getWorkingDir();
    

public:
    
    const std::string workingDir;
    QLabel* workingDirLabel;
    SimplicialMesh<3> mesh;
    SimplicialMeshActor* meshActor;
    DDconfigVtk* ddConfigVtk;
    
    
    DDqtVTKwidget(QWidget *parent);
    

};


} // namespace model

#endif







