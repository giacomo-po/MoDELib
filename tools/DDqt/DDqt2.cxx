#include <QApplication>

#include <vtkActor.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkAxes.h>
#include <vtkProperty.h>
#include <vtkCamera.h>

//#include <QVTKWidget.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <QVTKOpenGLNativeWidget.h>
#include <QSurfaceFormat>

//#define NO_QT
//#define OLD_WIDGET

int main(int argc, char** argv)
{
  // Create X,Y,Z axes at the origin
  vtkNew<vtkAxes> centerAxes;
  centerAxes->SetOrigin(0, 0, 0);
  centerAxes->SetSymmetric(1);
  centerAxes->SetComputeNormals(1);
  vtkNew<vtkPolyDataMapper> axesMapper;
  axesMapper->SetInputConnection(centerAxes->GetOutputPort());
  vtkSmartPointer<vtkActor> centerAxesActor = vtkSmartPointer<vtkActor>::New();
  centerAxesActor->SetMapper(axesMapper);
  centerAxesActor->GetProperty()->SetLighting(false);
  centerAxesActor->PickableOff();
  centerAxesActor->SetScale(0.4);

  vtkNew<vtkRenderer> renderer;
  renderer->AddActor(centerAxesActor);
  renderer->SetBackground(0.06, 0.2, 0.5);

  double pos[3] = { 1, 0.2, 1 };
  double focalPoint[3] = { 0, 0, 0 };
  double viewUp[3] = { 0, 1, 0 };
  renderer->GetActiveCamera()->SetPosition(pos);
  renderer->GetActiveCamera()->SetFocalPoint(focalPoint);
  renderer->GetActiveCamera()->SetViewUp(viewUp);

#ifdef NO_QT
  vtkNew<vtkRenderWindow> renderWindow;
  renderWindow->AddRenderer(renderer);
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindow->SetSize(512, 512);
  renderWindow->Render();
  renderWindowInteractor->Start();
#else
  QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
  QApplication app(argc, argv);
#  if defined OLD_WIDGET
    QVTKWidget widget;
#  else
    QVTKOpenGLNativeWidget widget;
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkNew<QVTKInteractor> renderWindowInteractor;
//    widget.installEventFilter((QWidget*)&widget);

//    vtkNew<vtkInteractorStyleTrackballCamera> style;
//    renderWindowInteractor->SetInteractorStyle( style );
//    renderWindowInteractor->SetRenderWindow(renderWindow);
    
    
//    widget.interactor()->Start();
    widget.setRenderWindow(renderWindow);
#  endif // OLD_WIDGET
    widget.resize(512, 512);
    widget.renderWindow()->AddRenderer(renderer);
    widget.show();
    app.exec();
#endif // NO_QT

  return EXIT_SUCCESS;
}
