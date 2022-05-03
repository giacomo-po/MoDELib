#include <stdio.h>
#include <QVTKOpenGLNativeWidget.h>
#include <QWidget>
#include <QApplication>
#include <QVBoxLayout>
#include <vtkSphereSource.h>
#include <vtkActor.h>
#include <vtkSphere.h>
#include <vtkSphereSource.h>
#include <vtkCone.h>
#include <vtkConeSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <QKeyEvent>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCamera.h>

class MyWidget : public QWidget {
public:
  MyWidget()
    :qvtk_widget_(new QVTKOpenGLNativeWidget(this))
  {
    QVBoxLayout* layout = new QVBoxLayout;
    setLayout(layout);
    layout->addWidget(qvtk_widget_);
    qvtk_widget_->setFixedSize(640,480);
    renderer_ = vtkSmartPointer<vtkRenderer>::New();
    renderer_->SetBackground(1.0, 1.0, 1.0);
    render_window_ = vtkSmartPointer<vtkRenderWindow>::New();
    render_window_->AddRenderer(renderer_);
    qvtk_widget_->SetRenderWindow(render_window_);

    qvtk_widget_->installEventFilter((QWidget*)this);

    tarckball_style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    auto qvtk_interactor = vtkSmartPointer<QVTKInteractor>::New();
    // Make sure that don't auto adjust camera clipping range, because of zoom in/out.
    // AutoAdjust will cause glitch because of wrong adjustment.
    tarckball_style->SetAutoAdjustCameraClippingRange(false);
    qvtk_interactor->SetInteractorStyle(tarckball_style);
    qvtk_interactor->SetRenderWindow(render_window_);

    Initialize();
  }

  vtkSmartPointer<vtkInteractorStyleTrackballCamera> tarckball_style;

  void Initialize(){
    auto sphereSource = vtkSmartPointer<vtkSphereSource>::New();
    sphereSource->SetCenter(1.0, 0.0, 0.0);
    sphereSource->Update();
    auto sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
    auto sphereActor = vtkSmartPointer<vtkActor>::New();
    sphereActor->SetMapper(sphereMapper);
    auto coneSource = vtkSmartPointer<vtkConeSource>::New();
    auto coneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    coneMapper->SetInputConnection(coneSource->GetOutputPort());
    vtkSmartPointer<vtkActor> coneActor = vtkSmartPointer<vtkActor>::New();
    coneActor->SetMapper(coneMapper);
    renderer_->AddActor(sphereActor);
    renderer_->AddActor(coneActor);
    ResetCamera();
    return;
  }

  void ResetCamera(){
    renderer_->ResetCamera();
    // Don't forget reset clipiing range after Reset camera.
    vtkCamera* camera = renderer_->GetActiveCamera();
    camera->SetClippingRange(0.001, 1e+8);
    render_window_->Render();
    return;
  }

  bool eventFilter(QObject *obj, QEvent *ev){
    if (ev->type() == QEvent::KeyPress ) {
      QKeyEvent* keyevent = dynamic_cast<QKeyEvent*>(ev);
      if(keyevent->key() == Qt::Key_Q){
        close();
        // eventFilter intercepts the event.
        return true;
      }
      else if(keyevent->key() == Qt::Key_R){
        ResetCamera();
        return true;
      }
      else{
        vtkCamera* camera = renderer_->GetActiveCamera();
        std::cout << tarckball_style->GetAutoAdjustCameraClippingRange()
          << " : " << camera->GetClippingRange()[0]
          << " : " << camera->GetClippingRange()[1] << std::endl;
      }
    }
    // eventFilter doesn't intercept the event.
    return false;
  }

private:
  QVTKOpenGLNativeWidget*const qvtk_widget_;
  vtkSmartPointer<vtkRenderer> renderer_;
  vtkSmartPointer<vtkRenderWindow> render_window_;
};

int main(int argc, char**argv){
  QApplication app(argc, argv);
  MyWidget mywidget;
  mywidget.show();
  app.exec();
  return 1;
}
