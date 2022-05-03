#include "widget.h"
#include "ui_widget.h"

#include <vtkGenericOpenGLRenderWindow.h>

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget),
    coneActor( vtkSmartPointer<vtkActor>::New() )
{
    ui->setupUi(this);

    vtkSmartPointer<vtkConeSource> cone =
            vtkSmartPointer<vtkConeSource>::New();

    vtkSmartPointer<vtkPolyDataMapper> coneMapper =
            vtkSmartPointer<vtkPolyDataMapper>::New();
    coneMapper->SetInputConnection( cone->GetOutputPort() );
    coneActor->SetMapper( coneMapper );
    vtkSmartPointer<vtkRenderer> backRenderer =
            vtkSmartPointer<vtkRenderer>::New();
    backRenderer->AddActor( coneActor );
    backRenderer->SetBackground( 1, 1, 1 );

    vtkSmartPointer<vtkGenericOpenGLRenderWindow> backWindow =
            vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    backWindow->AddRenderer( backRenderer );

    ui->qvtkWidget->SetRenderWindow( backWindow );
    //ui->qvtkWidget->GetRenderWindow()->AddRenderer( backRenderer );
}

Widget::~Widget()
{
    delete ui;
}

void Widget::on_pushButton_clicked()
{

}
