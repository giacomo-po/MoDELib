/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ChartActor_cpp_
#define model_ChartActor_cpp_

#include <fstream>
#include <string>
#include <vector>

#include <vtkTable.h>
#include <vtkFloatArray.h>
#include <vtkPlot.h>
#include <vtkPoints.h>
#include <vtkAxis.h>

#include <ChartActor.h>
#include <vtkTextProperty.h>


namespace model
{

        ChartActor::ChartActor(const DDtraitsIO& traitsIO_in,vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_in,vtkRenderer* const renderer_in) :
        /* init */ traitsIO(traitsIO_in)
        /* init */,mainLayout(new QGridLayout(this))
        /* init */,showChart(new QCheckBox(this))
        /* init */,chartUpdate(new QPushButton("Reload",this))
        /* init */,xComboBox(new QComboBox(this))
        /* init */,yComboBox(new QComboBox(this))
        /* init */,xAxisLog(new QCheckBox(this))
        /* init */,yAxisLog(new QCheckBox(this))
        /* init */,renderWindow(renderWindow_in)
        /* init */,renderer(renderer_in)
        /* init */,chart(vtkSmartPointer<vtkChartXY>::New())
        /* init */,chartScene(vtkSmartPointer<vtkContextScene>::New())
        /* init */,chartActor(vtkSmartPointer<vtkContextActor>::New())
        /* init */,table(vtkSmartPointer<vtkTable>::New())
        /* init */,currentTable(vtkSmartPointer<vtkTable>::New())
        /* init */,colors(vtkSmartPointer<vtkNamedColors>::New())
        /* init */,points(chart->AddPlot(vtkChart::LINE))
        /* init */,currentPoints(chart->AddPlot(vtkChart::LINE))
        {

            showChart->setChecked(false);
            showChart->setText("chart");
            chartUpdate->setEnabled(false);
            xComboBox->setEnabled(false);
            yComboBox->setEnabled(false);
            chartActor->SetVisibility(false);
            
            xAxisLog->setChecked(false);
            xAxisLog->setEnabled(false);
            xAxisLog->setText("log-scale");

            yAxisLog->setChecked(false);
            yAxisLog->setEnabled(false);
            yAxisLog->setText("log-scale");
            
            vtkColor3d color3d1 = colors->GetColor3d("black");
            points->SetColor(color3d1.GetRed(), color3d1.GetGreen(), color3d1.GetBlue());
            points->SetWidth(5.0);

            vtkColor3d color3d2 = colors->GetColor3d("magenta");
            currentPoints->SetColor(color3d2.GetRed(), color3d2.GetGreen(), color3d2.GetBlue());
            currentPoints->SetWidth(5.0);

            std::ifstream flabFile(traitsIO.flabFile);
            if(flabFile)
            {
                std::cout<<"Reading "<<traitsIO.flabFile<<" ... "<<std::flush;
                const auto t0= std::chrono::system_clock::now();

                std::string line;
                while (std::getline(flabFile, line))
                {
                    if(!line.empty())
                    {
//                        std::cout<<line<<std::endl;
                        xComboBox->addItem(QString::fromStdString(line));
                        yComboBox->addItem(QString::fromStdString(line));
                        
//                        vtkNew<vtkFloatArray> farr;
//                        farr->SetName(line.c_str());
//                        table->AddColumn(farr);
                    }
                }
                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            else
            {
                std::cout<<"Unable to read "<<traitsIO.flabFile<<std::endl;
            }

            mainLayout->addWidget(showChart,0,0,1,1);
            mainLayout->addWidget(chartUpdate,0,1,1,1);
            mainLayout->addWidget(xComboBox,1,0,1,1);
            mainLayout->addWidget(xAxisLog,1,1,1,1);
            mainLayout->addWidget(yComboBox,2,0,1,1);
            mainLayout->addWidget(yAxisLog,2,1,1,1);
            mainLayout->setColumnStretch(0, 1);
            mainLayout->setColumnStretch(1, 12);
            this->setLayout(mainLayout);

            connect(showChart,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(xAxisLog, SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(yAxisLog, SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(xComboBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));
            connect(yComboBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));
            connect(chartUpdate,SIGNAL(released()), this, SLOT(reload()));

            chartScene->AddItem(chart);
            chartActor->SetScene(chartScene);
            renderer->AddActor(chartActor);
            chartScene->SetRenderer(renderer);
            
        }
        
    void ChartActor::reload()
    {
        table->Initialize();

        const int cols(xComboBox->count());

        for(int c=0;c<cols;++c)
        {
            vtkNew<vtkFloatArray> farr;
            farr->SetName(xComboBox->itemText(c).toStdString().c_str());
            table->AddColumn(farr);
        }
        
        
        std::ifstream fFile(traitsIO.fFile);
        if(fFile)
        {
            std::cout<<"Reading "<<traitsIO.fFile<<" ... "<<std::flush;
            const auto t0= std::chrono::system_clock::now();
            double temp(0.0);
            std::string line;
            while (std::getline(fFile, line))
            {
                if(!line.empty())
                {
                    const auto row=table->InsertNextBlankRow(cols);
                    std::stringstream ss(line);
                    for(int col=0;col<cols;++col)
                    {
                        ss >> temp;
                        table->SetValue(row, col, temp);
                    }
                }
            }
            std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        else
        {
            std::cout<<"Unable to read "<<traitsIO.fFile<<std::endl;
        }

        table->Modified();
        updatePlots();
    }

    void ChartActor::updateConfiguration(const size_t& frameID)
    {
        
        if(frameID==0)
        {
            reload();
        }
        
        std::cout<<"Updating chart..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();

        
        currentTable->Initialize();
        
        for(long long k=0;k<table->GetNumberOfColumns();++k)
        {
            vtkNew<vtkFloatArray> farr;
            farr->SetName(table->GetColumnName(k));
            currentTable->AddColumn(farr);
        }
        
        for(long long k=0;k<table->GetNumberOfRows();++k)
        {
            if(table->GetValue(k,0)<=frameID)
            {
                const auto row=currentTable->InsertNextBlankRow(table->GetNumberOfColumns());
                for(int col=0;col<table->GetNumberOfColumns();++col)
                {
                    currentTable->SetValue(row, col, table->GetValue(row,col));
                }
            }
        }
        currentTable->Modified();
        updatePlots();
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;

    }

void ChartActor::updatePlots()
{
    if(xComboBox->count() && yComboBox->count())
    {
        if(   table->GetNumberOfColumns()>xComboBox->currentIndex()
           && table->GetNumberOfColumns()>yComboBox->currentIndex())
        {
            points->SetInputData(table, xComboBox->currentIndex(), yComboBox->currentIndex());
        }
        
        if(   currentTable->GetNumberOfColumns()>xComboBox->currentIndex()
           && currentTable->GetNumberOfColumns()>yComboBox->currentIndex())
        {
            currentPoints->SetInputData(currentTable, xComboBox->currentIndex(), yComboBox->currentIndex());
        }
        
        chart->RecalculateBounds();
        chart->Modified();
    }
}

 
        
        void ChartActor::modify()
        {
//            std::cout<<"modifying chart..."<<std::endl;

            xComboBox->setEnabled(showChart->isChecked());
            yComboBox->setEnabled(showChart->isChecked());
            xAxisLog->setEnabled(showChart->isChecked());
            yAxisLog->setEnabled(showChart->isChecked());
            chartUpdate->setEnabled(showChart->isChecked());


            chart->GetAxis(0)->SetTitle(yComboBox->currentText().toStdString());
            chart->GetAxis(1)->SetTitle(xComboBox->currentText().toStdString());

            chart->GetAxis(0)->SetLogScale(yAxisLog->isChecked());
            chart->GetAxis(1)->SetLogScale(xAxisLog->isChecked());

            chart->GetAxis(0)->SetGridVisible(false);
            chart->GetAxis(1)->SetGridVisible(false);

            chart->GetAxis(0)->SetNumberOfTicks(10);
            chart->GetAxis(1)->SetNumberOfTicks(10);

            chart->GetAxis(0)->GetLabelProperties()->SetFontSize(16);
            chart->GetAxis(1)->GetLabelProperties()->SetFontSize(16);

            chart->GetAxis(0)->GetTitleProperties()->SetFontSize(16);
            chart->GetAxis(1)->GetTitleProperties()->SetFontSize(16);
            
            
            updatePlots();

            
//            std::cout<<"xComboBox->currentIndex()="<<xComboBox->currentIndex()<<std::endl;
//            std::cout<<"yComboBox->currentIndex()="<<yComboBox->currentIndex()<<std::endl;
//            std::cout<<"table->GetNumberOfRows()="<<table->GetNumberOfRows()<<std::endl;
//            std::cout<<"table->GetNumberOfColumns()="<<table->GetNumberOfColumns()<<std::endl;
//            std::cout<<"currentTable->GetNumberOfRows()="<<currentTable->GetNumberOfRows()<<std::endl;
//            std::cout<<"currentTable->GetNumberOfColumns()="<<currentTable->GetNumberOfColumns()<<std::endl;

            
//            chart->ClearPlots();

//            if(xComboBox->count() && yComboBox->count())
//            {
//                if(   table->GetNumberOfColumns()>xComboBox->currentIndex()
//                   && table->GetNumberOfColumns()>yComboBox->currentIndex()
////                   && table->GetNumberOfRows()>2
//                   )
//                {
////                    vtkPlot* points = chart->AddPlot(vtkChart::LINE);
//                    points->SetInputData(table, xComboBox->currentIndex(), yComboBox->currentIndex());
////                    vtkColor3d color3d1 = colors->GetColor3d("black");
////                    points->SetColor(color3d1.GetRed(), color3d1.GetGreen(), color3d1.GetBlue());
////                    points->SetWidth(2.0);
//                }
//
//                if(   currentTable->GetNumberOfColumns()>xComboBox->currentIndex()
//                   && currentTable->GetNumberOfColumns()>yComboBox->currentIndex()
////                   && currentTable->GetNumberOfRows()>2
//                   )
//                {
////                    vtkPlot* currentPoints = chart->AddPlot(vtkChart::LINE);
//                    currentPoints->SetInputData(currentTable, xComboBox->currentIndex(), yComboBox->currentIndex());
////                    vtkColor3d color3d2 = colors->GetColor3d("magenta");
////                    currentPoints->SetColor(color3d2.GetRed(), color3d2.GetGreen(), color3d2.GetBlue());
////                    currentPoints->SetWidth(3.0);
//                }
//
//                chart->RecalculateBounds();
//                chart->Modified();
//
//            }
            
            chartActor->SetVisibility(showChart->isChecked());
            renderWindow->Render(); // this makes code crash if called before the app window is created
        }



} // namespace model
#endif







