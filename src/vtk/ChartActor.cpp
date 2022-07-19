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


namespace model
{

        ChartActor::ChartActor(const DDtraitsIO& traitsIO,vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_in,vtkRenderer* const renderer_in) :
        /* init */ mainLayout(new QGridLayout(this))
        /* init */,showChart(new QCheckBox(this))
        /* init */,xComboBox(new QComboBox(this))
        /* init */,yComboBox(new QComboBox(this))
        /* init */,xAxisLog(new QCheckBox(this))
        /* init */,yAxisLog(new QCheckBox(this))
        /* init */,renderWindow(renderWindow_in)
        /* init */,renderer(renderer_in)
        /* init */,chart(vtkSmartPointer<vtkChartXY>::New())
        /* init */,chartScene(vtkSmartPointer<vtkContextScene>::New())
        /* init */,chartActor(vtkSmartPointer<vtkContextActor>::New())
//        /* init */,points(chart->AddPlot(vtkChart::LINE))
        /* init */,table(vtkSmartPointer<vtkTable>::New())
        /* init */,colors(vtkSmartPointer<vtkNamedColors>::New())
        {

            showChart->setChecked(false);
            showChart->setText("chart");
            xComboBox->setEnabled(false);
            yComboBox->setEnabled(false);
            chartActor->SetVisibility(false);
            
            xAxisLog->setChecked(false);
            xAxisLog->setEnabled(false);
            xAxisLog->setText("log-scale");

            yAxisLog->setChecked(false);
            yAxisLog->setEnabled(false);
            yAxisLog->setText("log-scale");


//            xComboBox->setText("x-axis");
//            yComboBox->setText("y-axis");
//            vtkNew<vtkTable> table;

            std::ifstream flabFile(traitsIO.flabFile);
            std::string line;
            if(flabFile)
            {
                while (std::getline(flabFile, line))
                {
                    if(!line.empty())
                    {
                        xComboBox->addItem(QString::fromStdString(line));
                        yComboBox->addItem(QString::fromStdString(line));
                        
                        vtkNew<vtkFloatArray> farr;
                        farr->SetName(line.c_str());
                        table->AddColumn(farr);
                    }
                }
            }

            // 200064
            
            // https://stackoverflow.com/questions/116038/how-do-i-read-an-entire-file-into-a-stdstring-in-c/40903508#40903508
            
            const int cols(xComboBox->count());

            std::ifstream fFile(traitsIO.fFile);
            if(fFile)
            {
                float temp;
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
            }

            mainLayout->addWidget(showChart,0,0,1,1);
            mainLayout->addWidget(xComboBox,1,0,1,1);
            mainLayout->addWidget(xAxisLog,1,1,1,1);
            mainLayout->addWidget(yComboBox,2,0,1,1);
            mainLayout->addWidget(yAxisLog,2,1,1,1);
            mainLayout->setColumnStretch(0, 1);
            mainLayout->setColumnStretch(1, 12);
            this->setLayout(mainLayout);

            connect(showChart,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(xAxisLog,SIGNAL(stateChanged(int)), this, SLOT(modify()));
            connect(yAxisLog,SIGNAL(stateChanged(int)), this, SLOT(modify()));

            connect(xComboBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));
            connect(yComboBox,SIGNAL(currentIndexChanged(int)), this, SLOT(modify()));

            chartScene->AddItem(chart);
            chartActor->SetScene(chartScene);
            renderer->AddActor(chartActor);
            chartScene->SetRenderer(renderer);
            
            // Create a table with some points in it.
//            vtkNew<vtkNamedColors> colors;


//            vtkNew<vtkFloatArray> arrX;
//            arrX->SetName("X Axis");
//            table->AddColumn(arrX);
//
//              vtkNew<vtkFloatArray> arrC;
//              arrC->SetName("Cosine");
//              table->AddColumn(arrC);
//
//              vtkNew<vtkFloatArray> arrS;
//              arrS->SetName("Sine");
//              table->AddColumn(arrS);
//
//              vtkNew<vtkFloatArray> arrT;
//              arrT->SetName("Tan");
//              table->AddColumn(arrT);
//
//              // Test charting with a few more points...
//              int numPoints = 69;
//              float inc = 7.5 / (numPoints - 1.0);
//              table->SetNumberOfRows(numPoints);
//              for (int i = 0; i < numPoints; ++i)
//              {
//                table->SetValue(i, 0, i * inc);
//                table->SetValue(i, 1, cos(i * inc) + 0.0);
//                table->SetValue(i, 2, sin(i * inc) + 0.0);
//                table->SetValue(i, 3, tan(i * inc) + 0.5);
//              }

              // Add multiple line plots, setting the colors etc
            if(table->GetNumberOfRows()>0 && table->GetNumberOfColumns()>0)
            {
                vtkPlot* points = chart->AddPlot(vtkChart::LINE);
                points->SetInputData(table, 0, 0);
              vtkColor3d color3d = colors->GetColor3d("black");
                points->SetColor(color3d.GetRed(), color3d.GetGreen(), color3d.GetBlue());
                points->SetWidth(1.0);
            }
//              dynamic_cast<vtkPlotPoints*>(points)->SetMarkerStyle(vtkPlotPoints::CROSS);
//              points = chart->AddPlot(vtkChart::POINTS);
//              points->SetInputData(table, 0, 2);
//              points->SetColor(color3d.GetRed(), color3d.GetGreen(), color3d.GetBlue());
//              points->SetWidth(1.0);
//              dynamic_cast<vtkPlotPoints*>(points)->SetMarkerStyle(vtkPlotPoints::PLUS);
//              points = chart->AddPlot(vtkChart::POINTS);
//              points->SetInputData(table, 0, 3);
//              points->SetColor(color3d.GetRed(), color3d.GetGreen(), color3d.GetBlue());
//              points->SetWidth(1.0);
            
        }
        

        
 
        
        void ChartActor::modify()
        {
            xComboBox->setEnabled(showChart->isChecked());
            yComboBox->setEnabled(showChart->isChecked());
            xAxisLog->setEnabled(showChart->isChecked());
            yAxisLog->setEnabled(showChart->isChecked());

            chart->GetAxis(0)->SetTitle(yComboBox->currentText().toStdString());
            chart->GetAxis(1)->SetTitle(xComboBox->currentText().toStdString());

            chart->GetAxis(0)->SetLogScale(yAxisLog->isChecked());
            chart->GetAxis(1)->SetLogScale(xAxisLog->isChecked());

            chart->GetAxis(0)->SetGridVisible(false);
            chart->GetAxis(1)->SetGridVisible(false);
            
            chart->ClearPlots();

            if(xComboBox->count() && yComboBox->count())
            {
                if(table->GetNumberOfRows()>xComboBox->currentIndex() && table->GetNumberOfColumns()>yComboBox->currentIndex())
                {
                    vtkPlot* points = chart->AddPlot(vtkChart::LINE);
                    points->SetInputData(table, xComboBox->currentIndex(), yComboBox->currentIndex());
                                        
                    vtkColor3d color3d = colors->GetColor3d("magenta");
                    points->SetColor(color3d.GetRed(), color3d.GetGreen(), color3d.GetBlue());
                    points->SetWidth(3.0);
                    chart->RecalculateBounds();
                    chart->Modified();

                }
            }
            
            chartActor->SetVisibility(showChart->isChecked());
            renderWindow->Render();
        }



} // namespace model
#endif







