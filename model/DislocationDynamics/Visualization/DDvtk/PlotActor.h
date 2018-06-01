/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlotActor_H_
#define model_PlotActor_H_

#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkRenderer.h>
#include <vtkContextActor.h>
#include <vtkFloatArray.h>
#include <vtkPlot.h>
#include <vtkAxis.h>

#include <model/IO/IDreader.h>


namespace model
{
    
    struct PlotActor : public IDreader<'F',1,200,double>
    {
        vtkRenderer* const renderer;
        vtkSmartPointer<vtkTable> table;
        vtkSmartPointer<vtkChartXY> chart;
//        vtkSmartPointer<vtkTable> table2;
//        vtkSmartPointer<vtkChartXY> chart2;
        vtkSmartPointer<vtkContextActor> actor;
        
        int xCol;
        int yCol;

        
        /**********************************************************************/
        PlotActor(vtkRenderer* const ren,
                  const int& xCol,
                  const int& yCol,
                  const int& frameID,
                  const std::map<int,std::string>& FlabelsMap) :
        /* init */ renderer(ren),
        /* init */ table(vtkSmartPointer<vtkTable>::New()),
        chart(vtkSmartPointer<vtkChartXY>::New()),
//        /* init */ table2(vtkSmartPointer<vtkTable>::New()),
//        chart2(vtkSmartPointer<vtkChartXY>::New()),
        actor(vtkSmartPointer<vtkContextActor>::New())
        {
            
            //    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
            
            vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
            arrX->SetName("X");
            table->AddColumn(arrX);
            
            vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
            arrY->SetName("Y");
            table->AddColumn(arrY);
            
            vtkSmartPointer<vtkFloatArray> arrXc = vtkSmartPointer<vtkFloatArray>::New();
                        arrXc->SetName("Xc");
            table->AddColumn(arrXc);
            
            vtkSmartPointer<vtkFloatArray> arrYc = vtkSmartPointer<vtkFloatArray>::New();
                        arrYc->SetName("Yc");
            table->AddColumn(arrYc);
            
            //    vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
            //    arrS->SetName("Sine");
            //    table->AddColumn(arrS);
            
            
            //            model::IDreader<'F',1,200,double> vReader;
            //            Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
            
            
            if (this->isGood(0,true))
            {
                this->read(0,true);
                
//                int rowsC=0;
//                for (const auto& row : *this)
//                {
//                    if(row.first<=frameID)
//                    {
//                        rowsC++;
//                    }
//                }
                
                table->SetNumberOfRows(this->size());
//                table2->SetNumberOfRows(rowsC);

                int i=0;
                double x,y,xc,yc;
                for (const auto& row : *this)
                {
//                    if(row.first<=frameID)
//                    {
//                        rowsC++;
//                    }
                    
                    x=(xCol==0)? row.first : row.second[xCol-1];
                    y=(yCol==0)? row.first : row.second[yCol-1];
                    if(row.first<=frameID)
                    {
                        xc=x;
                        yc=y;
                    }
                    table->SetValue(i, 0, x);
                    table->SetValue(i, 1, y);
                    table->SetValue(i, 2, xc);
                    table->SetValue(i, 3, yc);
                    
//                    if(xCol==0)
//                    {
//                        table->SetValue(i, 0, row.first);
//                        if(row.first<=frameID)
//                        {
//                            table->SetValue(i,2, row.first);
//                        }
//                        else
//                        {
//                            table->SetValue(i,2, nan);
//                        }
//                    }
//                    else
//                    {
//                        table->SetValue(i, 0, row.second[xCol-1]);
//                        if(row.first<=frameID)
//                        {
//                            table->SetValue(i,2, row.second[xCol-1]);
//                        }
//                        else
//                        {
//                            table->SetValue(i,2, nan);
//                        }
//                    }
//                    
//                    if(yCol==0)
//                    {
//                        table->SetValue(i, 1, row.first);
//                        if(row.first<=frameID)
//                        {
//                            table->SetValue(i,3, row.first);
//                        }
//                        else
//                        {
//                            table->SetValue(i,3, nan);
//                        }
//
//
//                    }
//                    else
//                    {
//                        table->SetValue(i, 1, row.second[yCol-1]);
//                        if(row.first<=frameID)
//                        {
//                            table->SetValue(i,3, row.second[yCol-1]);
//                        }
//                        else
//                        {
//                            table->SetValue(i,3, nan);
//                        }
//
//
//
//                    }
//                    table->SetValue(i, 0, row.second[xCcol]);
//                    table->SetValue(i, 0, row.first);
//                    table->SetValue(i, 1, row.second[yCol]);
                    
                    i++;
                }
                
                

                
                vtkPlot *line = chart->AddPlot(vtkChart::LINE);
                line->SetInputData(table, 0, 1);
                line->SetColor(0, 0, 0, 255);
                line->SetWidth(1.0);
                
                line = chart->AddPlot(vtkChart::LINE);
                line->SetInputData(table, 2, 3);
                line->SetColor(255, 0, 255, 255);
                line->SetWidth(2.0);
        
                const auto iterX=FlabelsMap.find(xCol);
                
                if (iterX!=FlabelsMap.end())
                {
                    chart->GetAxis(1)->SetTitle(FlabelsMap.at(xCol));
                }
                else
                {
                    std::cout<<"column "<<xCol<<" not found"<<std::endl;
                }
                
                const auto iterY=FlabelsMap.find(yCol);
                if (iterY!=FlabelsMap.end())
                {
                    chart->GetAxis(0)->SetTitle(FlabelsMap.at(yCol));
                }
                else
                {
                    std::cout<<"column "<<yCol<<" not found"<<std::endl;
                }
        
                
                chart->GetAxis(1)->GetTitleProperties()->SetFontSize(20);
                chart->GetAxis(0)->GetTitleProperties()->SetFontSize(20);

                
//
//                //            VTK_CREATE(vtkContextActor, actor);
//                //            VTK_CREATE(APIDiagram, diagram);
                actor->GetScene()->AddItem(chart);
//                actor->GetScene()->AddItem(chart2);
                renderer->AddActor(actor);
            }
            else
            {
                model::cout<<"could not read runID from F/F_0.txt"<<std::endl;
                //        runID=0;
            }
            
            
            
            // Add multiple line plots, setting the colors etc
            //vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
            

            
            
            
        }
        
        /**********************************************************************/
        ~PlotActor()
        {
            renderer->RemoveActor(actor);
        }
    };
    
    
} // namespace model
#endif







