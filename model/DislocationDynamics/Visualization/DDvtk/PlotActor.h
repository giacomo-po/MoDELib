/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlotActor_H_
#define model_PlaotActor_H_

#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkRenderer.h>
#include <vtkContextActor.h>
#include <vtkFloatArray.h>
#include <vtkPlot.h>

#include <model/IO/IDreader.h>


namespace model
{
    
    struct PlotActor : public IDreader<'F',1,200,double>
    {
        vtkRenderer* const renderer;
        vtkSmartPointer<vtkTable> table;
        vtkSmartPointer<vtkChartXY> chart;
        vtkSmartPointer<vtkContextActor> actor;
        
        /**********************************************************************/
        PlotActor(vtkRenderer* const ren) :
        /* init */ renderer(ren),
        /* init */ table(vtkSmartPointer<vtkTable>::New()),
        chart(vtkSmartPointer<vtkChartXY>::New()),
        actor(vtkSmartPointer<vtkContextActor>::New())
        {
            
            //    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
            
            vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
            arrX->SetName("disp");
            table->AddColumn(arrX);
            
            vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
            arrC->SetName("stress");
            table->AddColumn(arrC);
            
            //    vtkSmartPointer<vtkFloatArray> arrS = vtkSmartPointer<vtkFloatArray>::New();
            //    arrS->SetName("Sine");
            //    table->AddColumn(arrS);
            
            
            //            model::IDreader<'F',1,200,double> vReader;
            //            Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
            
            
            if (this->isGood(0,true))
            {
                this->read(0,true);
                
                table->SetNumberOfRows(this->size());
                int i=0;
                int xCcol=0;
                int yCol=1;
                
                for (const auto& row : *this)
                {
                    table->SetValue(i, 0, row.second[xCcol]);
                    table->SetValue(i, 1, row.second[yCol]);
                    
                    i++;
                }
            }
            else
            {
                model::cout<<"could not read runID from F/F_0.txt"<<std::endl;
                //        runID=0;
            }
            
            
            
            // Add multiple line plots, setting the colors etc
            //vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
            
            vtkPlot *line = chart->AddPlot(vtkChart::LINE);
            line->SetInputData(table, 0, 1);
            line->SetColor(0, 0, 255, 255);
            line->SetWidth(1.0);
            
            //            VTK_CREATE(vtkContextActor, actor);
            //            VTK_CREATE(APIDiagram, diagram);
            actor->GetScene()->AddItem(chart);
            renderer->AddActor(actor);
            
            
            
        }
    };
    
    
} // namespace model
#endif







