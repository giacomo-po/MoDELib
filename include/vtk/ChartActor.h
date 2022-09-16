/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_ChartActor_H_
#define model_ChartActor_H_

#include <random>
#include <algorithm>

#include <QGroupBox>
#include <QGridLayout>
#include <QCheckBox>
#include <QLabel>
#include <QWidget>
#include <QComboBox>
#include <QPushButton>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkChartXY.h>
#include <vtkContextActor.h>
#include <vtkContextScene.h>
#include <vtkCubeSource.h>
#include <vtkFloatArray.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPlotPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTable.h>

#include <TextFileParser.h>
#include <DDtraitsIO.h>

namespace model
{

    struct ChartActor : public QWidget
    {
        
        Q_OBJECT
        
    private slots:
        
        void modify();
        void reload();
        void updatePlots();

        
    public:
        
        const DDtraitsIO& traitsIO;
        QGridLayout* mainLayout;
        QCheckBox* showChart;
        QPushButton* chartUpdate;
        QComboBox* xComboBox;
        QComboBox* yComboBox;
        QCheckBox* xAxisLog;
        QCheckBox* yAxisLog;

        vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow;
        vtkRenderer* const renderer;
        
        vtkSmartPointer<vtkChartXY> chart;
        vtkSmartPointer<vtkContextScene> chartScene;
        vtkSmartPointer<vtkContextActor> chartActor;
//        vtkPlot* points;
        vtkSmartPointer<vtkTable> table;
        vtkSmartPointer<vtkTable> currentTable;
        vtkSmartPointer<vtkNamedColors> colors;
        
        vtkPlot* points;
        vtkPlot* currentPoints;

        ChartActor(const DDtraitsIO& traitsIO_in,vtkSmartPointer<vtkGenericOpenGLRenderWindow>,vtkRenderer* const renderer_in);
        
        void updateConfiguration(const size_t& frameID);

        
    };
    
} // namespace model
#endif







