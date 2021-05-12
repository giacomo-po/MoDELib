/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopActor_H_
#define model_DislocationLoopActor_H_

#include <deque>
#include <string>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkMath.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPolyLine.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLabeledDataMapper.h>
#include <vtkFloatArray.h>

#include <IDreader.h>
#include <PlanarPolygon.h>
#include <DDconfigIO.h>
#include <MeshPlane.h>
#include <DDconfigVtkBase.h>

namespace model
{
    struct DislocationLoopActor : public DDconfigVtkBase
    {
        vtkSmartPointer<vtkPolyData> trianglePolyData;
        vtkSmartPointer<vtkPolyDataMapper> triangleMapper;
        vtkSmartPointer<vtkActor> triangleActor;
        
        /**********************************************************************/
        DislocationLoopActor(vtkRenderer* const renderer) :
        /* init */ trianglePolyData(vtkSmartPointer<vtkPolyData>::New()),
        /* init */ triangleMapper(vtkSmartPointer<vtkPolyDataMapper>::New()),
        /* init */ triangleActor(vtkSmartPointer<vtkActor>::New())
        {
            triangleMapper->SetInputData(trianglePolyData);
            triangleActor->SetMapper(triangleMapper);
            triangleActor->GetProperty()->SetColor(0.5,0.0,0.5);

            renderer->AddActor(triangleActor);
            modify();
        }
        
        /**********************************************************************/
        void updateConfiguration(const DDconfigIO<3>& configIO,vtkPolyData* const nodePolyData)
        {
            std::cout<<"Updating loops..."<<std::flush;
            const auto t1= std::chrono::system_clock::now();
            std::map<size_t,std::map<size_t,size_t>> loopEdgeMap;
            
            for(const auto& looplink : configIO.links())
            {
                loopEdgeMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
            }
            assert(loopEdgeMap.size()==configIO.links().size());
            
            vtkSmartPointer<vtkCellArray> triangles(vtkSmartPointer<vtkCellArray>::New());
            size_t loopLumber=1;
            for(const auto& loop : configIO.loops())
            {
                const auto loopFound=loopEdgeMap.find(loop.sID);
                assert(loopFound!=loopEdgeMap.end());
                
                std::vector<size_t> nodeIDs;
                nodeIDs.push_back(loopFound->second.begin()->first);
                std::deque<Eigen::Vector3d> nodePositions;
                auto vIter0(configIO.nodeMap().find(loopFound->second.begin()->first)); //source
                assert(vIter0!=configIO.nodeMap().end() && "VERTEX NOT FOUND IN V-FILE");
                nodePositions.push_back(vIter0->second->P);
                for(size_t k=0;k<loopFound->second.size();++k)
                {
                    const auto nodeFound=loopFound->second.find(*nodeIDs.rbegin());
                    if(k<loopFound->second.size()-1)
                    {
                        nodeIDs.push_back(nodeFound->second);
                        
                        auto vIter(configIO.nodeMap().find(nodeFound->second)); //source
                        assert(vIter!=configIO.nodeMap().end() && "VERTEX NOT FOUND IN V-FILE");
                        //Eigen::Map<const Eigen::Matrix<double,1,10>> vertexRow(vIter->second.data());
                        nodePositions.push_back(vIter->second->P);
                    }
                    else
                    {
                        assert(nodeFound->second==nodeIDs[0]);
                    }
                }
                
                if(fabs(loop.B.dot(loop.N))<FLT_EPSILON && loop.N.squaredNorm()>FLT_EPSILON)
                {
                    PlanarPolygon pp(loop.P,loop.N);
                    
                    pp.assignPoints(nodePositions);
                    std::deque<std::array<size_t, 3>> tri=pp.triangulate();
                    
                    for(const auto& triID : tri)
                    {
                        const size_t& nodeID0(nodeIDs[std::get<0>(triID)]);
                        const size_t& nodeID1(nodeIDs[std::get<1>(triID)]);
                        const size_t& nodeID2(nodeIDs[std::get<2>(triID)]);
                        
                        const size_t ptID0=std::distance(configIO.nodeMap().begin(),configIO.nodeMap().find(nodeID0));
                        const size_t ptID1=std::distance(configIO.nodeMap().begin(),configIO.nodeMap().find(nodeID1));
                        const size_t ptID2=std::distance(configIO.nodeMap().begin(),configIO.nodeMap().find(nodeID2));
                        
                        vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                        triangle->GetPointIds()->SetId ( 0, ptID0 );
                        triangle->GetPointIds()->SetId ( 1, ptID1 );
                        triangle->GetPointIds()->SetId ( 2, ptID2 );
                        triangles->InsertNextCell ( triangle );
                    }
                }
            }
            
            trianglePolyData->SetPoints (nodePolyData->GetPoints());
            trianglePolyData->SetPolys ( triangles );
            trianglePolyData->Modified();
            model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;
        }
        
        /**********************************************************************/
        void modify()
        {
            
            triangleActor->GetProperty()->SetOpacity(this->slippedAreaOpacity);
            
            if(showSlippedArea)
            {
                triangleActor->VisibilityOn();
            }
            else
            {
                triangleActor->VisibilityOff();
            }
        }
        
    };
    
} // namespace
#endif
