/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_InclusionActor_cpp_
#define model_InclusionActor_cpp_

#include <Eigen/Dense>

#include <QWidget>
#include <QGridLayout>
#include <QCheckBox>


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkLookupTable.h>
#include <vtkPolyhedron.h>
#include <vtkIdList.h>
//#include <vtkXMLUnstructuredGridWriter.h>
#include <InclusionActor.h>


// VTK documentation
// http://vtk.1045678.n5.nabble.com/VTK-slow-to-display-300-vtkSphereSource-in-real-time-td5740730.html
// http://www.vtk.org/Wiki/VTK/Examples/Cxx/Visualization/Scaleglyph3D
// http://hiraku0n.blogspot.com/2012/10/writing-scalar-and-vector-field-in-vtk.html
// http://vtk.1045678.n5.nabble.com/Glyphing-vtkImageData-scalars-3D-as-arrows-td3199837.html

namespace model
{

    const typename InclusionActor::PolyhedronInclusionNodeContainerType& InclusionActor::polyhedronInclusionNodes() const
    {
        return *this;
    }

typename InclusionActor::PolyhedronInclusionNodeContainerType& InclusionActor::polyhedronInclusionNodes()
{
    return *this;
}

    const typename InclusionActor::EshelbyInclusionContainerType& InclusionActor::eshelbyInclusions() const
    {
        return *this;
    }

typename InclusionActor::EshelbyInclusionContainerType& InclusionActor::eshelbyInclusions()
{
    return *this;
}


    void InclusionActor::updateConfiguration(const DDconfigIO<3>& configIO)
    {
        std::cout<<"Updating inclusions..."<<std::flush;
        const auto t0= std::chrono::system_clock::now();
        
        updateInclusions(configIO);
        
        vtkSmartPointer<vtkPoints> points(vtkSmartPointer<vtkPoints>::New());
        vtkSmartPointer<vtkFloatArray> diameters(vtkSmartPointer<vtkFloatArray>::New());
        vtkSmartPointer<vtkFloatArray> colors(vtkSmartPointer<vtkFloatArray>::New());
        
        diameters->SetName("diameters");
        colors->SetName("colors");
        
        for(const auto& inclusion : configIO.sphericalInclusions())
        {
            points->InsertNextPoint(inclusion.C(0), inclusion.C(1), inclusion.C(2));
            diameters->InsertNextValue(2.0*inclusion.a);  // origin of arrow
            colors->InsertNextValue(inclusion.phaseID);
        }
        
        grid->SetPoints(points);
        grid->GetPointData()->AddArray(diameters);
        grid->GetPointData()->SetActiveScalars("diameters"); // to set radius first
        grid->GetPointData()->AddArray(colors);
        
        vtkNew<vtkUnstructuredGrid> ugrid;
        vtkNew<vtkPoints> polyhedronPoints;
        vtkIdType polyhedronPointsIDs[configIO.polyhedronInclusionNodes().size()];
        long long k=0;
        for(const auto& node : configIO.polyhedronInclusionNodes())
        {
            polyhedronPoints->InsertNextPoint(node.P.data());
            polyhedronPointsIDs[k]=k;
            k++;
        }
        ugrid->SetPoints(polyhedronPoints);
        
        std::map<size_t,std::map<size_t,std::vector<size_t>>> faces;
        for(const auto& edge : configIO.polyhedronInclusionEdges())
        {
            const size_t& iID(edge.inclusionID);
            const size_t& fID(edge.faceID);
            const size_t& sourceID(edge.sourceID);
            faces[iID][fID].push_back(sourceID);
        }
        
        for(const auto& pair1 : faces)
        {
            vtkNew<vtkIdList> faces;
            for(const auto& pair2 : pair1.second)
            {
                faces->InsertNextId(pair2.second.size());
                for(const auto& nID : pair2.second)
                {
                    faces->InsertNextId(nID);
                }
            }
            ugrid->InsertNextCell(VTK_POLYHEDRON, configIO.polyhedronInclusionNodes().size(), polyhedronPointsIDs, pair1.second.size(), faces->GetPointer(0));
        }
        
        polyhedronMapper->SetInputData(ugrid);
        std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }

    InclusionActor::InclusionActor(vtkGenericOpenGLRenderWindow* const renWin,vtkRenderer* const renderer,const Polycrystal<3>& poly_in) :
    /* init */ renderWindow(renWin)
    /* init */,mainLayout(new QGridLayout(this))
    /* init */,showInclusions(new QCheckBox(this))
    /* init */,sliderInclusionOpacity(new QSlider(this))
    /* init */,grid(vtkSmartPointer<vtkUnstructuredGrid>::New())
    /* init */,sphereSource(vtkSmartPointer<vtkSphereSource>::New())
    /* init */,glyph3D(vtkSmartPointer<vtkGlyph3D>::New())
    /* init */,mapper(vtkSmartPointer<vtkPolyDataMapper>::New())
    /* init */,actor(vtkSmartPointer<vtkActor>::New())
    /* init */,lookUpColors(vtkSmartPointer<vtkLookupTable>::New())
    /* init */,polyhedronMapper(vtkSmartPointer<vtkDataSetMapper>::New())
    /* init */,polyhedronActor(vtkSmartPointer<vtkActor>::New())
    /* init */,poly(poly_in)
    {
        showInclusions->setChecked(true);
        showInclusions->setText("show inclusions");
        mainLayout->addWidget(showInclusions,0,0,1,1);
        this->setLayout(mainLayout);
        connect(showInclusions,SIGNAL(stateChanged(int)), this, SLOT(modify()));
        
        sliderInclusionOpacity->setMinimum(0);
        sliderInclusionOpacity->setMaximum(10);
        sliderInclusionOpacity->setValue(5);
        sliderInclusionOpacity->setOrientation(Qt::Horizontal);
        connect(sliderInclusionOpacity,SIGNAL(valueChanged(int)), this, SLOT(modify()));
        mainLayout->addWidget(sliderInclusionOpacity,1,0,1,1);

        
        sphereSource->SetPhiResolution(20);
        sphereSource->SetThetaResolution(20);
        glyph3D->SetInputData(grid);
        glyph3D->SetSourceConnection(sphereSource->GetOutputPort());
        mapper->SetInputConnection(glyph3D->GetOutputPort());
        mapper->ScalarVisibilityOff();
        actor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
        actor->GetProperty()->SetOpacity(0.4); // Give some color to the mesh. (1,1,1) is white
        actor->SetMapper(mapper);
        
        polyhedronActor->SetMapper(polyhedronMapper);
        polyhedronActor->GetProperty()->SetColor(0.5,0.5,0.5); // Give some color to the mesh. (1,1,1) is white
        polyhedronActor->GetProperty()->SetOpacity(0.4); // Give some color to the mesh. (1,1,1) is white
                
        renderer->AddActor(actor);
        renderer->AddActor(polyhedronActor);
    }

    void InclusionActor::modify()
    {
        actor->SetVisibility(showInclusions->isChecked());
        actor->GetProperty()->SetOpacity(sliderInclusionOpacity->value()/10.0);

        renderWindow->Render();
    }

void InclusionActor::updateInclusions(const DDconfigIO<3>& configIO)
{
    polyhedronInclusionNodes().clear();
    eshelbyInclusions().clear();
    EshelbyInclusionBase<dim>::force_count(0);
    for(const auto& inclusion : configIO.sphericalInclusions())
    {
//        std::cout<<"Creating spherical inclusion "<<inclusion.inclusionID<<std::endl;
        const std::pair<bool,const Simplex<dim,dim>*> searchPair(poly.mesh.search(inclusion.C));
        if(searchPair.first)
        {
            
            const auto& grain(poly.grain(searchPair.second->region->regionID));
            if(inclusion.phaseID<int(grain.singleCrystal->secondPhases().size()))
            {
            const auto secondPhase(grain.singleCrystal->secondPhases()[inclusion.phaseID]);
            EshelbyInclusionBase<dim>::set_count(inclusion.inclusionID);
            
            
            std::shared_ptr<EshelbyInclusionBase<dim>> iptr(new SphericalInclusion<dim>(inclusion.C,inclusion.a,inclusion.eT,poly.nu,poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase));
            
            eshelbyInclusions().emplace(inclusion.inclusionID,iptr);
            }
            else
            {
                throw std::runtime_error("phaseID does not exist in grain.");
            }
        }
    }
    
    for(const auto& piNode : configIO.polyhedronInclusionNodes())
    {
        polyhedronInclusionNodes().emplace(piNode.nodeID,piNode);
    }
    
    std::map<size_t,std::map<size_t,std::vector<size_t>>> faceMap;
    for(const auto& edge : configIO.polyhedronInclusionEdges())
    {
        const size_t& iID(edge.inclusionID);
        const size_t& fID(edge.faceID);
        const size_t& sourceID(edge.sourceID);
        faceMap[iID][fID].push_back(sourceID);
    }

    
    for(const auto& inclusion : configIO.polyhedronInclusions())
    {
//        std::cout<<"Creating polyhedron inclusion "<<inclusion.inclusionID<<std::endl;

        const auto faceIter(faceMap.find(inclusion.inclusionID));
        if(faceIter!=faceMap.end())
        {
            const auto& faces(faceIter->second);
            std::cout<<"    #faces= "<<faces.size()<<std::endl;
            std::set<const PolyhedronInclusionNodeIO<dim>*> uniquePolyNodes;
            for(const auto& pair : faces)
            {
                for(const auto& nodeID : pair.second)
                {
                    uniquePolyNodes.emplace(&polyhedronInclusionNodes().at(nodeID));
                }
            }
            std::cout<<"    #nodes= "<<uniquePolyNodes.size()<<std::endl;
            if(uniquePolyNodes.size()>=dim+1)
            {
                // Find grain
                std::set<size_t> grainIDs;
                for(const auto& nodePtr : uniquePolyNodes)
                {
                    const std::pair<bool,const Simplex<dim,dim>*> searchPair(poly.mesh.search(nodePtr->P));
                    if(searchPair.first)
                    {
                        grainIDs.insert(searchPair.second->region->regionID);
                    }
                    else
                    {
                        throw std::runtime_error("inclusion node outside mesh");
                    }
                }
                
                // Add inclusion
                if(grainIDs.size()==1)
                {
                    const auto& grain(poly.grain(*grainIDs.begin()));
                    if(inclusion.phaseID<int(grain.singleCrystal->secondPhases().size()))
                    {
                    const auto secondPhase(grain.singleCrystal->secondPhases()[inclusion.phaseID]);
                    EshelbyInclusionBase<dim>::set_count(inclusion.inclusionID);
                    std::shared_ptr<EshelbyInclusionBase<dim>> iptr(new PolyhedronInclusion<dim>( polyhedronInclusionNodes(),faces,inclusion.eT,poly.nu,poly.mu,inclusion.mobilityReduction,inclusion.phaseID,secondPhase));
                    eshelbyInclusions().emplace(inclusion.inclusionID,iptr);
                    }
                    else
                    {
                        throw std::runtime_error("phaseID does not exist in grain.");
                    }
                }
                else
                {
                    throw std::runtime_error("inclusion across grain boundaries");
                }
            }
            else
            {
                throw std::runtime_error("inclusion does not have enough nodes");
            }
        }
        else
        {
            throw std::runtime_error("inclusionID not found in faceMap");
        }
    }
    
}

}
#endif
