/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationActors_H_
#define model_DislocationActors_H_


#include <memory>

#include <vtkRenderer.h>

#include <model/IO/VertexReader.h>
#include <model/IO/EdgeReader.h>
#include <model/IO/IDreader.h>
#include <model/DislocationDynamics/Visualization/vtk/DislocationSegmentActor.h>
#include <model/DislocationDynamics/Visualization/vtk/DislocationNodeActor.h>
#include <model/DislocationDynamics/Visualization/vtk/PKActor.h>



namespace model
{
    struct DislocationActors :
    /* inherits from   */ public VertexReader<'V',10,double>,
    /* inherits from   */ public EdgeReader  <'K',15,double>,
    /* inherits from   */ public IDreader<'P',3,6,double>,
    //    /* inherits from   */ std::deque<DislocationSegmentActor>,
    /* inherits from   */ std::deque<DislocationNodeActor>,
    /* inherits from   */ std::deque<PKActor>
    {
        static constexpr int dim=3;
        typedef VertexReader<'V',10,double> VertexContainerType; // CHANGE THIS DOUBLE TO SCALARTYPE
        typedef EdgeReader  <'K',15,double>	EdgeContainerType; // CHANGE THIS DOUBLE TO SCALARTYPE
        typedef IDreader<'P',3,6,double> PKContainerType;
        
        
        //        long int currentFrameID;
        bool showTubes;
        bool showVertices;
        //		bool deformedConfig;
        bool showPlaneNormal;
        bool showBurgers;
        bool showVertexID;
        bool plotBoundarySegments;
        
        //        static int colorScheme;
        //        static bool plotBoundarySegments;
        //        static bool showQuadParticles;
        
        bool showSpecificVertex;
        int specificVertexID;
        bool showPK;
        double PKfactor;
        
        std::unique_ptr<DislocationSegmentActor> ddSegments;
        
        //        vtkSmartPointer<vtkRenderer> renderer;
        
        //        /**********************************************************************/
        //        DislocationActors(vtkSmartPointer<vtkRenderWindow>& renderWindow) :
        //        /* init */ renderer(vtkSmartPointer<vtkRenderer>::New())
        //        {
        //            renderWindow->AddRenderer(renderer);
        //        }
        
        /**********************************************************************/
        DislocationActors() :
        //        /* init list   */ currentFrameID(-1),
        /* init list   */ showTubes(false),
        /* init list   */ showVertices(false),
        //		/* init list   */ deformedConfig(false),
        /* init list   */ showPlaneNormal(false),
        /* init list   */ showBurgers(false),
        /* init list   */ showVertexID(false),
        /* init list   */ plotBoundarySegments(true),
        /* init list   */ showSpecificVertex(false),
        /* init list   */ specificVertexID(0),
        /* init list   */ showPK(false),
        /* init list   */ PKfactor(1000.0)
        {
            
        }
        
        /**********************************************************************/
        void clear()
        {
            //            segmentActors().clear();
            nodeActors().clear();
            pkActors().clear();
        }
        
        /**********************************************************************/
        std::deque<PKActor>& pkActors()
        {
            return *this;
        }
        
        /**********************************************************************/
        std::deque<DislocationNodeActor>& nodeActors()
        {
            return *this;
        }
        
        //        std::deque<DislocationSegmentActor>& segmentActors()
        //        {
        //            return *this;
        //        }
        
        /**********************************************************************/
        const VertexContainerType& vertexContainer() const
        {
            return *this;
        }
        
        VertexContainerType& vertexContainer()
        {
            return *this;
        }
        
        
        /**********************************************************************/
        const EdgeContainerType& edgeContainer() const
        {
            return *this;
        }
        
        EdgeContainerType& edgeContainer()
        {
            return *this;
        }
        
        /**********************************************************************/
        const PKContainerType& pkContainer() const
        {
            return *this;
        }
        
        PKContainerType& pkContainer()
        {
            return *this;
        }
        
        
        /*************************************************************************/
        bool isGood(const int& frameN, const bool& useTXT) const
        {
            return vertexContainer().isGood(frameN,useTXT) && edgeContainer().isGood(frameN,useTXT);
        }
        
        /*************************************************************************/
        void read(const int& frameN)
        {
            std::cout<<"Reading frame "<<frameN<<std::endl;
            if (isGood(frameN,false)) // bin format
            {
                vertexContainer().read(frameN,false);
                edgeContainer().read(frameN,false);
            }
            else // txt format
            {
                vertexContainer().read(frameN,true);
                edgeContainer().read(frameN,true);
            }
            
            //            if(showQuadParticles) // Show quadrature particles
            //            {
            //                QuadContainerType::read(frameN,true);
            //            }
            
            PKContainerType::read(frameN,true);
            
            // Create DislocationSegmentActor
            ddSegments.reset(new DislocationSegmentActor(vertexContainer(),edgeContainer()));
            
            for(const auto& node : vertexContainer())
            {
                nodeActors().emplace_back(node.second.segment<3>(0).template cast<float>());
            }
            
            if (showPK) // Show PK force
            {
                for(const auto& pk : pkContainer())
                {
                    Eigen::Map<const Eigen::Matrix<double,1,6>> val(pk.second.data());
                    pkActors().emplace_back(val.segment<3>(0).template cast<float>(),val.segment<3>(3).template cast<float>());
                }
            }
            
            std::cout<<"finished reading. "<<frameN<<std::endl;
        }
        
        /*************************************************************************/
        void update(const long int& frameID,vtkRenderer* renderer)
        
        {
            
            
            // Remove current actors from renderer
            //renderer->RemoveAllViewProps();
            //            for(auto& segment : segmentActors())
            //            {
            //                renderer->RemoveActor(segment.tubeActor());
            //                //                renderer->AddActor(actor.lineActor());
            //            }
            
            if(ddSegments.get()!=nullptr)
            {
                renderer->RemoveActor(ddSegments->tubeActor);
            }
            
            for(auto& node : nodeActors())
            {
                renderer->RemoveActor(node.actor);
                //                renderer->AddActor(actor.lineActor());
            }
            
            for(auto& pk : pkActors())
            {
                renderer->RemoveActor(pk.actor);
            }
            
            // Clear current actors
            clear();
            
            // Read new frame and store new actors
            read(frameID);
            
            // Add new actors to renderer
            //            for(auto& segment : segmentActors())
            //            {
            //                renderer->AddActor(segment.tubeActor());
            //                //                renderer->AddActor(actor.lineActor());
            //            }
            
            renderer->AddActor(ddSegments->tubeActor);
            
            
            for(auto& node : nodeActors())
            {
                renderer->AddActor(node.actor);
            }
            
            //            for(auto& pk : pkActors())
            //            {
            //                renderer->AddActor(pk.actor);
            //            }
            
            
        }
        
        /*************************************************************************/
        void modify()
        {
            //            ddSegments->modify();
            //            for(auto& segment : segmentActors())
            //            {
            //                segment.modify();
            //            }
            
            for(auto& node : nodeActors())
            {
                node.modify();
            }
        }
        
    };
    
    
} // namespace model
#endif







