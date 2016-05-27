/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationRenderer_H_
#define model_DislocationRenderer_H_

#include <vtkRenderer.h>

#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>
#include <model/DislocationDynamics/Visualization/vtk/DislocationSegmentActor.h>



namespace model
{
    struct DislocationRenderer :
    /* inherits from   */ public VertexReader<'V',9,double>,
                        /* inherits from   */ public EdgeReader  <'E',11,double>
    {
        static constexpr int dim=3;
        typedef VertexReader<'V',9,double> VertexContainerType; // CHANGE THIS DOUBLE TO SCALARTYPE
        typedef EdgeReader  <'E',11,double>	EdgeContainerType; // CHANGE THIS DOUBLE TO SCALARTYPE

        
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
        
//        vtkSmartPointer<vtkRenderer> renderer;
        
//        /**********************************************************************/
//        DislocationRenderer(vtkSmartPointer<vtkRenderWindow>& renderWindow) :
//        /* init */ renderer(vtkSmartPointer<vtkRenderer>::New())
//        {
//            renderWindow->AddRenderer(renderer);
//        }
        
        /**********************************************************************/
        DislocationRenderer() : showTubes(false),
        /* init list   */ showVertices(false),
        //		/* init list   */ deformedConfig(false),
        /* init list   */ showPlaneNormal(false),
        /* init list   */ showBurgers(false),
        /* init list   */ showVertexID(false),
        plotBoundarySegments(true),
        /* init list   */ showSpecificVertex(false),
        /* init list   */ specificVertexID(0),
        /* init list   */ showPK(false),
        /* init list   */ PKfactor(1000.0)
        {
            
        }

        
        const VertexReader<'V',9,double>& vertexContainer() const
        {
            return *this;
        }
        
        const EdgeReader<'E',11,double>& edgeContainer() const
        {
            return *this;
        }
        
        /* isGood *************************************************************/
        static bool isGood(const int& frameN, const bool& useTXT)
        {
            return VertexContainerType::isGood(frameN,useTXT) && EdgeContainerType::isGood(frameN,useTXT);
        }
        
        /* read ***************************************************************/
        void read(const int& frameN)
        {
            if (isGood(frameN,false)) // bin format
            {
                VertexContainerType::read(frameN,false);
                EdgeContainerType::read(frameN,false);
            }
            else // txt format
            {
                VertexContainerType::read(frameN,true);
                EdgeContainerType::read(frameN,true);
            }
            
//            if(showQuadParticles) // Show quadrature particles
//            {
//                QuadContainerType::read(frameN,true);
//            }
            
//            PKContainerType::read(frameN,true);
            
//            SingleSplinePlotterVectorType::clear(); // clear the current content of sspVector
//            SingleSplinePlotterVectorType::reserve(EdgeContainerType::size()); // reserve to speed-up push_back
//            for (EdgeContainerType::const_iterator itEdge=EdgeContainerType::begin(); itEdge !=EdgeContainerType::end(); ++itEdge)
            
        }
        
//        void addActors(vtkSmartPointer<vtkRenderer>& renderer) const
        void addActors(vtkRenderer* renderer) const

        {
            for (const auto& edge : edgeContainer())
                
            {
                VertexContainerType::const_iterator itSource(VertexContainerType::find(edge.first.first)); //source
                assert(itSource!=VertexContainerType::end() && "SOURCE VERTEX NOT FOUND IN V-FILE");
                VertexContainerType::const_iterator itSink(VertexContainerType::find(edge.first.second)); //sink
                assert(  itSink!=VertexContainerType::end() &&   "SINK VERTEX NOT FOUND IN V-FILE");
                
                
                Eigen::Matrix<float,dim,6> P0T0P1T1BN;
                
                const int sourceTfactor(edge.second(2*dim));
                const int   sinkTfactor(edge.second(2*dim+1));
                const int   snID(edge.second(2*dim+2));
                const bool sourceOnBoundary(itSource->second(2*dim+1));
                const bool   sinkOnBoundary(  itSink->second(2*dim+1));
                
                if(!(sourceOnBoundary && sinkOnBoundary) || plotBoundarySegments)
                {
                    
                    P0T0P1T1BN.col(0) = itSource->second.segment<dim>(0*dim).transpose().template cast<float>();	// source position
                    P0T0P1T1BN.col(2) =   itSink->second.segment<dim>(0*dim).transpose().template cast<float>();	// sink position
                    P0T0P1T1BN.col(1) = sourceTfactor*(itSource->second.segment<dim>(1*dim).transpose().template cast<float>());	// source tangent
                    P0T0P1T1BN.col(3) =  -sinkTfactor*(  itSink->second.segment<dim>(1*dim).transpose().template cast<float>());	// sink tangent
                    P0T0P1T1BN.col(4) = edge.second.segment<dim>(0*dim).transpose().template cast<float>();		// Burgers vector
                    P0T0P1T1BN.col(5) = edge.second.segment<dim>(1*dim).transpose().template cast<float>();		// plane normal
                    
                    //                    SingleSplinePlotterVectorType::emplace_back(P0T0P1T1BN,snID);
                    
                    //                    SingleSplinePlotterVectorType::push_back(new SingleSplinePlotterType(P0T0P1T1BN,snID));

                    model::DislocationSegmentActor sa(P0T0P1T1BN);

//                    renderer->AddActor(sa.tubeActor());
                    renderer->AddActor(sa.lineActor());

                }
            }
        }
        
    };
    
	
} // namespace model
#endif







