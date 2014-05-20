/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MESHPLOTTER_H_
#define model_MESHPLOTTER_H_


//#ifdef __APPLE__
//#include <OpenGL/OpenGL.h>
//#include <GLUT/glut.h> // for quad
//#else
//#include <GL/glut.h> // for quad
//#endif

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#include <GL/glx.h>
#endif

#include <float.h>
#include <vector>

#include <Eigen/Core>

#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>
#include <model/Mesh/SimplicialMesh.h>


namespace model
{
	
	class MeshPlotter :
    //    /*               */ public VertexReader<'N',5,float>,
    //	/*               */ public EdgeReader  <'T',3,float>,
    //    /*               */ public VertexReader<'Q',4,float> // quadrature points file
	/*               */ public VertexReader<'D',4,float>
    {
		
        //		typedef VertexReader<'N',5,float> NodeContainerType;
		typedef VertexReader<'D',4,float> DispContainerType;
        //		typedef VertexReader<'Q',4,float> QuadContainerType;
        //		typedef   EdgeReader<'T',3,float> EdgeContainerType;
		
        //		bool edgeFileIsGood;
        //		bool nodeFileIsGood;
		bool dispFileIsGood;
        //		bool quadFileIsGood;
		
		typedef std::vector<Eigen::Matrix<float,3,4>,Eigen::aligned_allocator<Eigen::Matrix<float,3,4> > > EdgeVectoType;
		EdgeVectoType edgeVector; // this is [P0 P1 D0 D1]
        
        //        SimplicialMesh<3> mesh;
        
        const SimplicialMesh<3>* const p_mesh;
		
	public:
		enum {showMeshStates=3};
		short unsigned int showMesh;
		
        //        bool showQuad;
        float dispScale;
		
        //        SimplicialMesh<3> mesh;
        
		/* Constructor ******************************************/
		MeshPlotter(const SimplicialMesh<3>* const p_mesh_in) :
        /* init list */ p_mesh(p_mesh_in),
        //        /* init list */ edgeFileIsGood(false),
        //		/* init list */ nodeFileIsGood(false),
		/* init list */ dispFileIsGood(false),
		/* init list */ showMesh(0),
        //        /* init list */ showQuad(false),
		/* init list */ dispScale(1.0f)
        //        /* init list */ mesh(1) // read N/N_0.txt and T/T_0.txt
        //        /* init list */ p_mesh(new SimplicialMesh<3>(1)) // read N/N_0.txt and T/T_0.txt
        {
            
            edgeVector.clear();
            edgeVector.reserve(SimplexObserver<3,1>::size()); // use reserve to speed-up push_back used later
            for (typename SimplexObserver<3,1>::const_iterator sIter=SimplexObserver<3,1>::simplexBegin();
                 /*                                         */ sIter!=SimplexObserver<3,1>::simplexEnd();++sIter)
            {
                if(sIter->second->isBoundarySimplex())
                {
                    Eigen::Matrix<float,3,4> temp(Eigen::Matrix<float,3,4>::Zero());
                    temp.col(0)=sIter->second->child(0).P0.cast<float>();
                    temp.col(1)=sIter->second->child(1).P0.cast<float>();
                    //					if (dispFileIsGood) {
                    //						VertexReader<'D',4,float>::iterator iterD1(DispContainerType::find(iterE->first.first));
                    //						VertexReader<'D',4,float>::iterator iterD2(DispContainerType::find(iterE->first.second));
                    //						assert(iterD1!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                    //						assert(iterD2!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                    //						temp.col(2)=iterD1->second.segment<3>(0);
                    //						temp.col(3)=iterD2->second.segment<3>(0);
                    //					}
                    edgeVector.push_back(temp);
                    
                }
            }
        }
		
		/* read *************************************************/
		void read(const int& frameN){
			// Read edge file T only if it exists, otherwise try to read 0
            //			edgeFileIsGood=EdgeContainerType::isGood(frameN,true);
            //			if (edgeFileIsGood){
            //				EdgeContainerType::read(frameN,true);
            //			}
            //			else{
            //				EdgeContainerType::read(0,true);
            //			}
			
			// Read node file N only if it exists, otherwise try to read 0
            //			nodeFileIsGood=NodeContainerType::isGood(frameN,true);
            //			if (nodeFileIsGood){
            //				NodeContainerType::read(frameN,true);
            //			}
            //			else{
            //				NodeContainerType::read(0,true);
            //			}
			
			dispFileIsGood=DispContainerType::isGood(frameN,true);
			if (dispFileIsGood){
				DispContainerType::read(frameN,true);
				
			}
			else{
				DispContainerType::read(0,true);
			}
            
            //            quadFileIsGood=QuadContainerType::isGood(frameN,true);
            //            if (quadFileIsGood){
            //				QuadContainerType::read(frameN,true);
            //
            //			}
            //			else{
            //				QuadContainerType::read(0,true);
            //			}
			//			assert(NodeContainerType==DispContainerType && "NUMBER OF NODES IN DISPLACEMENT FILE");
			
			
            //			edgeVector.clear();
            //			edgeVector.reserve(EdgeContainerType::size()); // use reserve to speed-up push_back used later
            //			for (EdgeContainerType::const_iterator iterE=EdgeContainerType::begin(); iterE!=EdgeContainerType::end();++iterE)
            //            {
            //				VertexReader<'N',5,float>::const_iterator iterN1 = NodeContainerType::find(iterE->first.first);
            //				VertexReader<'N',5,float>::const_iterator iterN2 = NodeContainerType::find(iterE->first.second);
            //				assert(iterN1!=NodeContainerType::end() && "MESH NODE NOT FOUND IN N FILE");
            //				assert(iterN2!=NodeContainerType::end() && "MESH NODE NOT FOUND IN N FILE");
            //				bool isBoundaryNode1(iterN1->second.operator()(3));
            //				bool isBoundaryNode2(iterN2->second.operator()(3));
            //
            //				if ( (isBoundaryNode1 && isBoundaryNode2) /*|| showInteriorBoundaryMesh*/ )
            //                {
            //					Eigen::Matrix<float,3,4> temp(Eigen::Matrix<float,3,4>::Zero());
            //					temp.col(0)=iterN1->second.segment<3>(0);
            //					temp.col(1)=iterN2->second.segment<3>(0);
            //					if (dispFileIsGood) {
            //						VertexReader<'D',4,float>::iterator iterD1(DispContainerType::find(iterE->first.first));
            //						VertexReader<'D',4,float>::iterator iterD2(DispContainerType::find(iterE->first.second));
            //						assert(iterD1!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
            //						assert(iterD2!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
            //						temp.col(2)=iterD1->second.segment<3>(0);
            //						temp.col(3)=iterD2->second.segment<3>(0);
            //					}
            //					edgeVector.push_back(temp);
            //				}
            //			}
		}
		
		/* plot *************************************************/
		void plot() const
        {
			if (showMesh>0)
            {
				float dispCorr(dispScale*(showMesh>1));
				//glDisable(GL_DEPTH_TEST);
				glEnable (GL_BLEND);
				glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				//glEnable(GL_ALPHA_TEST);
				//glAlphaFunc(GL_GREATER, 0.0f);
				//glEnable(GL_COLOR_MATERIAL);
				//glColor4f(1.0f, 1.0f, 1.0f, 0.3);
                glDisable(GL_COLOR_MATERIAL); // use glMaterialfv(...) to set material colors
                GLfloat materialColor[]={0.0, 0.0, 0.0, 0.1};
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialColor);      // ambient color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialColor);      // diffuse color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialColor);  // specular color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialColor);  // emission color for the material
                
                
				for (EdgeVectoType::const_iterator iterE=edgeVector.begin(); iterE!=edgeVector.end();++iterE)
                {
					glBegin(GL_LINES);
                    glVertex3f(iterE->operator()(0,0)+iterE->operator()(0,2)*dispCorr, iterE->operator()(1,0)+iterE->operator()(1,2)*dispCorr,iterE->operator()(2,0)+iterE->operator()(2,2)*dispCorr);
                    glVertex3f(iterE->operator()(0,1)+iterE->operator()(0,3)*dispCorr, iterE->operator()(1,1)+iterE->operator()(1,3)*dispCorr,iterE->operator()(2,1)+iterE->operator()(2,3)*dispCorr);
					glEnd();
				}
                
                glDisable(GL_BLEND);
				glEnable(GL_DEPTH_TEST); //Makes 3D drawing work when something is in front of something else
                
                //                if (showQuad){
                //                    GLUquadric* myQuad;
                //                    myQuad=gluNewQuadric();
                //                    float radius=10.0;
                //                    glColor4f(0.0f, 0.0f, 0.0f, 0.7);
                //                    for (QuadContainerType::const_iterator iterQ=QuadContainerType::begin(); iterQ!=QuadContainerType::end();++iterQ)
                //                    {
                //                        glTranslatef(  iterQ->second(0),  iterQ->second(1),  iterQ->second(2) );
                //                        gluSphere( myQuad , radius , 10 , 10 );
                //                        glTranslatef( -iterQ->second(0), -iterQ->second(1), -iterQ->second(2) );
                //
                //                    }
                //                    gluDeleteQuadric(myQuad); // free myQuad pointer
                //                }
                
			}
		}
		
	};
    /*********************************************************************/
    /*********************************************************************/	
} // namespace model
#endif

