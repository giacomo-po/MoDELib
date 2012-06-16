/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MESHPLOTTER_H_
#define model_MESHPLOTTER_H_


#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <float.h>
#include <vector>

#include <Eigen/Core>

#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>

namespace model {
	
	class MeshPlotter : public  VertexReader<'N',5,float>,
	/*               */ public VertexReader<'D',4,float>,
	/*               */ public EdgeReader  <'T',3,float>{
		
		typedef VertexReader<'N',5,float> NodeContainerType;
		typedef VertexReader<'D',4,float> DispContainerType;
		typedef EdgeReader<'T',3,float>   EdgeContainerType;
		
		bool edgeFileIsGood;
		bool nodeFileIsGood;
		bool dispFileIsGood;
		
		typedef std::vector<Eigen::Matrix<float,3,4>,Eigen::aligned_allocator<Eigen::Matrix<float,3,4> > > EdgeVectoType;
		EdgeVectoType edgeVector; // this is [P0 P1 D0 D1]
		
	public:
		enum {showMeshStates=3};
		short unsigned int showMesh;
		float dispScale;
		
		/* Constructor ******************************************/		
		MeshPlotter() : edgeFileIsGood(false),
		/* init list */ nodeFileIsGood(false),
		/* init list */ dispFileIsGood(false),
		/* init list */ showMesh(0),
		/* init list */ dispScale(1.0f){}
		
		/* read *************************************************/		
		void read(const int& frameN){
			// Read edge file T only if it exists, otherwise try to read 0
			edgeFileIsGood=EdgeContainerType::isGood<true>(frameN);
			if (edgeFileIsGood){
				EdgeContainerType::read<true>(frameN);
			}
			else{
				EdgeContainerType::read<true>(0);					
			}
			
			// Read node file N only if it exists, otherwise try to read 0
			nodeFileIsGood=NodeContainerType::isGood(frameN);
			if (nodeFileIsGood){
				NodeContainerType::read(frameN);
			}
			else{
				NodeContainerType::read(0);		
			}
			
			dispFileIsGood=DispContainerType::isGood(frameN);
			if (dispFileIsGood){
				DispContainerType::read(frameN);
				
			}
			else{
				DispContainerType::read(0);
			}
			
			//			assert(NodeContainerType==DispContainerType && "NUMBER OF NODES IN DISPLACEMENT FILE");
			
			
			edgeVector.clear();
			edgeVector.reserve(EdgeContainerType::size()); // use reserve to speed-up push_back used later
			for (EdgeContainerType::const_iterator iterE=EdgeContainerType::begin(); iterE!=EdgeContainerType::end();++iterE){
				VertexReader<'N',5,float>::const_iterator iterN1 = NodeContainerType::find(iterE->first.first);
				VertexReader<'N',5,float>::const_iterator iterN2 = NodeContainerType::find(iterE->first.second);
				assert(iterN1!=NodeContainerType::end() && "MESH NODE NOT FOUND IN N FILE");
				assert(iterN2!=NodeContainerType::end() && "MESH NODE NOT FOUND IN N FILE");
				bool isBoundaryNode1(iterN1->second.operator()(3));
				bool isBoundaryNode2(iterN2->second.operator()(3));
				
				if ( (isBoundaryNode1 && isBoundaryNode2) /*|| showInteriorBoundaryMesh*/ ) {
					Eigen::Matrix<float,3,4> temp(Eigen::Matrix<float,3,4>::Zero());
					temp.col(0)=iterN1->second.segment<3>(0);
					temp.col(1)=iterN2->second.segment<3>(0);
					if (dispFileIsGood) {
						VertexReader<'D',4,float>::iterator iterD1(DispContainerType::find(iterE->first.first));
						VertexReader<'D',4,float>::iterator iterD2(DispContainerType::find(iterE->first.second));
						assert(iterD1!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
						assert(iterD2!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
						temp.col(2)=iterD1->second.segment<3>(0);
						temp.col(3)=iterD2->second.segment<3>(0);
					}
					edgeVector.push_back(temp);
				}
			}
		}
		
		/* plot *************************************************/
		void plot() const {
			if (showMesh>0){
				float dispCorr(dispScale*(showMesh>1));
				glDisable(GL_DEPTH_TEST);
				glEnable (GL_BLEND);
				glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				glEnable(GL_ALPHA_TEST);
				glAlphaFunc(GL_GREATER, 0.0f);
				glEnable(GL_COLOR_MATERIAL);
				glColor4f(0.0f, 0.0f, 0.0f, 0.1);
				for (EdgeVectoType::const_iterator iterE=edgeVector.begin(); iterE!=edgeVector.end();++iterE){
					glBegin(GL_LINES); 
				glVertex3f(iterE->operator()(0,0)+iterE->operator()(0,2)*dispCorr, iterE->operator()(1,0)+iterE->operator()(1,2)*dispCorr,iterE->operator()(2,0)+iterE->operator()(2,2)*dispCorr); 
				glVertex3f(iterE->operator()(0,1)+iterE->operator()(0,3)*dispCorr, iterE->operator()(1,1)+iterE->operator()(1,3)*dispCorr,iterE->operator()(2,1)+iterE->operator()(2,3)*dispCorr); 
					glEnd();
				}
				glEnable(GL_DEPTH_TEST); //Makes 3D drawing work when something is in front of something else
			}
		}
		
	};
/*********************************************************************/
/*********************************************************************/	
} // namespace model
#endif

