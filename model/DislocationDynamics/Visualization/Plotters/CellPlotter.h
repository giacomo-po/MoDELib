/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CELLPLOTTER_H_
#define model_CELLPLOTTER_H_


#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <float.h>
#include <Eigen/Core>

#include <model/Network/Readers/VertexReader.h>

namespace model {
	
	
	class CellPlotter : public VertexReader<'C',14,double> {
		
	public:
		
		bool showCells;
		
		/* Constructor ***********************************************/
		CellPlotter() : showCells(false){}
		
		/* plot ******************************************************/
		void plot() const {
			if(showCells){
				glDisable(GL_DEPTH_TEST);
				glEnable (GL_BLEND);
				glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				
				glEnable(GL_ALPHA_TEST);
				glAlphaFunc(GL_GREATER, 0.0f);
				
				glEnable(GL_COLOR_MATERIAL);
				
				glColor4f(1.0f, 0.0f, 1.0f, 0.1);
				
				for (model::VertexReader<'C',14,double>::const_iterator cellIter=this->begin(); cellIter!=this->end(); ++cellIter){
					Eigen::Matrix<int,3,1> cellID(cellIter->second.segment<3>(0).cast<int>());
					float cellSize(cellIter->second(3));
					
					glBegin(GL_LINE_LOOP); 
					glVertex3f((cellID(0)-0.5)*cellSize,(cellID(1)-0.5)*cellSize,(cellID(2)-0.5)*cellSize);
					glVertex3f((cellID(0)-0.5)*cellSize,(cellID(1)-0.5)*cellSize,(cellID(2)-0.5)*cellSize); //increase x
					glVertex3f((cellID(0)+0.5)*cellSize,(cellID(1)+0.5)*cellSize,(cellID(2)-0.5)*cellSize); //increase y 
					glVertex3f((cellID(0)-0.5)*cellSize,(cellID(1)+0.5)*cellSize,(cellID(2)-0.5)*cellSize); //decrease x 
					glEnd();
					
					glBegin(GL_LINE_LOOP); 
					glVertex3f((cellID(0)-0.5)*cellSize,(cellID(1)-0.5)*cellSize,(cellID(2)+0.5)*cellSize); 
					glVertex3f((cellID(0)+0.5)*cellSize,(cellID(1)-0.5)*cellSize,(cellID(2)+0.5)*cellSize); //increase x 
					glVertex3f((cellID(0)+0.5)*cellSize,(cellID(1)+0.5)*cellSize,(cellID(2)+0.5)*cellSize); //increase y 
					glVertex3f((cellID(0)-0.5)*cellSize,(cellID(1)+0.5)*cellSize,(cellID(2)+0.5)*cellSize); //decrease x 
					glEnd();
					
					
					glBegin(GL_LINES); 
					glVertex3f((cellID(0)-0.5)*cellSize,(cellID(1)-0.5)*cellSize,(cellID(2)-0.5)*cellSize); 
					glVertex3f((cellID(0)-0.5)*cellSize,(cellID(1)-0.5)*cellSize,(cellID(2)+0.5)*cellSize); //increase x 
					
					glVertex3f((cellID(0)+0.5)*cellSize,(cellID(1)-0.5)*cellSize,(cellID(2)-0.5)*cellSize); 
					glVertex3f((cellID(0)+0.5)*cellSize,(cellID(1)-0.5)*cellSize,(cellID(2)+0.5)*cellSize); //increase x 
					
					glVertex3f((cellID(0)+0.5)*cellSize,(cellID(1)+0.5)*cellSize,(cellID(2)-0.5)*cellSize); 
					glVertex3f((cellID(0)+0.5)*cellSize,(cellID(1)+0.5)*cellSize,(cellID(2)+0.5)*cellSize); //increase x 
					
					glVertex3f((cellID(0)-0.5)*cellSize,(cellID(1)+0.5)*cellSize,(cellID(2)-0.5)*cellSize); 
					glVertex3f((cellID(0)-0.5)*cellSize,(cellID(1)+0.5)*cellSize,(cellID(2)+0.5)*cellSize); //increase x
					
					
					glEnd();
					
				}
				glEnable(GL_DEPTH_TEST); //Makes 3D drawing work when something is in front of something else
			}
		}
		
		
	};
	
}
#endif
/*********************************************************************/
/*********************************************************************/





