/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PLANEPLOTTER_H_
#define model_PLANEPLOTTER_H_


#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <float.h>
#include <Eigen/Core>

#include <model/Network/Readers/EdgeReader.h>

namespace model {
	
	
	class PlanePlotter : public EdgeReader<'G',8,float> {
		
	public:
		
		bool showGlidePlanes;
		
		/* Constructor ******************************************/
		PlanePlotter() : showGlidePlanes(false) {}
		
		/* plot *************************************************/
		void plot() const {
			if(showGlidePlanes){
				glDisable(GL_DEPTH_TEST);
				glEnable (GL_BLEND);
				glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				
				glEnable(GL_ALPHA_TEST);
				glAlphaFunc(GL_GREATER, 0.0f);
				
				glEnable(GL_COLOR_MATERIAL);
				glColor4f(1.0f, 0.0f, 0.0f, 0.1);
				for (model::EdgeReader<'G',8,float>::const_iterator iterG=this->begin(); iterG!=this->end();++iterG){
					glBegin(GL_LINES); 
					glVertex3f(iterG->second.operator()(0), iterG->second.operator()(1),iterG->second.operator()(2)); 
					glVertex3f(iterG->second.operator()(3), iterG->second.operator()(4),iterG->second.operator()(5)); 
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





