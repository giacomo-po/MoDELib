/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BITMAPPLOTTER_H_
#define model_BITMAPPLOTTER_H_


#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <Eigen/Core>


namespace model {
	
	class BitmapPlotter  {
		
        
       static void *font;
        
	public:
		
        /* plot ******************************************************/
        template <typename scalarType>
		static void renderString(const Eigen::Matrix<scalarType,3,1>& P, const std::string& string2render) {
  			glRasterPos3f(P(0), P(1), P(2));
			glColor3f(1.0f, 1.0f, 1.0f);
  			for (unsigned int k=0; k<string2render.length(); ++k) {
    			glutBitmapCharacter(font, string2render[k]);
  			}
		}
		
		
	};
	
    // Declare static data mebers
    void* BitmapPlotter::font=GLUT_BITMAP_HELVETICA_18;

    
}
#endif
/*********************************************************************/
/*********************************************************************/





