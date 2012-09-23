/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_BITMAPPLOTTER_H_
#define model_BITMAPPLOTTER_H_

#include <string>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

//#include <Eigen/Core>


namespace model {
	
	class BitmapPlotter  {
		
        
       static void *font;
        
	public:

        
        
        
        static void drawGLString(const GLfloat& x, const GLfloat& y, char *string)
        {
            /*! Renders 
             * \code
             * 		sprintf (outString, "Camera Position: (%0.1f, %0.1f, %0.1f)", gCamera.viewPos.x, gCamera.viewPos.y, gCamera.viewPos.z);
             * \endcode
             */
          glRasterPos2f(x, y);
          int len = (int) strlen(string);
          for (int i = 0; i < len; i++) {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, string[i]);
          }
        }
        
        
        
        /* renderString2 ******************************************************/
        template <typename scalarType>
		static void renderString(const scalarType& x, const scalarType& y, const std::string& string2render)
        {
  			glRasterPos2f(x, y);
			//glColor3f(1.0f, 1.0f, 1.0f);
  			for (unsigned int k=0; k<string2render.length(); ++k) {
    			glutBitmapCharacter(font, string2render[k]);
  			}
		}
        
        
        /* renderString3 ******************************************************/
        template <typename scalarType>
		static void renderString(const scalarType& x, const scalarType& y,const scalarType& z, const std::string& string2render)
        {
  			glRasterPos3f(x,y,z);
			//glColor3f(1.0f, 1.0f, 1.0f);
  			for (unsigned int k=0; k<string2render.length(); ++k) {
    			glutBitmapCharacter(font, string2render[k]);
  			}
		}
		
		
	};
	
    // Declare static data mebers
    void* BitmapPlotter::font=GLUT_BITMAP_HELVETICA_10;

    
}
#endif
/*********************************************************************/
/*********************************************************************/





