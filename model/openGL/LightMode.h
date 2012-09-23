/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LIGHTMODE_H_
#define model_LIGHTMODE_H_

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#include <GL/glx.h>
#endif

namespace model {
	
    struct LightMode{
    
        static void set(unsigned int mode)
        {
            GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
            GLfloat mat_shininess[] = {90.0};
            
            GLfloat position[4] = {7.0,-7.0,12.0,0.0};
            GLfloat ambient[4]  = {0.2,0.2,0.2,1.0};
            GLfloat diffuse[4]  = {1.0,1.0,1.0,1.0};
            GLfloat specular[4] = {1.0,1.0,1.0,1.0};
            
            glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
            glMaterialfv (GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
            
            glEnable(GL_COLOR_MATERIAL);
            glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
            
            switch (mode) {
                case 0:
                    break;
                case 1:
                    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
                    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_FALSE);
                    break;
                case 2:
                    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_FALSE);
                    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);
                    break;
                case 3:
                    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
                    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_FALSE);
                    break;
                case 4:
                    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
                    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);
                    break;
            }
            
            glLightfv(GL_LIGHT0,GL_POSITION,position);
            glLightfv(GL_LIGHT0,GL_AMBIENT,ambient);
            glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuse);
            glLightfv(GL_LIGHT0,GL_SPECULAR,specular);
            glEnable(GL_LIGHT0);
        }
        
    };
    

	
    /**************************************************************************/
} // namespace
#endif
