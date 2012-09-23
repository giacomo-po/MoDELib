/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CAMERA_H_
#define model_CAMERA_H_

#include <math.h>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#include <GL/glx.h>
#endif





typedef struct {
    GLdouble x,y,z;
} recVec;

namespace model {
    
    	
    struct Camera{
        
        static GLdouble aperture; // gCamera aperture
        static GLdouble focalLength; // Focal Length along view direction
        static recVec rotPoint; // Point to rotate about
        
    	static recVec viewPos; // View position
        static recVec viewDir; // View direction vector
        static recVec viewUp; // View up direction
        static GLdouble eyeSep; // Eye separation
        static  GLint screenWidth;
        static GLint screenHeight; // current window/screen height and width ?? why
        static GLfloat gWorldRotation[4];
        static GLfloat gShapeSize;

        /* reset **************************************************************/
        static void reset(){
            
            aperture = 40;
            focalLength = 15;
            
            rotPoint.x = 0.0;
            rotPoint.y =0.0;
            rotPoint.z =0.0;
            
            viewPos.x = 0.0;
            viewPos.y = 0.0;
            viewPos.z = -focalLength;
            
            viewDir.x = -viewPos.x;
            viewDir.y = -viewPos.y;
            viewDir.z = -viewPos.z;
            
            viewUp.x = 0;
            viewUp.y = 1;
            viewUp.z = 0;
        }
        
        

        
    };
    
    
    // static data
    
    GLdouble Camera::aperture=40; // gCamera aperture
    GLdouble Camera::focalLength=15; // Focal Length along view direction
    recVec   Camera::rotPoint={0.0,0.0,0.0}; // Point to rotate about
    recVec   Camera::viewPos={0.0,0.0,-Camera::focalLength}; // View position
    recVec   Camera::viewDir={0.0,0.0, Camera::focalLength}; // View direction vector
    recVec   Camera::viewUp={0.0,1.0,0.0}; // View up direction
    GLdouble Camera::eyeSep; // Eye separation
    GLint    Camera::screenWidth;
    GLint    Camera::screenHeight; // current window/screen height and width ?? why
    GLfloat  Camera::gWorldRotation[4] = {0.0, 0.0, 0.0, 0.0};
    GLfloat  Camera::gShapeSize = 1100.0f;
	
    /**************************************************************************/
} // namespace
#endif

