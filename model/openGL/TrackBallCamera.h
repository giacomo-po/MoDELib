/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_TRACKBALLCAMERA_H_
#define model_TRACKBALLCAMERA_H_


#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include <math.h>

#include <model/openGL/Camera.h>
#include <model/openGL/TrackBall.h>


#define DTOR 0.0174532925


#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

namespace model {
	
    

    
    class TrackBallCamera {

        /* mouseDolly *********************************************************/
        static void mouseDolly (int x, int y)
        {
            if (gDolly) {
                GLfloat dolly = (gDollyPanStartPoint[1] - y) * -Camera::viewPos.z / 200.0f;
                GLfloat eyeRelative = Camera::eyeSep / Camera::focalLength;
                Camera::focalLength += Camera::focalLength / Camera::viewPos.z * dolly;
                if (Camera::focalLength < 1.0)
                    Camera::focalLength = 1.0;
                Camera::eyeSep = Camera::focalLength * eyeRelative;
                Camera::viewPos.z += dolly;
                if (Camera::viewPos.z == 0.0) // do not let z = 0.0
                    Camera::viewPos.z = 0.0001;
                gDollyPanStartPoint[0] = x;
                gDollyPanStartPoint[1] = y;
                glutPostRedisplay();
            }
        }
        
        /* mousePan ***********************************************************/
        static void mousePan (int x, int y)
        {
            if (gPan) {
                GLfloat panX = (gDollyPanStartPoint[0] - x) / (900.0f / -Camera::viewPos.z);
                GLfloat panY = (gDollyPanStartPoint[1] - y) / (900.0f / -Camera::viewPos.z);
                Camera::viewPos.x -= panX;
                Camera::viewPos.y -= panY;
                gDollyPanStartPoint[0] = x;
                gDollyPanStartPoint[1] = y;
                glutPostRedisplay();
            }
        }
        
        /* mouseTrackball *****************************************************/
        static void mouseTrackball (int x, int y)
        {
            if (gTrackBall) {
                TrackBall::rollToTrackball (x, y, gTrackBallRotation);
                glutPostRedisplay();
            }
        }
        
        
    public:
        
        static GLboolean gDolly;
        static GLint gDollyPanStartPoint[2];
        static GLboolean gPan;
        static GLfloat gTrackBallRotation[4];
        static GLboolean gTrackBall;
       // static GLfloat Camera::gWorldRotation[4];

        
        /* mouse **************************************************************/
        static void mouse (int button, int state, int x, int y)
        {
            if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN)) {
                if (gDolly) { // if we are currently dollying, end dolly
                    mouseDolly (x, y);
                    gDolly = GL_FALSE;
                    glutMotionFunc (NULL);
                    gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
                    glutMotionFunc (NULL);
                } else if (gPan) {
                    mousePan (x, y);
                    gPan = GL_FALSE;
                    glutMotionFunc (NULL);
                    gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
                    glutMotionFunc (NULL);
                }
                TrackBall::start(x, y, 0, 0, Camera::screenWidth, Camera::screenHeight);
                glutMotionFunc (mouseTrackball);
                gTrackBall = GL_TRUE;
            } else if ((button == GLUT_LEFT_BUTTON) && (state == GLUT_UP)) {
                gTrackBall = GL_FALSE;
                glutMotionFunc (NULL);
                TrackBall::rollToTrackball (x, y, gTrackBallRotation);
                if (gTrackBallRotation[0] != 0.0)
                    TrackBall::addToRotation (gTrackBallRotation, Camera::gWorldRotation);
                gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
            }
            else if ((button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN)) {
                if (gTrackBall) {// if we are currently trackballing, end trackball
                    gTrackBall = GL_FALSE;
                    glutMotionFunc (NULL);
                    TrackBall::rollToTrackball (x, y, gTrackBallRotation);
                    if (gTrackBallRotation[0] != 0.0)
                        TrackBall::addToRotation (gTrackBallRotation, Camera::gWorldRotation);
                    gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
                } else if (gPan) {
                    mousePan (x, y);
                    gPan = GL_FALSE;
                    glutMotionFunc (NULL);
                    gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
                    glutMotionFunc (NULL);
                }
                gDollyPanStartPoint[0] = x;
                gDollyPanStartPoint[1] = y;
                glutMotionFunc (mouseDolly);
                gDolly = GL_TRUE;
            } else if ((button == GLUT_MIDDLE_BUTTON) && (state == GLUT_UP)) {
                mouseDolly (x, y);
                gDolly = GL_FALSE;
                glutMotionFunc (NULL);
                gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
                glutMotionFunc (NULL);
            }
            else if ((button == GLUT_RIGHT_BUTTON) && (state == GLUT_DOWN)) {
                if (gTrackBall) {// if we are currently trackballing, end trackball
                    gTrackBall = GL_FALSE;
                    glutMotionFunc (NULL);
                    TrackBall::rollToTrackball (x, y, gTrackBallRotation);
                    if (gTrackBallRotation[0] != 0.0)
                        TrackBall::addToRotation (gTrackBallRotation, Camera::gWorldRotation);
                    gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
                } else if (gDolly) {
                    mouseDolly (x, y);
                    gDolly = GL_FALSE;
                    glutMotionFunc (NULL);
                    gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
                    glutMotionFunc (NULL);
                }
                gDollyPanStartPoint[0] = x;
                gDollyPanStartPoint[1] = y;
                glutMotionFunc (mousePan);
                gPan = GL_TRUE;
            } else if ((button == GLUT_RIGHT_BUTTON) && (state == GLUT_UP)) {
                mousePan (x, y);
                gPan = GL_FALSE;
                glutMotionFunc (NULL);
                gTrackBallRotation [0] = gTrackBallRotation [1] = gTrackBallRotation [2] = gTrackBallRotation [3] = 0.0f;
                glutMotionFunc (NULL);
            }
        }
        
        /* reset **************************************************************/
        static void spaceballmotion (int x, int y, int z)
        {
            long deadZone = 105;
            float scale = -Camera::viewPos.z * 0.00000001f;
            if (abs (x) > deadZone) {
                GLfloat panX = abs (x) * x * scale;
                Camera::viewPos.x += panX;
            }
            if (abs (y) > deadZone) {
                GLfloat panY = abs (y) * y * scale;
                Camera::viewPos.y -= panY;
            }
            if (abs (z) > deadZone) {
                GLfloat dolly = abs (z) * z * scale;
                Camera::viewPos.z += dolly;
                if (Camera::viewPos.z == 0.0) // do not let z = 0.0
                    Camera::viewPos.z = 0.0001;
            }
            glutPostRedisplay();
        }
        
        /* reset **************************************************************/
        static void spaceballrotate (int rx, int ry, int rz)
        {
            long deadZone = 60;
            float rotation[4] = {0.0f, 0.0f, 0.0f, 0.0f};
            // handle rotations about each respective axis
            if (abs (rx) > deadZone) {
                rotation[0] = abs (rx) * -rx * 0.0000008f;
                rotation[1] = 1.0f;
                rotation[2] = 0.0f;
                rotation[3] = 0.0f;
                TrackBall::addToRotation (rotation, Camera::gWorldRotation);
            }
            if (abs (ry) > deadZone) {
                rotation[0] = abs (ry) * ry * 0.0000008f;
                rotation[1] = 0.0f;
                rotation[2] = 1.0f;
                rotation[3] = 0.0f;
                TrackBall::addToRotation (rotation, Camera::gWorldRotation);
            }
            if (abs(rz) > deadZone) {
                rotation[0] = abs (rz) * -rz * 0.0000008f;
                rotation[1] = 0.0f;
                rotation[2] = 0.0f;
                rotation[3] = 1.0f;
                TrackBall::addToRotation (rotation, Camera::gWorldRotation);
            }
            glutPostRedisplay();
        }
        
        
        /* reset **************************************************************/
        static void update(){
            GLdouble xmin, xmax, ymin, ymax;
            // far frustum plane
//            GLdouble zFar = -Camera::viewPos.z + Camera::gShapeSize * 0.5;
            GLdouble zFar =  Camera::gShapeSize * 0.5;

            // near frustum plane clamped at 1.0
            GLdouble zNear = MIN (-Camera::viewPos.z - Camera::gShapeSize * 0.5, 1.0);
            // window aspect ratio
            GLdouble aspect = Camera::screenWidth / (GLdouble)Camera::screenHeight;
            
            glMatrixMode (GL_PROJECTION);
            glLoadIdentity ();
            
            if (aspect > 1.0) {
                ymax = zNear * tan (Camera::aperture * 0.5 * DTOR);
                ymin = -ymax;
                xmin = ymin * aspect;
                xmax = ymax * aspect;
            } else {
                xmax = zNear * tan (Camera::aperture * 0.5 * DTOR);
                xmin = -xmax;
                ymin = xmin / aspect;
                ymax = xmax / aspect;
            }
            glFrustum(xmin, xmax, ymin, ymax, zNear, zFar);
            
            glMatrixMode (GL_MODELVIEW);
            glLoadIdentity ();
            gluLookAt (Camera::viewPos.x, Camera::viewPos.y, Camera::viewPos.z,
                       Camera::viewPos.x + Camera::viewDir.x,
                       Camera::viewPos.y + Camera::viewDir.y,
                       Camera::viewPos.z + Camera::viewDir.z,
                       Camera::viewUp.x, Camera::viewUp.y, Camera::viewUp.z);
			
            // track ball rotation
            glRotatef ( gTrackBallRotation[0],  gTrackBallRotation[1],  gTrackBallRotation[2],  gTrackBallRotation[3]);
            glRotatef (Camera::gWorldRotation[0], Camera::gWorldRotation[1], Camera::gWorldRotation[2], Camera::gWorldRotation[3]);
            
        }
        
        
    };
    
    //static data
    GLboolean TrackBallCamera::gDolly = GL_FALSE;
    GLint     TrackBallCamera::gDollyPanStartPoint[2] = {0, 0};
    GLboolean TrackBallCamera::gPan = GL_FALSE;
    GLfloat   TrackBallCamera::gTrackBallRotation[4] = {0.0, 0.0, 0.0, 0.0};
    GLboolean TrackBallCamera::gTrackBall = GL_FALSE;
	
    /**************************************************************************/
} // namespace
#endif

