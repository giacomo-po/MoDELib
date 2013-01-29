/*
 *  glutBasics.cpp
 *  GLUTTest
 *
 *  Created by GGS on June 17 2003.
 *  Copyright (c) 2003 Apple. All rights reserved.
 *
 
	Disclaimer:	IMPORTANT:  This Apple software is supplied to you by Apple Computer, Inc.
				("Apple") in consideration of your agreement to the following terms, and your
				use, installation, modification or redistribution of this Apple software
				constitutes acceptance of these terms.  If you do not agree with these terms,
				please do not use, install, modify or redistribute this Apple software.

				In consideration of your agreement to abide by the following terms, and subject
				to these terms, Apple grants you a personal, non-exclusive license, under Apple’s
				copyrights in this original Apple software (the "Apple Software"), to use,
				reproduce, modify and redistribute the Apple Software, with or without
				modifications, in source and/or binary forms; provided that if you redistribute
				the Apple Software in its entirety and without modifications, you must retain
				this notice and the following text and disclaimers in all such redistributions of
				the Apple Software.  Neither the name, trademarks, service marks or logos of
				Apple Computer, Inc. may be used to endorse or promote products derived from the
				Apple Software without specific prior written permission from Apple.  Except as
				expressly stated in this notice, no other rights or licenses, express or implied,
				are granted by Apple herein, including but not limited to any patent rights that
				may be infringed by your derivative works or by other works in which the Apple
				Software may be incorporated.

				The Apple Software is provided by Apple on an "AS IS" basis.  APPLE MAKES NO
				WARRANTIES, EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION THE IMPLIED
				WARRANTIES OF NON-INFRINGEMENT, MERCHANTABILITY AND FITNESS FOR A PARTICULAR
				PURPOSE, REGARDING THE APPLE SOFTWARE OR ITS USE AND OPERATION ALONE OR IN
				COMBINATION WITH YOUR PRODUCTS.

				IN NO EVENT SHALL APPLE BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL OR
				CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
				GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
				ARISING IN ANY WAY OUT OF THE USE, REPRODUCTION, MODIFICATION AND/OR DISTRIBUTION
				OF THE APPLE SOFTWARE, HOWEVER CAUSED AND WHETHER UNDER THEORY OF CONTRACT, TORT
				(INCLUDING NEGLIGENCE), STRICT LIABILITY OR OTHERWISE, EVEN IF APPLE HAS BEEN
				ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */
 
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
 
#include <GLUT/glut.h>
#include <OpenGL/glext.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>

//#include "trackball.h"
#include "SurfaceGeometry.h"


#include <model/openGL/LightMode.h>
#include <model/openGL/BitmapPlotter.h>
//#include <model/openGL/Camera.h>
#include <model/openGL/TrackBallCamera.h>

#include <model/DislocationDynamics/Visualization/DDreader.h>
#include <model/DislocationDynamics/Visualization/Plotters/SplinePlotter.h>


using namespace model;

DDreader ddr;
float alpha=0.5;
enum{dim=3};
SplinePlotter<dim,50,10,alpha> splinePlotter;


//#define DTOR 0.0174532925
////#define CROSSPROD(p1,p2,p3) \
////   p3.x = p1.y*p2.z - p1.z*p2.y; \
////   p3.y = p1.z*p2.x - p1.x*p2.z; \
////   p3.z = p1.x*p2.y - p1.y*p2.x
//   
//#ifndef MIN
//	#define MIN(a, b) ((a) < (b) ? (a) : (b))
//#endif

enum {
	kTextureSize = 256
};

//typedef struct {
//   GLdouble x,y,z;
//} recVec;

//typedef struct {
//	recVec viewPos; // View position
//	recVec viewDir; // View direction vector
//	recVec viewUp; // View up direction
//	recVec rotPoint; // Point to rotate about
//	GLdouble focalLength; // Focal Length along view direction
//	GLdouble aperture; // gCamera aperture
//	GLdouble eyeSep; // Eye separation
//	GLint screenWidth,screenHeight; // current window/screen height and width
//} recCamera;

//GLfloat gShapeSize = 11.0f;

//GLint gDollyPanStartPoint[2] = {0, 0};
//GLfloat gTrackBallRotation [4] = {0.0, 0.0, 0.0, 0.0};
//GLboolean gDolly = GL_FALSE;
//GLboolean gPan = GL_FALSE;
//GLboolean gTrackBall = GL_FALSE;
//GLfloat gWorldRotation [4] = {155.0, 0.0, -1.0, 0.0};

GLboolean gLines = GL_FALSE;
GLboolean gPolygons = GL_TRUE;

GLboolean gShowHelp = GL_TRUE;
GLboolean gShowInfo = GL_TRUE;

//recCamera gCamera;
recVec gOrigin = {0.0, 0.0, 0.0};

int gLastKey = ' ';

int gMainWindow = 0;

GLuint gPointList = NULL;
GLuint gWireList = NULL;
GLuint gSolidList = NULL;






void drawGLText (GLint window_width, GLint window_height)
{
	char outString [256] = "";
	GLint matrixMode;
	GLint vp[4];
	GLint lineSpacing = 13;
	GLint line = 0;
	GLint startOffest = 7;
	
	glGetIntegerv(GL_VIEWPORT, vp);
	glViewport(0, 0, window_width, window_height);
	
	glGetIntegerv(GL_MATRIX_MODE, &matrixMode);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glScalef(2.0f / window_width, -2.0f / window_height, 1.0f);
	glTranslatef(-window_width / 2.0f, -window_height / 2.0f, 0.0f);
	
	// draw 
	glDisable(GL_LIGHTING);
	glColor3f (1.0, 1.0, 1.0);
	if (gShowInfo) {
	
		sprintf (outString, "Camera Position: (%0.1f, %0.1f, %0.1f)", Camera::viewPos.x, Camera::viewPos.y, Camera::viewPos.z);
        BitmapPlotter::drawGLString (10, window_height - (lineSpacing * line++) - startOffest, outString);
		sprintf (outString, "Trackball Rotation: (%0.1f, %0.2f, %0.2f, %0.2f)", TrackBallCamera::gTrackBallRotation[0], TrackBallCamera::gTrackBallRotation[1], TrackBallCamera::gTrackBallRotation[2], TrackBallCamera::gTrackBallRotation[3]);
		BitmapPlotter::drawGLString (10, window_height - (lineSpacing * line++) - startOffest, outString);
		sprintf (outString, "World Rotation: (%0.1f, %0.2f, %0.2f, %0.2f)", Camera::gWorldRotation[0], Camera::gWorldRotation[1], Camera::gWorldRotation[2], Camera::gWorldRotation[3]);
		BitmapPlotter::drawGLString (10, window_height - (lineSpacing * line++) - startOffest, outString);
		sprintf (outString, "Aperture: %0.1f", Camera::aperture);
		BitmapPlotter::drawGLString (10, window_height - (lineSpacing * line++) - startOffest, outString);
		sprintf (outString, "Focus Distance: %0.1f", Camera::focalLength);
		BitmapPlotter::drawGLString (10, window_height - (lineSpacing * line++) - startOffest, outString);
        sprintf (outString, "Vertex Order: %u", splinePlotter.vertexOrder());
		BitmapPlotter::drawGLString (10, window_height - (lineSpacing * line++) - startOffest, outString);
	}
	
	if (gShowHelp) {
		line = 1;
		sprintf (outString, "Controls:\n");
		BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
		sprintf (outString, "left mouse button drag: rotate camera\n");
		BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
		sprintf (outString, "right mouse button drag: translate camera\n");
		BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
		sprintf (outString, "midle mouse button drag: dolly (zoom) camera\n");
		BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
		sprintf (outString, "arrows: aperture & focal length\n");
		BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
		sprintf (outString, "H: toggle help\n");
		BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
		sprintf (outString, "I: toggle info\n");
		BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
		sprintf (outString, "W: toggle wire frame\n");
		BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
	}
	
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(matrixMode);
	
	glViewport(vp[0], vp[1], vp[2], vp[3]);
}

/* special ********************************************************************/
void maindisplay(void)
{
	
    TrackBallCamera::update();
    
	glClearColor (0.2f, 0.2f, 0.4f, 1.0f);	// clear the surface
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if (gLines) {
		glColor3f (1.0, 1.0, 1.0); // no coloring
		glDisable(GL_LIGHTING);
		glCallList (gWireList);
	}
	if (gPolygons) {
		glEnable(GL_LIGHTING);
		glCallList (gSolidList);
	}
    
    splinePlotter.showTubes=true;
    splinePlotter.showVertices=true;
    splinePlotter.plot(10.0);

	
	drawGLText (Camera::screenWidth, Camera::screenHeight);
    
	glutSwapBuffers();
}

/* special ********************************************************************/
void special(int key, int px, int py)
{
  gLastKey = key;
  switch (key) {
	case GLUT_KEY_UP: // arrow forward, close in on world
		Camera::focalLength -= 0.5f;
		if (Camera::focalLength < 0.0f)
			Camera::focalLength = 0.0f;
		glutPostRedisplay();
		break;
	case GLUT_KEY_DOWN: // arrow back, back away from world
		Camera::focalLength += 0.5f;
		glutPostRedisplay();
		break;
	case GLUT_KEY_LEFT: // arrow left, smaller aperture
		Camera::aperture -= 0.5f;
		if (Camera::aperture < 0.0f)
			Camera::aperture = 0.0f;
		glutPostRedisplay();
		break;
	case GLUT_KEY_RIGHT: // arrow right, larger aperture
		Camera::aperture += 0.5f;
		glutPostRedisplay();
		break;
  }
}

/* key ************************************************************************/
void key(unsigned char inkey, int px, int py)
{
  gLastKey = inkey;
  switch (inkey) {
	case 27:
		exit(0);
		break;
	case 'h': // help
	case 'H':
		gShowHelp =  1 - gShowHelp;
		glutPostRedisplay();
		break;
	case 'i': // info
	case 'I':
		gShowInfo =  1 - gShowInfo;
		glutPostRedisplay();
		break;
	case 'w': // toggle wire
	case 'W':
		gLines = 1 - gLines;
		gPolygons = 1 - gPolygons;
		glutPostRedisplay();
		break;
  }
}

/* mouse **********************************************************************/
void mouse (int button, int state, int x, int y)
{
    TrackBallCamera::mouse(button,state,x,y);
}

/* spaceballmotion ************************************************************/
void spaceballmotion (int x, int y, int z)
{
    TrackBallCamera::spaceballmotion(x,y,z);
}

/* spaceballrotate ************************************************************/
void spaceballrotate (int rx, int ry, int rz)
{
    TrackBallCamera::spaceballrotate(rx,ry,rz);
}

/* reshape ********************************************************************/
void reshape (int w, int h)
{
	glViewport(0,0,(GLsizei)w,(GLsizei)h);
	Camera::screenWidth = w;
	Camera::screenHeight = h;
	glutPostRedisplay();
}

/* initGL *********************************************************************/
void initGL (void)
{
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glFrontFace(GL_CCW);
	glColor3f(1.0,1.0,1.0);
	BuildGeometry (kTranguloidTrefoil, 4, 64, 3, &gSolidList, &gWireList, &gPointList);
	glClearColor(0.0,0.0,0.0,0.0);         /* Background recColor */
    Camera::reset();
	glPolygonOffset (1.0, 1.0);
    LightMode::set(4);
	glEnable(GL_LIGHTING);
}

/* main ***********************************************************************/
int main (int argc, const char * argv[])
{
    glutInit(&argc, (char **)argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); // non-stereo for main window
	glutInitWindowPosition (300, 50);
	glutInitWindowSize (1024, 768);
	gMainWindow = glutCreateWindow("DDviewer");

    initGL(); // standard GL init

    glutReshapeFunc (reshape);
    glutDisplayFunc (maindisplay);
	glutKeyboardFunc (key);
	glutSpecialFunc (special);
	glutMouseFunc (mouse);
	glutSpaceballMotionFunc(spaceballmotion);
	glutSpaceballRotateFunc(spaceballrotate);
    
    
    ddr.list();
    splinePlotter.read(308);
    float xMin;
    float xMax;
    float yMin;
    float yMax;
    float zMin;
    float zMax;

    
//    splinePlotter.boundingBox(xMin,xMax,yMin,yMax,zMin,zMax);
//    
    //Camera::viewPos.x=0.5*(xMax+xMin);
    //Camera::viewPos.y=0.5*(yMax+yMin);
    //Camera::viewPos.z=0.5*(zMax+zMin);
//
    Camera::viewDir.x=0.5*(xMax+xMin);
    Camera::viewDir.y=0.5*(yMax+yMin);
   Camera::viewDir.z=0.5*(zMax+zMin);
    std::cout<<xMax<<" "<<xMin<<std::endl;

    
    glutMainLoop();
    return 0;
}





