/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <model/Dislocations/Visualization/DDgl.h>
#include <model/Geometry/Splines/SplineConsts.h>

float centripetalf=model::centripetal;

float chordalf=model::chordal;

model::DDgl<3,centripetalf> DDv;
//model::DDgl<3,chordalf> DDv;


namespace model {
	
	
	/*************************************************************/
	//Called when the window is resized
	void DDglut_handleResize(int w, int h) {
		//Tell opengl how to convert from coordinates to pixels values
		glViewport(0, 0, w, h);
		glMatrixMode(GL_PROJECTION);//Switch to set the camera perspective
		//set the camera perspective
		glLoadIdentity();//Resets the camera
		gluPerspective(45.0,// The camera angle which the eye swips
					   (double)w / (double)h,// The width-to-heigth ratio
					   1.0,// The near z clipping coordinate
					   200.0);// The far z clipping coordinate
		
		
	}
	
	/*************************************************************/
	//Draws the 3D scene
	void DDglut_drawScene(){
		DDv.drawScene();
	}
	
	
	/*************************************************************/
	//float _angle = -70.0f;
	
	void update(int value) {
		//	_angle += 1.5f;
		//	if (_angle > 360) {
		//		_angle -= 360;
		//	}
		
		glutPostRedisplay();
		glutTimerFunc(25, update, 0);
	}
	
	
	/*************************************************************/
	void DDglut_handleKeypress(unsigned char key, int x, int y) {// key and current mouse coordinates
		DDv.handleKeypress(key,x,y);
	}
	
	/*************************************************************/
	void DDglut_mouseButton(int button, int state, int x, int y){
		//! see DDviewer::DDglut_mouseButton
		DDv.mouseButton(button,state,x,y);
	}
	
	/*************************************************************/
	void DDglut_mouseMotion(int x, int y){
		DDv.mouseMotion(x,y);
	}
	
	/*************************************************************/
	void DDglut_specialKey(int key, int x, int y){
		DDv.specialKey(key,x,y);
	}
	
	void DDglut_keyUp(unsigned char key, int x, int y){
		DDv.keyUp(key,x,y);
	}

	
	
	/*************************************************************/
	int DDglut(int argc, char** argv) {
		
		//! \brief A Glut wrapper for model::DDviewer
		
		//Initialize GLUT
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
		glutInitWindowSize(DDv.windowWidth, DDv.windowHeight); //Set windows size
		
		//Create the window
		glutCreateWindow("DDviewer");
		DDv.initGL();
		//		initRendering();//Initialize rendering
		
		//Set handler functions for drawing, keypresses, and windows resizes
		glutDisplayFunc(DDglut_drawScene);
		glutKeyboardFunc(DDglut_handleKeypress);
		glutMouseFunc (DDglut_mouseButton);
		glutMotionFunc (DDglut_mouseMotion);	
		//	glutMouseWheelFunc (mouseWheel);
		glutSpecialFunc(DDglut_specialKey);
//		glutKeyboardFunc(DDglut_keyPressed); // Tell GLUT to use the method "keyPressed" for key presses  
		glutKeyboardUpFunc(DDglut_keyUp); // Tell GLUT to use the method "keyUp" for key up events
		glutReshapeFunc(DDglut_handleResize);
		
		glutTimerFunc(25, update, 0); //Add a timer
		
		glutMainLoop();//Start the main loop. glutMainLoop doesn't return.
		return 0;//This line is never reached, its purpose is to avoid compilation error
	}
	
} // namespace model

