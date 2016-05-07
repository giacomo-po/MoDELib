/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <model/Utilities/EigenDataReader.h>
#include <model/DislocationDynamics/Visualization/DDgl.h>
#include <model/Geometry/Splines/SplineConsts.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

//float centripetalf=model::centripetal;
//
//float chordalf=model::chordal;
//
//typedef model::DDgl<centripetalf> DDgl;


model::SimplicialMesh<3>* p_mesh;  // declaring as pointer is necessary for glut
model::DDgl* p_DDgl; // declaring as pointer is necessary for glut


namespace model
{
	
	
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
	void DDglut_DisplayFunc()
    {
		p_DDgl->displayFunc();
	}
	
	
	/*************************************************************/
	void update(int value)
    {
		//	_angle += 1.5f;
		//	if (_angle > 360) {
		//		_angle -= 360;
		//	}
		
		glutPostRedisplay();
		glutTimerFunc(25, update, 0);
	}
	
	
	/*************************************************************/
	void DDglut_handleKeypress(unsigned char key, int x, int y)
    {// key and current mouse coordinates
		p_DDgl->handleKeypress(key,x,y);
	}
	
	/*************************************************************/
	void DDglut_mouseButton(int button, int state, int x, int y)
    {
		//! see DDviewer::DDglut_mouseButton
		p_DDgl->mouseButton(button,state,x,y);
	}
	
	/*************************************************************/
	void DDglut_mouseMotion(int x, int y)
    {
		p_DDgl->mouseMotion(x,y);
	}
	
	/*************************************************************/
	void DDglut_specialKey(int key, int x, int y)
    {
		p_DDgl->specialKey(key,x,y);
	}
	
	void DDglut_keyUp(unsigned char key, int x, int y){
		p_DDgl->keyUp(key,x,y);
	}
    
	void DDglut_menu1(int item)
    {
        switch (item)
        {
            case 0:
                GL2tga::saveTGA=!GL2tga::saveTGA;
                break;
                
            case 1:
                GL2pdf<DDgl>::savePDF=!GL2pdf<DDgl>::savePDF;
                break;
                
            default:
                break;
        }
	}
    
    void DDglut_menu2(int item)
    {
        switch (item) {
            case DDgl::SingleSplinePlotterType::colorBurgers:
                DDgl::SplinePlotterType::colorScheme=DDgl::SingleSplinePlotterType::colorBurgers;
                break;
                
            case DDgl::SingleSplinePlotterType::colorSessile:
                DDgl::SplinePlotterType::colorScheme=DDgl::SingleSplinePlotterType::colorSessile;
                break;
                
            case DDgl::SingleSplinePlotterType::colorNormal:
                DDgl::SplinePlotterType::colorScheme=DDgl::SingleSplinePlotterType::colorNormal;
                break;
                
            case DDgl::SingleSplinePlotterType::colorEdgeScrew:
                DDgl::SplinePlotterType::colorScheme=DDgl::SingleSplinePlotterType::colorEdgeScrew;
                break;
                
            case DDgl::SingleSplinePlotterType::colorComponent:
                DDgl::SplinePlotterType::colorScheme=DDgl::SingleSplinePlotterType::colorComponent;
                break;
                
            default:
                DDgl::SplinePlotterType::colorScheme=DDgl::SingleSplinePlotterType::colorBurgers;
                break;
        }
	}
    
    void DDglut_menu3(int item)
    {
	//	p_DDgl->menu(item);
	}
	
	/*************************************************************/
	int DDglut(int argc, char** argv)
    {
        
        int meshID(0);
        EigenDataReader EDR;
        bool use_boundary=false;
        EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
        if (use_boundary)
        {
            EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
        }
        
        int windowHeight = 768;
        int windowWidth  = 1024;
        if (argc>1)
        {
            windowWidth=atoi(argv[1]);
            if (argc>2)
            {
                windowHeight=atoi(argv[2]);
            }
        }


        
        p_mesh = new model::SimplicialMesh<3>(meshID);
        p_DDgl = new model::DDgl(p_mesh,windowHeight,windowWidth);
		//! \brief A Glut wrapper for model::DDviewer
		
		//Initialize GLUT
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
		glutInitWindowSize(p_DDgl->windowWidth, p_DDgl->windowHeight); //Set windows size
		
		//Create the window
		glutCreateWindow("DDviewer");
		p_DDgl->initGL();
		//		initRendering();//Initialize rendering
		
		//Set handler functions for drawing, keypresses, and windows resizes
		glutDisplayFunc(DDglut_DisplayFunc);
		glutKeyboardFunc(DDglut_handleKeypress);
		glutMouseFunc (DDglut_mouseButton);
		glutMotionFunc (DDglut_mouseMotion);	
		//	glutMouseWheelFunc (mouseWheel);
		glutSpecialFunc(DDglut_specialKey);
//		glutKeyboardFunc(DDglut_keyPressed); // Tell GLUT to use the method "keyPressed" for key presses  
		glutKeyboardUpFunc(DDglut_keyUp); // Tell GLUT to use the method "keyUp" for key up events
		glutReshapeFunc(DDglut_handleResize);
        
        GLint m1_choice=glutCreateMenu(DDglut_menu1);
        glutAddMenuEntry("tga", 0);
        glutAddMenuEntry("pdf", 1);
        glutAddMenuEntry("eps", 2);

        
        GLint m2_choice=glutCreateMenu(DDglut_menu2);
        glutAddMenuEntry("Burgers vector", DDgl::SingleSplinePlotterType::colorBurgers);
        glutAddMenuEntry("Glissile-Sessile", DDgl::SingleSplinePlotterType::colorSessile);
        glutAddMenuEntry("Plane normal", DDgl::SingleSplinePlotterType::colorNormal);
        glutAddMenuEntry("Edge-Screw", DDgl::SingleSplinePlotterType::colorEdgeScrew);
        glutAddMenuEntry("Network Component", DDgl::SingleSplinePlotterType::colorComponent);

//        glutAddMenuEntry("Network Component", 3);
        
        glutCreateMenu(DDglut_menu3);
        glutAddSubMenu("Save image...",m1_choice);
        glutAddSubMenu("Coloring scheme",m2_choice);

        
        // Associate a mouse button with menu
        glutAttachMenu(GLUT_RIGHT_BUTTON);
		
		glutTimerFunc(25, update, 0); //Add a timer
		
		glutMainLoop();//Start the main loop. glutMainLoop doesn't return.
        
		return 0;//This line is never reached, its purpose is to avoid compilation error
	}
	
} // namespace model

