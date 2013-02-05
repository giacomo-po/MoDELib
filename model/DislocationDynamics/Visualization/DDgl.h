/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDGL_H_
#define model_DDGL_H_


#include <Eigen/Dense>
#include <Eigen/Geometry>


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <iterator> // std::advance
#include <time.h>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#include <GL/glx.h>
#endif



#include <model/DislocationDynamics/Visualization/DDreader.h>

#include <model/DislocationDynamics/Visualization/Plotters/SplinePlotter.h>
#include <model/DislocationDynamics/Visualization/Plotters/CellPlotter.h>
#include <model/DislocationDynamics/Visualization/Plotters/PlanePlotter.h>
#include <model/DislocationDynamics/Visualization/Plotters/MeshPlotter.h>
#include <model/DislocationDynamics/Visualization/Plotters/BitmapPlotter.h>

#include <model/openGL/gl2ps.h>




namespace model {
	
	std::string defaultColor    = "\033[0m";	  // the default color for the console
	std::string redBoldColor    = "\033[1;31m";   // a bold red color
	std::string greenBoldColor  = "\033[1;32m";   // a bold green color
	std::string blueBoldColor   = "\033[1;34m";   // a bold blue color
	std::string magentaColor    = "\033[0;35m";   // a magenta color
	
	
    //#include <model/DislocationDynamics/Visualization/eps/rendereps_out.h>
    
    
	
	/*************************************************************/
	/* DDgl **************************************************/
	/*************************************************************/
	template <short unsigned int dim, float & alpha>
	struct DDgl{
		
	public:
		DDgl(){
			assert(0 && "DDgl.h: template specialization not implemented");
		}
		
	};
	
	/*************************************************************/
	/* DDgl **************************************************/
	/*************************************************************/
	template <float & alpha>
	struct DDgl<3,alpha>{
		
        
        //  #include <model/DislocationDynamics/Visualization/eps/rendereps.h>
        
		enum {dim=3};
		int windowHeight;
		int windowWidth;
        
		DDreader ddr;
		DDreader::const_iterator ddrIter;
        
        
        
        float xMin;
        float xMax;
        float yMin;
        float yMax;
        float zMin;
        float zMax;
        
        
		int old_y;
		int old_x;
        
		float scale;
		float radius;
		
		float xMean;
		float yMean;
		float zMean;
		float boxSize;
		
		bool autoplay;
		bool rewind;
		
		
		bool g_bButton1Down;
		bool g_bButton2Down;
		
		bool saveTga;
		
		bool showAxes;
		
		//short unsigned int showMeshStates;
		//short unsigned int showMesh;
		
		int frameN;
		int savedframeN;
		
		//		int deltaframe;
        
        bool showControls;
        
		bool keyStates [256];
		
        unsigned int stepIncrement;
		
		
		typedef Eigen::Matrix<float,dim,1> VectorDimF;
		VectorDimF centerPoint;
		VectorDimF eyePoint;
		//VectorDimF center2eye;
		VectorDimF upVector;
		VectorDimF transVector;
		
		// The plotters
		SplinePlotter<dim,50,10,alpha> splinePlotter;
		CellPlotter  cellPlotter;
		PlanePlotter planePlotter;
		MeshPlotter  meshPlotter;
        

		
        /* tgaScreenShot ************************************************/
		void epsScreenShot(const std::string& filename) const {
        
            GLfloat viewport[4];
            glGetFloatv(GL_VIEWPORT, viewport);
            //char *creator
            
            std::fstream file(filename.c_str(), std::ios::out);
            //fputs("%!PS-Adobe-3.0 EPSF-3.0\n", file);
            file<<"%!PS-Adobe-3.0 EPSF-3.0\n";
//            fprintf(file,"%!PS-Adobe-3.0 EPSF-3.0\n");
            /* Notice %% for a single % in the fprintf calls. */
   //         fprintf(file, "%%%%Creator: %s (using OpenGL feedback)\n", file, "DDviewer");
 //           fprintf(file, "%%%%BoundingBox: %g %g %g %g\n",
   //                 viewport[0], viewport[1], viewport[2], viewport[3]);
          //  "%!PS-Adobe-3.0 EPSF-3.0\n"
            file<<"%%BoundingBox: "<<viewport[0]<<" "<<viewport[1]<<" "<<viewport[2]<<" "<<viewport[3]<<"\n";
     //       fputs("%%EndComments\n", file);
       //     fputs("\n", file);
         //   fputs("gsave\n", file);
           // fputs("\n", file);

            file.close();

        }

		
		/* tgaScreenShot ************************************************/
		void tgaScreenShot(const std::string& filename) const {
			
			
            //			std::stringstream filenameStream;
            //			filenameStream << "tga/image_" << frameN << ".tga";
            //
            //			std::cout<<"Saving file"<<filenameStream.str()<<std::endl;
			
			//const char* filename=filenameStream.str().c_str();
			int x=glutGet(GLUT_WINDOW_WIDTH);
			int y=glutGet(GLUT_WINDOW_HEIGHT);
			
			// get the image data
			long imageSize = x * y * 3;
			unsigned char* data = new unsigned char[imageSize];
			glReadPixels(0,0,x,y, GL_BGR,GL_UNSIGNED_BYTE,data);
			
			// split x and y sizes into bytes
			int xa= x % 256;
			int xb= (x-xa)/256;
			
			int ya= y % 256;
			int yb= (y-ya)/256;
			
			//assemble the header
			unsigned char header[18]={0,0,2,0,0,0,0,0,0,0,0,0,(char)xa,(char)xb,(char)ya,(char)yb,24,0};
			
			// write header and data to file
			std::fstream File(filename.c_str(), std::ios::out | std::ios::binary);
			File.write (reinterpret_cast<const char*>(header), sizeof (char)*18);
			File.write (reinterpret_cast<const char*>(data), sizeof (char)*imageSize);
			
			File.close();
			
			delete[] data;
			data=NULL;
		}
		
		
		/*************************************************************/
		/* initGL ****************************************************/
		void initGL() const {
            
            // Lights, see http://www.sjbaker.org/steve/omniv/opengl_lighting.html
            
			glEnable(GL_DEPTH_TEST);
			
            // LIGHTING
            
			glEnable(GL_LIGHTING); //Enable lighting
			
			//Add ambient light
//			GLfloat ambientColor[] = {0.2f, 0.2f, 0.2f, 1.0 }; //Color (0.2, 0.2, 0.2)
//			GLfloat ambientColor[] = {0.8f, 0.8f, 0.8f, 1.0 }; //Color (0.2, 0.2, 0.2)
//			glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);
			
			//Add positioned light
			glEnable(GL_LIGHT0); //Enable light #0
            GLfloat light0Pos[] =   {1.0f, 1.0f, 0.0f,1.0f};
            glLightfv(GL_LIGHT0, GL_POSITION, light0Pos);
            GLfloat light0Amb[] =   {0.0f, 0.0f, 0.0f, 1.0f}; // ambient color for GL_LIGHT0
            glLightfv(GL_LIGHT0, GL_AMBIENT, light0Amb);
            GLfloat light0Dif[] =   {1.0f, 1.0f, 1.0f, 1.0f}; // diffuse color for GL_LIGHT0
            glLightfv(GL_LIGHT0, GL_DIFFUSE, light0Dif);        
            GLfloat light0Spec[] =   {1.0f, 1.0f, 1.0f, 1.0f}; // specular color for GL_LIGHT0
            glLightfv(GL_LIGHT0, GL_SPECULAR, light0Spec);      
            
            //glLightfv(GL_LIGHT0, GL_SPECULAR, lightColor0);
            
            
            
			//GLfloat lightColor0[] = {1.0f, 1.0f, 1.0f, 1.0f}; //Color (0.5, 0.5, 0.5)
			//GLfloat lightPos0[] =   {0.5f, 0.5f, 0.5f, 1.0f}; //Positioned at (4, 0, 8)
            //	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
            //            			glLightfv(GL_LIGHT0, GL_SPECULAR, lightColor0);
            
						glEnable(GL_COLOR_MATERIAL);

			
			glEnable(GL_NORMALIZE); //Automatically normalize normals
			glShadeModel(GL_SMOOTH); //Enable smooth shading

			
            //Disable color materials, so that glMaterial calls work
//			glDisable(GL_COLOR_MATERIAL);
		}
		
		
		/* showControls *******************************************************/
        void drawGLText () {
            char outString [256] = "";
            GLint matrixMode;
            GLint vp[4];
            GLint lineSpacing = 13;
            GLint line = 0;
            GLint startOffest = 7;
            
            glGetIntegerv(GL_VIEWPORT, vp);
            glViewport(0, 0, windowWidth, windowHeight);
            
            glGetIntegerv(GL_MATRIX_MODE, &matrixMode);
            glMatrixMode(GL_PROJECTION);
            glPushMatrix();
            glLoadIdentity();
            
            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
            glScalef(2.0f / windowWidth, -2.0f / windowHeight, 1.0f);
            glTranslatef(-windowWidth / 2.0f, -windowHeight / 2.0f, 0.0f);
            
            // draw
            glDisable(GL_LIGHTING);
            glColor3f (0.0, 0.0, 1.0);
            if (showControls) {
                
//                sprintf (outString, "Camera Position: (%0.1f, %0.1f, %0.1f)", Camera::viewPos.x, Camera::viewPos.y, Camera::viewPos.z);
//                BitmapPlotter::drawGLString (10, window_height - (lineSpacing * line++) - startOffest, outString);
//                sprintf (outString, "Trackball Rotation: (%0.1f, %0.2f, %0.2f, %0.2f)", TrackBallCamera::gTrackBallRotation[0], TrackBallCamera::gTrackBallRotation[1], TrackBallCamera::gTrackBallRotation[2], TrackBallCamera::gTrackBallRotation[3]);
//                BitmapPlotter::drawGLString (10, window_height - (lineSpacing * line++) - startOffest, outString);
//                sprintf (outString, "World Rotation: (%0.1f, %0.2f, %0.2f, %0.2f)", Camera::gWorldRotation[0], Camera::gWorldRotation[1], Camera::gWorldRotation[2], Camera::gWorldRotation[3]);
//                BitmapPlotter::drawGLString (10, window_height - (lineSpacing * line++) - startOffest, outString);
//                sprintf (outString, "Aperture: %0.1f", Camera::aperture);
//                BitmapPlotter::drawGLString (10, window_height - (lineSpacing * line++) - startOffest, outString);
//                sprintf (outString, "Focus Distance: %0.1f", Camera::focalLength);
//                BitmapPlotter::drawGLString (10, windowHeight - (lineSpacing * line++) - startOffest, outString);
                sprintf (outString, "Vertex Order: %u", splinePlotter.vertexOrder());
                BitmapPlotter::drawGLString (10, windowHeight - (lineSpacing * line++) - startOffest, outString);
                sprintf (outString, "Edge Order: %u", splinePlotter.edgeOrder());
                BitmapPlotter::drawGLString (10, windowHeight - (lineSpacing * line++) - startOffest, outString);
                sprintf (outString, "Components: %u", splinePlotter.components());
                BitmapPlotter::drawGLString (10, windowHeight - (lineSpacing * line++) - startOffest, outString);
                
                
            }
            
            if (showControls) {
                line = 1;
                sprintf (outString, "Controls:\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "spacebar: hide/show controls\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "left mouse button drag: rotate camera\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "right mouse button drag: translate camera\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "a: hide/show axis\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "c: hide/show cells\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "e: hide/show edges as tubes\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "g: hide/show glide planes\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "i: set frame increment planes\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "k: locate specific vertex (by ID)\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "l: load specific frame\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "m: change mesh state (no mesh, undeformed, deformed)\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "n: hide/show glide plane normals\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "p: hide/show PK force\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "r (hold r and press -/+): decrese/increse edge and vertex radius\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "s: saves frames to tga folder\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "v: hide/show vertices\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "x: auto-zoom\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "z (hold z and press -/+): zoom\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
            }
            
            glPopMatrix();
            glMatrixMode(GL_PROJECTION);
            glPopMatrix();
            glMatrixMode(matrixMode);
            
            glViewport(vp[0], vp[1], vp[2], vp[3]);
            
            
            glEnable(GL_LIGHTING);

        }
		
		
		/*************************************************************/
		/* plotAxes *************************************************/
		void plotAxes(){
			if (showAxes){
				//glDisable(GL_DEPTH_TEST);
				//glDisable (GL_BLEND);
				//glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
				
                //				glDisable(GL_ALPHA_TEST);
                //	//			glAlphaFunc(GL_GREATER, 0.0f);
				
				//glEnable(GL_COLOR_MATERIAL);
                //	//			glColor4f(0.0f, 0.0f, 0.0f, 0.1);
                
 //               glColor3f(0.0f, 0.0f, 0.0f);
                
                glDisable(GL_COLOR_MATERIAL); // use glMaterialfv(...) to set material colors
                GLfloat materialColor[]={0.0, 0.0, 0.0, 0.1};
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialColor);      // ambient color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialColor);      // diffuse color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialColor);  // specular color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialColor);  // emission color for the material
                
				
                //				glColor4f(0.0f, 0.0f, 0.0f, 0.8);
				glBegin(GL_LINES);
				glVertex3f(xMin, 0.0f,0.0f);
				glVertex3f(xMax, 0.0f,0.0f);
				
				glVertex3f(0.0f, yMin,0.0f);
				glVertex3f(0.0f, yMax,0.0f);
				
				glVertex3f(0.0f, 0.0f,zMin);
				glVertex3f(0.0f, 0.0f,zMax);
				glEnd();
                
                
                
                
                
                //          VectorDimF P(VectorDimF::Zero());
                //        BitmapPlotter::renderString(P,"hello");
                
				
                //			glEnable(GL_DEPTH_TEST); //Makes 3D drawing work when something is in front of something else
                
                
                BitmapPlotter::renderString((VectorDimF()<<0.0f, 0.0f,zMax).finished(),
                                            /*                       */ static_cast<std::ostringstream*>( &(std::ostringstream() << zMax) )->str());
                
                BitmapPlotter::renderString((VectorDimF()<<0.0f, yMax,0.0f).finished(),
                                            /*                       */ static_cast<std::ostringstream*>( &(std::ostringstream() << yMax) )->str());
                
                
                BitmapPlotter::renderString((VectorDimF()<<xMax, 0.0f,0.0f).finished(),
                                            /*                       */ static_cast<std::ostringstream*>( &(std::ostringstream() << xMax) )->str());
                
                
                
                
			}
		}
		
		
		////////////////////////////////////////////////
		void autoCenter(){
			xMin=0.0f;
			xMax=0.0f;
			yMin=0.0f;
			yMax=0.0f;
			zMin=0.0f;
			zMax=0.0f;
			splinePlotter.boundingBox(xMin,xMax,yMin,yMax,zMin,zMax);
			xMean=0.5*(xMax+xMin);
			yMean=0.5*(yMax+yMin);
			zMean=0.5*(zMax+zMin);
			boxSize=xMax-xMin;
			boxSize= (boxSize > (yMax-yMin))? boxSize : (yMax-yMin);
			boxSize= (boxSize > (zMax-zMin))? boxSize : (zMax-zMin);
			//			scale=1.0;
			scale=boxSize;
			radius=boxSize/400.0;
			
			centerPoint<< xMean, yMean, zMean;
			eyePoint=centerPoint+VectorDimF::Ones().normalized()*scale;
			upVector   << 0.0,0.0,1.0;
			transVector<< 0.0,0.0,0.0;
		}
		
        
        //void
        //render(void)
        //{
        //
        //    gluLookAt(eyePoint(0),eyePoint(1),eyePoint(2), centerPoint(0), centerPoint(1), centerPoint(2), upVector(0), upVector(1), upVector(2)); //giacomo
        //
        //    glTranslatef(centerPoint(0),centerPoint(1),centerPoint(2));
        //    int object=1;
        //    //glPushMatrix();
        //    //glRotatef(angle, 0.0, 1.0, 0.0);
        //    switch (object) {
        //        case 0:
        //            glutSolidSphere(100.0, 10, 10);
        //            break;
        //        case 1:
        //            glutSolidTorus(0.5, 1.0, 15, 15);
        //            break;
        //        case 2:
        //            glutSolidTeapot(1.0);
        //            break;
        //    }
        //   // glPopMatrix();
        //}
        
        
		/* drawScene *************************************************/
		//template<short unsigned int dim, double & alpha>
		void drawScene() {
			
            
			
			float transparency(0.1);
			glClearColor(1.0f,1.0f,1.0f,transparency);
//			glClearColor(0.0f,0.0f,0.0f,transparency);

            //            glClearColor(0.3, 0.5, 0.8, 0.);

			
			// Clear information from the last draw
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			
			
			glMatrixMode(GL_PROJECTION); //erel
			glLoadIdentity(); //erel
			gluPerspective (60.0, (float)windowWidth/(float)windowHeight, 10.0, 10.0*boxSize); //erel
			glMatrixMode(GL_MODELVIEW);//Switch to the drawing perspective
			
			glLoadIdentity();
			
			//			gluLookAt(xMean+1.0*boxSize,yMean+1.0*boxSize,zMean+1.0*boxSize, xMean, yMean, zMean, 0.0, 1.0, 0.0); //erel
			//VectorDimF eyePoint = centerPoint +  scale* center2eye;
			gluLookAt(eyePoint(0),eyePoint(1),eyePoint(2), centerPoint(0), centerPoint(1), centerPoint(2), upVector(0), upVector(1), upVector(2)); //giacomo
			
			//GLfloat lightColor0[] = {1.0f, 1.0f, 1.0f, 1.0f}; //Color (0.5, 0.5, 0.5)
			GLfloat lightPos0[] = {eyePoint(0)-centerPoint(0), eyePoint(1)-centerPoint(1), eyePoint(2)-centerPoint(2), 0.0f}; //Positioned at (4, 0, 8)
			//glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
			glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
			
			
			
			//			glScalef(scale,scale,scale);
			
			glEnable(GL_MULTISAMPLE);
			splinePlotter.plot(radius);
			cellPlotter.plot();
			plotAxes();
			meshPlotter.plot();
			planePlotter.plot();
            
            
            //            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            //            glEnable(GL_LINE_STIPPLE);
            //            gl2psEnable(GL2PS_LINE_STIPPLE);
            //            glDisable(GL_LINE_STIPPLE);
            //            gl2psDisable(GL2PS_LINE_STIPPLE);
            //
            //
            //            glFlush();
			
			if (saveTga && savedframeN!=frameN) {
                
                std::stringstream filenameStream;
                filenameStream << "tga/image_" << frameN << ".tga";
                std::string filename=filenameStream.str();
                std::cout<<"Saving file"<<filename<<std::endl;
                
				tgaScreenShot(filename);

                std::stringstream filenameStream1;
                filenameStream1 << "eps/image_" << frameN << ".eps";
                std::string filename1=filenameStream1.str();
                std::cout<<"Saving file"<<filename1<<std::endl;

                epsScreenShot(filename1);

                
                savedframeN=frameN;
                
                // EPS
//                std::cout<<"I'm here 1"<<std::endl;
//                FILE *fp;
//                std::cout<<"I'm here 2"<<std::endl;
//                int state = GL2PS_OVERFLOW, buffsize = 0;
//                GLint* viewport;
//                std::cout<<"I'm here 3"<<std::endl;
//                glGetIntegerv(GL_VIEWPORT, viewport);
//                std::cout<<"I'm here 4"<<std::endl;
//                fp = fopen("out.eps", "wb");
//                printf("Writing 'out.eps'... ");
//                while(state == GL2PS_OVERFLOW){
//                    buffsize += 1024*1024;
//                    //                    gl2psBeginPage("test", "gl2psTestSimple", NULL, GL2PS_EPS, GL2PS_SIMPLE_SORT,
//                    //                                   GL2PS_DRAW_BACKGROUND | GL2PS_USE_CURRENT_VIEWPORT,
//                    //                                   GL_RGBA, 0, NULL, 0, 0, 0, buffsize, fp, "out.eps");
//                    gl2psBeginPage("test", "gl2psTestSimple", NULL, GL2PS_EPS, GL2PS_SIMPLE_SORT,
//                                   GL2PS_DRAW_BACKGROUND | GL2PS_USE_CURRENT_VIEWPORT,
//                                   GL_RGBA, 0, NULL, 0, 0, 0, buffsize, fp, "out.eps");
//                    //display();
//                    state = gl2psEndPage();
//                }
//                fclose(fp);
//                printf("Done!\n");
                
                
			}
			
            drawGLText();


			
			glutSwapBuffers(); //Send the 3D scene to the screen
            
			
			// Load next frame is "autoplay" is on
			if (autoplay){
				if (rewind){
					previousFrame();
				}
				else {
					nextFrame();
				}
			}
		}
		
		
		
		/*************************************************************/
		//Called when the window is resized
		void handleResize(int w, int h) {
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
		
		
	public:
		
		
		/* Constructor **************************************************/
		DDgl(){
			
			ddr.list();
			
			assert(ddr.size() && "NO FILES FOUND");
			ddrIter=ddr.begin();
			frameN=(ddrIter->first);
			rewind=false;
			
			windowHeight = 768;
			windowWidth  = 1024;
			
			old_y=0;
			old_x=0;
            
            stepIncrement=1;
			
			scale = 1.0f;
			
			radius=1.0;
			
			autoplay = false;
			
			g_bButton1Down = false;
			g_bButton2Down = false;
			
			showAxes=true;
			
            //			showMeshStates=3; // 0=no mesh, 1= mesh, 2=deformed mesh
            //			showMesh=0;
			
            showControls=true;
			
			saveTga=false;
			savedframeN=-1;
			
			xMean=0.0f;
			yMean=0.0f;
			zMean=0.0f;
			boxSize=0.0f;
			
			// Read first frame
			readFrame();
            
            autoCenter();
		}
		
		
		
		
		/*************************************************************/
		/* handleKeypress ********************************************/
		void zoom () {
			float mat[16];
			glGetFloatv(GL_MODELVIEW_MATRIX, mat);
			//mat[0], mat[4], mat[8] is the right vector
			//mat[1], mat[5], mat[9] is the up vector
			//mat[2], mat[6], mat[10] is the forward vector ( looking towards the screen )
			
			VectorDimF center2eye;
			center2eye<<mat[2], mat[6], mat[10];
			
			eyePoint=centerPoint+center2eye*scale;
		}
		
		
		
		
		/* keyUp *****************************************************/
		void keyUp (unsigned char key, int x, int y) {
			keyStates[key] = false; // Set the state of the current key to "not pressed"
		}
		
		
		/* handleKeypress ********************************************/
		void handleKeypress(unsigned char key, int x, int y) { // The key that was press and the current mouse coordinates
			
			keyStates[key]=true; // Set the state of the current key to pressed
			
			switch (key) {
				case 27: //Escape key
					exit(0);
					break;
					
					
				case '-':
					if (keyStates['z'] && keyStates['f']){
						scale*=1.005;
						zoom();
					}
					if (keyStates['z'] && !keyStates['f']){
						scale*=1.05;
						zoom();
					}
					if (keyStates['r']){
						radius*=.95;
					}
					if (keyStates['q']){
						meshPlotter.dispScale*=0.75f;;
					}
					if (keyStates['a']){
						splinePlotter.PKfactor*=0.75f;;
					}
					break;
					
				case '=':
					if (keyStates['z'] && keyStates['f']){
						scale*=0.995;
						zoom();
					}
					if (keyStates['z'] && !keyStates['f']){
						scale*=0.95;
						zoom();
					}
					if (keyStates['r']){
						radius*=1.05;
					}
					if (keyStates['q']){
						meshPlotter.dispScale*=1.5f;;
					}
					if (keyStates['a']){
						splinePlotter.PKfactor*=1.5f;;
					}
					break;
					
				case ' ':
					showControls=!showControls;
					break;
					
                    
				case 'a':
					showAxes=!showAxes;
					break;
					
					//				case 'b':
					//					showInteriorBoundaryMesh=!showInteriorBoundaryMesh;
					//					break;
					
				case 'c':
					cellPlotter.showCells=!cellPlotter.showCells;
					break;
					
				case 'd':
					splinePlotter.deformedConfig=!splinePlotter.deformedConfig;
					break;
					
				case 'e':
					splinePlotter.showTubes=!splinePlotter.showTubes;
					break;
					
					
				case 'g':
					planePlotter.showGlidePlanes=!planePlotter.showGlidePlanes;
					break;
					
					
                case 'i':
                    autoplay=false;
                    //if(splinePlotter.showSpecificVertex){
                    std::cout<<"Enter a stepIncrement: ";
                    std::cin>>stepIncrement;
                    //}
					break;
					
                case 'k':
                    autoplay=false;
					splinePlotter.showSpecificVertex=!splinePlotter.showSpecificVertex;
                    if(splinePlotter.showSpecificVertex){
                        std::cout<<"Enter a Vertex ID: ";
                        size_t temp;
                        std::cin>>temp;
                        splinePlotter.specificVertexID=temp;
                    }
					break;
                    
                    
                case 'l':
                    autoplay=false;
					std::cout<<"Enter a step number ane press Enter: ";
					int temp;
					std::cin>>temp;
					if (splinePlotter.isGood(temp)) {
						frameN=temp;
						ddrIter=ddr.begin();
						std::advance(ddrIter,temp);
						readFrame();
					}
					break;
                    
				case 'm':
					meshPlotter.showMesh=(meshPlotter.showMesh+1)%meshPlotter.showMeshStates;
					break;
					
				case 'n':
					splinePlotter.showPlaneNormal=!splinePlotter.showPlaneNormal;
					break;
                    
				case 'p':
					splinePlotter.showPK=!splinePlotter.showPK;
					break;
					
				case 's':
					saveTga=!saveTga;
					break;
					
				case 'v':
					splinePlotter.showVertices=!splinePlotter.showVertices;
					break;
					
				case 'x':
					autoCenter();
					break;
                    
				case 't':
					splinePlotter.showVertexID=!splinePlotter.showVertexID;
					break;
                    
                case '0':
					splinePlotter.colorScheme=0;
					break;
                    
                case '1':
					splinePlotter.colorScheme=1;
					break;
                    
                case '2':
					splinePlotter.colorScheme=2;
					break;
                    
                    
                    //                case 'b':
                    //                    int size=1;
                    //                    int doSort=1;
                    //					outputEPS(size,  doSort, "pippo.eps");
                    //					break;
                    
                    
			}
		}
		
		
        
		
		/* mouseMotion ***********************************************/
		void mouseMotion(int x, int y){
			
			float mat[16];
			glGetFloatv(GL_MODELVIEW_MATRIX, mat);
			//mat[0], mat[4], mat[8] is the right vector
			//mat[1], mat[5], mat[9] is the up vector
			//mat[2], mat[6], mat[10] is the forward vector ( looking towards the screen )
			
			VectorDimF center2eye;
			center2eye<<mat[2], mat[6], mat[10];
			
			
			VectorDimF rightVector(mat[0], mat[4], mat[8]);
			//			VectorDimF currentUpVector(mat[1], mat[5], mat[9]);
			//			currentUpVector<<mat[1], mat[5], mat[9];
			upVector<<mat[1], mat[5], mat[9];
			
			//			upVector=currentUpVector;
			
			
			float thetaX(0.0);
			float thetaY(0.0);
			float ytrans(0.0);
			float xtrans(0.0);
			
			
			float screeSize(0.5*(windowHeight+windowWidth));
			
			if (g_bButton1Down){
				center2eye<< mat[2], mat[6], mat[10];
				
				//				std::cout<<"center2eye was"<<center2eye.transpose()<<std::endl;
				thetaX = static_cast<float>(old_y-y)/screeSize*M_PI*2.0;
				thetaY = static_cast<float>(old_x-x)/screeSize *M_PI*2.0;
				Eigen::AngleAxis<float> Rx(thetaX,rightVector);
				//				Eigen::AngleAxis<float> Ry(thetaY,currentUpVector);
				Eigen::AngleAxis<float> Ry(thetaY,upVector);
				
				
				center2eye=(Ry*Rx*center2eye).normalized();
				
				upVector=(Ry*Rx*upVector).normalized(); // apply rotation about rightvector first
				
				//				std::cout<<"center2eye was"<<center2eye.transpose()<<std::endl;
				//				std::cout<<"center2eye was"<<center2eye.transpose()<<std::endl;
				
				
				//				Eigen::AngleAxis<float> Rx(thetaX,rightVector);
				////				Eigen::AngleAxis<float> Ry(thetaY,currentUpVector);
				//				Eigen::AngleAxis<float> Ry(thetaY,upVector);
				//				upVector=(Ry*Rx*upVector).normalized();
				//				center2eye=(Ry*Rx*center2eye).normalized();
				//				std::cout<<"center2eye is now"<<center2eye.transpose()<<std::endl;
				
				//				eyePoint = centerPoint +  center2eye;
			}
			else if(g_bButton2Down) {
				ytrans=static_cast<float>(old_y-y)/screeSize*boxSize*1.0;
				xtrans=static_cast<float>(old_x-x)/screeSize*boxSize*1.0;
				
				centerPoint+= (xtrans*rightVector - ytrans * upVector);
				
				
				//				thetaX=std::atan(ytrans/scale);
				//				scale=std::pow(std::pow(scale,2)+std::pow(ytrans,2),0.5f);
				
				//				thetaY=std::atan(-xtrans/scale);
				//				scale=std::pow(std::pow(scale,2)+std::pow(xtrans,2),0.5f);
				
				//				centerPoint+= (xtrans*rightVector - ytrans * currentUpVector);
				//				centerPoint+= (xtrans*rightVector - ytrans * upVector);
				//				centerPoint+= (xtrans*rightVector - ytrans * upVector);
				
			}
			
			//			Eigen::AngleAxis<float> Rx(thetaX,rightVector);
			//			//				Eigen::AngleAxis<float> Ry(thetaY,currentUpVector);
			//			Eigen::AngleAxis<float> Ry(thetaY,upVector);
			eyePoint=centerPoint+center2eye*scale;
			//std::pow(scale,2)+std::pow(xtrans,2)+std::pow(ytrans,2);
			//			scale=std::pow(std::pow(scale,2)+std::pow(xtrans,2)+std::pow(ytrans,2),0.5f);
			
			
			//			eyePoint = centerPoint +  center2eye*scale;
			
			old_y=y;
			old_x=x;
            
            
 //           GLfloat lightColor0[] = {1.0f, 1.0f, 1.0f, 1.0f}; //Color (0.5, 0.5, 0.5)
	//		GLfloat lightPos0[] =   {eyePoint(0), eyePoint(1), eyePoint(2), 1.0f}; //Positioned at (4, 0, 8)
	//		glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
	//		glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
            
            
            
		}
		
		
		
		/* readFrame **********************************************/
		void readFrame(){
			double t0(clock());
			std::cout<<greenBoldColor<<"Reading frame "<<frameN<<"..."<<defaultColor<<std::endl;
            
			splinePlotter.read(frameN);
			planePlotter.read<true>(frameN);
			cellPlotter.read(frameN);
			meshPlotter.read(frameN);
            
            
			std::cout<<" done ["<<(clock()-t0)/CLOCKS_PER_SEC<<defaultColor<<std::endl;
		}
		
		
		/* previousFrame **********************************************/
		void previousFrame(){
            DDreader::const_iterator nextIter(ddrIter);
            for (int i=0;i<stepIncrement;++i){
                --nextIter;
            }
			if (ddrIter!=ddr.begin()){
				//--ddrIter;
                ddrIter=nextIter;
				frameN=ddrIter->first;
				readFrame();
			}
		}
		
		/* nextFrame *************************************************/
		void nextFrame(){
			DDreader::const_iterator nextIter(ddrIter);
            for (int i=0;i<stepIncrement;++i){
                ++nextIter;
            }
			if (nextIter!=ddr.end()){
                //				++ddrIter;
                ddrIter=nextIter;
				frameN=ddrIter->first;
				readFrame();
			}
		}
		
		/* specialKey ************************************************/
		void specialKey(int key, int x, int y){
			switch(key){
				case GLUT_KEY_LEFT : 
					autoplay = true;
					rewind=true;
					break;
					
				case GLUT_KEY_RIGHT :
					autoplay = true;
					rewind=false;
					break;
				case GLUT_KEY_DOWN :		
					autoplay=false;
					previousFrame();
					break;
					
				case GLUT_KEY_UP :		
					autoplay=false;
					nextFrame();
					break;
			}
		}
		
		/* mouseButton ***********************************************/
		void mouseButton(int button, int state, int x, int y){
			// Respond to mouse button presses.
			// If button1 pressed, mark this state so we know in motion function.
			if (button == GLUT_LEFT_BUTTON){
				old_y=y;
				old_x=x;				
				g_bButton1Down = (state == GLUT_DOWN) ? true : false;
			}
			else if(button == GLUT_RIGHT_BUTTON){
				old_y=y;
				old_x=x;
				g_bButton2Down = (state == GLUT_DOWN) ? true : false;
			}
		}
		
	};
	
	/********************************************************************/
	/********************************************************************/
} // namespace model
#endif
