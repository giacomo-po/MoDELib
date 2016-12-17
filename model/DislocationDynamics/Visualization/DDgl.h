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



#include <model/Utilities/TerminalColors.h>
#include <model/DislocationDynamics/Visualization/DDreader.h>

#include <model/DislocationDynamics/Visualization/Plotters/SplinePlotter.h>
#include <model/DislocationDynamics/Visualization/Plotters/CellPlotter.h>
#include <model/DislocationDynamics/Visualization/Plotters/PlanePlotter.h>
#include <model/DislocationDynamics/Visualization/Plotters/MeshPlotter.h>

#include <model/openGL/BitmapPlotter.h>
#include <model/openGL/GL2tga.h>
#include <model/openGL/GL2pdf.h>

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#include <GL/glx.h>
#endif


namespace model {
	
	
	/*************************************************************/
	/* DDgl **************************************************/
	/*************************************************************/
	struct DDgl
    {
        
        static float alpha;

        
        //  #include <model/DislocationDynamics/Visualization/eps/rendereps.h>
        
        enum MouseActions {MOUSE_NONE,MOUSE_ROTATE,MOUSE_PAN};
        
        MouseActions mouseAction;
        
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
		
		
        //		bool g_bButton1Down;
        //		bool g_bButton2Down;
		
//		bool saveTga;
		//bool savePdf;
        
		
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
		typedef SplinePlotter<dim,50,10> SplinePlotterType;
		typedef SingleSplinePlotter<dim,50,10> SingleSplinePlotterType;

        SplinePlotterType splinePlotter;
		CellPlotter  cellPlotter;
		PlanePlotter planePlotter;
		MeshPlotter  meshPlotter;
        
        
		
        
		
		
		/*************************************************************/
		/* initGL ****************************************************/
		void initGL() const
        {
            
            
            
            
            
            // Lights, see http://www.sjbaker.org/steve/omniv/opengl_lighting.html
            
			glEnable(GL_DEPTH_TEST);
			
            // LIGHTING
            
			glEnable(GL_LIGHTING); //Enable lighting
			
			
			//Add positioned light
			glEnable(GL_LIGHT0); //Enable light #0
            GLfloat light0Pos[] =   {10000.0f, 10000.0f, 10000.0f,1.0f};
            glLightfv(GL_LIGHT0, GL_POSITION, light0Pos);
            GLfloat light0Amb[] =   {0.0f, 0.0f, 0.0f, 1.0f}; // ambient color for GL_LIGHT0
            glLightfv(GL_LIGHT0, GL_AMBIENT, light0Amb);
            GLfloat light0Dif[] =   {1.0f, 1.0f, 1.0f, 1.0f}; // diffuse color for GL_LIGHT0
            glLightfv(GL_LIGHT0, GL_DIFFUSE, light0Dif);
            GLfloat light0Spec[] =   {1.0f, 1.0f, 1.0f, 1.0f}; // specular color for GL_LIGHT0
            glLightfv(GL_LIGHT0, GL_SPECULAR, light0Spec);
            
            //glLightfv(GL_LIGHT0, GL_SPECULAR, lightColor0);
            
            
            
            
//			glEnable(GL_COLOR_MATERIAL);
            
			
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
                
                sprintf (outString, "Step: %u", frameN);
                BitmapPlotter::drawGLString (10, windowHeight - (lineSpacing * line++) - startOffest, outString);
                sprintf (outString, "Vertex Order: %u", splinePlotter.vertexOrder());
                BitmapPlotter::drawGLString (10, windowHeight - (lineSpacing * line++) - startOffest, outString);
                sprintf (outString, "Edge Order: %u", splinePlotter.edgeOrder());
                BitmapPlotter::drawGLString (10, windowHeight - (lineSpacing * line++) - startOffest, outString);
                sprintf (outString, "Components: %u", splinePlotter.components());
                BitmapPlotter::drawGLString (10, windowHeight - (lineSpacing * line++) - startOffest, outString);
                
                
            }
            
            if (showControls)
            {
                line = 1;
                sprintf (outString, "HELP:\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "spacebar: hide/show help\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "right mouse button: show menu\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "left mouse button drag: rotate camera\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "SHIFT + left mouse button drag: translate camera\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "a: hide/show axis\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "b: hide/show Burgers vector\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "c: hide/show cells\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "e: hide/show edges as tubes\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "g: hide/show glide planes\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "i: set frame increment\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "j: show mesh region boundaries\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "k: locate specific vertex (by ID)\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "l: load specific frame\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "m: change mesh state (no mesh, undeformed, deformed)\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "n: hide/show glide plane normals\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "b: hide/show Burgers vector\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "p: hide/show PK force\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "r (hold r and press -/+): decrese/increse edge and vertex radius\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "v: hide/show vertices\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "x: auto-zoom\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "z (hold z and press -/+): zoom\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
                sprintf (outString, "1-6: show boundary stress (components 11,22,33,12,23,13)\n");
                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);
//                sprintf (outString, "0: hide boundary stress\n");
//                BitmapPlotter::drawGLString (10, (lineSpacing * line++) + startOffest, outString);

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
		void plotAxes()
        {
			if (showAxes)
            {
                

                glDisable(GL_COLOR_MATERIAL); // use glMaterialfv(...) to set material colors
                glEnable(GL_DEPTH_TEST);
                GLfloat materialColor[]={0.0f, 0.0f, 0.0f, 1.0};
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialColor);      // ambient color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialColor);      // diffuse color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialColor);  // specular color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialColor);  // emission color for the material

                // Plot axes
				glBegin(GL_LINES);
				glVertex3f(xMin, 0.0f,0.0f);
				glVertex3f(xMax, 0.0f,0.0f);
				
				glVertex3f(0.0f, yMin,0.0f);
				glVertex3f(0.0f, yMax,0.0f);
				
				glVertex3f(0.0f, 0.0f,zMin);
				glVertex3f(0.0f, 0.0f,zMax);
				glEnd();
                
                // Plot labels x1,x2,x3
                BitmapPlotter::renderString((VectorDimF()<<1.05f*xMax, 0.0f,0.0f).finished(),"x1");
                BitmapPlotter::renderString((VectorDimF()<<0.0f, 1.05f*yMax,0.0f).finished(),"x2");
                BitmapPlotter::renderString((VectorDimF()<<0.0f, 0.0f,1.05f*zMax).finished(),"x3");

                // Plot box limits
                BitmapPlotter::renderString((VectorDimF()<<xMax, 0.0f,0.0f).finished(),
                                            /*                       */ static_cast<std::ostringstream*>( &(std::ostringstream() << xMax) )->str());
                BitmapPlotter::renderString((VectorDimF()<<0.0f, yMax,0.0f).finished(),
                                            /*                       */ static_cast<std::ostringstream*>( &(std::ostringstream() << yMax) )->str());
                BitmapPlotter::renderString((VectorDimF()<<0.0f, 0.0f,zMax).finished(),
                                            /*                       */ static_cast<std::ostringstream*>( &(std::ostringstream() << zMax) )->str());
			}
		}
		
		
        /**********************************************************************/
		void autoCenter()
        {
			xMin=0.0f;
			xMax=0.0f;
			yMin=0.0f;
			yMax=0.0f;
			zMin=0.0f;
			zMax=0.0f;
            if(meshPlotter.p_mesh->simplices().size())
            {
                
                std::cout<<"mesh xMin="<<meshPlotter.p_mesh->xMin().transpose()<<std::endl;
                std::cout<<"mesh xMax="<<meshPlotter.p_mesh->xMax().transpose()<<std::endl;
                
                xMin=meshPlotter.p_mesh->xMin()(0);
                xMax=meshPlotter.p_mesh->xMax()(0);
                yMin=meshPlotter.p_mesh->xMin()(1);
                yMax=meshPlotter.p_mesh->xMax()(1);
                zMin=meshPlotter.p_mesh->xMin()(2);
                zMax=meshPlotter.p_mesh->xMax()(2);
            }
            else
            {
                splinePlotter.boundingBox(xMin,xMax,yMin,yMax,zMin,zMax);
            }
			xMean=0.5*(xMax+xMin);
			yMean=0.5*(yMax+yMin);
			zMean=0.5*(zMax+zMin);
			boxSize=xMax-xMin;
			boxSize= (boxSize > (yMax-yMin))? boxSize : (yMax-yMin);
			boxSize= (boxSize > (zMax-zMin))? boxSize : (zMax-zMin);
			//			scale=1.0;
            if (boxSize==0.0)
            {
                boxSize=1000;
            }
			scale=boxSize;
			radius=boxSize/400.0;
			
			centerPoint<< xMean, yMean, zMean;
			eyePoint=centerPoint+VectorDimF::Ones().normalized()*scale;
			upVector   << 0.0,0.0,1.0;
			transVector<< 0.0,0.0,0.0;
            std::cout<<"eyePoint="<<eyePoint.transpose()<<std::endl;
            std::cout<<"centerPoint="<<centerPoint.transpose()<<std::endl;

		}

        /**********************************************************************/
        void displayFunc()
        {
            
            
            drawScene();
            
            // Save screenshot
            if(savedframeN!=frameN)
            {
                // tga files
                if (GL2tga::saveTGA)
                {
                    std::stringstream filenameStream;
                    filenameStream << "tga/image_" << frameN;
                    GL2tga::saveAs(filenameStream.str());
                }
                
                // pdf files
                if (GL2pdf<DDgl>::savePDF)
                {
                    std::stringstream filenameStream;
                    filenameStream << "pdf/image_" << frameN;
                    GL2pdf<DDgl>(*this).saveAs(filenameStream.str());
                }
                savedframeN=frameN; // update savedframeN
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
        
        /**********************************************************************/
		void drawScene()
        {
			
            float transparency(0.1);
			glClearColor(1.0f,1.0f,1.0f,transparency);
            //			glClearColor(0.0f,0.0f,0.0f,transparency);
            
            //            glClearColor(0.3, 0.5, 0.8, 0.);
            
			
			// Clear information from the last draw
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			
			
			glMatrixMode(GL_PROJECTION); //erel
			glLoadIdentity(); //erel
			gluPerspective (60.0, (float)windowWidth/(float)windowHeight, 0.01*boxSize, 10.0*boxSize); //erel
			glMatrixMode(GL_MODELVIEW);//Switch to the drawing perspective
			
			glLoadIdentity();
			
			//			gluLookAt(xMean+1.0*boxSize,yMean+1.0*boxSize,zMean+1.0*boxSize, xMean, yMean, zMean, 0.0, 1.0, 0.0); //erel
			//VectorDimF eyePoint = centerPoint +  scale* center2eye;
			gluLookAt(eyePoint(0),eyePoint(1),eyePoint(2), centerPoint(0), centerPoint(1), centerPoint(2), upVector(0), upVector(1), upVector(2)); //giacomo
			
//			GLfloat lightPos0[] = {eyePoint(0)-centerPoint(0), eyePoint(1)-centerPoint(1), eyePoint(2)-centerPoint(2), 0.0f}; //Positioned at (4, 0, 8)
			GLfloat lightPos0[] = {eyePoint(0), eyePoint(1), eyePoint(2), 0.0f};

			glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
            
            
            meshPlotter.plot();

			splinePlotter.plot(radius);
			cellPlotter.plot();
            plotAxes();
			planePlotter.plot();
            
			
			
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
        DDgl(const SimplicialMesh<3>* const p_mesh,	int windowH, int windowW) :
//        /* init list */ p_mesh(p_mesh_in),
        /* init list */ stepIncrement(1),
        /* init list */ meshPlotter(p_mesh),
         windowHeight(windowH),
         windowWidth(windowW)
        {
			
 //               mesh.readMesh(0);
            
            EigenDataReader EDR;
            EDR.readScalarInFile("./DDinput.txt","parametrizationExponent",SingleSplinePlotterType::alpha);


            
			ddr.list();
			
			assert(ddr.size() && "NO FILES FOUND");
			ddrIter=ddr.begin();
			frameN=(ddrIter->first);
			rewind=false;
			
//			windowHeight = 768;
//			windowWidth  = 1024;
			
			old_y=0;
			old_x=0;
            
            
			
			scale = 1.0f;
			
			radius=1.0;
			
			autoplay = false;
			
            //			g_bButton1Down = false;
            //			g_bButton2Down = false;
			
            mouseAction=MOUSE_NONE;
            
			showAxes=false;
			
            //			showMeshStates=3; // 0=no mesh, 1= mesh, 2=deformed mesh
            //			showMesh=0;
			
            showControls=true;
			
			//saveTga=false;
			GL2tga::saveTGA=false;
			GL2pdf<DDgl>::savePDF=false;
            
			savedframeN=-1;
			
			xMean=0.0f;
			yMean=0.0f;
			zMean=0.0f;
			boxSize=0.0f;
			
            
            
            
            
			// Read first frame
			readFrame();
            
            autoCenter();
		}
		
		/**********************************************************************/
		void zoom ()
        {
			float mat[16];
			glGetFloatv(GL_MODELVIEW_MATRIX, mat);
			//mat[0], mat[4], mat[8] is the right vector
			//mat[1], mat[5], mat[9] is the up vector
			//mat[2], mat[6], mat[10] is the forward vector ( looking towards the screen )
			
			VectorDimF center2eye;
			center2eye<<mat[2], mat[6], mat[10];
			
			eyePoint=centerPoint+center2eye*scale;
		}
		

        /**********************************************************************/
		void keyUp (unsigned char key, int x, int y) {
			keyStates[key] = false; // Set the state of the current key to "not pressed"
		}
		
		
        /**********************************************************************/
		void handleKeypress(unsigned char key, int x, int y)
        { // The key that was press and the current mouse coordinates
			
			keyStates[key]=true; // Set the state of the current key to pressed
			
			switch (key) {
				case 27: //Escape key
					exit(0);
					break;
					
					
				case '-':
					if (keyStates['z'] && keyStates['f']){
						scale*=1.005;
//                        scale+=0.001*boxSize;
                        zoom();
					}
					if (keyStates['z'] && !keyStates['f']){
						scale*=1.05;
//                        scale+=0.01*boxSize;
                        zoom();
					}
					if (keyStates['r']){
						radius*=.95;
					}
					if (keyStates['d']){
						meshPlotter.dispScale*=0.75f;;
					}
					if (keyStates['p']){
						splinePlotter.PKfactor*=0.75f;;
					}
					break;
					
				case '=':
					if (keyStates['z'] && keyStates['f']){
						scale*=0.995;
//                        scale-=0.001*boxSize;
                        zoom();
					}
					if (keyStates['z'] && !keyStates['f']){
						scale*=0.95;
//                        scale-=0.01*boxSize;
                        zoom();
					}
					if (keyStates['r']){
						radius*=1.05;
					}
					if (keyStates['d']){
						meshPlotter.dispScale*=1.5f;;
					}
					if (keyStates['p']){
						splinePlotter.PKfactor*=1.5f;;
					}
					break;
					
				case ' ':
					showControls=!showControls;
					break;
					
                    
				case 'a':
					showAxes=!showAxes;
					break;
                    
                case 'b':
					splinePlotter.showBurgers=!splinePlotter.showBurgers;
					break;
                    
				case 'c':
					cellPlotter.showCells=!cellPlotter.showCells;
					break;
					
//				case 'd':
//					splinePlotter.deformedConfig=!splinePlotter.deformedConfig;
//					break;
					
				case 'e':
					splinePlotter.showTubes=!splinePlotter.showTubes;
					break;
					
					
				case 'g':
					planePlotter.showGlidePlanes=!planePlotter.showGlidePlanes;
					break;
					
					
                case 'i':
                    autoplay=false;
                    //if(splinePlotter.showSpecificVertex){
                    std::cout<<"Enter a step increment and press Enter: ";
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
					if (splinePlotter.isGood(temp,true) || splinePlotter.isGood(temp,false))
                    {
						frameN=temp;
						ddrIter=ddr.begin();
						std::advance(ddrIter,temp);
						readFrame();
					}
					break;
                    
				case 'm':
					meshPlotter.showMesh=(meshPlotter.showMesh+1)%meshPlotter.showMeshStates;
					break;
					
                case 'j':
                    meshPlotter.showRegionBoundaries=!meshPlotter.showRegionBoundaries;
                    break;

				case 'n':
					splinePlotter.showPlaneNormal=!splinePlotter.showPlaneNormal;
					break;
                    
				case 'p':
					splinePlotter.showPK=!splinePlotter.showPK;
					break;
                    
                case 'q':
                    splinePlotter.showQuadParticles=!splinePlotter.showQuadParticles;
                    break;
                    
                case 'r':
                    SingleSplinePlotterType::use_BurgersNorm=!SingleSplinePlotterType::use_BurgersNorm;
                    break;
                    
//				case 's':
//					saveTga=!saveTga;
//					break;
					
                case 's':
                    autoplay=false;
                    meshPlotter.showSpecificSimplex=!meshPlotter.showSpecificSimplex;
                    if(meshPlotter.showSpecificSimplex)
                    {
                        std::cout<<"Enter a Simplex ID: ";
                        size_t a;
                        size_t b;
                        size_t c;
                        size_t d;
                        std::cin>>a>>b>>c>>d;
                        Eigen::Matrix<size_t,4,1> xID;
                        xID<<a,b,c,d;
                        meshPlotter.specificSimplexID=SimplexTraits<3,3>::sortID(xID);
                    }
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
					splinePlotter.plotBoundarySegments=!splinePlotter.plotBoundarySegments;
					break;
                    
                case '1':
                    if (meshPlotter.stressCol==0)
                    {
                        meshPlotter.plotBndStress=!meshPlotter.plotBndStress;
                    }
                    else
                    {
                        meshPlotter.stressCol=0;
                        meshPlotter.plotBndStress=true;
                    }
					break;
                    
                case '2':
                    if (meshPlotter.stressCol==1)
                    {
                        meshPlotter.plotBndStress=!meshPlotter.plotBndStress;
                    }
                    else
                    {
                        meshPlotter.stressCol=1;
                        meshPlotter.plotBndStress=true;
                    }
					break;
                    
                case '3':
                    if (meshPlotter.stressCol==2)
                    {
                        meshPlotter.plotBndStress=!meshPlotter.plotBndStress;
                    }
                    else
                    {
                        meshPlotter.stressCol=2;
                        meshPlotter.plotBndStress=true;
                    }
					break;
                    
                case '4':
                    if (meshPlotter.stressCol==3)
                    {
                        meshPlotter.plotBndStress=!meshPlotter.plotBndStress;
                    }
                    else
                    {
                        meshPlotter.stressCol=3;
                        meshPlotter.plotBndStress=true;
                    }
					break;
                    
                case '5':
                    if (meshPlotter.stressCol==4)
                    {
                        meshPlotter.plotBndStress=!meshPlotter.plotBndStress;
                    }
                    else
                    {
                        meshPlotter.stressCol=4;
                        meshPlotter.plotBndStress=true;
                    }
					break;
                    
                case '6':
                    if (meshPlotter.stressCol==5)
                    {
                        meshPlotter.plotBndStress=!meshPlotter.plotBndStress;
                    }
                    else
                    {
                        meshPlotter.stressCol=5;
                        meshPlotter.plotBndStress=true;
                    }
					break;
                    
//                case '0':
//					splinePlotter.colorScheme=0;
//					break;
//                    
//                case '1':
//					splinePlotter.colorScheme=1;
//					break;
//                    
//                case '2':
//					splinePlotter.colorScheme=2;
//					break;
                    
                    
                    //                case 'b':
                    //                    int size=1;
                    //                    int doSort=1;
                    //					outputEPS(size,  doSort, "pippo.eps");
                    //					break;
                    
                    
			}
            
            if(meshPlotter.plotBndStress)
            {
                meshPlotter.read(frameN);
            }
		}
		
		
        
		
		/* mouseMotion ***********************************************/
		void mouseMotion(int x, int y){
			
            
            if (mouseAction != MOUSE_NONE) // ctrl key is pressed
            {
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
                
                
                if (mouseAction == MOUSE_PAN) // ctrl key is pressed
                {// panning
                    ytrans=static_cast<float>(old_y-y)/screeSize*boxSize*1.0;
                    xtrans=static_cast<float>(old_x-x)/screeSize*boxSize*1.0;
                    
                    centerPoint+= (xtrans*rightVector - ytrans * upVector);
                    
                }
                if (mouseAction == MOUSE_ROTATE) // ctrl key is pressed
                {// rotation
                    center2eye<< mat[2], mat[6], mat[10];
                    
                    //				std::cout<<"center2eye was"<<center2eye.transpose()<<std::endl;
                    thetaX = static_cast<float>(old_y-y)/screeSize*M_PI*2.0;
                    thetaY = static_cast<float>(old_x-x)/screeSize *M_PI*2.0;
                    Eigen::AngleAxis<float> Rx(thetaX,rightVector);
                    //				Eigen::AngleAxis<float> Ry(thetaY,currentUpVector);
                    Eigen::AngleAxis<float> Ry(thetaY,upVector);
                    
                    
                    center2eye=(Ry*Rx*center2eye).normalized();
                    
                    upVector=(Ry*Rx*upVector).normalized(); // apply rotation about rightvector first
                }
                
                
                
                
                eyePoint=centerPoint+center2eye*scale;
            }
            
			
            
			
			old_y=y;
			old_x=x;
            
            
            
            
            
		}
		
        /* mouseButton ***********************************************/
		void mouseButton(const int& button, const int& state, const int& x, const int& y)
        {/*! Callback function executed when a mouse button is pressed.
          */
            if (button == GLUT_LEFT_BUTTON) // mouse button is GLUT_LEFT_BUTTON
            {
                // Store mouse position when left button is pressed or released
				old_y=y;
				old_x=x;
                
                if (state == GLUT_DOWN) // mouse button is pressed down (as opposed to released)
                {
                    if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) // SHIFT is pressed
                    {
                        mouseAction=MOUSE_PAN;
                    }
                    else // SHIFT is not pressed
                    {
                        mouseAction=MOUSE_ROTATE;
                    }
                }
                else // mouse button is released
                {
                    mouseAction=MOUSE_NONE;
                }
			}
            
		}
		
		
		/* readFrame **********************************************/
		void readFrame()
        {
			double t0(clock());
			std::cout<<greenBoldColor<<"Reading frame "<<frameN<<"..."<<defaultColor<<std::endl;
            
			splinePlotter.read(frameN);
			planePlotter.read(frameN,true);
			cellPlotter.read(frameN,true);
			meshPlotter.read(frameN);
            
            
//			std::cout<<" done ["<<(clock()-t0)/CLOCKS_PER_SEC<<defaultColor<<std::endl;
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
		
        
		
	};
    
    	float DDgl::alpha=0.5;
	
	/********************************************************************/
	/********************************************************************/
} // namespace model
#endif
