/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPLINEPLOTTER_H_
#define model_SPLINEPLOTTER_H_


#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <string.h>
#include <sstream>
#include <float.h>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>


#include <Eigen/Core>
#include <Eigen/Geometry>

#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>

#include <model/DislocationDynamics/Visualization/Plotters/BitmapPlotter.h>
#include <model/DislocationDynamics/Visualization/RGBmap.h>



namespace model {
	
	
	/*********************************************************************/
	/*********************************************************************/
	template <int dim, int Np, int Nc, float & alpha>
	class SingleSplinePlotter{
		
	public:
		typedef float scalarType;
		typedef Eigen::Matrix<scalarType,dim,Np> MatrixDimNp;
		typedef Eigen::Matrix<scalarType,dim,Nc> MatrixDimNc;
		typedef Eigen::Matrix<scalarType,dim,1>  VectorDim;
		
	private:
		//using HermiteCubicSplineAxis<dim,Np,alpha>::tubeAxis;
		//using HermiteCubicSplineAxis<dim,Np,alpha>::tubeTangents;
		
		MatrixDimNp tubeAxis;
		MatrixDimNp tubeTangents;
		std::vector<MatrixDimNc> tubeCircles;
		
		
		//		MatrixDimNc circle;
		VectorDim planeNormal;
		VectorDim burgers;
		VectorDim chord;
		
		//		scalarType radius;
		
		scalarType specularity;
		scalarType emissivity;
		scalarType shininess;
		scalarType transparency;
		
//		VectorDim colorVector;
        
        const bool isSessile;
		
		//	unsigned int edgeTypes;	// 0 = full axis only, 1 = full tubes, 2 =  tubes with direction
		
		
		/* getCircle **************************************************************/
		MatrixDimNc getCircle(const int& k) const {
			MatrixDimNc circle;
			for (int c=0;c<Nc;++c){
				circle.col(c)=  Eigen::AngleAxisf(2.0f * M_PI/Nc*c,tubeTangents.col(k)) * tubeTangents.col(k).cross(planeNormal);
			}
			return circle;
		}
		
		//	enum {edgeTypes=3};
		
	public:
		
		/* Constructor ************************************************************/
		SingleSplinePlotter(const Eigen::Matrix<scalarType,dim,6>& P0T0P1T1BN) :
		/* init list */ planeNormal(P0T0P1T1BN.col(5).normalized()),
		/* init list */ burgers(P0T0P1T1BN.col(4).normalized()),
		/* init list */ specularity(0.5f),
		/* init list */ emissivity(0.0f),
		/* init list */ shininess(1.0f),
		/* init list */ transparency(0.1f),
        /* init list */ isSessile(std::fabs(planeNormal.dot(burgers))>FLT_EPSILON)
        {
            
			tubeCircles.reserve(Np); // set the capacity of the vector to speed-up push_back
			chord = P0T0P1T1BN.col(2)-P0T0P1T1BN.col(0);
			scalarType g = std::pow(chord.norm(),alpha);
			
			for (int k=0;k<Np;++k){
				scalarType u1=k*1.0/(Np-1);
				scalarType u2=u1*u1;
				scalarType u3=u2*u1;
				
				// Compute positions along axis
				tubeAxis.col(k) =   ( 2.0f*u3-3.0f*u2+1.0f) * P0T0P1T1BN.col(0)
				/*************/ + g*(      u3-2.0f*u2+u1)   * P0T0P1T1BN.col(1)
				/*************/ +   (-2.0f*u3+3.0f*u2)      * P0T0P1T1BN.col(2)
				/*************/ + g*(      u3-u2)           * P0T0P1T1BN.col(3);
				
				// Compute tangents along axis
				tubeTangents.col(k)= (     ( 6.0f*u2-6.0f*u1)      * P0T0P1T1BN.col(0)
									  /*                  */ + g*( 3.0f*u2-4.0f*u1+1.0f) * P0T0P1T1BN.col(1)
									  /*                  */ +   (-6.0f*u2+6.0f*u1)      * P0T0P1T1BN.col(2)
									  /*                  */ + g*( 3.0f*u2-2.0f*u1)      * P0T0P1T1BN.col(3) ).normalized();
				
				// Compute unit vectors in radial direction orthogonal to axis
				tubeCircles.push_back(getCircle(k));
			}			
		}
		
		
        void flipColor(VectorDim& colorVector) const
        {
            if(colorVector(0)<0.0){
				colorVector*=-1;
			}
			else if(colorVector(0)==0.0){
				if(colorVector(1)<0.0){
					colorVector*=-1;
				}
				else if(colorVector(1)==0.0){
					if(colorVector(2)<0.0){
						colorVector*=-1;
					}
				}
			}
			
			//			VectorDim colorVector = burgers + VectorDim::Ones(dim) * burgers.norm();
			colorVector = (colorVector + VectorDim::Ones(dim) * colorVector.norm()).eval();
			
            //		colorVector << 0.0f,0.6f,0.4f;
			colorVector.normalize();
        
        }
        
		
		/*********************************************************************/
		void plot(const scalarType& radius, const bool& showTubes, const bool& showPlaneNormal, const bool& showBurgers, const int& colorScheme) const {
			
            // 1- Define the color
            VectorDim colorVector;
            switch (colorScheme) {
                case 1:
                    colorVector(0)= isSessile? 1.0 : 0.1;
                    colorVector(1)= isSessile? 0.5 : 0.4;
                    colorVector(2)= isSessile? 0.0 : 0.9;
                    break;
                    
                case 2:
                    colorVector = planeNormal;
                    flipColor(colorVector);
                    break;
                    
                default:
                    colorVector = burgers;
                    flipColor(colorVector);
                    break;
            }
            
			glDisable(GL_COLOR_MATERIAL); // use glMaterialfv(...) to set material colors
			glEnable(GL_DEPTH_TEST);
            GLfloat materialColor[]={colorVector(0), colorVector(1), colorVector(2), 1.0};
			GLfloat materialSpecular[] = {specularity, specularity, specularity, 1.0f};
			GLfloat materialEmission[] = {emissivity, emissivity, emissivity, 1.0f};
            // note that glMaterialfv(...) works when GL_COLOR_MATERIAL is disabled
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialColor);      // ambient color for the material
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialColor);      // diffuse color for the material 
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialColor);  // specular color for the material
			glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialEmission);  // emission color for the material
			glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess); //The shininess parameter
			
			if (showTubes)
            {
				for (int k=1;k<Np;++k){
					//					tubeCircles[k]=getCircle(k);
					glBegin(GL_TRIANGLE_STRIP);
					//glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
					for (int c=0;c<Nc;++c){
						glNormal3f(tubeCircles[k-1](0,c),tubeCircles[k-1](1,c),tubeCircles[k-1](2,c));
						glVertex3f(tubeAxis(0,k-1)+radius*tubeCircles[k-1](0,c),tubeAxis(1,k-1)+radius*tubeCircles[k-1](1,c),tubeAxis(2,k-1)+radius*tubeCircles[k-1](2,c));
						glNormal3f(tubeCircles[k](0,c),tubeCircles[k](1,c),tubeCircles[k](2,c));
						//if(showTubes==1){
						glVertex3f(tubeAxis(0,k)+radius*tubeCircles[k](0,c),tubeAxis(1,k)+radius*tubeCircles[k](1,c),tubeAxis(2,k)+radius*tubeCircles[k](2,c));
						//}
						//else{
						//	glVertex3f(tubeAxis(0,k),tubeAxis(1,k),tubeAxis(2,k));
						//}
					}
					glNormal3f(tubeCircles[k-1](0,0),tubeCircles[k-1](1,0),tubeCircles[k-1](2,0));
					glVertex3f(tubeAxis(0,k-1)+radius*tubeCircles[k-1](0,0),tubeAxis(1,k-1)+radius*tubeCircles[k-1](1,0),tubeAxis(2,k-1)+radius*tubeCircles[k-1](2,0));
					glNormal3f(tubeCircles[k](0,0),tubeCircles[k](1,0),tubeCircles[k](2,0));
					//if(showTubes==1){
					glVertex3f(tubeAxis(0,k)+radius*tubeCircles[k](0,0),tubeAxis(1,k)+radius*tubeCircles[k](1,0),tubeAxis(2,k)+radius*tubeCircles[k](2,0));
					//}
					//else{
					glVertex3f(tubeAxis(0,k),tubeAxis(1,k),tubeAxis(2,k));
					//}
					glEnd();
					//tubeCircles[k-1]=tubeCircles[k];
				}
				
			}
			else
            {
				glBegin(GL_LINE_STRIP);
				for (int k=0;k<Np;++k)
                {
					glVertex3f(tubeAxis(0,k),tubeAxis(1,k),tubeAxis(2,k));
				}
				glEnd();
			}
			
			if(showPlaneNormal)
            {
				glBegin(GL_LINES);
				int k(Np*0.5);
				glVertex3f(tubeAxis(0,k),tubeAxis(1,k),tubeAxis(2,k));
				glVertex3f(tubeAxis(0,k)+planeNormal(0)*10.0*radius,tubeAxis(1,k)+planeNormal(1)*10.0*radius,tubeAxis(2,k)+planeNormal(2)*10.0*radius);
				glEnd();
			}
            if(showBurgers)
            {
				glBegin(GL_LINES);
				int k(Np*0.5);
				glVertex3f(tubeAxis(0,k),tubeAxis(1,k),tubeAxis(2,k));
				glVertex3f(tubeAxis(0,k)+burgers(0)*10.0*radius,tubeAxis(1,k)+burgers(1)*10.0*radius,tubeAxis(2,k)+burgers(2)*10.0*radius);
				glEnd();
			}
			
		}

	};
	
	
	
	/*************************************************************/
	/*************************************************************/
	template <int dim, int Np, int Nc, float & alpha>
	class SplinePlotter :
	/* inherits from   */ public VertexReader<'V',11,double>, // CHANGE THIS DOUBLE TO SCALARTYPE
	/* inherits from   */ public EdgeReader  <'E',11,double>,
	/*                 */ public VertexReader<'P',7,double>,
	//	/* inherits from   */ private std::vector<SingleSplinePlotter<dim,Np,Nc,alpha> >{ // CHANGE THIS DOUBLE TO SCALARTYPE
	/* inherits from   */ private boost::ptr_vector<SingleSplinePlotter<dim,Np,Nc,alpha> >{ // ptr_vector::push_back doesn't use copy constructor so creation of SingleSplinePlotter will be faster // CHANGE THIS DOUBLE TO SCALARTYPE
		
		typedef float scalarType;
		typedef model::VertexReader<'V',11,double> VertexContainerType; // CHANGE THIS DOUBLE TO SCALARTYPE
		typedef model::EdgeReader  <'E',11,double>	EdgeContainerType; // CHANGE THIS DOUBLE TO SCALARTYPE
		typedef SingleSplinePlotter<dim,Np,Nc,alpha> SingleSplinePlotterType;
        //		typedef std::vector<SingleSplinePlotterType> SingleSplinePlotterVectorType;
		typedef boost::ptr_vector<SingleSplinePlotterType> SingleSplinePlotterVectorType;
		typedef typename SingleSplinePlotterType::VectorDim VectorDim;
        
        
        std::set<int> SIDs; // use std::set to automatically sort sID's

        
	public:
		
		bool showTubes;
		bool showVertices;
		bool deformedConfig;
		bool showPlaneNormal;
		bool showBurgers;
        bool showVertexID;
        int colorScheme;
        bool showSpecificVertex;
        int specificVertexID;
        bool showPK;
        double PKfactor;
		
        /* Constructor ***********************************************/
		SplinePlotter() : showTubes(false),
		/* init list   */ showVertices(false),
		/* init list   */ deformedConfig(false),
		/* init list   */ showPlaneNormal(false),
		/* init list   */ showBurgers(false),
        /* init list   */ showVertexID(false),
        /* init list   */ colorScheme(0),
        /* init list   */ showSpecificVertex(false),
        /* init list   */ specificVertexID(0),
        /* init list   */ showPK(false),
        /* init list   */ PKfactor(1000.0){}
		
		/* isGood ***************************************************/
		static bool isGood(const int& frameN)
        {
			return VertexContainerType::isGood(frameN) && EdgeContainerType::isGood<true>(frameN);
		}
		
		/* read ***************************************************/
		void read(const int& frameN)
        {
			VertexContainerType::read(frameN);
			EdgeContainerType::read<true>(frameN);
			VertexReader<'P',7,double>::read(frameN);
            
			SingleSplinePlotterVectorType::clear(); // clear the current content of sspVector
			SingleSplinePlotterVectorType::reserve(EdgeContainerType::size()); // reserve to speed-up push_back
			for (EdgeContainerType::const_iterator itEdge=EdgeContainerType::begin(); itEdge !=EdgeContainerType::end(); ++itEdge)
            {
				VertexContainerType::const_iterator itSource(VertexContainerType::find(itEdge->first.first)); //source
				assert(itSource!=VertexContainerType::end() && "SOURCE VERTEX NOT FOUND IN V-FILE");
				VertexContainerType::const_iterator itSink(VertexContainerType::find(itEdge->first.second)); //sink
				assert(  itSink!=VertexContainerType::end() &&   "SINK VERTEX NOT FOUND IN V-FILE");
				
				
				Eigen::Matrix<scalarType,dim,6> P0T0P1T1BN;
				
				int sourceTfactor(itEdge->second(2*dim));
				int   sinkTfactor(itEdge->second(2*dim+1));
				
				if(deformedConfig)
                {
					P0T0P1T1BN.col(0) = itSource->second.segment<dim>(7).transpose().template cast<float>();	// source position
					P0T0P1T1BN.col(2) =   itSink->second.segment<dim>(7).transpose().template cast<float>();	// sink position
				}
				else
                {
					P0T0P1T1BN.col(0) = itSource->second.segment<dim>(0*dim).transpose().template cast<float>();	// source position
					P0T0P1T1BN.col(2) =   itSink->second.segment<dim>(0*dim).transpose().template cast<float>();	// sink position
				}
				P0T0P1T1BN.col(1) = sourceTfactor*(itSource->second.segment<dim>(1*dim).transpose().template cast<float>());	// source tangent
				P0T0P1T1BN.col(3) =  -sinkTfactor*(  itSink->second.segment<dim>(1*dim).transpose().template cast<float>());	// sink tangent
				P0T0P1T1BN.col(4) = itEdge->second.segment<dim>(0*dim).transpose().template cast<float>();		// Burgers vector
				P0T0P1T1BN.col(5) = itEdge->second.segment<dim>(1*dim).transpose().template cast<float>();		// plane normal
				
				SingleSplinePlotterVectorType::push_back(new SingleSplinePlotterType(P0T0P1T1BN));
			}
		}
		
		
        //		/* plot ******************************************************/
        //		void renderBitmapString(const VectorDim& P,	void *font, const std::string& string2render) const {
        //  			glRasterPos3f(P(0), P(1), P(2));
        //			glColor3f(1.0f, 1.0f, 1.0f);
        //  			for (unsigned int k=0; k<string2render.length(); ++k) {
        //    			glutBitmapCharacter(font, string2render[k]);
        //  			}
        //		}
        
		/* plot ******************************************************/
		void plot(const scalarType& radius)
        {
			
			for (typename SingleSplinePlotterVectorType::const_iterator itEdge=SingleSplinePlotterVectorType::begin(); itEdge!=SingleSplinePlotterVectorType::end(); ++itEdge)
            {
				itEdge->plot(radius,showTubes,showPlaneNormal,showBurgers,colorScheme);
			}
            
            
            SIDs.clear();
            for (typename VertexContainerType::const_iterator vIter=VertexContainerType::begin();vIter!=VertexContainerType::end();++vIter)
            {
                SIDs.insert(vIter->second(6));
            }

			
			if(showVertices)
            {
				// Find color range based on sID
				//std::set<int> SIDs; // use std::set to automatically sort sID's
				float sIDmax = *SIDs.rbegin();
				float sIDmin = *SIDs.begin();
				
				// Loop and plot spheres
				GLUquadric* myQuad;
				myQuad=gluNewQuadric();
				for (typename VertexContainerType::const_iterator vIter=VertexContainerType::begin();vIter!=VertexContainerType::end();++vIter)
                {
                    //				float snID=(vIter->second(6)-sIDmin)/(sIDmax-sIDmin+1.0);
                    RGBcolor rgb(RGBmap::getColor(vIter->second(6),sIDmin,sIDmax));
                    GLfloat materialColor[] = {rgb.r, rgb.g, rgb.b, 1.0};
                    
                    //					GLfloat materialColor[] = {0.8*snID, 0.6*snID, 0.2*snID, 1.0};
					glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, materialColor);
					glTranslatef(  vIter->second(0),  vIter->second(1),  vIter->second(2) );
					gluSphere( myQuad , radius*1.2 , 10 , 10 );
					glTranslatef( -vIter->second(0), -vIter->second(1), -vIter->second(2) );
                    //					renderBitmapString(vIter->second.template cast<float>().template segment<dim>(0),
                    //               /*              */ GLUT_BITMAP_HELVETICA_18,
                    //					/*              */ static_cast<std::ostringstream*>( &(std::ostringstream() << vIter->first) )->str());
                    if (showVertexID || (showSpecificVertex && specificVertexID==vIter->first))
                    {
                        VectorDim PT(vIter->second.template cast<float>().template segment<dim>(0));
                        BitmapPlotter::renderString(PT,
                                                    /*                       */ static_cast<std::ostringstream*>( &(std::ostringstream() << vIter->first) )->str());
                    }
                    
                    
				}
				gluDeleteQuadric(myQuad); // free myQuad pointer
			}
            
            if (showPK){
                for (typename VertexReader<'P',7,double>::const_iterator vIter=VertexReader<'P',7,double>::begin();vIter!=VertexReader<'P',7,double>::end();++vIter){
                    //std::cout<<"I'm here"<< 	PKfactor<<std::endl;
                    glBegin(GL_LINES);
                    glVertex3f(vIter->second(0),vIter->second(1),vIter->second(2));
                    glVertex3f(vIter->second(0)+vIter->second(3)*PKfactor,vIter->second(1)+vIter->second(4)*PKfactor,vIter->second(2)+vIter->second(5)*PKfactor);
                    glEnd();
                }
                
            }
            
            
		}
		
		/* boundingBox ******************************************************/
		void boundingBox(float& xMin, float& xMax, float& yMin, float& yMax, float& zMin, float& zMax) const
        {
			xMin=0.0f;
			xMax=0.0f;
			yMin=0.0f;
			yMax=0.0f;
			zMin=0.0f;
			zMax=0.0f;
			for (typename VertexContainerType::const_iterator vIter=VertexContainerType::begin();vIter!=VertexContainerType::end();++vIter)
            {
				xMax= (xMax>vIter->second(0))? xMax : vIter->second(0);
				xMin= (xMin<vIter->second(0))? xMin : vIter->second(0);
				yMax= (yMax>vIter->second(1))? yMax : vIter->second(1);
				yMin= (yMin<vIter->second(1))? yMin : vIter->second(1);
				zMax= (zMax>vIter->second(2))? zMax : vIter->second(2);
				zMin= (zMin<vIter->second(2))? zMin : vIter->second(2);
			}
		}
        
        
        /* vertexOrder ******************************************************/
        unsigned int vertexOrder() const
        {
            return VertexReader<'V',11,double>::size();
        }
        
        /* vertexOrder ******************************************************/
        unsigned int edgeOrder() const
        {
            return EdgeReader<'E',11,double>::size();
        }
        
        /* vertexOrder ******************************************************/
        unsigned int components() const
        {
            return SIDs.size();
        }
		
	};
	
}
#endif
/*********************************************************************/
/*********************************************************************/







