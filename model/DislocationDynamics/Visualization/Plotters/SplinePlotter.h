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
#include <set>
#include <iterator> // std::distance

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>
#include <model/Utilities/IDreader.h>


#include <model/openGL/BitmapPlotter.h>
#include <model/openGL/RGBmap.h>



namespace model
{
	
	/*********************************************************************/
	/*********************************************************************/
	template <int dim, int Np, int Nc>
	class SingleSplinePlotter
    {
		
	public:
        
        static float alpha;
        static bool use_BurgersNorm;
        
		typedef float scalarType;
		typedef Eigen::Matrix<scalarType,dim,Np> MatrixDimNp;
		typedef Eigen::Matrix<scalarType,dim,Nc> MatrixDimNc;
		typedef Eigen::Matrix<scalarType,dim,1>  VectorDim;
        
        
        const int snID;
        
		
	private:
		
		MatrixDimNp tubeAxis;
		MatrixDimNp tubeTangents;
		std::vector<MatrixDimNc> tubeCircles;
		
        VectorDim planeNormal;
		VectorDim burgers;
        scalarType burgersNorm;
		VectorDim chord;
        
        const bool isSessile;
		
		
		/* getCircle **************************************************************/
		MatrixDimNc getCircle(const int& k) const {
			MatrixDimNc circle;
			for (int c=0;c<Nc;++c){
				circle.col(c)=  Eigen::AngleAxisf(2.0f * M_PI/Nc*c,tubeTangents.col(k)) * tubeTangents.col(k).cross(planeNormal);
			}
			return circle;
		}
		
	public:
		
        enum{colorBurgers=0,colorSessile=1,colorNormal=2,colorEdgeScrew=3,colorComponent=4};
        
        
		/* Constructor ************************************************************/
		SingleSplinePlotter(const Eigen::Matrix<scalarType,dim,6>& P0T0P1T1BN, const int& snID_in) :
        /* init list */ snID(snID_in),
		/* init list */ planeNormal(P0T0P1T1BN.col(5).normalized()),
		/* init list */ burgers(P0T0P1T1BN.col(4)),
        /* init list */ burgersNorm(use_BurgersNorm? burgers.norm() : 1.0),
        /* init list */ isSessile(std::fabs(planeNormal.dot(burgers.normalized()))>FLT_EPSILON)
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
		
        /*********************************************************************/
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
		void plot(const scalarType& radius, const bool& showTubes, const bool& showPlaneNormal, const bool& showBurgers, const int& colorScheme,const int& ids, const int& sIDmin,const int& sIDmax) const {
			
            VectorDim screwColor=VectorDim::UnitY();
            VectorDim edgeColor=VectorDim::UnitZ();
            
            // 1- Define the color
            VectorDim colorVector;
//            switch (colorScheme)
//            {
//                case colorSessile:
//                    colorVector(0)= isSessile? 1.0 : 0.1;
//                    colorVector(1)= isSessile? 0.5 : 0.4;
//                    colorVector(2)= isSessile? 0.0 : 0.9;
//                    break;
//                    
//                case colorNormal:
//                    colorVector = planeNormal;
//                    flipColor(colorVector);
//                    break;
//                    
//                case colorEdgeScrew:
//                    colorVector = planeNormal;
//                    flipColor(colorVector);
//                    break;
//                    
//                case colorComponent:
//                {
//                    RGBcolor rgb(RGBmap::getColor(ids,sIDmin,sIDmax));
//                    colorVector << rgb.r, rgb.g, rgb.b;
//                }
//                    break;
//                    
//                default:
//                    colorVector = burgers.normalized();
//                    flipColor(colorVector);
//                    break;
//            }
//            
//            //			glDisable(GL_COLOR_MATERIAL); // use glMaterialfv(...) to set material colors
//            //			glEnable(GL_DEPTH_TEST);
//            GLfloat materialAmbient[]={colorVector(0), colorVector(1), colorVector(2), 1.0};
//            GLfloat materialEmission[]={colorVector(0)*0.1f, colorVector(1)*0.1f, colorVector(2)*0.1f, 1.0};
//
//            //			GLfloat materialSpecular[] = {specularity, specularity, specularity, 1.0f};
//            //			GLfloat materialEmission[] = {emissivity, emissivity, emissivity, 1.0f};
//            // note that glMaterialfv(...) works when GL_COLOR_MATERIAL is disabled
//			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialAmbient);      // ambient color for the material
//			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialAmbient);      // diffuse color for the material
//			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialAmbient);  // specular color for the material
//            glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialEmission);  // emission color for the material

            //			glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialEmission);  // emission color for the material
            //			glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess); //The shininess parameter
			
			if (showTubes)
            {
				for (int k=1;k<Np;++k)
                {
                    switch (colorScheme)
                    {
                        case colorSessile:
                            colorVector(0)= isSessile? 1.0 : 0.1;
                            colorVector(1)= isSessile? 0.5 : 0.4;
                            colorVector(2)= isSessile? 0.0 : 0.9;
                            break;
                            
                        case colorNormal:
                            colorVector = planeNormal;
                            flipColor(colorVector);
                            break;
                            
                        case colorComponent:
                        {
                            RGBcolor rgb(RGBmap::getColor(ids,sIDmin,sIDmax));
                            colorVector << rgb.r, rgb.g, rgb.b;
                        }
                            break;
                            
                        case colorEdgeScrew:
                        {
                            const float u = std::fabs(tubeTangents.col(k).normalized().dot(burgers.normalized()));
//                            RGBcolor rgb(RGBmap::getColor(std::fabs(tubeTangents.col(k).normalized().dot(burgers.normalized())),0,1));
//                            colorVector << rgb.r, rgb.g, rgb.b;
                            colorVector=screwColor*u+edgeColor*(1-u);
                        }
                            break;
                            
                        default:
                            colorVector = burgers.normalized();
                            flipColor(colorVector);
                            break;
                    }
                    
                    //			glDisable(GL_COLOR_MATERIAL); // use glMaterialfv(...) to set material colors
                    //			glEnable(GL_DEPTH_TEST);
                    GLfloat materialAmbient[]={colorVector(0), colorVector(1), colorVector(2), 1.0};
                    GLfloat materialEmission[]={colorVector(0)*0.1f, colorVector(1)*0.1f, colorVector(2)*0.1f, 1.0};
                    
                    //			GLfloat materialSpecular[] = {specularity, specularity, specularity, 1.0f};
                    //			GLfloat materialEmission[] = {emissivity, emissivity, emissivity, 1.0f};
                    // note that glMaterialfv(...) works when GL_COLOR_MATERIAL is disabled
                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialAmbient);      // ambient color for the material
                    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialAmbient);      // diffuse color for the material
                    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialAmbient);  // specular color for the material
                    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialEmission);  // emission color for the material
                    
					//					tubeCircles[k]=getCircle(k);
					glBegin(GL_TRIANGLE_STRIP);
					//glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
					for (int c=0;c<Nc;++c){
						glNormal3f(tubeCircles[k-1](0,c),tubeCircles[k-1](1,c),tubeCircles[k-1](2,c));
						glVertex3f(tubeAxis(0,k-1)+burgersNorm*radius*tubeCircles[k-1](0,c),tubeAxis(1,k-1)+burgersNorm*radius*tubeCircles[k-1](1,c),tubeAxis(2,k-1)+burgersNorm*radius*tubeCircles[k-1](2,c));
						glNormal3f(tubeCircles[k](0,c),tubeCircles[k](1,c),tubeCircles[k](2,c));
						//if(showTubes==1){
						glVertex3f(tubeAxis(0,k)+burgersNorm*radius*tubeCircles[k](0,c),tubeAxis(1,k)+burgersNorm*radius*tubeCircles[k](1,c),tubeAxis(2,k)+burgersNorm*radius*tubeCircles[k](2,c));
						//}
						//else{
						//	glVertex3f(tubeAxis(0,k),tubeAxis(1,k),tubeAxis(2,k));
						//}
					}
					glNormal3f(tubeCircles[k-1](0,0),tubeCircles[k-1](1,0),tubeCircles[k-1](2,0));
					glVertex3f(tubeAxis(0,k-1)+burgersNorm*radius*tubeCircles[k-1](0,0),tubeAxis(1,k-1)+burgersNorm*radius*tubeCircles[k-1](1,0),tubeAxis(2,k-1)+burgersNorm*radius*tubeCircles[k-1](2,0));
					glNormal3f(tubeCircles[k](0,0),tubeCircles[k](1,0),tubeCircles[k](2,0));
					//if(showTubes==1){
					glVertex3f(tubeAxis(0,k)+burgersNorm*radius*tubeCircles[k](0,0),tubeAxis(1,k)+burgersNorm*radius*tubeCircles[k](1,0),tubeAxis(2,k)+burgersNorm*radius*tubeCircles[k](2,0));
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
                    switch (colorScheme)
                    {
                        case colorSessile:
                            colorVector(0)= isSessile? 1.0 : 0.1;
                            colorVector(1)= isSessile? 0.5 : 0.4;
                            colorVector(2)= isSessile? 0.0 : 0.9;
                            break;
                            
                        case colorNormal:
                            colorVector = planeNormal;
                            flipColor(colorVector);
                            break;
                            
                        case colorComponent:
                        {
                            RGBcolor rgb(RGBmap::getColor(ids,sIDmin,sIDmax));
                            colorVector << rgb.r, rgb.g, rgb.b;
                        }
                            break;
                            
                        case colorEdgeScrew:
                        {
                            RGBcolor rgb(RGBmap::getColor(std::fabs(tubeTangents.col(k).normalized().dot(burgers.normalized())),0,1));
                            colorVector << rgb.r, rgb.g, rgb.b;
                        }
                            break;
                            
                        default:
                            colorVector = burgers.normalized();
                            flipColor(colorVector);
                            break;
                    }
                    
                    //			glDisable(GL_COLOR_MATERIAL); // use glMaterialfv(...) to set material colors
                    //			glEnable(GL_DEPTH_TEST);
                    GLfloat materialAmbient[]={colorVector(0), colorVector(1), colorVector(2), 1.0};
                    GLfloat materialEmission[]={colorVector(0)*0.1f, colorVector(1)*0.1f, colorVector(2)*0.1f, 1.0};
                    
                    //			GLfloat materialSpecular[] = {specularity, specularity, specularity, 1.0f};
                    //			GLfloat materialEmission[] = {emissivity, emissivity, emissivity, 1.0f};
                    // note that glMaterialfv(...) works when GL_COLOR_MATERIAL is disabled
                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialAmbient);      // ambient color for the material
                    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialAmbient);      // diffuse color for the material
                    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialAmbient);  // specular color for the material
                    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialEmission);  // emission color for the material
                    
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
    
    template <int dim, int Np, int Nc>
	float SingleSplinePlotter<dim,Np,Nc>::alpha=0.5;

    template <int dim, int Np, int Nc>
    bool SingleSplinePlotter<dim,Np,Nc>::use_BurgersNorm=false;
	
	
	/*************************************************************/
	/*************************************************************/
	template <int dim, int Np, int Nc>
	class SplinePlotter :
	/* inherits from   */ public VertexReader<'V',10,double>, // CHANGE THIS DOUBLE TO SCALARTYPE
	/* inherits from   */ public EdgeReader  <'E',11,double>,
	/*                 */ public IDreader<'P',3,6,double>,
    /*                 */ public IDreader<'Y',1,3,double>,
    /*                 */ public VertexReader<'L',10,double>,
    /* inherits from   */ public IDreader<'Q',3,13,double>,
    /* inherits from   */ private std::vector<SingleSplinePlotter<dim,Np,Nc> >
    { // ptr_vector::push_back doesn't use copy constructor so creation of SingleSplinePlotter will be faster // CHANGE THIS DOUBLE TO SCALARTYPE
		
		typedef float scalarType;
		typedef VertexReader<'V',10,double> VertexContainerType; // CHANGE THIS DOUBLE TO SCALARTYPE
		typedef EdgeReader  <'E',11,double>	EdgeContainerType; // CHANGE THIS DOUBLE TO SCALARTYPE
        typedef IDreader<'P',3,6,double> PKContainerType;
        typedef IDreader<'Q',3,13,double> QuadContainerType;

		typedef SingleSplinePlotter<dim,Np,Nc> SingleSplinePlotterType;
        typedef std::vector<SingleSplinePlotterType> SingleSplinePlotterVectorType;

        typedef typename SingleSplinePlotterType::VectorDim VectorDim;
        
        typedef VertexReader<'L',10,double> GBDreaderType;
        
        typedef IDreader<'Y',1,3,double> VelocityReader;
        
        std::set<int> SIDs; // use std::set to automatically sort sID's
        
        
	public:
        
		
		bool showTubes;
		bool showVertices;
        //		bool deformedConfig;
		bool showPlaneNormal;
		bool showBurgers;
        bool showVertexID;
        static int colorScheme;
        static bool plotBoundarySegments;
        static bool showQuadParticles;
        static bool showVelocity;

        bool showSpecificVertex;
        int specificVertexID;
        bool showPK;
        double PKfactor;
        
        const QuadContainerType& quadContainer() const
        {
            return *this;
        }
		
        /* Constructor ********************************************************/
		SplinePlotter() : showTubes(false),
		/* init list   */ showVertices(false),
        //		/* init list   */ deformedConfig(false),
		/* init list   */ showPlaneNormal(false),
		/* init list   */ showBurgers(false),
        /* init list   */ showVertexID(false),
        /* init list   */ showSpecificVertex(false),
        /* init list   */ specificVertexID(0),
        /* init list   */ showPK(false),
        /* init list   */ PKfactor(1000.0)
        {
        
            GBDreaderType::read(0,true);
            
        }
		
		/* isGood *************************************************************/
		static bool isGood(const int& frameN, const bool& useTXT)
        {
			return VertexContainerType::isGood(frameN,useTXT) && EdgeContainerType::isGood(frameN,useTXT);
		}
		
		/* read ***************************************************************/
		void read(const int& frameN)
        {
            if (isGood(frameN,false)) // bin format
            {
                VertexContainerType::read(frameN,false);
                EdgeContainerType::read(frameN,false);
            }
            else // txt format
            {
                VertexContainerType::read(frameN,true);
                EdgeContainerType::read(frameN,true);
            }
            
            if(showQuadParticles) // Show quadrature particles
            {
                QuadContainerType::read(frameN,true);
            }
            
            if(showVelocity)
            {
                VelocityReader::read(frameN,true);
            }
            
			PKContainerType::read(frameN,true);
            
			SingleSplinePlotterVectorType::clear(); // clear the current content of sspVector
			SingleSplinePlotterVectorType::reserve(EdgeContainerType::size()); // reserve to speed-up push_back
			for (EdgeContainerType::const_iterator itEdge=EdgeContainerType::begin(); itEdge !=EdgeContainerType::end(); ++itEdge)
            {
				VertexContainerType::const_iterator itSource(VertexContainerType::find(itEdge->first.first)); //source
				assert(itSource!=VertexContainerType::end() && "SOURCE VERTEX NOT FOUND IN V-FILE");
				VertexContainerType::const_iterator itSink(VertexContainerType::find(itEdge->first.second)); //sink
				assert(  itSink!=VertexContainerType::end() &&   "SINK VERTEX NOT FOUND IN V-FILE");
				
				
				Eigen::Matrix<scalarType,dim,6> P0T0P1T1BN;
				
				const int sourceTfactor(itEdge->second(2*dim));
				const int   sinkTfactor(itEdge->second(2*dim+1));
				const int   snID(itEdge->second(2*dim+2));
                const bool sourceOnBoundary(itSource->second(2*dim+1));
                const bool   sinkOnBoundary(  itSink->second(2*dim+1));
                
                if(!(sourceOnBoundary && sinkOnBoundary) || plotBoundarySegments)
                {
                    
                    P0T0P1T1BN.col(0) = itSource->second.segment<dim>(0*dim).transpose().template cast<float>();	// source position
                    P0T0P1T1BN.col(2) =   itSink->second.segment<dim>(0*dim).transpose().template cast<float>();	// sink position
                    P0T0P1T1BN.col(1) = sourceTfactor*(itSource->second.segment<dim>(1*dim).transpose().template cast<float>());	// source tangent
                    P0T0P1T1BN.col(3) =  -sinkTfactor*(  itSink->second.segment<dim>(1*dim).transpose().template cast<float>());	// sink tangent
                    P0T0P1T1BN.col(4) = itEdge->second.segment<dim>(0*dim).transpose().template cast<float>();		// Burgers vector
                    P0T0P1T1BN.col(5) = itEdge->second.segment<dim>(1*dim).transpose().template cast<float>();		// plane normal

                    SingleSplinePlotterVectorType::emplace_back(P0T0P1T1BN,snID);
                }
			}
            
            Eigen::Matrix<scalarType,dim,6> P0T0P1T1BN(Eigen::Matrix<scalarType,dim,6>::Zero());

            for (const auto& gbd : GBdislocations())
            {
                //const Eigen::Matrix<double,1,3*dim> temp=gbd.second;
                
                P0T0P1T1BN.col(0) = gbd.second.template segment<dim>(0*dim).transpose().template cast<float>();	// source position
                P0T0P1T1BN.col(2) = gbd.second.template segment<dim>(1*dim).transpose().template cast<float>();	// sink position
//                P0T0P1T1BN.col(1) = P0T0P1T1BN.col(2)-P0T0P1T1BN.col(0);	// source tangent
//                P0T0P1T1BN.col(3) = P0T0P1T1BN.col(1);	// sink tangent
                P0T0P1T1BN.col(4) = gbd.second.template segment<dim>(2*dim).transpose().template cast<float>();		// Burgers vector
                P0T0P1T1BN.col(5) = P0T0P1T1BN.col(4);		// plane normal
                
                SingleSplinePlotterVectorType::emplace_back(P0T0P1T1BN,0);
            }
		}
		
        /**********************************************************************/
        const GBDreaderType& GBdislocations() const
        {
            return *this;
        }
        
		/* plot ***************************************************************/
		void plot(const scalarType& radius)
        {
            
            // Collect IDs of Network Components
            SIDs.clear();
            for (typename VertexContainerType::const_iterator vIter=VertexContainerType::begin();vIter!=VertexContainerType::end();++vIter)
            {
                SIDs.insert(vIter->second(6));
            }
            float sIDmax(*SIDs.rbegin());
            float sIDmin(*SIDs.begin());
            
            // Plot Segments
            glDisable(GL_COLOR_MATERIAL); // use glMaterialfv(...) to set material colors
            glEnable(GL_DEPTH_TEST);
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100.0f); //The shininess parameter
			for (typename SingleSplinePlotterVectorType::const_iterator itEdge=SingleSplinePlotterVectorType::begin(); itEdge!=SingleSplinePlotterVectorType::end(); ++itEdge)
            {
				itEdge->plot(radius,showTubes,showPlaneNormal,showBurgers,colorScheme,std::distance(SIDs.begin(),SIDs.find(itEdge->snID)),0,SIDs.size());
			}
			
            // Plot Nodes
			if(showVertices) // Show vertices
            {
                GLfloat materialColor[] = {0.0, 0.0, 0.0, 1.0};
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialColor);      // ambient color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialColor);      // diffuse color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialColor);  // specular color for the material

                for (typename VertexContainerType::const_iterator vIter=VertexContainerType::begin();vIter!=VertexContainerType::end();++vIter)
                {
                if(showVelocity)
                {
                    float vf=10.0;
                    const auto velIter=VelocityReader::find(std::array<int, 1>{ {vIter->first}});
                    if(velIter!=VelocityReader::end())
                    {
//                        std::cout<<vIter->first<<" "<<vIter->second[0]<<" "<<vIter->second[1]<<" "<<vIter->second[2]<<std::endl;
//                        
//                        std::cout<<vIter->first<<" "<<velIter->second[0]<<" "<<velIter->second[1]<<" "<<velIter->second[2]<<std::endl;
                        
                        glColor3f(1.0f, 0.0f, 0.0f);
                        
                        glBegin(GL_LINES);
                        glVertex3f(vIter->second[0],vIter->second[1],vIter->second[2]);
                        glVertex3f(vIter->second[0]+velIter->second[0]*vf,vIter->second[1]+velIter->second[1]*vf,vIter->second[2]+velIter->second[2]*vf);
                        glEnd();
                    }
                  //  else
                  //  {
                        //std::cout<<vIter->first<<" not found"<<std::endl;
                   // }
                }
                }
                
				// Loop and plot spheres
				GLUquadric* myQuad;
				myQuad=gluNewQuadric();
				for (typename VertexContainerType::const_iterator vIter=VertexContainerType::begin();vIter!=VertexContainerType::end();++vIter)
                {
                    

                    
                    
                    GLfloat materialColor[] = {0.0, 0.0, 0.0, 1.0};
                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialColor);      // ambient color for the material
                    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialColor);      // diffuse color for the material
                    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialColor);  // specular color for the material
                    const bool vertexOnBoundary(vIter->second(2*dim+1));
                    
                    if(!vertexOnBoundary || plotBoundarySegments)
                    {
                        glColor3f(0.0f, 0.0f, 0.0f);
                        glTranslatef(  vIter->second(0),  vIter->second(1),  vIter->second(2) );
                        gluSphere( myQuad , radius*1.2 , 10 , 10 );
                        glTranslatef( -vIter->second(0), -vIter->second(1), -vIter->second(2) );
                        
                        
                        if (showVertexID || (showSpecificVertex && specificVertexID==vIter->first))
                        {
                            VectorDim PT(vIter->second.template cast<float>().template segment<dim>(0));
                            BitmapPlotter::renderString(PT,
                                                        /*                       */ static_cast<std::ostringstream*>( &(std::ostringstream() << vIter->first) )->str());
                        }
                        
                    }
                    

                    

                    
                    
				}
				gluDeleteQuadric(myQuad); // free myQuad pointer
			}
            
            if(showQuadParticles) // Show QuadParticles
            {


                // Loop and plot spheres
                GLUquadric* myQuad;
                myQuad=gluNewQuadric();
                for(const auto& quad : quadContainer())
                {
                    GLfloat materialColor[] = {0.0, 1.0, quad.second[8], 1.0};
                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialColor);      // ambient color for the material
                    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialColor);      // diffuse color for the material
                    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialColor);  // specular color for the material
                    glTranslatef( quad.second[0], quad.second[1], quad.second[2] );
                    gluSphere( myQuad , radius*1.2 , 10 , 10 );
                    glTranslatef(-quad.second[0],-quad.second[1],-quad.second[2] );
                }
                gluDeleteQuadric(myQuad); // free myQuad pointer
            }
            
            if (showPK) // Show PK force
            {
                for (typename PKContainerType::const_iterator vIter=PKContainerType::begin();vIter!=PKContainerType::end();++vIter)
                {
                    glBegin(GL_LINES);
                    glVertex3f(vIter->second[0],vIter->second[1],vIter->second[2]);
                    glVertex3f(vIter->second[0]+vIter->second[3]*PKfactor,vIter->second[1]+vIter->second[4]*PKfactor,vIter->second[2]+vIter->second[5]*PKfactor);
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
            return VertexContainerType::size();
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
	
    
    //static data
    template <int dim, int Np, int Nc>
	int SplinePlotter<dim,Np,Nc>::colorScheme=0;

    template <int dim, int Np, int Nc>
    bool SplinePlotter<dim,Np,Nc>::plotBoundarySegments=false;
    
    template <int dim, int Np, int Nc>
    bool SplinePlotter<dim,Np,Nc>::showQuadParticles=false;

    
    template <int dim, int Np, int Nc>
    bool SplinePlotter<dim,Np,Nc>::showVelocity=true;


}
#endif
