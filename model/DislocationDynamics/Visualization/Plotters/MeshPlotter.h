/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MESHPLOTTER_H_
#define model_MESHPLOTTER_H_


//#ifdef __APPLE__
//#include <OpenGL/OpenGL.h>
//#include <GLUT/glut.h> // for quad
//#else
//#include <GL/glut.h> // for quad
//#endif

#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#include <GL/glx.h>
#endif

#include <float.h>
#include <vector>
#include <deque>
#include <map>

#include <Eigen/Core>

#include <model/Network/Readers/VertexReader.h>
#include <model/Network/Readers/EdgeReader.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/Utilities/IDreader.h>


namespace model
{
    
    class MeshPlotter :
    /*               */ public VertexReader<'D',4,float>,
    /*               */ public IDreader<'I',1,3,int>
    {
        
        //		typedef VertexReader<'N',5,float> NodeContainerType;
        typedef VertexReader<'D',4,float> DispContainerType;
        typedef IDreader<'I',1,3,int> SpecialTrianglesContainerType;

        //		typedef VertexReader<'Q',4,float> QuadContainerType;
        //		typedef   EdgeReader<'T',3,float> EdgeContainerType;
        
        //		bool edgeFileIsGood;
        //		bool nodeFileIsGood;
        bool dispFileIsGood;
        bool specialTrianglesFileIsGood;

        //		bool quadFileIsGood;
        
        //        typedef std::vector<Eigen::Matrix<float,3,4>,Eigen::aligned_allocator<Eigen::Matrix<float,3,4> > > EdgeVectoType;
        
        typedef std::map<std::pair<size_t,size_t>,Eigen::Matrix<float,3,2> > EdgeVectoType;
        EdgeVectoType edgeVector; // this is [P0 P1 D0 D1]
        
        typedef Simplex<3,3>::SimplexIDType SimplexIDType;
        
        //        SimplicialMesh<3> mesh;
        
        
        std::deque<std::pair<Eigen::Matrix<float,3,1>,Eigen::Matrix<float,6,1>>> deq;
        
        std::deque<Eigen::Matrix<float,3,3> > regionsBndDeq;
        std::deque<std::pair<int,int>> regionsBndIDDeq;
        
        
        std::deque<Eigen::Matrix<float,3,3> > specialTrianglesDeq;

        
        Eigen::Matrix<float,6,1> sMin;
        Eigen::Matrix<float,6,1> sMax;
        
        int rMax;
        
        const SpecialTrianglesContainerType& specialTriangles() const
        {
            return *this;
        }
        
    public:
        
        const SimplicialMesh<3>* const p_mesh;
        
        
        enum {showMeshStates=3};
        short unsigned int showMesh;
        
        //        bool showQuad;
        float dispScale;
        
        static bool plotBndStress;
        static unsigned int stressCol;
        static bool showRegionBoundaries;
        static bool showSpecificSimplex;
        static SimplexIDType specificSimplexID;
        //        SimplicialMesh<3> mesh;
        
        /* Constructor ******************************************/
        MeshPlotter(const SimplicialMesh<3>* const p_mesh_in) :
        /* init list */ dispFileIsGood(false),
        /* init list */ specialTrianglesFileIsGood(false),
        /* init list */ p_mesh(p_mesh_in),
        //        /* init list */ edgeFileIsGood(false),
        //		/* init list */ nodeFileIsGood(false),
        rMax(0),
        /* init list */ showMesh(1),
        //        /* init list */ showQuad(false),
        /* init list */ dispScale(1.0f)
        //        /* init list */ mesh(1) // read N/N_0.txt and T/T_0.txt
        //        /* init list */ p_mesh(new SimplicialMesh<3>(1)) // read N/N_0.txt and T/T_0.txt
        {
            
            edgeVector.clear();
            //            edgeVector.reserve(SimplexObserver<3,1>::size()); // use reserve to speed-up push_back used later
//            for (typename SimplexObserver<3,1>::const_iterator sIter=SimplexObserver<3,1>::simplexBegin();
//                 /*                                         */ sIter!=SimplexObserver<3,1>::simplexEnd();++sIter)
//            {
            for (const auto& sIter : p_mesh->template observer<1>())
            {
            if(sIter.second->isBoundarySimplex())
                {
                    std::pair<size_t,size_t> key=std::make_pair(sIter.second->child(0).xID(0),
                                                                sIter.second->child(1).xID(0));
                    
                    
                    Eigen::Matrix<float,3,2> temp(Eigen::Matrix<float,3,2>::Zero());
                    temp.col(0)=sIter.second->child(0).P0.cast<float>();
                    temp.col(1)=sIter.second->child(1).P0.cast<float>();
                    
                    edgeVector.emplace(key,temp);
                    
                    //					if (dispFileIsGood)
                    //                    {
                    //						VertexReader<'D',4,float>::iterator iterD1(DispContainerType::find(sIter.second->child(0).xID(0)));
                    //						VertexReader<'D',4,float>::iterator iterD2(DispContainerType::find(sIter.second->child(1).xID(0)));
                    //						assert(iterD1!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                    //						assert(iterD2!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                    //						temp.col(2)=iterD1->second.segment<3>(0);
                    //						temp.col(3)=iterD2->second.segment<3>(0);
                    //					}
                    //                    edgeVector.push_back(temp);
                    
                }
            }
            
//            for (typename SimplexObserver<3,2>::const_iterator sIter=SimplexObserver<3,2>::simplexBegin();
//                 /*                                         */ sIter!=SimplexObserver<3,2>::simplexEnd();++sIter)
            for (const auto& sIter : p_mesh->template observer<2>())
            {
                if(sIter.second->isRegionBoundarySimplex())
                {
                    regionsBndDeq.emplace_back(sIter.second->vertexPositionMatrix().template cast<float>());
                    (*(sIter.second->parents().begin()))->region->regionID;
                    
                    int rID1=(*(sIter.second->parents(). begin()))->region->regionID;
                    int rID2=(*(sIter.second->parents().rbegin()))->region->regionID;
                    
                    if(rID1>rMax)
                    {
                        rMax=rID1;
                    }
                    if(rID2>rMax)
                    {
                        rMax=rID2;
                    }
                    
                    regionsBndIDDeq.emplace_back(rID1,rID2);
                }
            }
            
        }
        
        /* read *************************************************/
        void read(const int& frameN)
        {
            // Read edge file T only if it exists, otherwise try to read 0
            //			edgeFileIsGood=EdgeContainerType::isGood(frameN,true);
            //			if (edgeFileIsGood){
            //				EdgeContainerType::read(frameN,true);
            //			}
            //			else{
            //				EdgeContainerType::read(0,true);
            //			}
            
            // Read node file N only if it exists, otherwise try to read 0
            //			nodeFileIsGood=NodeContainerType::isGood(frameN,true);
            //			if (nodeFileIsGood){
            //				NodeContainerType::read(frameN,true);
            //			}
            //			else{
            //				NodeContainerType::read(0,true);
            //			}
            
            dispFileIsGood=DispContainerType::isGood(frameN,true);
            if (dispFileIsGood)
            {
                DispContainerType::read(frameN,true);
            }
            else
            {
                DispContainerType::read(0,true);
            }
            
            specialTrianglesDeq.clear();
            specialTrianglesFileIsGood=SpecialTrianglesContainerType::isGood(frameN,true);
            if(specialTrianglesFileIsGood)
            {
                float dispCorr(dispScale*(showMesh>1));

                
                SpecialTrianglesContainerType::read(frameN,true);
                
                for(const auto& triangleID : specialTriangles())
                {
                    const Eigen::Matrix<size_t,3,1> xID=(Eigen::Matrix<size_t,3,1>()<<triangleID.second[0],triangleID.second[1],triangleID.second[2]).finished();
                    
                    
                    Eigen::Matrix<float,3,3> vertexPositions=p_mesh->template observer<2>().simplex(xID).vertexPositionMatrix().template cast<float>();
                    
                    if(dispFileIsGood)
                    {
                        VertexReader<'D',4,float>::const_iterator iterD1(DispContainerType::find(triangleID.second[0]));
                        VertexReader<'D',4,float>::const_iterator iterD2(DispContainerType::find(triangleID.second[1]));
                        VertexReader<'D',4,float>::const_iterator iterD3(DispContainerType::find(triangleID.second[2]));

                        assert(iterD1!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                        assert(iterD2!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                        assert(iterD3!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");

                        vertexPositions.col(0)+=dispCorr*iterD1->second.segment<3>(0);
                        vertexPositions.col(1)+=dispCorr*iterD2->second.segment<3>(0);
                        vertexPositions.col(2)+=dispCorr*iterD3->second.segment<3>(0);
                    }

                    
                    
                    specialTrianglesDeq.emplace_back(vertexPositions);
//                    
//                    
//                    glBegin(GL_TRIANGLES);
//                    glVertex3f(mat.col(0)(0),mat.col(0)(1),mat.col(0)(2));
//                    glVertex3f(mat.col(1)(0),mat.col(1)(1),mat.col(1)(2));
//                    glVertex3f(mat.col(2)(0),mat.col(2)(1),mat.col(2)(2));
//                    glEnd();
                    
                }
                
                
            }
//            else
//            {
//                specialTrianglesDeq.clear();
//            }
            
            //            quadFileIsGood=QuadContainerType::isGood(frameN,true);
            //            if (quadFileIsGood){
            //				QuadContainerType::read(frameN,true);
            //
            //			}
            //			else{
            //				QuadContainerType::read(0,true);
            //			}
            //			assert(NodeContainerType==DispContainerType && "NUMBER OF NODES IN DISPLACEMENT FILE");
            
            
            //			edgeVector.clear();
            //			edgeVector.reserve(EdgeContainerType::size()); // use reserve to speed-up push_back used later
            //			for (EdgeContainerType::const_iterator iterE=EdgeContainerType::begin(); iterE!=EdgeContainerType::end();++iterE)
            //            {
            //				VertexReader<'N',5,float>::const_iterator iterN1 = NodeContainerType::find(iterE->first.first);
            //				VertexReader<'N',5,float>::const_iterator iterN2 = NodeContainerType::find(iterE->first.second);
            //				assert(iterN1!=NodeContainerType::end() && "MESH NODE NOT FOUND IN N FILE");
            //				assert(iterN2!=NodeContainerType::end() && "MESH NODE NOT FOUND IN N FILE");
            //				bool isBoundaryNode1(iterN1->second.operator()(3));
            //				bool isBoundaryNode2(iterN2->second.operator()(3));
            //
            //				if ( (isBoundaryNode1 && isBoundaryNode2) /*|| showInteriorBoundaryMesh*/ )
            //                {
            //					Eigen::Matrix<float,3,4> temp(Eigen::Matrix<float,3,4>::Zero());
            //					temp.col(0)=iterN1->second.segment<3>(0);
            //					temp.col(1)=iterN2->second.segment<3>(0);
            //					if (dispFileIsGood) {
            //						VertexReader<'D',4,float>::iterator iterD1(DispContainerType::find(iterE->first.first));
            //						VertexReader<'D',4,float>::iterator iterD2(DispContainerType::find(iterE->first.second));
            //						assert(iterD1!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
            //						assert(iterD2!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
            //						temp.col(2)=iterD1->second.segment<3>(0);
            //						temp.col(3)=iterD2->second.segment<3>(0);
            //					}
            //					edgeVector.push_back(temp);
            //				}
            //			}
            
            
            if (plotBndStress)
            {
                deq.clear();
                float x,y,z,s11,s22,s33,s12,s23,s13;
                
                std::ostringstream fullName;
                fullName<<"S/S_"<<frameN<<".txt";
                FILE *fp =fopen(fullName.str().c_str(), "r");
                
                if (fp!=NULL)
                {
                    
                    Eigen::Matrix<float,3,1> P;
                    Eigen::Matrix<float,6,1> S;
                    
                    while (fscanf (fp, "%f%f%f%f%f%f%f%f%f", &x,&y,&z,&s11,&s22,&s33,&s12,&s23,&s13)==9)
                    {
                        P<<x,y,z;
                        S<<s11,s22,s33,s12,s23,s13;
                        deq.emplace_back(P,S);
                    }
                    fclose(fp);
                    
                    sMin=deq[0].second;
                    sMax=deq[0].second;
                    
                    
                    for (int k=0;k<deq.size();++k)
                    {
                        for(int c=0;c<6;++c)
                        {
                            if(deq[k].second(c)>sMax(c))
                            {
                                sMax(c)=deq[k].second(c);
                            }
                            
                            if(deq[k].second(c)<sMin(c))
                            {
                                sMin(c)=deq[k].second(c);
                            }
                        }
                    }
                    
                    std::cout<<"sigma_min="<<sMin<<std::endl;
                    std::cout<<"sigma_max="<<sMax<<std::endl;
                    
                }
                
                std::cout<<fullName.str()<<" has "<<deq.size()<<"rows"<<std::endl;
                
                //                float sMin=FLT_MAX;
                //                float sMax=FLT_MIN;
                
                //                sMin=Eigen::Matrix<float,1,6>::Constant(FLT_MAX);
                //                sMax=Eigen::Matrix<float,1,6>::Constant(FLT_MIN);
                
            }
            
        }
        
        /* plot *************************************************/
        void plot() const
        {
            
            if (plotBndStress)
            {
                glEnable(GL_DEPTH_TEST);
                
                glEnable(GL_COLOR_MATERIAL); // use glColorMaterial(...) to set material colors
                //                glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
                
                glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT);
                glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
                glColorMaterial(GL_FRONT_AND_BACK, GL_SPECULAR);
                glColorMaterial(GL_FRONT_AND_BACK, GL_EMISSION);
                
                for (int k=0;k<deq.size()/3;++k)
                {
                    RGBcolor clr0=RGBmap::getColor(deq[3*k+0].second(stressCol),sMin(stressCol),sMax(stressCol));
                    RGBcolor clr1=RGBmap::getColor(deq[3*k+1].second(stressCol),sMin(stressCol),sMax(stressCol));
                    RGBcolor clr2=RGBmap::getColor(deq[3*k+2].second(stressCol),sMin(stressCol),sMax(stressCol));
                    //                   glShadeModel(GL_SMOOTH);
                    glBegin(GL_TRIANGLES);
                    glColor4f(clr0.r, clr0.g, clr0.b,0.5);
                    glVertex3f(deq[3*k+0].first(0),deq[3*k+0].first(1),deq[3*k+0].first(2));
                    glColor4f(clr1.r, clr1.g, clr1.b,0.5);
                    glVertex3f(deq[3*k+1].first(0),deq[3*k+1].first(1),deq[3*k+1].first(2));
                    glColor4f(clr2.r, clr2.g, clr2.b,0.5);
                    glVertex3f(deq[3*k+2].first(0),deq[3*k+2].first(1),deq[3*k+2].first(2));
                    glEnd();
                }
                
            }
            
            glDisable(GL_DEPTH_TEST);
            glEnable (GL_BLEND);
            glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            
            if (showMesh>0)
            {
                float dispCorr(dispScale*(showMesh>1));
                                
                glDisable(GL_COLOR_MATERIAL); // use glMaterialfv(...) to set material colors
                GLfloat materialColor[]={0.0, 0.0, 0.0, 0.1};
                glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialColor);      // ambient color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialColor);      // diffuse color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialColor);  // specular color for the material
                glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialColor);  // emission color for the material
                
                for (auto& edge : edgeVector)
                {
                    Eigen::Matrix<float,3,2> disp=Eigen::Matrix<float,3,2>::Zero();
                    
                    if(dispFileIsGood)
                    {
                        VertexReader<'D',4,float>::const_iterator iterD1(DispContainerType::find(edge.first.first));
                        VertexReader<'D',4,float>::const_iterator iterD2(DispContainerType::find(edge.first.second));
                        assert(iterD1!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                        assert(iterD2!=DispContainerType::end() && "MESH NODE NOT FOUND IN D FILE");
                        disp.col(0)=iterD1->second.segment<3>(0);
                        disp.col(1)=iterD2->second.segment<3>(0);
                    }
                    
                    glBegin(GL_LINES);
                    glVertex3f(edge.second(0,0)+disp(0,0)*dispCorr, edge.second(1,0)+disp(1,0)*dispCorr,edge.second(2,0)+disp(2,0)*dispCorr);
                    glVertex3f(edge.second(0,1)+disp(0,1)*dispCorr, edge.second(1,1)+disp(1,1)*dispCorr,edge.second(2,1)+disp(2,1)*dispCorr);
                    glEnd();
                }
                
//                // test boundaryNormal
//                for (typename SimplexObserver<3,0>::const_iterator sIter=SimplexObserver<3,0>::simplexBegin();
//                     /*                                         */ sIter!=SimplexObserver<3,0>::simplexEnd();++sIter)
//                {
//                    Eigen::Matrix<double,3,1> outNormal=sIter.second->outNormal()*100.0;
//                    glBegin(GL_LINES);
//                    glVertex3f(sIter.second->P0(0),sIter.second->P0(1),sIter.second->P0(2));
//                    glVertex3f(sIter.second->P0(0)+outNormal(0),sIter.second->P0(1)+outNormal(1),sIter.second->P0(2)+outNormal(2));
//                    glEnd();
//                }
                
                if(showSpecificSimplex)
                {
                    
                    glEnable(GL_COLOR_MATERIAL); // use glColor4f to set color
                    
                    //                    std::cout<<specificSimplexID<<std::endl;
                    
                    const auto iter=p_mesh->template observer<3>().find(specificSimplexID);
                    if(iter!=p_mesh->template observer<3>().end())
                    {
                        Eigen::Matrix<float,3,4> mat=iter->second->vertexPositionMatrix().template cast<float>();
                        
                        
                        glColor4f(0.5, 0.0, 0.5,0.1);
                        
                        glBegin(GL_TRIANGLES);
                        glVertex3f(mat.col(0)(0),mat.col(0)(1),mat.col(0)(2));
                        glVertex3f(mat.col(1)(0),mat.col(1)(1),mat.col(1)(2));
                        glVertex3f(mat.col(2)(0),mat.col(2)(1),mat.col(2)(2));
                        glEnd();
                        
                        glBegin(GL_TRIANGLES);
                        glVertex3f(mat.col(0)(0),mat.col(0)(1),mat.col(0)(2));
                        glVertex3f(mat.col(1)(0),mat.col(1)(1),mat.col(1)(2));
                        glVertex3f(mat.col(3)(0),mat.col(3)(1),mat.col(3)(2));
                        glEnd();
                        
                        glBegin(GL_TRIANGLES);
                        glVertex3f(mat.col(1)(0),mat.col(1)(1),mat.col(1)(2));
                        glVertex3f(mat.col(2)(0),mat.col(2)(1),mat.col(2)(2));
                        glVertex3f(mat.col(3)(0),mat.col(3)(1),mat.col(3)(2));
                        glEnd();
                        
                        glBegin(GL_TRIANGLES);
                        glVertex3f(mat.col(0)(0),mat.col(0)(1),mat.col(0)(2));
                        glVertex3f(mat.col(2)(0),mat.col(2)(1),mat.col(2)(2));
                        glVertex3f(mat.col(3)(0),mat.col(3)(1),mat.col(3)(2));
                        glEnd();
                        
                        glColor4f(0.0, 0.0, 0.0,0.8);
                        
                        
                        glBegin(GL_LINES);
                        glVertex3f(mat.col(0)(0),mat.col(0)(1),mat.col(0)(2));
                        glVertex3f(mat.col(1)(0),mat.col(1)(1),mat.col(1)(2));
                        glEnd();
                        
                        glBegin(GL_LINES);
                        glVertex3f(mat.col(0)(0),mat.col(0)(1),mat.col(0)(2));
                        glVertex3f(mat.col(2)(0),mat.col(2)(1),mat.col(2)(2));
                        glEnd();
                        
                        glBegin(GL_LINES);
                        glVertex3f(mat.col(0)(0),mat.col(0)(1),mat.col(0)(2));
                        glVertex3f(mat.col(3)(0),mat.col(3)(1),mat.col(3)(2));
                        glEnd();
                        
                        glBegin(GL_LINES);
                        glVertex3f(mat.col(1)(0),mat.col(1)(1),mat.col(1)(2));
                        glVertex3f(mat.col(2)(0),mat.col(2)(1),mat.col(2)(2));
                        glEnd();
                        
                        glBegin(GL_LINES);
                        glVertex3f(mat.col(1)(0),mat.col(1)(1),mat.col(1)(2));
                        glVertex3f(mat.col(3)(0),mat.col(3)(1),mat.col(3)(2));
                        glEnd();
                        
                        glBegin(GL_LINES);
                        glVertex3f(mat.col(2)(0),mat.col(2)(1),mat.col(2)(2));
                        glVertex3f(mat.col(3)(0),mat.col(3)(1),mat.col(3)(2));
                        glEnd();
                        
                        for(int c=0;c<4;++c)
                        {
                            Eigen::Matrix<float,3,1> P=mat.col(c);
                            BitmapPlotter::renderString(P, static_cast<std::ostringstream*>( &(std::ostringstream() << specificSimplexID(c)) )->str());
                        }
                        
                    }
                    else
                    {
                        std::cout<<"Simplex not found"<<std::endl;
                    }
                }
                
            }
            
            if(showRegionBoundaries)
            {
                glEnable(GL_COLOR_MATERIAL); // use glColor4f to set color
                
                
                for(int k=0; k<regionsBndDeq.size();++k)
                {
                    
                    auto& mat=regionsBndDeq[k];
                    
                    RGBcolor clr0=RGBmap::getColor(float(regionsBndIDDeq[k].first+regionsBndIDDeq[k].second)*0.5,0,rMax);
                    
                    glColor4f(clr0.r, clr0.g, clr0.b,0.1);
                    
                    
                    //                    GLfloat materialColor1[]={0.0, 0.8, 0.0, 0.1};
                    //                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialColor1);      // ambient color for the material
                    //                    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialColor1);      // diffuse color for the material
                    //                    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialColor1);  // specular color for the material
                    //                    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialColor1);  // emission color for the material
                    glBegin(GL_TRIANGLES);
                    glVertex3f(mat.col(0)(0),mat.col(0)(1),mat.col(0)(2));
                    glVertex3f(mat.col(1)(0),mat.col(1)(1),mat.col(1)(2));
                    glVertex3f(mat.col(2)(0),mat.col(2)(1),mat.col(2)(2));
                    glEnd();
                    
                    //                    GLfloat materialColor2[]={0.0, 0.2, 0.0, 0.1};
                    //                    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, materialColor2);      // ambient color for the material
                    //                    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, materialColor2);      // diffuse color for the material
                    //                    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, materialColor2);  // specular color for the material
                    //                    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, materialColor2);  // emission color for the material
                    glBegin(GL_LINE_LOOP);
                    glVertex3f(mat.col(0)(0),mat.col(0)(1),mat.col(0)(2));
                    glVertex3f(mat.col(1)(0),mat.col(1)(1),mat.col(1)(2));
                    glVertex3f(mat.col(2)(0),mat.col(2)(1),mat.col(2)(2));
                    glEnd();
                }
            }
            
            if(specialTrianglesFileIsGood)
            {
                glEnable(GL_COLOR_MATERIAL); // use glColor4f to set color

                glColor4f(0.0, 1.0, 0.5,0.5);

                
                for(const auto& mat : specialTrianglesDeq)
                {
                    glBegin(GL_TRIANGLES);
                    glVertex3f(mat.col(0)(0),mat.col(0)(1),mat.col(0)(2));
                    glVertex3f(mat.col(1)(0),mat.col(1)(1),mat.col(1)(2));
                    glVertex3f(mat.col(2)(0),mat.col(2)(1),mat.col(2)(2));
                    glEnd();
                
                }
            }
            
            glDisable(GL_BLEND);
            
            glEnable(GL_DEPTH_TEST); //Makes 3D drawing work when something is in front of something else
            
            
            
        }
        
        
        
    };
    
    // Declare static data
    bool MeshPlotter::plotBndStress=false;
    unsigned int MeshPlotter::stressCol=0;
    bool MeshPlotter::showRegionBoundaries=false;
    bool MeshPlotter::showSpecificSimplex=false;
    Simplex<3,3>::SimplexIDType MeshPlotter::specificSimplexID;
    
} // namespace model
#endif

