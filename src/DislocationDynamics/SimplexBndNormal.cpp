/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexBndNormal_cpp_
#define model_SimplexBndNormal_cpp_

#include <SimplexBndNormal.h>

namespace model
{
    
    
    
        
        /**********************************************************************/
     typename SimplexBndNormal::VectorDim SimplexBndNormal::get_boundaryNormal(const VectorDim& P,
                                     const Simplex<dim,dim>& simplex,
                                     const double& dmax)
        {
            
            //std::cout<<simplex.xID<<std::endl;
            VectorDim temp(VectorDim::Zero());
            bool found=false;
            
            // 1) check if P is close to a boundary vertex
            for(int f=0;f<Simplex<dim,dim>::nFaces;++f) // loop over faces of current Simplex in the path
            {
                for(int e=0;e<3;++e)
                {
                    for(int v=0;v<2;++v)
                    {
                        if((simplex.child(f).child(e).child(v).P0-P).norm()<dmax && simplex.child(f).child(e).child(v).isBoundarySimplex())
                        {
                            temp=simplex.child(f).child(e).child(v).outNormal();
                            found=true;
                        }
                        if(found)
                        {
                            break;
                        }
                    }
                    if(found)
                    {
                        break;
                    }
                }
                if(found)
                {
                    break;
                }
            }
            
            if(!found)
            {// No boundary vertex was found. Check if P is close to a boundary edge
                for(int f=0;f<Simplex<dim,dim>::nFaces;++f) // loop over faces of current Simplex in the path
                {
                    for(int e=0;e<3;++e)
                    {
                        const VectorDim& x1=simplex.child(f).child(e).child(0).P0;
                        const VectorDim& x2=simplex.child(f).child(e).child(1).P0;
                        const double d=(P-x1).cross(P-x2).norm()/(x2-x1).norm(); // distance to edge
                        if(d<dmax && simplex.child(f).child(e).isBoundarySimplex()) // point close to edge, and edge is on boundary
                        {
                            temp=simplex.child(f).child(e).outNormal();
                            found=true;
                        }
                        if(found)
                        {
                            break;
                        }
                    }
                    if(found)
                    {
                        break;
                    }
                }
            }
            
            // 3) check if P is close to a boundary face
            if(!found) // no boundary face was found
            {
                const Eigen::Matrix<double,dim+1,1> bary(simplex.pos2bary(P));
                
                for(int f=0;f<Simplex<dim,dim>::nFaces;++f) // loop over faces of current Simplex in the path
                {// volume of tetrahedron V=1/3*A*h. Volume = V0*bary. So h=3*V0*bary/A
                    const double h=3.0*simplex.vol0*bary(f)/simplex.child(f).vol0; // distance from f-th face
//                    std::cout<<"h="<<h<<" "<<simplex.child(f).isBoundarySimplex()<<std::endl;
                    
                    if(simplex.child(f).isBoundarySimplex() && std::fabs(h)<dmax)
                    {
//                        temp=simplex.nda.col(f).normalized();
                        temp=simplex.child(f).outNormal();
                        found=true;
                    }
                    if(found)
                    {
                        break;
                    }
                }
            }
            
            return temp;
        }
        
} // close namespace
#endif

