/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplexBndNormal_H_
#define model_SimplexBndNormal_H_

#include <Eigen/Dense>
#include <model/Mesh/Simplex.h>

namespace model
{
    
    class SimplexBndNormal
    {
        
    public:
        
        // make dim available outside class
        constexpr static int dim=3;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        /**********************************************************************/
        static VectorDim get_boundaryNormal(const VectorDim& P,
                                     const Simplex<dim,dim>& simplex,
                                     const double& dmax)
        {
            
            
            VectorDim temp(VectorDim::Zero());
            
                const Eigen::Matrix<double,dim+1,1> bary(simplex.pos2bary(P));
            
            
                // 1) check if P is close to a boundary face
                for(int f=0;f<Simplex<dim,dim>::nFaces;++f) // loop over faces of current Simplex in the path
                {
                    
                    const double h=3.0*simplex.vol0/simplex.child(f).vol0*bary(f); // distance from f-th face
                    
                    if(simplex.child(f).isBoundarySimplex() && std::fabs(h)<dmax)
                    {
                        temp=simplex.nda.col(f).normalized();
                        break;
                    }
                }
                
                
                // 2) check if P is close to a boundary edge
                if(temp.squaredNorm()==0.0) // no boundary face was found
                {
                    for(int f=0;f<Simplex<dim,dim>::nFaces;++f) // loop over faces of current Simplex in the path
                    {
                        
                        for(int e=0;e<3;++e)
                        {
                            
                            const VectorDim& x1=simplex.child(f).child(e).child(0).P0;
                            const VectorDim& x2=simplex.child(f).child(e).child(1).P0;
                            const double d=(P-x1).cross(P-x2).norm()/(x2-x1).norm(); // distance to edge
                                                        
                            if(d<dmax && simplex.child(f).child(e).isBoundarySimplex()) // point close to edge, and edge is on boundary
                            {
                                for(const auto& parent : simplex.child(f).child(e).parents())
                                {
                                    if(parent->isBoundarySimplex())
                                    {
                                        temp+=simplex.nda.col(f).normalized();
                                    }
                                }
                            }
                        }
                    }
                    const double tempNorm(temp.norm());
                    if(tempNorm)
                    {
                        temp.normalize();
                    }
                    
                }
                
                // check if P is close to a boundary vertex
                if(temp.squaredNorm()==0.0) // no boundary face was found
                {
                    for(int f=0;f<Simplex<dim,dim>::nFaces;++f) // loop over faces of current Simplex in the path
                    {
                        for(int e=0;e<3;++e)
                        {
                            for(int v=0;v<2;++v)
                            {
                                if((simplex.child(f).child(e).child(v).P0-P).norm()<dmax && simplex.child(f).child(e).child(v).isBoundarySimplex())
                                {
                                    temp+=simplex.nda.col(f).normalized();

                                }
                            }
                        }
                    }
                    
                    const double tempNorm(temp.norm());
                    if(tempNorm)
                    {
                        temp.normalize();
                    }
                }
                
//            }
            return temp;
        }
        
        
    };
    
    
    
    
} // close namespace
#endif

