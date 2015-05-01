/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_LineMeshIntersection_h_
#define model_LineMeshIntersection_h_

#include <model/LatticeMath/LatticeLine.h>
#include <model/Mesh/SimplicialMesh.h>

namespace model
{
    class LineMeshIntersection
    {
        static constexpr int dim=3;
        
        typedef LatticeVector<dim> LatticeVectorType;
        
        /**********************************************************************/
        LatticeVectorType get_intersection(const LatticeLine& line,
                                           LatticeVectorType L1,
                                           const SimplicialMesh<dim>& mesh,
                                           const Simplex<dim,dim>* guess) const
        {
            const int d2=line.d.squaredNorm(); // minimum (squared) distance along line
            assert(d2>0 && "LINE HAS ZERO DIRECTION");
            
            LatticeVectorType L0=line.P;
            std::pair<bool,const Simplex<dim,dim>*> temp=mesh.searchWithGuess(L0.cartesian(),guess);
            assert(temp.first && "LINE-MESH-INTERSECTION, STARTING POINT NOT INSIDE MESH");
            
            //  Bring L1 to line
            L1=LatticeVectorType(line.snapToLattice(L1.cartesian()));
            
            assert((L1-L0).squaredNorm()>0 && "L0 and L1 are the same");
            
            int n=(L1-L0).dot(line.d)/d2;
            
            temp=mesh.searchWithGuess(L1.cartesian(),guess);
            while (temp.first) // repeat untile L1 is found outside
            {
                L0=L1;
                n*=2;
                L1=line.P+n*line.d;
                temp=mesh.searchWithGuess(L1.cartesian(),temp.second);
            }
            
            // now L1 is ouside mesh, and L0 is inside
            
            while((L0-L1).squaredNorm()>d2)
            {
                LatticeVectorType L3(line.snapToLattice(0.5*(L0.cartesian()+L1.cartesian())));
                temp=mesh.searchWithGuess(L3.cartesian(),temp.second);
                if (temp.first) // midpoint inside
                {
                    L0=L3;
                }
                else
                {
                    L1=L3;
                }
            }
            
            return L0;
        }
        
    public:
        const LatticeVectorType L;
        const std::pair<bool,const Simplex<dim,dim>*> search;
        
        /**********************************************************************/
        LineMeshIntersection(const LatticeLine& line,
                             LatticeVectorType L1,
                             const SimplicialMesh<dim>& mesh,
                             const Simplex<dim,dim>* guess) :
        /* init */ L(get_intersection(line,L1,mesh,guess)),
        /* init */ search(mesh.searchWithGuess(L.cartesian(),guess))
        {
        
        }
        
        
    };
    
} // end namespace
#endif
