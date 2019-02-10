/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SimplicialMeshFace_H_
#define model_SimplicialMeshFace_H_

#include <chrono>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <utility>      // std::pair, std::make_pair
#include <set>

#include <SimplexTraits.h>
#include <Simplex.h>
#include <ConvexHull.h>



namespace model
{
    
    
    /**************************************************************************/
    /**************************************************************************/
    template<int dim>
    class SimplicialMeshFace : public std::set<const Simplex<dim,dim-1>*>
    {
     
        Eigen::Matrix<double,dim,1> n;
        Eigen::Matrix<double,dim,1> c;
        std::vector<const Simplex<dim,0>*> _hull;
        
    public:
        
        /**********************************************************************/
        SimplicialMeshFace(const Simplex<dim,dim-1>* const pS) :
        /* init */ n(pS->outNormal())
        /* init */,c(pS->center())
        {
            assert((n.norm()-1.0)<FLT_EPSILON);
            this->insert(pS);
        }
        
        /**********************************************************************/
        const std::vector<const Simplex<dim,0>*>& convexHull() const
        {
            return _hull;
        }
        
        /**********************************************************************/
        void finalize()
        {
            // Recompute c and n to minimize floating point errors
            n.setZero();
            c.setZero();
            for(const auto& simplex : *this)
            {
                n+=simplex->outNormal();
                c+=simplex->center();
            }
            n/=this->size();
            c/=this->size();
            assert(fabs(n.norm()-1.0)<FLT_EPSILON);
            
            std::set<const Simplex<dim,0>*> vertices;
            for(const auto& simplex : *this)
            {// Collect vertices on this face
                for(const auto& v : simplex->vertices())
                {
                    vertices.insert(v);
                }
            }
            
            if(vertices.size())
            {
                
                // Compute local rotation matrix
                Eigen::Matrix<double,dim,dim> R;
                R.col(0)=((*vertices.begin())->P0-center()).normalized();
                R.col(2)=outNormal();
                R.col(1)=R.col(2).cross(R.col(0));

                ConvexHull<2,Simplex<dim,0>> hull;
                for(const auto& v : vertices)
                {// compute local position in the plane as insert in ConvexHull
                    const Eigen::Matrix<double,dim,1> x(R.transpose()*(v->P0-center()));
                    HullPoint<2,Simplex<dim,0>> hp({x[0],x[1]},v);
                    hull.push_back(hp);
                }
                
                
                const auto hullPts=hull.getHull();
                _hull.clear();
                for(const auto& hp : hullPts)
                {
                    _hull.push_back(hp.t);
                }
                
            }
            
            
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& outNormal() const
        {
            return n;
        }
        
        /**********************************************************************/
        const Eigen::Matrix<double,dim,1>& center() const
        {
            return c;
        }
        
    };
    
}
#endif
