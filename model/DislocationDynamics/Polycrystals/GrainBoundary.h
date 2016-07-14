/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GrainBoundary_H_
#define model_GrainBoundary_H_

#include <cfloat> // FLT_EPSILON
#include <Eigen/Core>
#include <model/Math/RoundEigen.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/Mesh/MeshRegionObserver.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/LatticeMath/LatticePlane.h>

namespace model
{
    
    
    
    template <int dim>
    class GrainBoundary : private std::map<int,LatticePlaneBase>
    {
        
        typedef MeshRegionBoundary<Simplex<dim,dim-1> > MeshRegionBoundaryType;

        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;

        //        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
//        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        
        
        void storeLatticePlane(const Grain<dim>& grain,
                               const VectorDimD& normal)
        {
            
            ReciprocalLatticeDirectionType R=grain.reciprocalLatticeDirection(normal);
            
            model::cout<<"Grain boundary for grain"<< grain.grainID<<" is "<<R.transpose()<<std::endl;
            model::cout<<"Cartesian Grain boundary for grain"<< grain.grainID<<" is "<<R.cartesian().transpose()<<std::endl;

            VectorDimD v1;
            VectorDimD v2;
            if(grain.grainID==1)
            {
                v1<<-1.0,0.0,-3.0;
                v2<< 0.0,1.0, 0.0;
            }
            else
            {
                v1<<-1.0,0.0, 3.0;
                v2<< 0.0,1.0, 0.0;
            }
            v1=grain.get_C2G()*v1*sqrt(2.0);
            v2=grain.get_C2G()*v2*sqrt(2.0);
            
            this->emplace(std::piecewise_construct,
             std::forward_as_tuple(grain.grainID),
             std::forward_as_tuple(LatticeVectorType(v1,grain.covBasis(),grain.contraBasis()),
                                   LatticeVectorType(v2,grain.covBasis(),grain.contraBasis())
                                   )
                          );
            
        }

        
    public:
        
        
        const MeshRegionBoundaryType& regionBoundary;
        
        const Grain<dim>& grainFirst;
        const Grain<dim>& grainSecond;
        
        
        
        /**********************************************************************/
        GrainBoundary(const MeshRegionBoundaryType& regionbnd_in,
                      const Grain<dim>& grainFirst_in,
                      const Grain<dim>& grainSecond_in) :
        /* init */ regionBoundary(regionbnd_in),
        /* init */ grainFirst(grainFirst_in),
        /* init */ grainSecond(grainSecond_in)
        {
            model::cout<<"Creating GrainBoundary ("<<regionBoundary.regionBndID.first<<" "<<regionBoundary.regionBndID.second<<")"<<std::endl;
            
        }
        
        /**********************************************************************/
        void createLatticePlanes()
        {
            this->clear();
            
            const MatrixDimD vertexMatrix=(*regionBoundary.begin())->vertexPositionMatrix();
            const VectorDimD normal((vertexMatrix.col(1)-vertexMatrix.col(0)).cross(vertexMatrix.col(2)-vertexMatrix.col(0)));
            
            std::cout<<"cartesian BG normal="<<normal.transpose()<<std::endl;
            
            storeLatticePlane(grainFirst,normal);
            storeLatticePlane(grainSecond,normal);
            
        }
        
        const LatticePlaneBase& latticePlane(const int& k) const
        {
            return this->at(k);
        }
        
    };
    
    
} // end namespace
#endif

