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
    class GrainBoundary :
    /* base */ private std::map<int,const Grain<dim>* const>,
    /* base */ private std::map<int,LatticePlane>
    {
        typedef MeshRegionBoundary<Simplex<dim,dim-1> > MeshRegionBoundaryType;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeVector<dim> ReciprocalLatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        typedef Grain<dim> GrainType;
        typedef std::map<int,const GrainType* const> GrainContainerType;
        typedef std::map<int,LatticePlane> LatticePlaneContainerType;
        
        /**********************************************************************/
        void storeLatticePlane(const Grain<dim>& grain,
                               const VectorDimD& normal)
        {
//            std::cout<<std::setprecision(15)<<std::scientific<<normal<<std::endl;
            ReciprocalLatticeDirectionType R=grain.reciprocalLatticeDirection(normal);
            
            model::cout<<"Grain boundary normal for grain"<< grain.grainID<<std::endl;
            model::cout<<"cartesian components="<<R.cartesian().transpose()<<std::endl;
            model::cout<<"lattice components="<<R.transpose()<<std::endl;

            //            model::cout<<"Cartesian Grain boundary for grain"<< grain.grainID<<" is "<<R.cartesian().transpose()<<std::endl;
            
            
            // CHENAGE THIS PART READING FROM TABLE OF GB CASES
            VectorDimD v1;
            VectorDimD v2;

            v1<< 3.0,0.0,-1.0;
            v1*=sqrt(2.0);
            v2<< 0.0,1.0, 0.0;
            v2*=sqrt(2.0);
            
            if(grain.grainID==1)
            {
            }
            else
            {
                v2*=-1.0; // when implementing GB table, normals must be opposite sign
            }
//            v1=grain.get_C2G()*v1*sqrt(2.0);
//            v2=grain.get_C2G()*v2*sqrt(2.0);
            
            LatticeVectorType O(grain.covBasis(),grain.contraBasis());
            
            // Here is ok
            LatticePlaneBase pb(LatticeVectorType(v1,grain.covBasis(),grain.contraBasis()),
                                LatticeVectorType(v2,grain.covBasis(),grain.contraBasis()));
            
            
            const auto temp = LatticePlaneContainerType::emplace(std::piecewise_construct,
                                               std::forward_as_tuple(grain.grainID),
                                               std::forward_as_tuple(O,pb));
            
            std::cout<<"GB plane normal="<<temp.first->second.n.cartesian().transpose()<<std::endl;

            
        }
        
        
    public:
        
        
        const MeshRegionBoundaryType& regionBoundary;
        
        /**********************************************************************/
        GrainBoundary(const MeshRegionBoundaryType& regionbnd_in,
                      const Grain<dim>& grainFirst,
                      const Grain<dim>& grainSecond) :
        /* init */ regionBoundary(regionbnd_in)
        {
            model::cout<<"Creating GrainBoundary ("<<regionBoundary.regionBndID.first<<" "<<regionBoundary.regionBndID.second<<")"<<std::endl;
            GrainContainerType::emplace(grainFirst.grainID,&grainFirst);
            GrainContainerType::emplace(grainSecond.grainID,&grainSecond);
            
        }
        
        /**********************************************************************/
        const GrainContainerType& grains() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const Grain<dim>& grain(const int& k) const
        {
            return *(grains().at(k));
        }
        
        /**********************************************************************/
        const std::map<int,LatticePlane>& latticePlanes() const
        {
            return *this;
        }
        
        /**********************************************************************/
        const LatticePlane& latticePlane(const int& k) const
        {
            return latticePlanes().at(k);
        }
        
        /**********************************************************************/
        void createLatticePlanes()
        {
            LatticePlaneContainerType::clear();
            const Simplex<dim,dim-1>& triangle(**regionBoundary.begin());
            for(const auto& tet : triangle.parents())
            {
                const size_t faceID=tet->childOrder(triangle.xID);
                const VectorDimD outNormal=tet->nda.col(faceID);
                storeLatticePlane(grain(tet->region->regionID),outNormal);
            }
            
            // Check that all tringle vertices are contained by both GB planes
            for(const auto& triangle : regionBoundary)
            {
                const auto Ps=triangle->vertexPositionMatrix();
                for(int j=0;j<Ps.cols();++j)
                {
                    for(const auto& grain : grains())
                    {
                        assert(latticePlane(grain.second->grainID).contains(Ps.col(j)) && "TRIANGLE VERTEX NOT CONTAINED IN GBPLANE");
                    }
                }                
            }
        }
        
        
        
    };
    
    
} // end namespace
#endif

