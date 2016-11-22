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
#include <vector>
#include <deque>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <model/Math/RoundEigen.h>
#include <model/Mesh/SimplicialMesh.h>
#include <model/Mesh/MeshRegionObserver.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundaryType.h>
#include <model/LatticeMath/LatticePlane.h>
#include <model/DislocationDynamics/StressStraight.h>


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
        typedef Eigen::Matrix<double,2*dim+1,1> BGkeyType;
        
        /**********************************************************************/
        void storeLatticePlane(const Grain<dim>& grain,
                               const VectorDimD& normal)
        {
            const ReciprocalLatticeDirectionType R=grain.reciprocalLatticeDirection(normal);
            
            model::cout<<"GB normal for grain "<< grain.grainID<<":"<<std::endl;
            model::cout<<"  cartesian components="<<R.cartesian().transpose()<<std::endl;
            model::cout<<"  crystallographic components="<<(grain.get_C2G().transpose()*R.cartesian()).transpose()<<std::endl;
            
            LatticeVectorType L0(grain.covBasis(),grain.contraBasis());
            bool latticePointFound=false;
            
            
            for(const auto& simplexPtr : regionBoundary.simplices())
            {
                for(size_t d=0;d<dim;++d)
                {
                    const VectorDimD P0=simplexPtr->vertexPositionMatrix().col(d);
                    L0=grain.snapToLattice(P0);
                    
                    if(fabs((L0.cartesian()-P0).dot(normal))<FLT_EPSILON) // L0 belongs to mesh plane
                    {
                        latticePointFound=true;
                        break; //break inner loop
                    }
                }
                
                if(latticePointFound==true) // L0 belongs to mesh plane
                {
                    break; //break outer loop
                }
            }
            
            assert(latticePointFound && "None of the GB mesh vertices, snapped to the lattice, belongs to GB plane.");
            
            LatticePlaneBase pb(R);
            
            const auto temp = LatticePlaneContainerType::emplace(std::piecewise_construct,
                                                                 std::forward_as_tuple(grain.grainID),
                                                                 std::forward_as_tuple(L0,pb));
            
            assert(temp.first->second.n.cartesian().normalized().cross(normal.normalized()).norm()<FLT_EPSILON && "LatticePlane normal and triangle normal are not the same.");
            
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
                const double outNorm(outNormal.norm());
                assert(outNorm>FLT_EPSILON && "Simplex normal has zero norm.");
                storeLatticePlane(grain(tet->region->regionID),outNormal/outNorm);
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
        
        /**********************************************************************/
        void computeCrystallographicRotationAxis()
        {
            
            // Compute the relative rotation matrix between grains
            const MatrixDimD R(grain(grainBndID.first).get_C2G().transpose()*grain(grainBndID.second).get_C2G());
            
            // Eigen-decompose R
            Eigen::EigenSolver<MatrixDimD> es(R);
            
            // Determine _crystallographicRotationAxis as the eigenvector of R corresponding to eigenvalue=1
            for(unsigned int j=0;j<dim;++j)
            {
                if(fabs(es.eigenvalues()(j).real()-1.0)<FLT_EPSILON && fabs(es.eigenvalues()(j).imag())<FLT_EPSILON)
                {
                    for(unsigned int d=0;d<dim;++d)
                    {
                        _crystallographicRotationAxis=es.eigenvectors().col(j).real();
                    }
                    break;
                }
            }
            const double axisNorm(_crystallographicRotationAxis.norm());
            assert(axisNorm>FLT_EPSILON);
            _crystallographicRotationAxis/=axisNorm;
            
            _rotationAxis=grain(grainBndID.first).get_C2G()*_crystallographicRotationAxis;
            assert((_rotationAxis-grain(grainBndID.second).get_C2G()*_crystallographicRotationAxis).norm()<FLT_EPSILON && "rotationAxis inconsistency.");
            
            cosTheta=0.5*(R.trace()-1.0);
            
            std::cout<<yellowColor<<"   Rotation axis="<<_crystallographicRotationAxis.transpose()<<std::endl;
            std::cout<<yellowColor<<"   Rotation angle="<<acos(cosTheta)*180.0/M_PI<<" deg"<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        void findGrainBoundaryType(const std::deque<GrainBoundaryType<dim>>& bgTypes)
        {
            
            const VectorDimD n1=(grain(grainBndID.first).get_C2G().transpose()*latticePlane(grainBndID.first).n.cartesian()).normalized();
            const VectorDimD n2=(grain(grainBndID.second).get_C2G().transpose()*latticePlane(grainBndID.second).n.cartesian()).normalized();

            std::cout<<"n1="<<n1.transpose()<<std::endl;
            std::cout<<"n2="<<n2.transpose()<<std::endl;
            
            for (const auto& gbt : bgTypes)
            {
                bool axisFound=false;
                for(const auto& a : gbt.axisPermutations)
                {
                    axisFound=( (_crystallographicRotationAxis-a.normalized()).norm()<FLT_EPSILON ||
                               (_crystallographicRotationAxis+a.normalized()).norm()<FLT_EPSILON );
                    if(axisFound)
                    {
                        break;
                    }
                }
                
                bool n1Found=false;
                for(const auto& n : gbt.n1Permutations)
                {
                    n1Found=( (n1-n.normalized()).norm()<FLT_EPSILON || (n1+n.normalized()).norm()<FLT_EPSILON );
                    if(n1Found)
                    {
                        break;
                    }
                }
                
                bool n2Found=false;
                for(const auto& n : gbt.n2Permutations)
                {
                    n2Found=( (n2-n.normalized()).norm()<FLT_EPSILON || (n2+n.normalized()).norm()<FLT_EPSILON );
                    if(n2Found)
                    {
                        break;
                    }
                }
                
                if(axisFound && n1Found && n2Found)
                {
                    p_gbType=&gbt;
                    break;
                }
                else
                {
                    std::cout<<"axisFound"<<axisFound<<std::endl;
                    std::cout<<"n1Found"<<n1Found<<std::endl;
                    std::cout<<"n2Found"<<n2Found<<std::endl;
                }

                
            }
            
            if(p_gbType!=NULL)
            {
                std::cout<<yellowColor<<"   GB type is "<<p_gbType->name<<std::endl;
            }
            else
            {
                assert(0 && "GRAIN BONUDARY TYPE NOT FOUND.");
            }

            
        }
        
        /**********************************************************************/
        void populateGBdislocations(std::vector<StressStraight<dim>>& vD) const
        {
            
            const VectorDimD dir=rotationAxis().normalized();
            const VectorDimD p=dir.cross(latticePlane(grainBndID.first).n.cartesian()).normalized()*grainBoundaryType().dislocationSpacing;
            
            
            VectorDimD P0(VectorDimD::Zero());
            VectorDimD P1(VectorDimD::Zero());
            VectorDimD b(VectorDimD::Zero());
            
            for(int k=0;k<10;++k)
            {
                P0=latticePlane(grainBndID.first).P.cartesian()-100.0*dir+k*p;
                P1=latticePlane(grainBndID.first).P.cartesian()+100.0*dir+k*p;
                b=latticePlane(grainBndID.first).n.cartesian();
                vD.emplace_back(P0,P1,b);
            }
            
        }
        
        VectorDimD _crystallographicRotationAxis;
        VectorDimD _rotationAxis;
        
        double cosTheta; // cosine of relative rotation angle between grains
        const GrainBoundaryType<dim>* p_gbType;
        
    public:
        
        const MeshRegionBoundaryType& regionBoundary;
        const std::pair<size_t,size_t>& grainBndID;
        
        /**********************************************************************/
        GrainBoundary(const MeshRegionBoundaryType& regionbnd_in,
                      const Grain<dim>& grainFirst,
                      const Grain<dim>& grainSecond) :
        /* init */ _crystallographicRotationAxis(VectorDimD::Zero()),
        /* init */ _rotationAxis(_crystallographicRotationAxis),
        /* init */ cosTheta(1.0),
        /* init */ p_gbType(NULL),
        /* init */ regionBoundary(regionbnd_in),
        /* init */ grainBndID(regionBoundary.regionBndID)
        {
            model::cout<<"Creating GrainBoundary ("<<grainBndID.first<<" "<<grainBndID.second<<")"<<std::endl;
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
        template<typename PolycrystalType>
        void initializeGrainBoundary(PolycrystalType& poly)
        {
            model::cout<<yellowColor<<"GrainBoundary ("<<grainBndID.first<<" "<<grainBndID.second<<")"<<defaultColor<<std::endl;
            computeCrystallographicRotationAxis();
            createLatticePlanes();
            findGrainBoundaryType(poly.grainBoundaryTypes());
            populateGBdislocations(poly.grainBoundaryDislocations());
        }
        
        /**********************************************************************/
        const VectorDimD& crystallographicRotationAxis() const
        {
            return _crystallographicRotationAxis;
        }
        
        /**********************************************************************/
        const VectorDimD& rotationAxis() const
        {
            return _rotationAxis;
        }
        
        /**********************************************************************/
        double rotationAngle() const
        {
            return acos(cosTheta);
        }
        
        /**********************************************************************/
        const GrainBoundaryType<dim>& grainBoundaryType() const
        {
            return *p_gbType;
        }
    };
    
} // end namespace
#endif

