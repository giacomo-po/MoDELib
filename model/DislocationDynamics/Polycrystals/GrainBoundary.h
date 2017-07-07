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
#include <model/LatticeMath/LineMeshIntersection.h>
#include <model/LatticeMath/CSL.h>
#include <model/LatticeMath/DSCL.h>
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
            
            LatticeVectorType L0(grain.lattice());
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
                
                if(!n1Found && ! n2Found) // search n1 in n2Permutations and viceversa
                {
                    for(const auto& n : gbt.n2Permutations)
                    {
                        n1Found=( (n1-n.normalized()).norm()<FLT_EPSILON || (n1+n.normalized()).norm()<FLT_EPSILON );
                        if(n1Found)
                        {
                            break;
                        }
                    }
                    
                    for(const auto& n : gbt.n1Permutations)
                    {
                        n2Found=( (n2-n.normalized()).norm()<FLT_EPSILON || (n2+n.normalized()).norm()<FLT_EPSILON );
                        if(n2Found)
                        {
                            break;
                        }
                    }
                }
                
                if(axisFound && n1Found && n2Found)
                {
                    p_gbType=&gbt;
                    break;
                }
                
            }
            
            if(p_gbType!=NULL)
            {
                std::cout<<yellowColor<<"   GB type is "<<p_gbType->name<<defaultColor<<std::endl;
                std::cout<<yellowColor<<"   GB energy= "<<p_gbType->energyDensity<<defaultColor<<std::endl;
                std::cout<<yellowColor<<"   GB dislocation spacing= "<<p_gbType->dislocationSpacing<<defaultColor<<std::endl;
                std::cout<<yellowColor<<"   Frank-Bilby dislocation spacing= "<<p_gbType->FrankBilby_dislocationSpacing<<defaultColor<<std::endl;
                std::cout<<yellowColor<<"   Read-Shockley_energyDensity= "<<p_gbType->ReadShockley_energyDensity<<defaultColor<<std::endl;
                
                
                //
                //                const double energyDensity;
                //                const double dislocationSpacing;
            }
            else
            {
                assert(0 && "GRAIN BONUDARY TYPE NOT FOUND.");
            }
            
            
        }
        

        
        /**********************************************************************/
        void populateGBdislocations(std::vector<StressStraight<dim>>& vD, const SimplicialMesh<dim>& mesh) const __attribute__ ((deprecated))
        {
            if(use_GBdislocations)
            {
                const int gbRegionID=grainBndID.first;
                
                const VectorDimD normal=latticePlane(gbRegionID).n.cartesian().normalized();
                const VectorDimD dir=rotationAxis().normalized();
                const VectorDimD p=dir.cross(latticePlane(gbRegionID).n.cartesian()).normalized()*grainBoundaryType().dislocationSpacing;
                
                std::cout<<"Grain boundary dislocations being added with spacing of "<<grainBoundaryType().dislocationSpacing<<" and along direction "<<p.normalized().transpose()<<std::endl;
                VectorDimD P0(VectorDimD::Zero());
                VectorDimD P1(VectorDimD::Zero());
                assert(mesh.search(P0).first    &&  mesh.search(P1).first && "Another method is needed to initialize the stress-straight segments if the mesh does not intersect the origin");
                bool positionOK=true;
                int k=-1;
                while(positionOK)
                {
                    k++;
                    P0=latticePlane(gbRegionID).P.cartesian()*0-100.0*dir+k*p;
                    P1=latticePlane(gbRegionID).P.cartesian()*0+100.0*dir+k*p;
                    LatticeVectorType L0=latticePlane(gbRegionID).snapToLattice(P0);
                    LatticeVectorType L1=latticePlane(gbRegionID).snapToLattice(P1);
                    LatticeDirection<dim> Dir(grain(gbRegionID).latticeDirection(dir));
                    //The original points have to first be within the mesh!!
                    if(  mesh.search(L0.cartesian()+Dir.cartesian()).first    &&
                       mesh.search(L0.cartesian()-Dir.cartesian()).first    &&
                       mesh.search(L1.cartesian()+Dir.cartesian()).first    &&
                       mesh.search(L1.cartesian()-Dir.cartesian()).first      )
                    {
                        //Determine the dislocation line - mesh intersection at either end of the projected dislocation position, and snap positions to that lattice point.
                        LatticeVectorType L0=latticePlane(gbRegionID).snapToLattice(P0);
                        LatticeVectorType L1=latticePlane(gbRegionID).snapToLattice(P1);
                        
                        LatticeLine line(L0,Dir);
                        LineMeshIntersection lmi(line,L0+Dir,mesh,mesh.search(P0).second);
                        LatticeDirection<dim> negativeDir(grain(gbRegionID).latticeDirection((-1.0)*dir));
                        LatticeLine line2(L0,negativeDir);
                        LineMeshIntersection lmi2(line2,L0-Dir,mesh,mesh.search(P0).second);
                        //
                        if(  mesh.search(lmi.L.cartesian()).second->region->regionID==gbRegionID  &&
                           mesh.search(lmi2.L.cartesian()).second->region->regionID==gbRegionID)
                        {
                            vD.emplace_back(lmi.L.cartesian(),lmi2.L.cartesian(),grainBoundaryType().Burgers*normal);
                        }
                        else
                        {
                            positionOK=false;
                        }
                    }
                }
                k=0;
                while(positionOK)
                {
                    k++;
                    P0=latticePlane(gbRegionID).P.cartesian()*0-1.0*dir-k*p;
                    P1=latticePlane(gbRegionID).P.cartesian()*0+1.0*dir-k*p;
                    LatticeVectorType L0=latticePlane(gbRegionID).snapToLattice(P0);
                    LatticeVectorType L1=latticePlane(gbRegionID).snapToLattice(P1);
                    LatticeDirection<dim> Dir(grain(gbRegionID).latticeDirection(dir));
                    //The original points have to first be within the mesh!!
                    if(  mesh.search(L0.cartesian()+Dir.cartesian()).first    &&
                       mesh.search(L0.cartesian()-Dir.cartesian()).first    &&
                       mesh.search(L1.cartesian()+Dir.cartesian()).first    &&
                       mesh.search(L1.cartesian()-Dir.cartesian()).first      )
                    {
                        //Determine the dislocation line - mesh intersection at either end of the projected dislocation position, and snap positions to that lattice point.
                        LatticeVectorType L0=latticePlane(gbRegionID).snapToLattice(P0);
                        LatticeVectorType L1=latticePlane(gbRegionID).snapToLattice(P1);
                        
                        LatticeLine line(L0,Dir);
                        LineMeshIntersection lmi(line,L0+Dir,mesh,mesh.search(P0).second);
                        LatticeDirection<dim> negativeDir(grain(gbRegionID).latticeDirection((-1.0)*dir));
                        LatticeLine line2(L0,negativeDir);
                        LineMeshIntersection lmi2(line2,L0-Dir,mesh,mesh.search(P0).second);
                        //
                        if(  mesh.search(lmi.L.cartesian()).second->region->regionID==gbRegionID  &&
                           mesh.search(lmi2.L.cartesian()).second->region->regionID==gbRegionID)
                        {
                            vD.emplace_back(lmi.L.cartesian(),lmi2.L.cartesian(),grainBoundaryType().Burgers*normal);
                        }
                        else
                        {
                            
                            positionOK=false;
                        }
                    }
                }
            }
        }
        
        
        CSL<dim> _csl;
        DSCL<dim> _dscl;
        VectorDimD _crystallographicRotationAxis;
        VectorDimD _rotationAxis;
        
        double cosTheta; // cosine of relative rotation angle between grains
        const GrainBoundaryType<dim>* p_gbType;
        
    public:
        
        static bool use_GBdislocations;
        const MeshRegionBoundaryType& regionBoundary;
        const std::pair<int,int>& grainBndID;
        
        
        /**********************************************************************/
        GrainBoundary(const MeshRegionBoundaryType& regionbnd_in,
                      const Grain<dim>& grainFirst,
                      const Grain<dim>& grainSecond) :
        /* init */ _csl(grainFirst,grainSecond),
        /* init */ _dscl(grainFirst,grainSecond),
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
            _csl.update();
            _dscl.update();
            computeCrystallographicRotationAxis();
            createLatticePlanes();
            findGrainBoundaryType(poly.grainBoundaryTypes());
            populateGBdislocations(poly.grainBoundaryDislocations(),poly.mesh);
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
        
        /**********************************************************************/
        const CSL<dim>& csl() const
        {
            return _csl;
        }
        
        /**********************************************************************/
        const DSCL<dim>& dscl() const
        {
            return _dscl;
        }
    };
    
    template <int dim>
    bool GrainBoundary<dim>::use_GBdislocations=false;
    
} // end namespace
#endif

                
                //        /**********************************************************************/
                //        void populateGBdislocations(std::vector<StressStraight<dim>>& vD) const
                //        {
                //            if(use_GBdislocations)
                //            {
                //            const VectorDimD dir=rotationAxis().normalized();
                //            const VectorDimD p=dir.cross(latticePlane(grainBndID.first).n.cartesian()).normalized()*grainBoundaryType().dislocationSpacing;
                //
                //
                //            VectorDimD P0(VectorDimD::Zero());
                //            VectorDimD P1(VectorDimD::Zero());
                //
                //            VectorDimD Q(latticePlane(grainBndID.first).P.cartesian());
                //            Q(1)=0.0;
                //
                //            for(int k=0;k<10;++k)
                //            {
                //                P0=Q-100.0*dir+k*p;
                //                P1=Q+100.0*dir+k*p;
                //                std::cout<<"GB DISLOCATIONS MAY BE WRONG SIGN"<<std::endl;
                //                vD.emplace_back(P0,P1,grainBoundaryType().Burgers*latticePlane(grainBndID.first).n.cartesian().normalized());
                //            }
                //            }
                //
                //        }
