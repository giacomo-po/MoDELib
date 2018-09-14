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
//#include <model/LatticeMath/LineMeshIntersection.h>
#include <model/LatticeMath/CSL.h>
#include <model/LatticeMath/DSCL.h>
#include <model/DislocationDynamics/ElasticFields/StressStraight.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlane.h>
#include <model/Mesh/MeshPlane.h>


namespace model
{
    
    
    
    template <int dim>
    class GrainBoundary : public MeshPlane<dim>,
    /* base            */ private std::map<int,const Grain<dim>* const>,
    /* base            */ private std::map<int,const Eigen::Matrix<double,dim,1>>
//    /* base */ private std::map<int,std::shared_ptr<GlidePlane<dim>>>
    {
//        static constexpr int dim=NetworkType::dim;
        
        enum AxisPlaneRelation{TILT=0,TWIST=1,MIXED=2};
        
        typedef MeshRegionBoundary<Simplex<dim,dim-1> > MeshRegionBoundaryType;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        typedef Grain<dim> GrainType;
        typedef std::map<int,const GrainType* const> GrainContainerType;
        typedef std::map<int,const Eigen::Matrix<double,dim,1>> OutNormalContainerType;
        typedef std::map<int,std::shared_ptr<GlidePlane<dim>>> GlidePlaneContainerType;
        typedef Eigen::Matrix<double,2*dim+1,1> BGkeyType;
        typedef GlidePlane<dim> GlidePlaneType;
        
        
//        static_assert(false,"NEED TO CREATE GLIDE PLANES FOR DSCL AND CSL, AND ADD THEM TO NODES WHICH ARE ON THE GB. THEN ADD THEM TO THE NODES. CAREFUL THAT REMOVE LOOP LINK DO NOT REMOVE THESE GLIDE PLANES");
//        static_assert(false,"IN DDVTK SHOW THE MESH FOR THE GBs, and let the GP SHADE THE PLANE");
        
        
//        /**********************************************************************/
////        template<typename NetworkType>
//        void storeLatticePlane(GlidePlaneObserver<dim>& dn,
//                               const SimplicialMesh<dim>& mesh,
//                               const Grain<dim>& grain,
//                               const VectorDimD& normal)
//        {
//            //std::cout<<"storeLatticePlane"<<std::endl;
//            //std::cout<<"normal="<<normal.transpose()<<std::endl;
//            
//            const ReciprocalLatticeDirectionType R=grain.reciprocalLatticeDirection(normal);
//            
//            model::cout<<"   GB plane normal for grain "<< grain.grainID<<":"<<defaultColor<<std::endl;
//            model::cout<<"      cartesian components="<<R.cartesian().normalized().transpose()<<defaultColor<<std::endl;
//            model::cout<<"      crystallographic components="<<grain.rationalApproximation((grain.C2G.transpose()*R.cartesian())).transpose()<<defaultColor<<std::endl;
////            model::cout<<"   crystallographic components="<<(grain.C2G.transpose()*R.cartesian()).transpose()<<defaultColor<<std::endl;
//            model::cout<<"      interplanar spacing="<<1.0/R.cartesian().norm()<<defaultColor<<std::endl;
//            
////            LatticeVectorType L0(grain.lattice());
////            bool latticePointFound=false;
////            
//////            //std::cout<<"here 0"<<std::endl;
////            
////            for(const auto& simplexPtr : regionBoundary.simplices())
////            {
//////                            //std::cout<<"here 1"<<std::endl;
////                for(size_t d=0;d<dim;++d)
////                {
////                    const VectorDimD P0=simplexPtr->vertexPositionMatrix().col(d);
////                    L0=grain.snapToLattice(P0);
////                    
////                    if(fabs((L0.cartesian()-P0).dot(normal))<FLT_EPSILON) // L0 belongs to mesh plane
////                    {
////                        std::cout<<"GrainBonudary: using origin "<<P0.transpose()<<std::endl;
////                        latticePointFound=true;
////                        break; //break inner loop
////                    }
////                }
////                
////                if(latticePointFound==true) // L0 belongs to mesh plane
////                {
////                    break; //break outer loop
////                }
////            }
////            
////            
////            
////            assert(latticePointFound && "None of the GB mesh vertices, snapped to the lattice, belongs to GB plane.");
////            
////            LatticePlaneBase pb(R);
////            const auto temp = GlidePlaneContainerType::emplace(std::piecewise_construct,
////                                                                 std::forward_as_tuple(grain.grainID),
////                                                                 std::forward_as_tuple(&dn,dn.mesh,grain,L0.cartesian(),pb.cartesian()));
////            assert(temp.first->second.unitNormal.cross(normal.normalized()).norm()<FLT_EPSILON && "LatticePlane normal and triangle normal are not the same.");
//            
//            const VectorDimD P0=(*regionBoundary.simplices().begin())->vertexPositionMatrix().col(0); // a point on the plane
////            const auto temp = GlidePlaneContainerType::emplace(std::piecewise_construct,
////                                                               std::forward_as_tuple(grain.grainID),
////                                                               std::forward_as_tuple(&dn,mesh,grain,P0,normal));
////            assert(temp.first->second.unitNormal.cross(normal.normalized()).norm()<FLT_EPSILON && "LatticePlane normal and triangle normal are not the same.");
//
//            const auto temp = GlidePlaneContainerType::emplace(std::piecewise_construct,
//                                                               std::forward_as_tuple(grain.grainID),
//                                                               std::forward_as_tuple(new GlidePlaneType(&dn,mesh,grain.lattice(),grain.grainID,grain.grainID,P0,normal)));
//
//            assert(temp.second);
//            assert(temp.first->second->unitNormal.cross(normal.normalized()).norm()<FLT_EPSILON && "LatticePlane normal and triangle normal are not the same.");
//            temp.first->second->addParentSharedPtr(&temp.first->second,false,-1); // add shared pointer to GlidePlane
//            
//        }
        
        /**********************************************************************/
//        template<typename NetworkType>

//        static MeshPlane<dim> getMeshPlane(const MeshRegionBoundaryType& regionBnd)
//        {
//            //std::cout<<"createLatticePlanes"<<std::endl;
//            
//            
////            GlidePlaneContainerType::clear();
//            const Simplex<dim,dim-1>& triangle(**regionBnd.begin());
//            const Simplex<dim,dim>& tet(**triangle.parents().begin());
//            const size_t faceID=tet.childOrder(triangle.xID);
//            const VectorDimD N=tet.nda.col(faceID);
//            const VectorDimD P=triangle.vertexPositionMatrix().col(0);
//            
//            MeshPlane<dim> mp(P,N);
//
////            for(const auto& tet : triangle.parents())
////            {
////                const size_t faceID=tet->childOrder(triangle.xID);
////                const VectorDimD outNormal=tet->nda.col(faceID);
////                const double outNorm(outNormal.norm());
////                assert(outNorm>FLT_EPSILON && "Simplex normal has zero norm.");
////                storeLatticePlane(dn,mesh,grain(tet->region->regionID),outNormal/outNorm);
////            }
////            
////            assert(glidePlanes().size()==2);
////            assert((glidePlanes().begin()->second->unitNormal+glidePlanes().rbegin()->second->unitNormal).norm()<FLT_EPSILON && "plane normals are not opposite to each other.");
//            
//            // Check that all tringle vertices are contained by both GB planes
//            for(const auto& triangle : regionBnd)
//            {
//                const auto Ps=triangle->vertexPositionMatrix();
//                for(int j=0;j<Ps.cols();++j)
//                {
//                                            assert(mp.contains(Ps.col(j)) && "TRIANGLE VERTEX NOT CONTAINED IN GBPLANE");
////                    for(const auto& grain : grains())
////                    {
////                        assert(glidePlane(grain.second->grainID).contains(Ps.col(j)) && "TRIANGLE VERTEX NOT CONTAINED IN GBPLANE");
////                    }
//                }
//            }
//            
//            return mp;
//        }
        
//        void createLatticePlanes(GlidePlaneObserver<dim>& dn,
//                                 const SimplicialMesh<dim>& mesh)
//        {
//            //std::cout<<"createLatticePlanes"<<std::endl;
//
//            
//            GlidePlaneContainerType::clear();
//            const Simplex<dim,dim-1>& triangle(**regionBoundary.begin());
//            for(const auto& tet : triangle.parents())
//            {
//                const size_t faceID=tet->childOrder(triangle.xID);
//                const VectorDimD outNormal=tet->nda.col(faceID);
//                const double outNorm(outNormal.norm());
//                assert(outNorm>FLT_EPSILON && "Simplex normal has zero norm.");
//                storeLatticePlane(dn,mesh,grain(tet->region->regionID),outNormal/outNorm);
//            }
//            
//            assert(glidePlanes().size()==2);
//            assert((glidePlanes().begin()->second->unitNormal+glidePlanes().rbegin()->second->unitNormal).norm()<FLT_EPSILON && "plane normals are not opposite to each other.");
//            
//            // Check that all tringle vertices are contained by both GB planes
//            for(const auto& triangle : regionBoundary)
//            {
//                const auto Ps=triangle->vertexPositionMatrix();
//                for(int j=0;j<Ps.cols();++j)
//                {
//                    for(const auto& grain : grains())
//                    {
//                        assert(glidePlane(grain.second->grainID).contains(Ps.col(j)) && "TRIANGLE VERTEX NOT CONTAINED IN GBPLANE");
//                    }
//                }
//            }
//            
//            if(fabs(_rotationAxis.dot(glidePlanes().begin()->second->unitNormal))<FLT_EPSILON)
//            {
//                model::cout<<"TILT BOUNDARY"<<std::endl;
//                // here check mirror symm equation r=d−2(d⋅n)n to check if symm tilt or asymm tilt
//            }
//            else if(_rotationAxis.cross(glidePlanes().begin()->second->unitNormal).norm()<FLT_EPSILON)
//            {
//                model::cout<<"TWIST BOUNDARY"<<std::endl;
//
//            }
//            else
//            {
//                model::cout<<"MIXED BOUNDARY"<<std::endl;
//
//            }
//            
//            
//        }
        
        /**********************************************************************/
        void computeCrystallographicRotationAxis()
        {
            // Compute the relative rotation matrix R such that C2G2=R*C2G1
            const MatrixDimD R(grain(grainBndID.first).C2G.transpose()*grain(grainBndID.second).C2G);
            
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
            
            _rotationAxis=grain(grainBndID.first).C2G*_crystallographicRotationAxis;
            assert((_rotationAxis-grain(grainBndID.second).C2G*_crystallographicRotationAxis).norm()<FLT_EPSILON && "rotationAxis inconsistency.");
            
            cosTheta=0.5*(R.trace()-1.0);
            
            model::cout<<"   Rotation axis="<<_crystallographicRotationAxis.transpose()<<std::endl;
            model::cout<<"   Rotation angle="<<acos(cosTheta)*180.0/M_PI<<" deg"<<defaultColor<<std::endl;
        }
        
//        /**********************************************************************/
//        void findGrainBoundaryType(const std::deque<GrainBoundaryType<dim>>& bgTypes)
//        {
//            //std::cout<<"findGrainBoundaryType"<<std::endl;
//            
//            const VectorDimD n1=(grain(grainBndID.first).C2G.transpose()*glidePlane(grainBndID.first).unitNormal).normalized();
//            const VectorDimD n2=(grain(grainBndID.second).C2G.transpose()*glidePlane(grainBndID.second).unitNormal).normalized();
//            
//            
//            for (const auto& gbt : bgTypes)
//            {
//                bool axisFound=false;
//                for(const auto& a : gbt.axisPermutations)
//                {
//                    axisFound=( (_crystallographicRotationAxis-a.normalized()).norm()<FLT_EPSILON ||
//                               (_crystallographicRotationAxis+a.normalized()).norm()<FLT_EPSILON );
//                    if(axisFound)
//                    {
//                        break;
//                    }
//                }
//                
//                bool n1Found=false;
//                for(const auto& n : gbt.n1Permutations)
//                {
//                    n1Found=( (n1-n.normalized()).norm()<FLT_EPSILON || (n1+n.normalized()).norm()<FLT_EPSILON );
//                    if(n1Found)
//                    {
//                        break;
//                    }
//                }
//                
//                bool n2Found=false;
//                for(const auto& n : gbt.n2Permutations)
//                {
//                    n2Found=( (n2-n.normalized()).norm()<FLT_EPSILON || (n2+n.normalized()).norm()<FLT_EPSILON );
//                    if(n2Found)
//                    {
//                        break;
//                    }
//                }
//                
//                if(!n1Found && ! n2Found) // search n1 in n2Permutations and viceversa
//                {
//                    for(const auto& n : gbt.n2Permutations)
//                    {
//                        n1Found=( (n1-n.normalized()).norm()<FLT_EPSILON || (n1+n.normalized()).norm()<FLT_EPSILON );
//                        if(n1Found)
//                        {
//                            break;
//                        }
//                    }
//                    
//                    for(const auto& n : gbt.n1Permutations)
//                    {
//                        n2Found=( (n2-n.normalized()).norm()<FLT_EPSILON || (n2+n.normalized()).norm()<FLT_EPSILON );
//                        if(n2Found)
//                        {
//                            break;
//                        }
//                    }
//                }
//                
//                if(axisFound && n1Found && n2Found)
//                {
//                    p_gbType=&gbt;
//                    break;
//                }
//                
//            }
//            
//            if(p_gbType!=NULL)
//            {
//                model::cout<<yellowColor<<"   GB type is "<<p_gbType->name<<defaultColor<<std::endl;
//                model::cout<<yellowColor<<"   GB energy= "<<p_gbType->energyDensity<<defaultColor<<std::endl;
//                model::cout<<yellowColor<<"   GB dislocation spacing= "<<p_gbType->dislocationSpacing<<defaultColor<<std::endl;
//                model::cout<<yellowColor<<"   Frank-Bilby dislocation spacing= "<<p_gbType->FrankBilby_dislocationSpacing<<defaultColor<<std::endl;
//                model::cout<<yellowColor<<"   Read-Shockley_energyDensity= "<<p_gbType->ReadShockley_energyDensity<<defaultColor<<std::endl;
//                
//                
//                //
//                //                const double energyDensity;
//                //                const double dislocationSpacing;
//            }
//            else
//            {
//                assert(0 && "GRAIN BONUDARY TYPE NOT FOUND.");
//            }
//            
//            
//        }
        

        
        /**********************************************************************/
        void populateGBdislocations(std::vector<StressStraight<dim>>& vD, const SimplicialMesh<dim>& mesh) const __attribute__ ((deprecated))
        {
            //std::cout<<"populateGBdislocations"<<std::endl;

            
            if(use_GBdislocations)
            {
                assert(0 && "Re-Enable this");
                
//                const int gbRegionID=grainBndID.first;
//                
//                const VectorDimD normal=latticePlane(gbRegionID).unitNormal.normalized();
//                const VectorDimD dir=rotationAxis().normalized();
//                const VectorDimD p=dir.cross(latticePlane(gbRegionID).unitNormal).normalized()*grainBoundaryType().dislocationSpacing;
//                
//                model::cout<<"Grain boundary dislocations being added with spacing of "<<grainBoundaryType().dislocationSpacing<<" and along direction "<<p.normalized().transpose()<<std::endl;
//                VectorDimD P0(VectorDimD::Zero());
//                VectorDimD P1(VectorDimD::Zero());
//                assert(mesh.search(P0).first    &&  mesh.search(P1).first && "Another method is needed to initialize the stress-straight segments if the mesh does not intersect the origin");
//                bool positionOK=true;
//                int k=-1;
//                while(positionOK)
//                {
//                    k++;
//                    P0=latticePlane(gbRegionID).P.cartesian()*0-100.0*dir+k*p;
//                    P1=latticePlane(gbRegionID).P.cartesian()*0+100.0*dir+k*p;
//                    LatticeVectorType L0=latticePlane(gbRegionID).snapToLattice(P0);
//                    LatticeVectorType L1=latticePlane(gbRegionID).snapToLattice(P1);
//                    LatticeDirection<dim> Dir(grain(gbRegionID).latticeDirection(dir));
//                    //The original points have to first be within the mesh!!
//                    if(  mesh.search(L0.cartesian()+Dir.cartesian()).first    &&
//                       mesh.search(L0.cartesian()-Dir.cartesian()).first    &&
//                       mesh.search(L1.cartesian()+Dir.cartesian()).first    &&
//                       mesh.search(L1.cartesian()-Dir.cartesian()).first      )
//                    {
//                        //Determine the dislocation line - mesh intersection at either end of the projected dislocation position, and snap positions to that lattice point.
//                        LatticeVectorType L0=latticePlane(gbRegionID).snapToLattice(P0);
//                        LatticeVectorType L1=latticePlane(gbRegionID).snapToLattice(P1);
//                        
//                        LatticeLine line(L0,Dir);
//                        LineMeshIntersection lmi(line,L0+Dir,mesh,mesh.search(P0).second);
//                        LatticeDirection<dim> negativeDir(grain(gbRegionID).latticeDirection((-1.0)*dir));
//                        LatticeLine line2(L0,negativeDir);
//                        LineMeshIntersection lmi2(line2,L0-Dir,mesh,mesh.search(P0).second);
//                        //
//                        if(  mesh.search(lmi.L.cartesian()).second->region->regionID==gbRegionID  &&
//                           mesh.search(lmi2.L.cartesian()).second->region->regionID==gbRegionID)
//                        {
//                            vD.emplace_back(lmi.L.cartesian(),lmi2.L.cartesian(),grainBoundaryType().Burgers*normal);
//                        }
//                        else
//                        {
//                            positionOK=false;
//                        }
//                    }
//                }
//                k=0;
//                while(positionOK)
//                {
//                    k++;
//                    P0=latticePlane(gbRegionID).P.cartesian()*0-1.0*dir-k*p;
//                    P1=latticePlane(gbRegionID).P.cartesian()*0+1.0*dir-k*p;
//                    LatticeVectorType L0=latticePlane(gbRegionID).snapToLattice(P0);
//                    LatticeVectorType L1=latticePlane(gbRegionID).snapToLattice(P1);
//                    LatticeDirection<dim> Dir(grain(gbRegionID).latticeDirection(dir));
//                    //The original points have to first be within the mesh!!
//                    if(  mesh.search(L0.cartesian()+Dir.cartesian()).first    &&
//                       mesh.search(L0.cartesian()-Dir.cartesian()).first    &&
//                       mesh.search(L1.cartesian()+Dir.cartesian()).first    &&
//                       mesh.search(L1.cartesian()-Dir.cartesian()).first      )
//                    {
//                        //Determine the dislocation line - mesh intersection at either end of the projected dislocation position, and snap positions to that lattice point.
//                        LatticeVectorType L0=latticePlane(gbRegionID).snapToLattice(P0);
//                        LatticeVectorType L1=latticePlane(gbRegionID).snapToLattice(P1);
//                        
//                        LatticeLine line(L0,Dir);
//                        LineMeshIntersection lmi(line,L0+Dir,mesh,mesh.search(P0).second);
//                        LatticeDirection<dim> negativeDir(grain(gbRegionID).latticeDirection((-1.0)*dir));
//                        LatticeLine line2(L0,negativeDir);
//                        LineMeshIntersection lmi2(line2,L0-Dir,mesh,mesh.search(P0).second);
//                        //
//                        if(  mesh.search(lmi.L.cartesian()).second->region->regionID==gbRegionID  &&
//                           mesh.search(lmi2.L.cartesian()).second->region->regionID==gbRegionID)
//                        {
//                            vD.emplace_back(lmi.L.cartesian(),lmi2.L.cartesian(),grainBoundaryType().Burgers*normal);
//                        }
//                        else
//                        {
//                            
//                            positionOK=false;
//                        }
//                    }
//                }
            }
        }
        
        
//        CSL<dim> _csl;
//        DSCL<dim> _dscl;
        VectorDimD _crystallographicRotationAxis;
        VectorDimD _rotationAxis;
        
        double cosTheta; // cosine of relative rotation angle between grains
        const GrainBoundaryType<dim>* p_gbType;
        
    public:
        
        static bool use_GBdislocations;
        const MeshRegionBoundaryType& regionBoundary;
        const std::pair<int,int>& grainBndID;
        
        
        /**********************************************************************/
//        template<typename NetworkType>
        GrainBoundary(const MeshRegionBoundaryType& regionbnd_in,
                      Grain<dim>& grainFirst,
                      Grain<dim>& grainSecond,
                      GlidePlaneObserver<dim>& dn,
                      const SimplicialMesh<dim>& mesh) :
//        /* init */ MeshPlane<dim>(getMeshPlane(regionbnd_in)),
        /* init */ MeshPlane<dim>(mesh,grainFirst.grainID,grainSecond.grainID),
//        /* init */ _csl(grainFirst,grainSecond),
//        /* init */ _dscl(grainFirst,grainSecond),
        /* init */ _crystallographicRotationAxis(VectorDimD::Zero()),
        /* init */ _rotationAxis(_crystallographicRotationAxis),
        /* init */ cosTheta(1.0),
        /* init */ p_gbType(NULL),
        /* init */ regionBoundary(regionbnd_in),
        /* init */ grainBndID(regionBoundary.regionBndID)
        {
            model::cout<<greenBoldColor<<"Creating GrainBoundary ("<<grainBndID.first<<" "<<grainBndID.second<<")"<<defaultColor<<std::endl;
//            model::cout<<defaultColor<<"   CSL:  sigma= "<<_csl.sigma()<<std::endl;
//            model::cout<<defaultColor<<"   DSCL: sigma= "<<_dscl.sigma()<<std::endl;
            
            // Populate pointers to parent Grains
            GrainContainerType::emplace(grainFirst.grainID,&grainFirst);
            GrainContainerType::emplace(grainSecond.grainID,&grainSecond);

            // Populate outNormals to grains
            const Simplex<dim,dim-1>& triangle(**regionBoundary.begin());
            for(const auto& tet : triangle.parents())
            {
                const size_t faceID=tet->childOrder(triangle.xID);
                const VectorDimD outNormal=tet->nda.col(faceID);
                const double outNorm(outNormal.norm());
                const int rID=tet->region->regionID;
                assert(outNorm>FLT_EPSILON && "Simplex normal has zero norm.");
                assert(rID==grainFirst.grainID || rID==grainSecond.grainID);
                OutNormalContainerType::emplace(rID,outNormal/outNorm);
            }

            // Initialize
            initializeGrainBoundary(dn,mesh);
            grainFirst.grainBoundaries().emplace(grainBndID,this);
            grainSecond.grainBoundaries().emplace(grainBndID,this);

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
        const OutNormalContainerType& outNormals() const
        {
            return *this;
        }
        
        const VectorDimD& outNormal(const int& k) const
        {
            return outNormals().at(k);
        }
        
//        /**********************************************************************/
//        const GlidePlaneContainerType& glidePlanes() const
//        {
//            return *this;
//        }
        
//        /**********************************************************************/
//        const GlidePlaneType& glidePlane(const int& k) const
//        {
//            return *glidePlanes().at(k).get();
//        }
        
        /**********************************************************************/
//        template<typename NetworkType>
        void initializeGrainBoundary(GlidePlaneObserver<dim>& dn,
                                     const SimplicialMesh<dim>& mesh)
        {
//            model::cout<<yellowColor<<"GrainBoundary ("<<grainBndID.first<<" "<<grainBndID.second<<")"<<defaultColor<<std::endl;
//            _csl.update(true);
//            _dscl.update(true);
            computeCrystallographicRotationAxis();
//            createLatticePlanes(dn,mesh);
//            findGrainBoundaryType(poly.grainBoundaryTypes());
//            populateGBdislocations(dn.poly.grainBoundaryDislocations(),dn.poly.mesh);
            
            //                grain(gb.first.second).emplace(gb.first,&gb.second);
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
        
//        /**********************************************************************/
//        const CSL<dim>& csl() const
//        {
//            return _csl;
//        }
//        
//        /**********************************************************************/
//        const DSCL<dim>& dscl() const
//        {
//            return _dscl;
//        }
        
        /**********************************************************************/
        std::string tag() const
        {
            return "(" + std::to_string(grainBndID.first)+"," +std::to_string(grainBndID.second)+")";
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
                //            const VectorDimD p=dir.cross(latticePlane(grainBndID.first).unitNormal).normalized()*grainBoundaryType().dislocationSpacing;
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
                //                //std::cout<<"GB DISLOCATIONS MAY BE WRONG SIGN"<<std::endl;
                //                vD.emplace_back(P0,P1,grainBoundaryType().Burgers*latticePlane(grainBndID.first).unitNormal.normalized());
                //            }
                //            }
                //
                //        }
