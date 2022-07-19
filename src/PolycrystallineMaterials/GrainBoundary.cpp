/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GrainBoundary_cpp_
#define model_GrainBoundary_cpp_

#include <numbers>
#include <cfloat> // FLT_EPSILON
#include <vector>
#include <deque>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
//#include <RoundEigen.h>
#include <SimplicialMesh.h>
#include <MeshRegionObserver.h>
#include <Grain.h>
//#include <GrainBoundaryType.h>
#include <LatticePlane.h>
//#include <LineMeshIntersection.h>
//#include <CSL.h>
//#include <DSCL.h>
//#include <StressStraight.h>
//#include <GlidePlane.h>
#include <MeshPlane.h>
#include <GrainBoundary.h>


namespace model
{

    template <int dim>
    void GrainBoundary<dim>::computeCrystallographicRotationAxis()
    {
        // Compute the relative rotation matrix R such that C2G2=R*C2G1
        const MatrixDimD R(grain(grainBndID.first).singleCrystal->C2G.transpose()*grain(grainBndID.second).singleCrystal->C2G);
        
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
        
        _rotationAxis=grain(grainBndID.first).singleCrystal->C2G*_crystallographicRotationAxis;
        assert((_rotationAxis-grain(grainBndID.second).singleCrystal->C2G*_crystallographicRotationAxis).norm()<FLT_EPSILON && "rotationAxis inconsistency.");
        
        cosTheta=0.5*(R.trace()-1.0);
        
        std::cout<<"   Rotation axis="<<_crystallographicRotationAxis.transpose()<<std::endl;
        std::cout<<"   Rotation angle="<<acos(cosTheta)*180.0/std::numbers::pi<<" deg"<<defaultColor<<std::endl;
    }

    template <int dim>
    GrainBoundary<dim>::GrainBoundary(const MeshRegionBoundaryType& regionbnd_in,
                                      const std::shared_ptr<PlanarMeshFace<dim>>& face_in,
                                      const Grain<dim>& grainFirst,
                                      const Grain<dim>& grainSecond
                                      ) :
    //        /* init */ MeshPlane<dim>(getMeshPlane(regionbnd_in)),
    //        /* init */ MeshPlane<dim>(mesh,grainFirst.grainID,grainSecond.grainID),
    /* init */ MeshPlane<dim>(*face_in,grainFirst.grainID,grainSecond.grainID),
    //        /* init */ _csl(grainFirst,grainSecond),
    //        /* init */ _dscl(grainFirst,grainSecond),
    /* init */ _crystallographicRotationAxis(VectorDimD::Zero()),
    /* init */ _rotationAxis(_crystallographicRotationAxis),
    /* init */ cosTheta(1.0),
    //        /* init */ p_gbType(NULL),
    /* init */ regionBoundary(regionbnd_in),
    /* init */ face(face_in),
    /* init */ grainBndID(regionBoundary.regionBndID)
    {
        std::cout<<greenBoldColor<<"Creating GrainBoundary ("<<grainBndID.first<<" "<<grainBndID.second<<")"<<defaultColor<<std::endl;
        //            std::cout<<defaultColor<<"   CSL:  sigma= "<<_csl.sigma()<<std::endl;
        //            std::cout<<defaultColor<<"   DSCL: sigma= "<<_dscl.sigma()<<std::endl;
        
        // Populate pointers to parent Grains
        GrainContainerType::emplace(grainFirst.grainID,&grainFirst);
        GrainContainerType::emplace(grainSecond.grainID,&grainSecond);
        
        // Populate outNormals to grains
        const Simplex<dim,dim-1>& triangle(**regionBoundary.simplices().begin());
        for(const auto& tet : triangle.parents())
        {
            const size_t faceID=tet.second->childOrder(triangle.xID);
            const VectorDimD outNormal=tet.second->nda.col(faceID);
            const double outNorm(outNormal.norm());
            const size_t rID=tet.second->region->regionID;
            assert(outNorm>FLT_EPSILON && "Simplex normal has zero norm.");
            assert(rID==grainFirst.grainID || rID==grainSecond.grainID);
            OutNormalContainerType::emplace(rID,outNormal/outNorm);
        }
        
        // Initialize
        //            initializeGrainBoundary(dn,mesh);
        computeCrystallographicRotationAxis();
//        grainFirst.grainBoundaries().emplace(grainBndID,this);
//        grainSecond.grainBoundaries().emplace(grainBndID,this);
        
    }

    template <int dim>
    const typename GrainBoundary<dim>::GrainContainerType& GrainBoundary<dim>::grains() const
    {
        return *this;
    }

    template <int dim>
    const Grain<dim>& GrainBoundary<dim>::grain(const size_t& k) const
    {
        return *(grains().at(k));
    }

    template <int dim>
    const typename GrainBoundary<dim>::OutNormalContainerType& GrainBoundary<dim>::outNormals() const
    {
        return *this;
    }

    template <int dim>
    const typename GrainBoundary<dim>::VectorDimD& GrainBoundary<dim>::outNormal(const size_t& k) const
    {
        return outNormals().at(k);
    }

    template <int dim>
    const typename GrainBoundary<dim>::VectorDimD& GrainBoundary<dim>::crystallographicRotationAxis() const
    {
        return _crystallographicRotationAxis;
    }

    template <int dim>
    const typename GrainBoundary<dim>::VectorDimD& GrainBoundary<dim>::rotationAxis() const
    {
        return _rotationAxis;
    }

    template <int dim>
    double GrainBoundary<dim>::rotationAngle() const
    {
        return acos(cosTheta);
    }

    template <int dim>
    std::string GrainBoundary<dim>::tag() const
    {
        return "(" + std::to_string(grainBndID.first)+"," +std::to_string(grainBndID.second)+")";
    }

    template <int dim>
    bool GrainBoundary<dim>::use_GBdislocations=false;

    template class GrainBoundary<3>;

}
#endif


