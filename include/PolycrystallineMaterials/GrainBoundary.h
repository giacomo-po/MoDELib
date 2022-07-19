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


namespace model
{

    template <int dim>
    class GrainBoundary : public MeshPlane<dim>,
    /* base            */ private std::map<size_t,const Grain<dim>* const>,
    /* base            */ private std::map<size_t,const Eigen::Matrix<double,dim,1>>
    {
        enum AxisPlaneRelation{TILT=0,TWIST=1,MIXED=2};

        typedef MeshRegionBoundary<dim> MeshRegionBoundaryType;
        typedef Eigen::Matrix<long int,dim,1> VectorDimI;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef ReciprocalLatticeDirection<dim> ReciprocalLatticeDirectionType;
        typedef Grain<dim> GrainType;
        typedef std::map<size_t,const GrainType* const> GrainContainerType;
        typedef std::map<size_t,const Eigen::Matrix<double,dim,1>> OutNormalContainerType;
        typedef Eigen::Matrix<double,2*dim+1,1> BGkeyType;

        void computeCrystallographicRotationAxis();

        VectorDimD _crystallographicRotationAxis;
        VectorDimD _rotationAxis;
        double cosTheta; // cosine of relative rotation angle between grains

    public:

        static bool use_GBdislocations;
        const MeshRegionBoundaryType& regionBoundary;
        const std::shared_ptr<PlanarMeshFace<dim>>& face;
        const std::pair<int,int>& grainBndID;

        GrainBoundary(const MeshRegionBoundaryType& regionbnd_in,
                      const std::shared_ptr<PlanarMeshFace<dim>>& face_in,
                      const Grain<dim>& grainFirst,
                      const Grain<dim>& grainSecond
                      );
        const GrainContainerType& grains() const;
        const Grain<dim>& grain(const size_t& k) const;
        const OutNormalContainerType& outNormals() const;
        const VectorDimD& outNormal(const size_t& k) const;
        const VectorDimD& crystallographicRotationAxis() const;
        const VectorDimD& rotationAxis() const;
        double rotationAngle() const;
        std::string tag() const;
    };

}
#endif


