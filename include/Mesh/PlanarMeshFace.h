/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarMeshFace_H_
#define model_PlanarMeshFace_H_

#include <chrono>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <utility>      // std::pair, std::make_pair
#include <set>

#include <Eigen/Dense>
#include <SimplexTraits.h>
#include <Simplex.h>
#include <ConvexHull.h>
#include <Plane.h>



namespace model
{
    template<int dim>
    struct PlanarMeshFace : public StaticID<PlanarMeshFace<dim>>
    /*                   */,public std::set<const Simplex<dim,dim-1>*>
    {
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Simplex<dim,dim-1> SimplexType;
        typedef std::set<const SimplexType*> SimplexContainerType;
        const std::pair<int,int> regionIDs;

    private:
        
        Eigen::Matrix<double,dim,1> n;
        Eigen::Matrix<double,dim,1> c;
        std::vector<const Simplex<dim,0>*> _hull;
        
    public:
        
        std::pair<VectorDim,const PlanarMeshFace<dim>* > periodicFacePair;

        PlanarMeshFace(const Simplex<dim,dim-1>* const pS);
        const std::set<const Simplex<dim,dim-1>*>& internalSimplices() const;
        const std::vector<const Simplex<dim,0>*>& convexHull() const;
        void finalize();
        const SimplexContainerType& simplices() const;
        SimplexContainerType& simplices();
        bool isExternal() const;
        Eigen::Matrix<double,dim,1> outNormal() const;
        Eigen::Matrix<double,dim,1> outNormal(const int& k) const;
        const Eigen::Matrix<double,dim,1>& center() const;
        Plane<dim> asPlane() const;
    };
    
}
#endif
