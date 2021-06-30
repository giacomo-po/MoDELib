/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarPolygon_CPP_
#define model_PlanarPolygon_CPP_

#include <cfloat>
#include <tuple>
#include <map>
#include <deque>
#include <Eigen/Dense>
#include <assert.h>
#include <earcut.hpp>

#include <PlanarPolygon.h>

namespace model
{
    /**********************************************************************/
    Eigen::Matrix3d PlanarPolygon::getR(const Eigen::Vector3d& x,
                                const Eigen::Vector3d& z)
    {
        const double xNorm(x.norm());
        const double zNorm(z.norm());
        assert(xNorm>FLT_EPSILON);
        assert(zNorm>FLT_EPSILON);
        assert(fabs(x.dot(z)<FLT_EPSILON*xNorm*zNorm));
        Eigen::Matrix3d temp(Eigen::Matrix3d::Identity());
        temp.col(2)=z/zNorm;
        temp.col(0)=x/xNorm;
        temp.col(1)=temp.col(2).cross(temp.col(0));
        return temp;
    }
    
    /**********************************************************************/
    std::array<double, 2> PlanarPolygon::projectRotate(const Eigen::Vector3d& v)
    {
        const Eigen::Vector3d vpr=R.transpose()*(v-v.dot(n)*n); // projected point
        return std::array<double, 2>({vpr(0),vpr(1)});
    }
    
    /**********************************************************************/
    PlanarPolygon::PlanarPolygon(const Eigen::Vector3d& x,        // global coordinate of the x axis on the plane
                  const Eigen::Vector3d& z) :      // global z coordinate
    R(getR(x,z)),
    n(R.col(2))
    {
        
    }
    
    /**********************************************************************/
    void PlanarPolygon::assignPoints(const std::deque<Eigen::Vector3d>& points)
    {
        //assert(this->size()==0 && "Polygon must be empty before adding external points. Call clear first.");
        this->clear();
        this->resize(1); // first vector is external polyline, further vectors define holes
        for(const auto& v3 : points)
        {
            this->operator[](0).push_back(projectRotate(v3));
        }
    }
    
    /**********************************************************************/
    void PlanarPolygon::addHole(const std::deque<Eigen::Vector3d>& hole)
    {
        assert(this->size()==1 && "Call assignPoints before adding holes.");
        this->push_back(std::deque<std::array<double, 2>>()); // first vector is external polyline, further vectors define holes
        for(const auto& v3 : hole)
        {
            this->rbegin()->push_back(projectRotate(v3));
        }
    }
    
    /**********************************************************************/
    std::deque<std::array<size_t, 3>> PlanarPolygon::triangulate() const
    {
        std::vector<size_t> indices=mapbox::earcut<size_t>(*this);
        assert((indices.size()%3)==0);
        const size_t nTri=indices.size()/3;
        std::deque<std::array<size_t, 3>> temp;
        for(size_t k=0;k<nTri;++k)
        {
            temp.push_back(std::array<size_t, 3>({indices[3*k+0],indices[3*k+1],indices[3*k+2]}));
        }
        return temp;
    }

    
} /* namespace model */
#endif
