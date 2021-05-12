/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlanarPolygon_H_
#define model_PlanarPolygon_H_

#include <cfloat>
#include <tuple>
#include <map>
#include <Eigen/Dense>
#include <assert.h>
#include <earcut.hpp>

namespace model
{
    
    class PlanarPolygon : public std::deque<std::deque<std::array<double, 2>>>
    {
        
        using Point2d = std::array<double, 2>;
        
    public:
        
        const Plane<3> plane;
        
        PlanarPolygon(const Eigen::Vector3d& p,const Eigen::Vector3d& n) :
        /* init */ plane(p,n)
        {
            
        }
        
        /**********************************************************************/
        void assignPoints(const std::deque<Eigen::Vector3d>& points)
        {
            //assert(this->size()==0 && "Polygon must be empty before adding external points. Call clear first.");
            this->clear();
            this->resize(1); // first vector is external polyline, further vectors define holes
            for(const auto& v3 : points)
            {
                //this->operator[](0).push_back(projectRotate(v3));
                const auto v2(plane.localPosition(v3));
                this->operator[](0).push_back(Point2d{v2(0),v2(1)});
            }
        }
        
        /**********************************************************************/
        void addHole(const std::deque<Eigen::Vector3d>& hole)
        {
            assert(this->size()==1 && "Call assignPoints before adding holes.");
            this->push_back(std::deque<std::array<double, 2>>()); // first vector is external polyline, further vectors define holes
            for(const auto& v3 : hole)
            {
                const auto v2(plane.localPosition(v3));
                this->rbegin()->push_back(Point2d{v2(0),v2(1)});
            }
        }
        
        /**********************************************************************/
        std::deque<std::array<size_t, 3>> triangulate() const
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
        
        
    };
    
}
#endif
