/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PlaneLineIntersection_h_
#define model_PlaneLineIntersection_h_

#include <model/LatticeMath/LatticeVector.h>
#include <model/LatticeMath/LatticePlane.h>
#include <model/LatticeMath/LatticeLine.h>

namespace model
{
    
    struct PlaneLineIntersection
    {
        typedef Eigen::Matrix<long int,3,1> VectorDimI;

        enum IntersectionType
        {
            parallel=0,
            coincident=1,
            intersecting=2,
            offLattice=3
        };
        
        const long int num;
        const long int den;
        const IntersectionType intersectionType;
        const LatticeVector<3> P;

        /**********************************************************************/
        PlaneLineIntersection(const LatticePlane& plane, const LatticeLine& line) :
        /* init */ num(plane.P.dot(plane.n)-line.P.dot(plane.n)),
        /* init */ den(line.d.dot(plane.n)),
        /* init */ intersectionType(den!=0? ((num%den)==0? intersecting : offLattice) : (num==0? coincident : parallel)),
        /* init */ P( (intersectionType==intersecting || intersectionType==offLattice )? (line.P+num/den*line.d).eval() : (intersectionType==coincident? line.P : VectorDimI::Zero() ))
        {
//            std::cout<<"num="<<num<<std::endl;
//            std::cout<<"den="<<den<<std::endl;
        }
        
        /**********************************************************************/
        friend std::ostream& operator << (std::ostream& os, const PlaneLineIntersection& pli)
        {
            os  <<pli.num<<"\n"
            /**/<<pli.den<<"\n"
            /**/<<pli.intersectionType<<"\n"
            /**/<<pli.P.cartesian().transpose();
            return os;
        }
        
    };
    
} // end namespace
#endif
