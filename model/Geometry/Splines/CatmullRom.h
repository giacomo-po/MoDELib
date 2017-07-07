/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CatmullRom_H_
#define model_CatmullRom_H_

#include <set>
#include <Eigen/Dense>
#include <model/Geometry/Splines/SplineSegment.h>
#include <model/LoopNetwork/LoopLink.h>

namespace model
{
    
//    template <int dim>
    struct CatmullRom
    {
        
//        typedef  VectorDim;
        
        /**********************************************************************/
        template <typename T>
        static Eigen::Matrix<double,T::dim,1> loopTangent(const std::set<LoopLink<T>*>& segmentSet)
        {
            
            double cT=0.0;
            for(const auto& link : segmentSet)
            {
                const double cL=link->pLink->parametricChordLength();
                cT+=cL;
            }
            
            Eigen::Matrix<double,T::dim,1> temp=Eigen::Matrix<double,T::dim,1>::Zero();
            for(const auto& link : segmentSet)
            {
                const double cL=link->pLink->parametricChordLength();
                temp+=(link->sink()->get_P()-link->source()->get_P())/cL*(cT-cL)/cT;
            }
            
            return temp;
        }
        
    };
    
}
#endif
