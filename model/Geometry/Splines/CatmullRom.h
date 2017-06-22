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

#include <Eigen/Dense>

namespace model
{
    
    template <int dim>
    struct CatmullRom
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        /**********************************************************************/
        template <typename T>
        static VectorDim loopTangent(const T& pair)
        {
            
            double cT=0.0;
            for(const auto& link : pair.second)
            {
                const double cL=link->pLink->parametricChordLength();
                cT+=cL;
            }
            
            VectorDim temp=VectorDim::Zero();
            for(const auto& link : pair.second)
            {
                const double cL=link->pLink->parametricChordLength();
                temp+=(link->sink()->get_P()-link->source()->get_P())/cL*(cT-cL)/cT;
            }
            
            return temp;
        }
        
    };
    
}
#endif
