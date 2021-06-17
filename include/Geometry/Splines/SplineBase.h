/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPLINEBASE_H_
#define model_SPLINEBASE_H_

namespace model
{
    
    template<int dim,int corder>
    struct SplineBase
    {
        static constexpr double chordal=1.0;
        static constexpr double centripetal=0.5;
        static constexpr double uniform=0.0;
        
        static constexpr int Ncoeff= 2*(corder+1);
        static constexpr int  pOrder= 2*corder+1;
        static constexpr int  Ndof  = dim*Ncoeff;
        //        static constexpr int  eigenSize=pOrder*pOrder;
        
    };
        
}

#endif
