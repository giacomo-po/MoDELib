/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoop_cpp_
#define model_DislocationLoop_cpp_


#include <DislocationLoop.h>

namespace model
{
    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationLoop<dim,corder,InterpolationType>::DislocationLoop(LoopNetworkType* const net,
                                     const FlowType& flow) :
    /* init */ Loop<DislocationLoop>(net,flow)
    {
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    std::shared_ptr<typename TypeTraits<DislocationLoop<dim,corder,InterpolationType>>::LoopType> DislocationLoop<dim,corder,InterpolationType>::clone() const
    {
        return std::shared_ptr<typename TypeTraits<DislocationLoop<dim,corder,InterpolationType>>::LoopType>(new DislocationLoop(this->p_network(),this->flow()));
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationLoop<dim,corder,InterpolationType>::initFromFile(const std::string& fileName)
    {
        verboseDislocationLoop=TextFileParser(fileName).readScalar<int>("verboseDislocationLoop",true);
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    int DislocationLoop<dim,corder,InterpolationType>::verboseDislocationLoop=0;

    
}
#endif
