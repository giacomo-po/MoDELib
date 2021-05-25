/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2019 by Yash Pachaury      <ypachaur@purdue.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopLink_cpp_
#define model_DislocationLoopLink_cpp_

#include <DislocationLoopLink.h>

namespace model
{
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationLoopLink<dim,corder,InterpolationType>::DislocationLoopLink(const std::shared_ptr<LoopNodeType>& so,
                                             const std::shared_ptr<LoopNodeType>& si,
                                             const std::shared_ptr<LoopType>& pL) :
    /* init */ LoopLink<LoopLinkType>(so,si,pL)
    {
        
        
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationLoopLink<dim,corder,InterpolationType>::hasNetworkLink() const
    {
        return !((this->source->get_P()-this->sink->get_P()).norm()<FLT_EPSILON
        /*    */ && this->source->periodicPlaneEdge
        /*    */ &&   this->sink->periodicPlaneEdge
        /*    */ && this->source->loop()->sID==this->sink->loop()->sID);
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    std::shared_ptr<PeriodicPlanePatch<dim>> DislocationLoopLink<dim,corder,InterpolationType>::periodicPlanePatch() const
    {
        
        return (this->source->periodicPlanePatch()==this->sink->periodicPlanePatch())? this->source->periodicPlanePatch() : std::shared_ptr<PeriodicPlanePatch<dim>>(nullptr);
        
    }

    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationLoopLink<dim,corder,InterpolationType>::initFromFile(const std::string& fileName)
    {
        verboseDislocationLoopLink=TextFileParser(fileName).readScalar<int>("verboseDislocationLoopLink",true);
    }
    
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    int DislocationLoopLink<dim,corder,InterpolationType>::verboseDislocationLoopLink=0;

    
    
}
#endif
