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
    
    template <int dim, short unsigned int corder>
    DislocationLoopLink<dim,corder>::DislocationLoopLink(const std::shared_ptr<LoopNodeType>& so,
                                             const std::shared_ptr<LoopNodeType>& si,
                                             const std::shared_ptr<LoopType>& pL) :
    /* init */ LoopLink<LoopLinkType>(so,si,pL)
    {
        
        
    }
    
    template <int dim, short unsigned int corder>
    bool DislocationLoopLink<dim,corder>::hasNetworkLink() const
    {
        // return !((this->source->get_P()-this->sink->get_P()).squaredNorm()<FLT_EPSILON 
        // /*    */ && this->source->periodicPlaneEdge
        // /*    */ &&   this->sink->periodicPlaneEdge
        // /*    */ && this->source->loop()->sID==this->sink->loop()->sID);
        if (this->source->periodicPlaneEdge && this->sink->periodicPlaneEdge)
        {
            // bool temp((this->source->loop()->sID==this->sink->loop()->sID));
            // if (temp)
            // {
                // std::cout<<"Checking Twin"<<std::endl;
                // std::cout<<this->source->periodicPlaneEdge->twin->tag()<<"=>= " <<this->sink->periodicPlaneEdge->tag()<<std::endl;
                // std::cout<<this->sink->periodicPlaneEdge->twin->tag()<<"=>= " <<this->source->periodicPlaneEdge->tag()<<std::endl;
                if (this->source->periodicPlaneEdge->twin && this->sink->periodicPlaneEdge->twin)
                {
                    // temp *= ((this->source->periodicPlaneEdge->twin == this->sink->periodicPlaneEdge.get()) &&
                    //          (this->sink->periodicPlaneEdge->twin == this->source->periodicPlaneEdge.get()));
                    return !((this->source->periodicPlaneEdge->twin == this->sink->periodicPlaneEdge.get()) &&
                             (this->sink->periodicPlaneEdge->twin == this->source->periodicPlaneEdge.get()));

                }
                else
                {
                    // temp=true;
                    return true;
                }
            // }
            // std::cout<<"Return temo "<<temp<<std::endl;
            // std::cout<<" This->tag() "<<this->tag()<<" has networklink "<<temp<<std::endl;
            // return !temp;
        }
        else
        {
            return true;
        }
        
    }
    
    template <int dim, short unsigned int corder>
    std::shared_ptr<PeriodicPlanePatch<dim>> DislocationLoopLink<dim,corder>::periodicPlanePatch() const
    {
        
        return (this->source->periodicPlanePatch()==this->sink->periodicPlanePatch())? this->source->periodicPlanePatch() : std::shared_ptr<PeriodicPlanePatch<dim>>(nullptr);
        
    }

    
    template <int dim, short unsigned int corder>
    void DislocationLoopLink<dim,corder>::initFromFile(const std::string& fileName)
    {
        verboseDislocationLoopLink=TextFileParser(fileName).readScalar<int>("verboseDislocationLoopLink",true);
    }
    
    template <int dim, short unsigned int corder>
    const typename DislocationLoopLink<dim,corder>::LoopLinkType * DislocationLoopLink<dim,corder>::twinnedLink() const
    {
        if (hasNetworkLink())
        {
            if (this->networkLink()->loopLinks().size()==2 && this->networkLink()->hasZeroBurgers())
            {
                if (*this->networkLink()->loopLinks().begin()==this)
                {
                    return (*this->networkLink()->loopLinks().rbegin());
                }
                else
                {
                    return (*this->networkLink()->loopLinks().begin());
                }
                
            }
            else
            {
                return nullptr;
            }
            
        }
        else
        {
            return nullptr;
        }
        
    }

    template <int dim, short unsigned int corder>
    int DislocationLoopLink<dim,corder>::verboseDislocationLoopLink=0;

    template class DislocationLoopLink<3,0>;
    
}
#endif
