/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNode_cpp_
#define model_DislocationNode_cpp_

#include <DislocationNode.h>



namespace model
{
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationNode<dim,corder,InterpolationType>::DislocationNode(LoopNetworkType* const net) :
    /* init */ NetworkNode<DislocationNode>(net)
    {

    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    std::shared_ptr<DislocationNode<dim,corder,InterpolationType>> DislocationNode<dim,corder,InterpolationType>::clone() const
    {
        return std::shared_ptr<DislocationNode<dim,corder,InterpolationType>>(new DislocationNode(this->p_network()));
    }

    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNode<dim,corder,InterpolationType>::initFromFile(const std::string& fileName)
    {
        use_velocityFilter=TextFileParser(fileName).readScalar<double>("use_velocityFilter",true);
        velocityReductionFactor=TextFileParser(fileName).readScalar<double>("velocityReductionFactor",true);
        assert(velocityReductionFactor>0.0 && velocityReductionFactor<=1.0);
        verboseDislocationNode=TextFileParser(fileName).readScalar<int>("verboseDislocationNode",true);
    }
    
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    int DislocationNode<dim,corder,InterpolationType>::verboseDislocationNode=0;

    template <int dim, short unsigned int corder, typename InterpolationType>
    bool DislocationNode<dim,corder,InterpolationType>::use_velocityFilter=true;

    template <int dim, short unsigned int corder, typename InterpolationType>
    double DislocationNode<dim,corder,InterpolationType>::velocityReductionFactor=0.75;
    
}
#endif
