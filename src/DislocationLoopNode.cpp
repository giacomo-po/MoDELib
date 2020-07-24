/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po         <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez   <ramirezbrf@gmail.com>.
 * Copyright (C) 2011 by Tamer Crsoby       <tamercrosby@gmail.com>,
 * Copyright (C) 2011 by Can Erel           <canerel55@gmail.com>,
 * Copyright (C) 2011 by Mamdouh Mohamed    <msm07d@fsu.edu>
 * Copyright (C) 2019 by Yash Pachaury      <ypachaur@purdue.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationLoopNode_cpp_
#define model_DislocationLoopNode_cpp_

#include <DislocationLoopNode.h>


namespace model
{
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationLoopNode<dim,corder,InterpolationType>::DislocationLoopNode(typename DislocationLoopNode<dim,corder,InterpolationType>::LoopNetworkType* const net,
                                 const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::LoopType>& loop,
                                 const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::NetworkNodeType>& networkNode) :
    /* init */ LoopNode<DislocationLoopNode>(net,loop,networkNode)
    {
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    std::shared_ptr<DislocationLoopNode<dim,corder,InterpolationType>> DislocationLoopNode<dim,corder,InterpolationType>::clone(const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::LoopType>& otherLoop,
                                                        const std::shared_ptr<typename DislocationLoopNode<dim,corder,InterpolationType>::NetworkNodeType>& otherNetworkNode) const
    {
        return std::shared_ptr<DislocationLoopNode<dim,corder,InterpolationType>>(new DislocationLoopNode(this->p_network(),otherLoop,otherNetworkNode));
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationLoopNode<dim,corder,InterpolationType>::initFromFile(const std::string& fileName)
    {
        verboseDislocationLoopNode=TextFileParser(fileName).readScalar<int>("verboseDislocationLoopNode",true);
    }


    template <int dim, short unsigned int corder, typename InterpolationType>
    int DislocationLoopNode<dim,corder,InterpolationType>::verboseDislocationLoopNode=0;
}
#endif
