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

#ifndef model_DislocationLoopNode_h_
#define model_DislocationLoopNode_h_

#include <TypeTraits.h>
#include <LoopNode.h>

#ifndef NDEBUG
#define VerboseDislocationLoopNode(N,x) if(verboseDislocationLoopNode>=N){std::cout<<x;}
#else
#define VerboseDislocationLoopNode(N,x)
#endif

namespace model
{
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    class DislocationLoopNode : public LoopNode<DislocationLoopNode<dim,corder,InterpolationType>>
    {
        
        
        
        public:
        
        typedef TypeTraits<DislocationLoop<dim,corder,InterpolationType>> TraitsType;
        typedef typename TraitsType::LoopNetworkType LoopNetworkType;
        typedef typename TraitsType::LoopType LoopType;
        typedef typename TraitsType::LoopNodeType LoopNodeType;
        typedef typename TraitsType::LoopLinkType LoopLinkType;
        typedef typename TraitsType::NetworkNodeType NetworkNodeType;
        typedef typename TraitsType::NetworkLinkType NetworkLinkType;
        typedef typename TraitsType::FlowType FlowType;
        
        
        static int verboseDislocationLoopNode;

        
        DislocationLoopNode(LoopNetworkType* const,
                      const std::shared_ptr<LoopType>&,
                      const std::shared_ptr<NetworkNodeType>&);
        
        std::shared_ptr<DislocationLoopNode> clone(const std::shared_ptr<LoopType>&,
                                             const std::shared_ptr<NetworkNodeType>&) const;
        
        static void initFromFile(const std::string&);

    };
    

    
}
#endif
