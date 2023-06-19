/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui  <cuiyinan@g.ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNodeContraction_H_
#define model_DislocationNodeContraction_H_

#include <memory>
#include <tuple>
#include <Eigen/Dense>
#include <TypeTraits.h>
#include <GlidePlane.h>
//#include <BoundingMeshSegments.h>
#include <ConfinedDislocationObject.h>
#include <Grain.h>
#include <FiniteLineSegment.h>
#include <PlaneLineIntersection.h>

#ifndef NDEBUG
#define VerboseNodeContraction(N,x) if(verboseNodeContraction>=N){std::cout<<x;}
#else
#define VerboseNodeContraction(N,x)
#endif

namespace model
{
    template <typename DislocationNetworkType>
    class DislocationNodeContraction
    {
        
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;
        typedef typename TypeTraits<DislocationNetworkType>::NetworkNodeType NetworkNodeType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        
        DislocationNetworkType& DN;

    public:
        
        const int verboseNodeContraction;

        /**********************************************************************/
        DislocationNodeContraction(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        /* init */,verboseNodeContraction(TextFileParser(DN.simulationParameters.traitsIO.ddFile).readScalar<int>("verboseNodeContraction",true))
        {
            
        }
        /**********************************************************************/
        //Working version of contract boundary and contract second and virtual
//         void contractBoundary(std::shared_ptr<NetworkNodeType> nA,
//                               std::shared_ptr<NetworkNodeType> nB)
//         {
//             if (DN.simulationParameters.isPeriodicSimulation())
//             {
//                 VerboseNodeContraction(1,"Contracting boundary nodes"<<std::endl;);

//                 std::set<std::pair<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>>> bndNodesToContract;
//                 for (const auto &nALN : nA->loopNodes())
//                 {
//                     const auto nApPrev(nALN->periodicPrev());
//                     const auto nApNext(nALN->periodicNext());
//                     for (const auto &nBLN : nB->loopNodes())
//                     {
//                         const auto nBpPrev(nBLN->periodicPrev());
//                         const auto nBpNext(nBLN->periodicNext());
//                         if ((nApPrev && nApNext) && (nBpPrev && nBpNext))
//                         {
//                             if (nApPrev->networkNode==nBpPrev->networkNode)
//                             {
//                                 VerboseNodeContraction(1,"Case a"<<std::endl;);
//                                 // std::cout<<"Prev networkNOde "<<nApPrev->networkNode->tag()<<std::endl;
//                                 // std::cout<<" bnd next prev size "<<nALN->boundaryPrev().size()<<" "<<nBLN->boundaryPrev().size()<<std::endl;
//                                 //bnd prev needs to be contracted
//                                 for (const auto& bndPrevA : nALN->boundaryPrev())
//                                 {
//                                     for (const auto& bndPrevB : nBLN->boundaryPrev())
//                                     {
//                                         // std::cout<<" bndPrevA->networkNode->tag()=>bndPrevB->networkNode->tag() "<<bndPrevA->networkNode->tag()
//                                         // <<bndPrevB->networkNode->tag()<<std::endl;

//                                         // if (bndPrevA->networkNode->isMovableTo(bndPrevB->networkNode->get_P()))
//                                         if (bndPrevA->networkNode->meshFaces()==bndPrevB->networkNode->meshFaces())
//                                         {
//                                             const auto minNode(std::min(bndPrevA->networkNode->sID,bndPrevB->networkNode->sID)==bndPrevA->networkNode->sID?
//                                             bndPrevA->networkNode : bndPrevB->networkNode);
//                                             const auto maxNode(std::max(bndPrevA->networkNode->sID,bndPrevB->networkNode->sID)==bndPrevA->networkNode->sID?
//                                             bndPrevA->networkNode : bndPrevB->networkNode);
//                                             bndNodesToContract.insert(std::make_pair(minNode,maxNode));
//                                         }
//                                     }
//                                 }
//                             }
//                             if (nApNext->networkNode==nBpNext->networkNode)
//                             {
//                                 VerboseNodeContraction(1,"Case b"<<std::endl;);

//                                 //bnd next needs to be contracted
//                                 for (const auto &bndNextA : nALN->boundaryNext())
//                                 {
//                                     for (const auto &bndNextB : nBLN->boundaryNext())
//                                     {
//                                         if (bndNextA->networkNode->meshFaces()==bndNextB->networkNode->meshFaces())
//                                         {
//                                             const auto minNode (std::min(bndNextA->networkNode->sID,bndNextB->networkNode->sID)==bndNextA->networkNode->sID?
//                                             bndNextA->networkNode : bndNextB->networkNode);
//                                             const auto maxNode (std::max(bndNextA->networkNode->sID,bndNextB->networkNode->sID)==bndNextA->networkNode->sID?
//                                             bndNextA->networkNode : bndNextB->networkNode);
//                                             bndNodesToContract.insert(std::make_pair(minNode,maxNode));
//                                         }
//                                     }
//                                 }
//                             }
//                             if (nApPrev->networkNode==nBpNext->networkNode)
//                             {
//                                 VerboseNodeContraction(1,"Case c"<<std::endl;);

//                                 //bnd prev of A needs to be contracted to bndNext of B
//                                 for (const auto &bndPrevA : nALN->boundaryPrev())
//                                 {
//                                     for (const auto &bndNextB : nBLN->boundaryNext())
//                                     {
//                                         if (bndPrevA->networkNode->meshFaces()==bndNextB->networkNode->meshFaces())
//                                         {
//                                             const auto minNode(std::min(bndPrevA->networkNode->sID,bndNextB->networkNode->sID)==bndPrevA->networkNode->sID?
//                                             bndPrevA->networkNode : bndNextB->networkNode);
//                                             const auto maxNode(std::max(bndPrevA->networkNode->sID,bndNextB->networkNode->sID)==bndPrevA->networkNode->sID?
//                                             bndPrevA->networkNode : bndNextB->networkNode);
//                                             bndNodesToContract.insert(std::make_pair(minNode,maxNode));
//                                         }
//                                     }
//                                 }
//                             }
//                             if (nApNext->networkNode==nBpPrev->networkNode)
//                             {
//                                 VerboseNodeContraction(1,"Case d"<<std::endl;);

//                                 //bnd next of A needs to be contracted to bndPrev of B
//                                 for (const auto &bndNextA : nALN->boundaryNext())
//                                 {
//                                     for (const auto &bndPrevB : nBLN->boundaryPrev())
//                                     {
//                                         if (bndNextA->networkNode->meshFaces()==bndPrevB->networkNode->meshFaces())
//                                         {
//                                             const auto minNode(std::min(bndNextA->networkNode->sID,bndPrevB->networkNode->sID)==bndNextA->networkNode->sID?
//                                             bndNextA->networkNode : bndPrevB->networkNode);
//                                             const auto maxNode(std::max(bndNextA->networkNode->sID,bndPrevB->networkNode->sID)==bndNextA->networkNode->sID?
//                                             bndNextA->networkNode : bndPrevB->networkNode);
//                                             bndNodesToContract.insert(std::make_pair(minNode,maxNode));
//                                         }
//                                     }
//                                 }
//                             }
//                         }
//                     }
//                 }
//                for (const auto& pairIter : bndNodesToContract)
//                {
//                    VerboseNodeContraction(1, "Contracting "<<pairIter.first->tag()<<", "<<pairIter.second->tag() << std::endl;);

//                    DN.contractNetworkNodes(pairIter.first,pairIter.second);
//                } 
//             }

//         }
//         /**********************************************************************/
//         bool contractSecondAndVirtual(std::shared_ptr<NetworkNodeType> nA,
//                                       std::shared_ptr<NetworkNodeType> nB)
//         {

//             if ((nB->get_P()-nA->get_P()).squaredNorm()>FLT_EPSILON)
//             {
//                 // std::cout<<"Na and Nb are "<<nA->tag()<<"=>"<<nB->tag()<<std::endl;
//                 assert(nB->isMovableTo(nA->get_P()));
//                 bool movedB(nB->trySet_P(nA->get_P()));
//                 assert(movedB);
//             }
//             switch (DN.simulationParameters.simulationType)
//             {
//                 case DefectiveCrystalParameters::FINITE_FEM:
//                 {
//                     if(nB->virtualBoundaryNode())
//                     {
//                         assert(nA->virtualBoundaryNode());
//                         DN.contractNetworkNodes(nA->virtualBoundaryNode(),nB->virtualBoundaryNode());
//                     }
//                     break;
//                 }
                    
// //                case DefectiveCrystalParameters::PERIODIC:
// //                {
// //                 // FINISH HERE
// //                    break;
// //                }
//             }

//             bool tempContract(DN.contractNetworkNodes(nA,nB));  //First contract the internal nodes then the boundary (Essential for cutLoop)
//             contractBoundary(nA,nB);
//             return tempContract;
//         }
        
//Populates the boundary nodes to be contracted to be contracted later by contract second and virtual
        //Populates only if the positions of the nodes are at the same location (Same Nodal location is a necessary requirement (May need to not allow
        // contraction if all the nodes cannot be contracted but this condition give above to only contract those nodes which are at the same location
        // found successful)) ( bug : Contraction of network node that may already have been contracted)
        // std::set<std::pair<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>>> contractBoundary(std::shared_ptr<NetworkNodeType> nA,
        //                                                                                                         std::shared_ptr<NetworkNodeType> nB)
        // {
        //         std::set<std::pair<std::shared_ptr<NetworkNodeType>,std::shared_ptr<NetworkNodeType>>> bndNodesToContract;

        //     if (DN.simulationParameters.isPeriodicSimulation())
        //     {
        //         VerboseNodeContraction(1,"Contracting boundary nodes"<<std::endl;);

        //         for (const auto &nALN : nA->loopNodes())
        //         {
        //             const auto nApPrev(nALN->periodicPrev());
        //             const auto nApNext(nALN->periodicNext());
        //             for (const auto &nBLN : nB->loopNodes())
        //             {
        //                 const auto nBpPrev(nBLN->periodicPrev());
        //                 const auto nBpNext(nBLN->periodicNext());
        //                 if ((nApPrev && nApNext) && (nBpPrev && nBpNext))
        //                 {
        //                     if (nApPrev->networkNode==nBpPrev->networkNode)
        //                     {
        //                         VerboseNodeContraction(1,"Case a"<<std::endl;);
        //                         // std::cout<<"Prev networkNOde "<<nApPrev->networkNode->tag()<<std::endl;
        //                         // std::cout<<" bnd next prev size "<<nALN->boundaryPrev().size()<<" "<<nBLN->boundaryPrev().size()<<std::endl;
        //                         //bnd prev needs to be contracted
        //                         for (const auto& bndPrevA : nALN->boundaryPrev())
        //                         {
        //                             for (const auto& bndPrevB : nBLN->boundaryPrev())
        //                             {
        //                                 // std::cout<<" bndPrevA->networkNode->tag()=>bndPrevB->networkNode->tag() "<<bndPrevA->networkNode->tag()
        //                                 // <<bndPrevB->networkNode->tag()<<std::endl;

        //                                 // if (bndPrevA->networkNode->isMovableTo(bndPrevB->networkNode->get_P()))
        //                                 if (bndPrevA->networkNode->meshFaces()==bndPrevB->networkNode->meshFaces() && 
        //                                 (bndPrevA->networkNode->get_P()-bndPrevB->networkNode->get_P()).norm()<FLT_EPSILON)
        //                                 {
        //                                     const auto minNode(std::min(bndPrevA->networkNode->sID,bndPrevB->networkNode->sID)==bndPrevA->networkNode->sID?
        //                                     bndPrevA->networkNode : bndPrevB->networkNode);
        //                                     const auto maxNode(std::max(bndPrevA->networkNode->sID,bndPrevB->networkNode->sID)==bndPrevA->networkNode->sID?
        //                                     bndPrevA->networkNode : bndPrevB->networkNode);
        //                                     bndNodesToContract.insert(std::make_pair(minNode,maxNode));
        //                                 }
        //                             }
        //                         }
        //                     }
        //                     if (nApNext->networkNode==nBpNext->networkNode)
        //                     {
        //                         VerboseNodeContraction(1,"Case b"<<std::endl;);

        //                         //bnd next needs to be contracted
        //                         for (const auto &bndNextA : nALN->boundaryNext())
        //                         {
        //                             for (const auto &bndNextB : nBLN->boundaryNext())
        //                             {
        //                                 if (bndNextA->networkNode->meshFaces()==bndNextB->networkNode->meshFaces() && 
        //                                 (bndNextA->networkNode->get_P()-bndNextB->networkNode->get_P()).norm()<FLT_EPSILON)
        //                                 {
        //                                     const auto minNode (std::min(bndNextA->networkNode->sID,bndNextB->networkNode->sID)==bndNextA->networkNode->sID?
        //                                     bndNextA->networkNode : bndNextB->networkNode);
        //                                     const auto maxNode (std::max(bndNextA->networkNode->sID,bndNextB->networkNode->sID)==bndNextA->networkNode->sID?
        //                                     bndNextA->networkNode : bndNextB->networkNode);
        //                                     bndNodesToContract.insert(std::make_pair(minNode,maxNode));
        //                                 }
        //                             }
        //                         }
        //                     }
        //                     if (nApPrev->networkNode==nBpNext->networkNode)
        //                     {
        //                         VerboseNodeContraction(1,"Case c"<<std::endl;);

        //                         //bnd prev of A needs to be contracted to bndNext of B
        //                         for (const auto &bndPrevA : nALN->boundaryPrev())
        //                         {
        //                             for (const auto &bndNextB : nBLN->boundaryNext())
        //                             {
        //                                 if (bndPrevA->networkNode->meshFaces()==bndNextB->networkNode->meshFaces() && 
        //                                 (bndPrevA->networkNode->get_P()-bndNextB->networkNode->get_P()).norm()<FLT_EPSILON)
        //                                 {
        //                                     const auto minNode(std::min(bndPrevA->networkNode->sID,bndNextB->networkNode->sID)==bndPrevA->networkNode->sID?
        //                                     bndPrevA->networkNode : bndNextB->networkNode);
        //                                     const auto maxNode(std::max(bndPrevA->networkNode->sID,bndNextB->networkNode->sID)==bndPrevA->networkNode->sID?
        //                                     bndPrevA->networkNode : bndNextB->networkNode);
        //                                     bndNodesToContract.insert(std::make_pair(minNode,maxNode));
        //                                 }
        //                             }
        //                         }
        //                     }
        //                     if (nApNext->networkNode==nBpPrev->networkNode)
        //                     {
        //                         VerboseNodeContraction(1,"Case d"<<std::endl;);

        //                         //bnd next of A needs to be contracted to bndPrev of B
        //                         for (const auto &bndNextA : nALN->boundaryNext())
        //                         {
        //                             for (const auto &bndPrevB : nBLN->boundaryPrev())
        //                             {
        //                                 if (bndNextA->networkNode->meshFaces()==bndPrevB->networkNode->meshFaces() && 
        //                                 (bndNextA->networkNode->get_P()-bndPrevB->networkNode->get_P()).norm()<FLT_EPSILON)
        //                                 {
        //                                     const auto minNode(std::min(bndNextA->networkNode->sID,bndPrevB->networkNode->sID)==bndNextA->networkNode->sID?
        //                                     bndNextA->networkNode : bndPrevB->networkNode);
        //                                     const auto maxNode(std::max(bndNextA->networkNode->sID,bndPrevB->networkNode->sID)==bndNextA->networkNode->sID?
        //                                     bndNextA->networkNode : bndPrevB->networkNode);
        //                                     bndNodesToContract.insert(std::make_pair(minNode,maxNode));
        //                                 }
        //                             }
        //                         }
        //                     }
        //                 }
        //             }
        //         }
        //     //    for (const auto& pairIter : bndNodesToContract)
        //     //    {
        //     //        VerboseNodeContraction(1, "Contracting "<<pairIter.first->tag()<<", "<<pairIter.second->tag() << std::endl;);

        //     //        DN.contractNetworkNodes(pairIter.first,pairIter.second);
        //     //    } 
        //     }
        //     return bndNodesToContract;

        // }
        //Populates the boundary nodes to be contracted to be contracted later by contract second and virtual
        //Populates only if the positions of the nodes are at the same location (Same Nodal location is a necessary requirement (May need to not allow
        // contraction if all the nodes cannot be contracted but this condition give above to only contract those nodes which are at the same location
        // found successful))
        std::vector<std::set<std::shared_ptr<NetworkNodeType>>> contractBoundary(std::shared_ptr<NetworkNodeType> nA,
                                                                                 std::shared_ptr<NetworkNodeType> nB)
        {
            std::vector<std::set<std::shared_ptr<NetworkNodeType>>> bndNodesToContract;

//See if you want to implement contract network node based on sIDs .. this will greatly help with contracting younger nodes in contract boundary nodes

            if (DN.simulationParameters.isPeriodicSimulation())
            {
                VerboseNodeContraction(1,"Populating boundary nodes to be contracted"<<std::endl;);

                for (const auto &nALN : nA->loopNodes())
                {
                    const auto nApPrev(nALN->periodicPrev());
                    const auto nApNext(nALN->periodicNext());
                    for (const auto &nBLN : nB->loopNodes())
                    {
                        const auto nBpPrev(nBLN->periodicPrev());
                        const auto nBpNext(nBLN->periodicNext());
                        if ((nApPrev && nApNext) && (nBpPrev && nBpNext))
                        {
                            if (nApPrev->networkNode==nBpPrev->networkNode)
                            {
                                VerboseNodeContraction(1,"Case a"<<std::endl;);
                                // std::cout<<"Prev networkNOde "<<nApPrev->networkNode->tag()<<std::endl;
                                // std::cout<<" bnd next prev size "<<nALN->boundaryPrev().size()<<" "<<nBLN->boundaryPrev().size()<<std::endl;
                                //bnd prev needs to be contracted
                                for (const auto& bndPrevA : nALN->boundaryPrev())
                                {
                                    for (const auto& bndPrevB : nBLN->boundaryPrev())
                                    {
                                        // std::cout<<" bndPrevA->networkNode->tag()=>bndPrevB->networkNode->tag() "<<bndPrevA->networkNode->tag()
                                        // <<bndPrevB->networkNode->tag()<<std::endl;

                                        // if (bndPrevA->networkNode->isMovableTo(bndPrevB->networkNode->get_P()))
                                        if (bndPrevA->networkNode->meshFaces()==bndPrevB->networkNode->meshFaces() && 
                                        (bndPrevA->networkNode->get_P()-bndPrevB->networkNode->get_P()).norm()<FLT_EPSILON)
                                        {
                                            const auto minNode(std::min(bndPrevA->networkNode->sID,bndPrevB->networkNode->sID)==bndPrevA->networkNode->sID?
                                            bndPrevA->networkNode : bndPrevB->networkNode);
                                            const auto maxNode(std::max(bndPrevA->networkNode->sID,bndPrevB->networkNode->sID)==bndPrevA->networkNode->sID?
                                            bndPrevA->networkNode : bndPrevB->networkNode);
                                            //search in the container if the nodes are already present.. if yes place in that location
                                            bool nodeAlreadyPresentInSet(false);
                                            if (bndNodesToContract.size()==0)
                                            {
                                                std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                                nodeSet.insert(minNode);
                                                nodeSet.insert(maxNode);
                                                bndNodesToContract.push_back(nodeSet);
                                            }
                                            for (auto& setIter : bndNodesToContract)
                                            {
                                                if (setIter.find(minNode)!=setIter.end() || setIter.find(maxNode)!=setIter.end())
                                                {
                                                    setIter.insert(minNode);
                                                    setIter.insert(maxNode);
                                                    nodeAlreadyPresentInSet=true;
                                                    break;
                                                }
                                            }
                                            if (!nodeAlreadyPresentInSet)
                                            {
                                                std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                                nodeSet.insert(minNode);
                                                nodeSet.insert(maxNode);
                                                bndNodesToContract.push_back(nodeSet);
                                            }
                                        }
                                    }
                                }
                            }
                            if (nApNext->networkNode==nBpNext->networkNode)
                            {
                                VerboseNodeContraction(1,"Case b"<<std::endl;);

                                //bnd next needs to be contracted
                                for (const auto &bndNextA : nALN->boundaryNext())
                                {
                                    for (const auto &bndNextB : nBLN->boundaryNext())
                                    {
                                        if (bndNextA->networkNode->meshFaces()==bndNextB->networkNode->meshFaces() && 
                                        (bndNextA->networkNode->get_P()-bndNextB->networkNode->get_P()).norm()<FLT_EPSILON)
                                        {
                                            const auto minNode (std::min(bndNextA->networkNode->sID,bndNextB->networkNode->sID)==bndNextA->networkNode->sID?
                                            bndNextA->networkNode : bndNextB->networkNode);
                                            const auto maxNode (std::max(bndNextA->networkNode->sID,bndNextB->networkNode->sID)==bndNextA->networkNode->sID?
                                            bndNextA->networkNode : bndNextB->networkNode);
                                            //search in the container if the nodes are already present.. if yes place in that location
                                            bool nodeAlreadyPresentInSet(false);
                                            if (bndNodesToContract.size()==0)
                                            {
                                                std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                                nodeSet.insert(minNode);
                                                nodeSet.insert(maxNode);
                                                bndNodesToContract.push_back(nodeSet);
                                            }
                                            for (auto& setIter : bndNodesToContract)
                                            {
                                                if (setIter.find(minNode)!=setIter.end() || setIter.find(maxNode)!=setIter.end())
                                                {
                                                    setIter.insert(minNode);
                                                    setIter.insert(maxNode);
                                                    nodeAlreadyPresentInSet=true;
                                                    break;
                                                }
                                            }
                                            if (!nodeAlreadyPresentInSet)
                                            {
                                                std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                                nodeSet.insert(minNode);
                                                nodeSet.insert(maxNode);
                                                bndNodesToContract.push_back(nodeSet);
                                            }
                                            // bndNodesToContract.insert(std::make_pair(minNode,maxNode));
                                        }
                                    }
                                }
                            }
                            if (nApPrev->networkNode==nBpNext->networkNode)
                            {
                                VerboseNodeContraction(1,"Case c"<<std::endl;);

                                //bnd prev of A needs to be contracted to bndNext of B
                                for (const auto &bndPrevA : nALN->boundaryPrev())
                                {
                                    for (const auto &bndNextB : nBLN->boundaryNext())
                                    {
                                        if (bndPrevA->networkNode->meshFaces()==bndNextB->networkNode->meshFaces() && 
                                        (bndPrevA->networkNode->get_P()-bndNextB->networkNode->get_P()).norm()<FLT_EPSILON)
                                        {
                                            const auto minNode(std::min(bndPrevA->networkNode->sID,bndNextB->networkNode->sID)==bndPrevA->networkNode->sID?
                                            bndPrevA->networkNode : bndNextB->networkNode);
                                            const auto maxNode(std::max(bndPrevA->networkNode->sID,bndNextB->networkNode->sID)==bndPrevA->networkNode->sID?
                                            bndPrevA->networkNode : bndNextB->networkNode);
                                            //search in the container if the nodes are already present.. if yes place in that location
                                            bool nodeAlreadyPresentInSet(false);
                                            if (bndNodesToContract.size()==0)
                                            {
                                                std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                                nodeSet.insert(minNode);
                                                nodeSet.insert(maxNode);
                                                bndNodesToContract.push_back(nodeSet);
                                            }
                                            for (auto& setIter : bndNodesToContract)
                                            {
                                                if (setIter.find(minNode)!=setIter.end() || setIter.find(maxNode)!=setIter.end())
                                                {
                                                    setIter.insert(minNode);
                                                    setIter.insert(maxNode);
                                                    nodeAlreadyPresentInSet=true;
                                                }
                                            }
                                            if (!nodeAlreadyPresentInSet)
                                            {
                                                std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                                nodeSet.insert(minNode);
                                                nodeSet.insert(maxNode);
                                                bndNodesToContract.push_back(nodeSet);
                                            }
                                            // bndNodesToContract.insert(std::make_pair(minNode,maxNode));
                                        }
                                    }
                                }
                            }
                            if (nApNext->networkNode==nBpPrev->networkNode)
                            {
                                VerboseNodeContraction(1,"Case d"<<std::endl;);

                                //bnd next of A needs to be contracted to bndPrev of B
                                for (const auto &bndNextA : nALN->boundaryNext())
                                {
                                    for (const auto &bndPrevB : nBLN->boundaryPrev())
                                    {
                                        if (bndNextA->networkNode->meshFaces()==bndPrevB->networkNode->meshFaces() && 
                                        (bndNextA->networkNode->get_P()-bndPrevB->networkNode->get_P()).norm()<FLT_EPSILON)
                                        {
                                            const auto minNode(std::min(bndNextA->networkNode->sID,bndPrevB->networkNode->sID)==bndNextA->networkNode->sID?
                                            bndNextA->networkNode : bndPrevB->networkNode);
                                            const auto maxNode(std::max(bndNextA->networkNode->sID,bndPrevB->networkNode->sID)==bndNextA->networkNode->sID?
                                            bndNextA->networkNode : bndPrevB->networkNode);
                                            //search in the container if the nodes are already present.. if yes place in that location
                                            bool nodeAlreadyPresentInSet(false);
                                            if (bndNodesToContract.size()==0)
                                            {
                                                std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                                nodeSet.insert(minNode);
                                                nodeSet.insert(maxNode);
                                                bndNodesToContract.push_back(nodeSet);
                                            }
                                            for (auto& setIter : bndNodesToContract)
                                            {
                                                if (setIter.find(minNode)!=setIter.end() || setIter.find(maxNode)!=setIter.end())
                                                {
                                                    setIter.insert(minNode);
                                                    setIter.insert(maxNode);
                                                    nodeAlreadyPresentInSet=true;
                                                }
                                            }
                                            if (!nodeAlreadyPresentInSet)
                                            {
                                                std::set<std::shared_ptr<NetworkNodeType>> nodeSet;
                                                nodeSet.insert(minNode);
                                                nodeSet.insert(maxNode);
                                                bndNodesToContract.push_back(nodeSet);
                                            }
                                            // bndNodesToContract.insert(std::make_pair(minNode,maxNode));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            //    for (const auto& vecIter : bndNodesToContract)
            //    {
            //        std::cout<<"Contracting  : "<<std::flush<<std::endl;
            //         for (const auto& setIter : vecIter)
            //         {
            //             std::cout<<setIter->sID<<", ";
            //         }
            //         std::cout<<std::endl;
            //     //    DN.contractNetworkNodes(pairIter.first,pairIter.second);
            //    } 
            }
            return bndNodesToContract;

        }
        /**********************************************************************/
        bool contractSecondAndVirtual(std::shared_ptr<NetworkNodeType> nA,
                                      std::shared_ptr<NetworkNodeType> nB)
        {

            if ((nB->get_P()-nA->get_P()).squaredNorm()>FLT_EPSILON)
            {
                // std::cout<<"Na and Nb are "<<nA->tag()<<"=>"<<nB->tag()<<std::endl;
                assert(nB->isMovableTo(nA->get_P()));
                bool movedB(nB->trySet_P(nA->get_P()));
                if(!movedB)
                {
                    throw std::runtime_error("Node B could not be moved to position.");
                }
//                assert(movedB);
            }
//            switch (DN.simulationParameters.simulationType)
//            {
//                case DefectiveCrystalParameters::FINITE_FEM:
//                {
//                    if(nB->virtualBoundaryNode())
//                    {
//                        assert(nA->virtualBoundaryNode());
//                        DN.contractNetworkNodes(nA->virtualBoundaryNode(),nB->virtualBoundaryNode());
//                    }
//                    break;
//                }
//                    
////                case DefectiveCrystalParameters::PERIODIC:
////                {
////                 // FINISH HERE
////                    break;
////                }
//            }

            auto bndNodesToContract(contractBoundary(nA,nB));

            bool tempContract(DN.contractNetworkNodes(nA,nB));  //First contract the internal nodes then the boundary (Essential for cutLoop)
            //Now contract the boundary nodes
            for (const auto& vecIter : bndNodesToContract)
            {
                for (const auto &setIter : vecIter)
                {
                    if (setIter!= *vecIter.begin())
                    {
                       VerboseNodeContraction(1, "Contracting Boundary"<<(*vecIter.begin())->sID<<", "<<setIter->sID << std::endl;);
                       DN.contractNetworkNodes(*vecIter.begin(),setIter);
                    }
                }
            }
            return tempContract;
        }
        /**********************************************************************/
        bool contractYoungest(std::shared_ptr<NetworkNodeType> nA,
                              std::shared_ptr<NetworkNodeType> nB)
        {
            return nA->sID<nB->sID? contractSecondAndVirtual(nA,nB) : contractSecondAndVirtual(nB,nA);
            
        }
        
        /**********************************************************************/
        bool contractToPosition(std::shared_ptr<NetworkNodeType> nA,
                                std::shared_ptr<NetworkNodeType> nB,
                                const VectorDim& X,
                                const double& maxRange)
        {
            
            bool movedA=false;
            bool movedB=false;
            
            if(   nA->isMovableTo(X)
               && nB->isMovableTo(X)
               && (nA->get_P()-X).norm()+(nB->get_P()-X).norm()<maxRange)
            {
                // std::cout<<" Before motion positions are "<<nA->get_P().transpose()<<"=>"<<nB->get_P().transpose()<<std::endl;
                movedA=nA->trySet_P(X);
                movedB=nB->trySet_P(X);
                // std::cout<<" After motion positions are "<<nA->get_P().transpose()<<"=>"<<nB->get_P().transpose()<<std::endl;
                // std::cout<<" X is "<<X.transpose()<<std::endl;
                // movedA=nA->set_P(X);
                // movedB=nB->set_P(X);
                VerboseNodeContraction(2,"contractToPosition"<<std::endl;);
                VerboseNodeContraction(2,"movedA="<<movedA<<std::endl;);
                VerboseNodeContraction(2,"movedB="<<movedB<<std::endl;);
                assert(movedA && movedB && "COULD NOT MOVE NODES");
            }
            
            return (movedA && movedB)? contractYoungest(nA,nB) : false;
        }
        
        

//Contract function for not allowing the node contraction at the boundary or outside the boundary (Nodes should not be placed at the boundary)
        bool contract(std::shared_ptr<NetworkNodeType> nA,
                      std::shared_ptr<NetworkNodeType> nB)
        {
            if (!(DN.simulationParameters.isPeriodicSimulation() && (nA->isBoundaryNode() || nB->isBoundaryNode()))) //Added for avoiding node contraction with boundary node
            {
                VerboseNodeContraction(1, "DislocationNodeContraction::contract " << nA->sID << " " << nB->sID << std::endl;);

                const bool nAisMovable = nA->isMovableTo(nB->get_P());
                const bool nBisMovable = nB->isMovableTo(nA->get_P());

                if (nAisMovable && nBisMovable)
                {
                    VerboseNodeContraction(1, "DislocationNodeContraction case 1a" << std::endl;);
                    return contractYoungest(nA, nB);
                }
                else if (nAisMovable && !nBisMovable)
                {
                    VerboseNodeContraction(1, "DislocationNodeContraction case 1b" << std::endl;);
                    return contractSecondAndVirtual(nB, nA);
                }
                else if (!nAisMovable && nBisMovable)
                {
                    VerboseNodeContraction(1, "DislocationNodeContraction case 1c" << std::endl;);
                    return contractSecondAndVirtual(nA, nB);
                }
                else
                { // nA and nB cannot be moved to each other. The calculation of a third point is necessary
                    const ConfinedDislocationObject<dim> cdoA(*nA);
                    const ConfinedDislocationObject<dim> cdoB(*nB);

                    VerboseNodeContraction(1, "DislocationNodeContraction case 1d" << std::endl;);

                    const double maxRange = 4.0 * (nA->get_P() - nB->get_P()).norm();

                    if (nA->isOnBoundary() || nB->isOnBoundary())
                    { // either one of the nodes is a boundary node. Therefore the contraction point must be a boundary node

                        //BoundingMeshSegments<dim> temp(nA->boundingBoxSegments(),nB->boundingBoxSegments());
                        const ConfinedDislocationObject<dim> cdo(*nA, *nB);
                        const BoundingMeshSegments<dim> &temp(cdo.boundingBoxSegments());

                        VerboseNodeContraction(1, "temp.size=" << temp.size() << std::endl;);

                        switch (temp.size())
                        {

                        case 0:
                        { // no intersection of the bounding boxes
                            VerboseNodeContraction(1, "DislocationNodeContraction case 1f" << std::endl;);
                            return false;
                            break;
                        }

                        case 1:
                        {
                            const std::shared_ptr<MeshBoundarySegment<dim>> &mbs(temp.front());
                            if ((mbs->P0 - mbs->P1).norm() < FLT_EPSILON)
                            { // a unique intersection point of the bounding boxes exist
                                VerboseNodeContraction(1, "DislocationNodeContraction case 1d" << std::endl;);
                                const VectorDim X = 0.5 * (mbs->P0 + mbs->P1);
                                return contractToPosition(nA, nB, X, maxRange);
                            }
                            else
                            { // two possible intersection points of the bounding boxes exist
                                const bool firstIsCloser = (nA->get_P() - mbs->P0).norm() + (nB->get_P() - mbs->P0).norm() < (nA->get_P() - mbs->P1).norm() + (nB->get_P() - mbs->P1).norm();
                                const VectorDim X = firstIsCloser ? mbs->P0 : mbs->P1;
                                const VectorDim Y = firstIsCloser ? mbs->P1 : mbs->P0;

                                VerboseNodeContraction(1, "DislocationNodeContraction case 1dX" << std::endl;);
                                const bool Xcontracted = contractToPosition(nA, nB, X, maxRange);
                                if (Xcontracted)
                                {
                                    return true;
                                }
                                else
                                {
                                    VerboseNodeContraction(1, "DislocationNodeContraction case 1dY" << std::endl;);
                                    const bool Ycontracted = contractToPosition(nA, nB, Y, maxRange);
                                    if (Ycontracted)
                                    {
                                        return true;
                                    }
                                    else
                                    {
                                        return false;
                                    }
                                }
                            }
                            break;
                        }

                        default:
                        { // bounding boxes intersect in more than one line

                            std::map<double, VectorDim> vertexMap;
                            for (const auto &seg : temp)
                            {
                                const double firstRange = (nA->get_P() - seg->P0).norm() + (nB->get_P() - seg->P0).norm();
                                if (firstRange < maxRange)
                                {
                                    vertexMap.insert(std::make_pair(firstRange, seg->P0));
                                }

                                const double secondRange = (nA->get_P() - seg->P1).norm() + (nB->get_P() - seg->P1).norm();
                                if (secondRange < maxRange)
                                {
                                    vertexMap.insert(std::make_pair(secondRange, seg->P1));
                                }
                            }

                            bool success = false;
                            for (const auto &pair : vertexMap)
                            {
                                success = contractToPosition(nA, nB, pair.second, maxRange);
                                if (success)
                                {
                                    break;
                                }
                            }
                            return success;

                            break;
                        }
                        }
                    }
                    else
                    { // neither nA nor nB are on bounding box
                        VerboseNodeContraction(1, "DislocationNodeContraction case 5" << std::endl;);
                        if (cdoA.glidePlaneIntersections() && cdoB.glidePlaneIntersections())
                        { // both nodes confined by more then one plane
                            SegmentSegmentDistance<dim> ssd(cdoA.glidePlaneIntersections()->P0, cdoA.glidePlaneIntersections()->P1,
                                                            cdoB.glidePlaneIntersections()->P0, cdoB.glidePlaneIntersections()->P1);

                            const auto iSeg = ssd.intersectionSegment();
                            if (iSeg.size() == 1)
                            { // incident intersection
                                VerboseNodeContraction(1, "DislocationNodeContraction case 5a" << std::endl;);
                                if (DN.simulationParameters.isPeriodicSimulation())
                                {
                                    if (((cdoA.glidePlaneIntersections()->snap(std::get<0>(iSeg[0]))-cdoA.glidePlaneIntersections()->snapToInfiniteLine(std::get<0>(iSeg[0]))).squaredNorm()<FLT_EPSILON)
                                && ((cdoB.glidePlaneIntersections()->snap(std::get<0>(iSeg[0]))-cdoB.glidePlaneIntersections()->snapToInfiniteLine(std::get<0>(iSeg[0]))).squaredNorm()<FLT_EPSILON))
                                    {
                                        if ((std::get<1>(iSeg[0])>FLT_EPSILON && std::get<1>(iSeg[0])<(1.0-FLT_EPSILON)) &&
                                        (std::get<2>(iSeg[0])>FLT_EPSILON && std::get<2>(iSeg[0])<(1.0-FLT_EPSILON)) )
                                        {
                                            return contractToPosition(nA, nB, std::get<0>(iSeg[0]), maxRange);
                                        }
                                        else
                                        {
                                            return false;
                                        }
                                        
                                    }
                                    else
                                    {
                                        return false;
                                    }
                                }
                                else
                                {
                                    return contractToPosition(nA, nB, std::get<0>(iSeg[0]), maxRange);
                                }
                                
                            }
                            else if (iSeg.size() == 2)
                            { // coincident intersection
                                VerboseNodeContraction(1, "DislocationNodeContraction case 5b" << std::endl;);
                                FiniteLineSegment<dim> ls(std::get<0>(iSeg[0]), std::get<0>(iSeg[1]));
                                
                                if (DN.simulationParameters.isPeriodicSimulation())
                                {
                                    if ((ls.snap(0.5 * (nA->get_P() + nB->get_P())) - ls.snapToInfiniteLine(0.5 * (nA->get_P() + nB->get_P()))).squaredNorm() < FLT_EPSILON)
                                    {

                                        if ((std::get<1>(iSeg[0]) > FLT_EPSILON && std::get<1>(iSeg[0]) < (1.0 - FLT_EPSILON)) &&
                                            (std::get<2>(iSeg[0]) > FLT_EPSILON && std::get<2>(iSeg[0]) < (1.0 - FLT_EPSILON)) &&
                                            (std::get<1>(iSeg[1]) > FLT_EPSILON && std::get<1>(iSeg[1]) < (1.0 - FLT_EPSILON)) &&
                                            (std::get<2>(iSeg[1]) > FLT_EPSILON && std::get<2>(iSeg[1]) < (1.0 - FLT_EPSILON))) //May not be necessary
                                        {
                                            return contractToPosition(nA, nB, ls.snap(0.5 * (nA->get_P() + nB->get_P())), maxRange);
                                        }
                                        else
                                        {
                                            return false;
                                        }
                                        
                                    }
                                    else
                                    {
                                        return false;
                                    }
                                }
                                else
                                {
                                     return contractToPosition(nA, nB, ls.snap(0.5 * (nA->get_P() + nB->get_P())), maxRange);
                                }
                            }
                            else
                            { // parallel or skew intersection
                                VerboseNodeContraction(1, "DislocationNodeContraction case 5c" << std::endl;);
                                return false;
                            }
                        }
                        else if (cdoA.glidePlaneIntersections() && !cdoB.glidePlaneIntersections())
                        { // nA confined by more then one plane, nB confined by zero or one plane
                            if (nB->glidePlanes().size() == 0)
                            { // nB confined by zero planes
                                VerboseNodeContraction(1, "DislocationNodeContraction case 6AA" << std::endl;);
                                if (cdoA.glidePlaneIntersections()->contains(nB->get_P()))
                                {
                                    VerboseNodeContraction(1, "DislocationNodeContraction case 6AA1" << std::endl;);
                                    return contractToPosition(nA, nB, nB->get_P(), maxRange);
                                }
                                else
                                {
                                    VerboseNodeContraction(1, "DislocationNodeContraction case 6AA2" << std::endl;);
                                    return false;
                                }
                            }
                            else if (nB->glidePlanes().size() == 1)
                            { // nB confined by one plane
                                PlaneLineIntersection<dim> pli((*nB->glidePlanes().begin())->P,
                                                               (*nB->glidePlanes().begin())->unitNormal,
                                                               cdoA.glidePlaneIntersections()->P0,                                    // origin of line
                                                               cdoA.glidePlaneIntersections()->P1 - cdoA.glidePlaneIntersections()->P0 // line direction
                                );

                                // THERE SHOULD BE A PlaneSegmentIntersection class, which intersects the plane with a finite segment

                                if (pli.type == PlaneLineIntersection<dim>::COINCIDENT)
                                { // nothing to do, _glidePlaneIntersections remains unchanged
                                    VerboseNodeContraction(1, "DislocationNodeContraction case 6a" << std::endl;);
                                    if (DN.simulationParameters.isPeriodicSimulation())
                                    {
                                                //There should be no need to check for the intersection here for internal nodes being placed at the boundary

                                        if ((cdoA.glidePlaneIntersections()->snap(0.5 * (nA->get_P() + nB->get_P()))-cdoA.glidePlaneIntersections()->snapToInfiniteLine(0.5 * (nA->get_P() + nB->get_P()))).squaredNorm()<FLT_EPSILON)
                                        {
                                            return contractToPosition(nA, nB, cdoA.glidePlaneIntersections()->snap(0.5 * (nA->get_P() + nB->get_P())), maxRange);
                                        }
                                        else
                                        {
                                            return false;
                                        }
                                    }
                                    else
                                    {
                                        return contractToPosition(nA, nB, cdoA.glidePlaneIntersections()->snap(0.5 * (nA->get_P() + nB->get_P())), maxRange);
                                    }
                                }
                                else if (pli.type == PlaneLineIntersection<dim>::INCIDENT)
                                { // _glidePlaneIntersections becomes a singular point
                                    VerboseNodeContraction(1, "DislocationNodeContraction case 6b" << std::endl;);
                                    FiniteLineSegment<dim> cutLine(cdoA.glidePlaneIntersections()->P0, cdoA.glidePlaneIntersections()->P1);
                                    if ((pli.P - cutLine.snap(pli.P)).squaredNorm() < FLT_EPSILON)
                                    { // intersection point is inside mesh
                                        VerboseNodeContraction(1, "DislocationNodeContraction case 6b1" << std::endl;);
                                        if (DN.simulationParameters.isPeriodicSimulation())
                                        {
                                            if ((cutLine.snap(pli.P)-cutLine.snapToInfiniteLine(pli.P)).squaredNorm()<FLT_EPSILON)
                                            {
                                                if ((cutLine.P0-pli.P).norm()>FLT_EPSILON && (cutLine.P1-pli.P).norm()>FLT_EPSILON) //Indicates position at the boundary
                                                {
                                                    return contractToPosition(nA, nB, pli.P, maxRange);
                                                }
                                                else
                                                {
                                                    return false;
                                                } 
                                            }
                                            else
                                            {
                                                return false;
                                            }
                                        }
                                        else
                                        {
                                            return contractToPosition(nA, nB, pli.P, maxRange);
                                        }
                                    }
                                    else
                                    {
                                        return false;
                                    }
                                }
                                else
                                { // parallel planes, cannot contract
                                    VerboseNodeContraction(1, "DislocationNodeContraction case 6c" << std::endl;);
                                    return false;
                                }
                            }
                            else
                            {
                                assert(false && " no cdoB.glidePlaneIntersections() with more than one glide planes");
                            }
                        }
                        else if (!cdoA.glidePlaneIntersections() && cdoB.glidePlaneIntersections())
                        { // same as previous case, so switch nA and nB
                            VerboseNodeContraction(1, "DislocationNodeContraction case 7a" << std::endl;);
                            return contract(nB, nA);
                        }
                        else
                        { // both nodes confined by only one plane

                            if (nA->glidePlanes().size() == 1 && nB->glidePlanes().size() == 1)
                            {
                                const PlanePlaneIntersection<dim> &ppi(DN.glidePlaneFactory.glidePlaneIntersection(*nA->glidePlanes().begin(), *nB->glidePlanes().begin()));

                                if (ppi.type == PlanePlaneIntersection<dim>::COINCIDENT)
                                { // the contraction point can be the averago of nA and nB, which should be internal for convex domains
                                    VerboseNodeContraction(1, "DislocationNodeContraction case 8a" << std::endl;);
                                    return contractToPosition(nA, nB, (*nA->glidePlanes().begin())->snapToPlane(0.5 * (nA->get_P() + nB->get_P())), maxRange);
                                    //                                return contractToPosition(nA,nB,cdoA.glidePlaneIntersections()->snap(0.5*(nA->get_P()+nB->get_P())),maxRange);
                                }
                                else if (ppi.type == PlanePlaneIntersection<dim>::INCIDENT)
                                {
                                    VerboseNodeContraction(1, "DislocationNodeContraction case 8b" << std::endl;);
                                    const ConfinedDislocationObject<dim> cdo(*nA, *nB);
                                    const BoundingMeshSegments<dim> &temp(cdo.boundingBoxSegments());
                                    switch (temp.size())
                                    {
                                        case 2:
                                        {
                                            const std::shared_ptr<MeshBoundarySegment<dim>> &seg0(temp.front());
                                            const std::shared_ptr<MeshBoundarySegment<dim>> &seg1(temp.back());

                                            FiniteLineSegment<dim> cutLine(0.5 * (seg0->P0 + seg0->P1), 0.5 * (seg1->P0 + seg1->P1));
                                            if (DN.simulationParameters.isPeriodicSimulation())
                                            {
                                                //There should be no need to check for the intersection here for internal nodes being placed at the boundary

                                                if ((cutLine.snap(0.5 * (nA->get_P() + nB->get_P())) - cutLine.snapToInfiniteLine(0.5 * (nA->get_P() + nB->get_P()))).squaredNorm() < FLT_EPSILON)
                                                {
                                                    return contractToPosition(nA, nB, cutLine.snap(0.5 * (nA->get_P() + nB->get_P())), maxRange);
                                                }
                                                else
                                                {
                                                    return false;
                                                }
                                            }
                                            else
                                            {
                                                return contractToPosition(nA, nB, cutLine.snap(0.5 * (nA->get_P() + nB->get_P())), maxRange);
                                            }

                                            break;
                                        }
                                        case 1:
                                        {
                                            const std::shared_ptr<MeshBoundarySegment<dim>> &seg0(temp.front());

                                            FiniteLineSegment<dim> cutLine(seg0->P0, seg0->P1);
                                            if (DN.simulationParameters.isPeriodicSimulation())
                                            {
                                                //There should be no need to check for the intersection here for internal nodes being placed at the boundary

                                                if ((cutLine.snap(0.5 * (nA->get_P() + nB->get_P())) - cutLine.snapToInfiniteLine(0.5 * (nA->get_P() + nB->get_P()))).squaredNorm() < FLT_EPSILON)
                                                {
                                                    return contractToPosition(nA, nB, cutLine.snap(0.5 * (nA->get_P() + nB->get_P())), maxRange);
                                                }
                                                else
                                                {
                                                    return false;
                                                }
                                            }
                                            else
                                            {
                                                return contractToPosition(nA, nB, cutLine.snap(0.5 * (nA->get_P() + nB->get_P())), maxRange);
                                            }
                                            break;
                                        }
                                        default:
                                        { // intersection line outside mesh
                                            return false;
                                        }
                                    }
                                }
                                else
                                {
                                    VerboseNodeContraction(1, "DislocationNodeContraction case 8c" << std::endl;);
                                    return false;
                                }
                            }
                            else if (nA->glidePlanes().size() == 0 && nB->glidePlanes().size() == 1)
                            { // nA has no GlidePlane
                                const VectorDim dir(nA->invariantDirectionOfMotion());
                                if (dir.norm() > FLT_EPSILON)
                                { // nB can be moved along dir
                                    PlaneLineIntersection<dim> pli((*nB->glidePlanes().begin())->P,
                                                                   (*nB->glidePlanes().begin())->unitNormal,
                                                                   nA->get_P(), // origin of line
                                                                   dir          // line direction
                                    );

                                    if (pli.type == PlaneLineIntersection<dim>::COINCIDENT)
                                    { // nothing to do, _glidePlaneIntersections remains unchanged
                                        VerboseNodeContraction(1, "DislocationNodeContraction case 9a" << std::endl;);
                                        return contractToPosition(nA, nB, nA->get_P(), maxRange);
                                    }
                                    else if (pli.type == PlaneLineIntersection<dim>::INCIDENT)
                                    { // _glidePlaneIntersections becomes a singular point
                                        VerboseNodeContraction(1, "DislocationNodeContraction case 9b" << std::endl;);
                                        // This is a wrong case
                                        FiniteLineSegment<dim> cutLine(cdoA.glidePlaneIntersections()->P0, cdoA.glidePlaneIntersections()->P1);
                                        return contractToPosition(nA, nB, pli.P, maxRange);
                                        //
                                        //                                    if((pli.P-cutLine.snap(pli.P)).squaredNorm()<FLT_EPSILON)
                                        //                                    {// intersection point is inside mesh
                                        //                                        VerboseNodeContraction(1,"DislocationNodeContraction case96b1"<<std::endl;);
                                        //                                    }
                                        //                                    else
                                        //                                    {
                                        //                                        return false;
                                        //                                    }
                                    }
                                    else
                                    { // parallel planes, cannot contract
                                        VerboseNodeContraction(1, "DislocationNodeContraction case 9c" << std::endl;);
                                        return false;
                                    }
                                }
                                else
                                {
                                    return false;
                                }
                            }
                            else if (nA->glidePlanes().size() == 0 && nB->glidePlanes().size() == 0)
                            { // nB has no GlidePlane
                                return false;
                            }
                            else
                            { // neither nA nor nB have a GlidePlane
                                assert(nA->glidePlanes().size() <= 1);
                                assert(nB->glidePlanes().size() <= 1);
                                return false;
                            }
                        }
                    }
                }
            }

            else
            {
                VerboseNodeContraction(1, "DislocationNodeContraction boundary node contraction case " << std::endl;);

                if ((nA->get_P() - nB->get_P()).squaredNorm() < FLT_EPSILON)
                {
                    return contractYoungest(nA, nB);
                }
                else
                {
                    return false;
                }
            }
        }
    };
    
    // Static data
//    template <typename DislocationNetworkType>
//    int DislocationNodeContraction<DislocationNetworkType>::verboseNodeContraction=0;
    
}
#endif

