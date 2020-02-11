/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PolygonClipper_H_
#define model_PolygonClipper_H_
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <math.h>       /* fmod */

#include <memory>
#include <vector>
#include <Eigen/Dense>
#include <StaticID.h>
#include <CompareVectorsByComponent.h>
#include <Polygon2D.h>
#include <SegmentSegmentDistance.h>
#include <stdio.h>



namespace model
{

    struct LoopPathClipperNode;

    /******************************************************************************/
    struct ClipperEdge
    {
        const std::shared_ptr<LoopPathClipperNode> source;
        const std::shared_ptr<LoopPathClipperNode> sink;
        const int wn;
        
        ClipperEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
                const std::shared_ptr<LoopPathClipperNode>& sink_in,
                const int& wn_in);
        std::string tag() const;
        const ClipperEdge* next(std::set<const ClipperEdge*>&) const;
        const ClipperEdge* twin() const;
        bool isLoop() const;
        bool isPatch() const;
        bool isReversePatch() const;
        ClipperEdge(const ClipperEdge& other) =delete;
        ClipperEdge(ClipperEdge&& other) =delete;
        ClipperEdge& operator=(const ClipperEdge& other)=delete;
        ClipperEdge& operator=(ClipperEdge&& other)=delete;
    };

    /******************************************************************************/
    struct LoopEdge : public ClipperEdge
    {
        LoopEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
                const std::shared_ptr<LoopPathClipperNode>& sink_in,
                const int& wn_in);
    };

    struct PatchEdge : public ClipperEdge
    {
        PatchEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
                const std::shared_ptr<LoopPathClipperNode>& sink_in,
                const int& wn_in);
    };

    //Added by Yash
    struct ReversePatchEdge : public ClipperEdge
    {
        ReversePatchEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
                        const std::shared_ptr<LoopPathClipperNode>& sink_in,
                        const int& wn_in);
    };


    /******************************************************************************/
    struct LoopPathClipperNode : public Eigen::Matrix<double,2,1>,
                                public StaticID<LoopPathClipperNode>
    {
        const LoopEdge* loopPrev;
        const LoopEdge* loopNext;
        const PatchEdge* patchPrev;
        const PatchEdge* patchNext;
        const ReversePatchEdge *reversepatchPrev;
        const ReversePatchEdge *reversepatchNext;

        
        LoopPathClipperNode(const Eigen::Matrix<double,2,1>& pos) :
        /* init */ Eigen::Matrix<double,2,1>(pos)
        /* init */,loopPrev(nullptr)
        /* init */,loopNext(nullptr)
        /* init */,patchPrev(nullptr)
        /* init */,patchNext(nullptr)
        /* init */,reversepatchPrev(nullptr)  //Added by Yash
        /* init */,reversepatchNext(nullptr)  //Added by Yash
        {
            
        }
        
        LoopPathClipperNode(const LoopPathClipperNode& other) =delete;
        LoopPathClipperNode(LoopPathClipperNode&& other) =delete;
        LoopPathClipperNode& operator=(const LoopPathClipperNode& other)=delete;
        LoopPathClipperNode& operator=(LoopPathClipperNode&& other)=delete;
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const LoopPathClipperNode& node)
        {
            os<<node.sID<<" "<<node.transpose();
            return os;
        }
        
    };

    /******************************************************************************/
    ClipperEdge::ClipperEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
                    const std::shared_ptr<LoopPathClipperNode>& sink_in,
                    const int& wn_in) :
    /* init */ source(source_in)
    /* init */,sink(sink_in)
    /* init */,wn(wn_in)
    {
    }

    bool ClipperEdge::isLoop() const
    {
        return source->loopNext==this && sink->loopPrev==this;
    }

    bool ClipperEdge::isPatch() const
    {
        // std::cout<<" Source--->Sink "<<source->sID<<"--->"<<sink->sID<<std::endl;
        // std::cout<<" SourcePN--->SinkPN "<<source->patchNext->source->sID<<"("<<source->patchNext->sink->sID<<")"<<std::endl;

        return source->patchNext==this && sink->patchPrev==this;
    }

    //Added by Yash
    bool ClipperEdge::isReversePatch() const
    {
        // std::cout << " Source--->Sink " << source->sID << "--->" << sink->sID << std::endl;
        // std::cout<<" SourceRPN--->SinkRPN "<<source->reversepatchNext->source->sID<<"("<<source->reversepatchNext->sink->sID<<")"<<std::endl;

        return source->reversepatchNext==this && sink->reversepatchPrev==this;
    }

    const ClipperEdge* ClipperEdge::next(std::set<const ClipperEdge*>& validStarts) const
    {
        const bool isLoopEdg(isLoop());
        const bool isPatchEdg(isPatch());
        const bool isReversePatchEdg(isReversePatch());
        assert(( (!isLoopEdg && (isPatchEdg!=isReversePatchEdg)) || (isLoopEdg && (!isPatchEdg && !isReversePatchEdg)))  && "edge must be either loop or patch or reverse patch");
        if (isLoopEdg)
        {
            //Currently on loop
            assert(sink->loopNext && "loopNext must exist");
            if (sink->patchNext)
            {
                if ((sink->patchNext->wn > 0 && sink->loopNext->wn > 0) ||
                    (sink->patchNext->wn < 0 && sink->loopNext->wn < 0))
                {
                    //Bifurcation exists add into validstart
                    if (!sink->patchNext->twin())
                    {
                        validStarts.insert(sink->patchNext);
                    }

                    if (abs(sink->patchNext->wn) > abs(sink->loopNext->wn))
                    {
                        if (!sink->patchNext->twin())
                        {
                            return sink->patchNext;
                        }
                        else
                        {
                            if (validStarts.find(sink->patchNext->twin()) != validStarts.end())
                            {
                                return sink->patchNext->twin();
                            }
                        }
                        
                    }
                    else
                    {
                        
                        if (validStarts.find(sink->loopNext) != validStarts.end())
                        {
                            return sink->loopNext;
                        }
                    }
                }
                else if ((sink->patchNext->wn>0 && sink->loopNext->wn<0) ||
                (sink->patchNext->wn<0 && sink->loopNext->wn>0))
                {
                    //Bifurcation exists add into validstart
                    if (!sink->reversepatchNext->twin())
                    {
                        validStarts.insert(sink->reversepatchNext);
                    }
                    if (abs(sink->reversepatchNext->wn) > abs(sink->loopNext->wn))
                    {
                        if (!sink->reversepatchNext->twin())
                        {
                            return sink->reversepatchNext;
                        }
                        else 
                        {
                            if (validStarts.find(sink->reversepatchNext->twin()) != validStarts.end())
                            {
                                return sink->reversepatchNext->twin();
                            }
                        }
                    }
                    else
                    {

                        if (validStarts.find(sink->loopNext) != validStarts.end())
                        {
                            return sink->loopNext;
                        }
                    }
                }
                else if (sink->loopNext->wn==0 && (sink->patchNext->wn!=0 || sink->reversepatchNext->wn!=0))
                {
                    if (!sink->loopNext->twin())
                    {
                        if (sink->patchNext->wn != 0 && sink->reversepatchNext->wn == 0)
                        {
                                validStarts.insert(sink->patchNext);
                                return sink->patchNext;
                        }
                        else if (sink->patchNext->wn == 0 && sink->reversepatchNext->wn != 0)
                        {
                            validStarts.insert(sink->reversepatchNext);
                            return sink->reversepatchNext;
                        }
                        else if (sink->patchNext->wn != 0 && sink->reversepatchNext->wn != 0)
                        {
                            // std::cout<<"Coming here \n";
                            // std::cout<<"Coming here \n";

                            if (sink->loopPrev->wn != 0)
                            {
                                if ((sink->loopPrev->wn < 0 && sink->patchNext->wn < 0 && sink->reversepatchNext->wn > 0) || (sink->loopPrev->wn > 0 && sink->patchNext->wn > 0 && sink->reversepatchNext->wn < 0))
                                {
                                    validStarts.insert(sink->patchNext);
                                    return sink->patchNext;
                                }
                                else if ((sink->loopPrev->wn < 0 && sink->patchNext->wn > 0 && sink->reversepatchNext->wn < 0) || (sink->loopPrev->wn > 0 && sink->patchNext->wn < 0 && sink->reversepatchNext->wn > 0))
                                {
                                    validStarts.insert(sink->reversepatchNext);
                                    return sink->reversepatchNext;
                                }
                                else if ((sink->loopPrev->wn < 0 && sink->patchNext->wn < 0 && sink->reversepatchNext->wn < 0) || (sink->loopPrev->wn > 0 && sink->patchNext->wn > 0 && sink->reversepatchNext->wn > 0))
                                {
                                    assert(false && "Both patch next and reverse patch next can be traversed in this case. Need to look for this case");
                                }
                            }
                        }
                    }
                    else
                    {
                        if (validStarts.find(sink->loopNext) != validStarts.end())
                        {
                            return sink->loopNext;
                        }
                    }
                    
                }
                else if (sink->loopNext->wn!=0 && (sink->patchNext->wn==0 || sink->reversepatchNext->wn==0))
                {
                    if (!sink->loopNext->twin())
                    {
                        if (sink->loopNext->wn != 0 && ((sink->patchNext->wn == 0 && sink->reversepatchNext->wn == 0)))
                        {
                            if (validStarts.find(sink->loopNext) != validStarts.end())
                            {
                                return sink->loopNext;
                            }
                        }
                        else if ((sink->loopNext->wn > 0 && sink->patchNext->wn > 0 && sink->reversepatchNext->wn == 0) ||
                                (sink->loopNext->wn < 0 && sink->patchNext->wn < 0 && sink->reversepatchNext->wn == 0))
                        {
                            validStarts.insert(sink->patchNext);
                            return sink->patchNext;
                        }
                        else if ((sink->loopNext->wn > 0 && sink->patchNext->wn == 0 && sink->reversepatchNext->wn > 0) ||
                                (sink->loopNext->wn < 0 && sink->patchNext->wn == 0 && sink->reversepatchNext->wn < 0))
                        {
                            validStarts.insert(sink->reversepatchNext);
                            return sink->reversepatchNext;
                        }
                    }
                    else
                    {
                        if (validStarts.find(sink->loopNext) != validStarts.end())
                        {
                            return sink->loopNext;
                        }
                    }
                }
                else if (sink->loopNext->wn==0 && sink->patchNext->wn==0 && sink->reversepatchNext->wn==0)
                {
                    //either the patchnext or the reverse patch next should be the twin of loopnext
                    bool looptwin((sink->loopNext->twin()==sink->patchNext)||(sink->loopNext->twin()==sink->reversepatchNext));
                    assert(looptwin);
                    if (validStarts.find(sink->loopNext) != validStarts.end())
                        {
                            return sink->loopNext;
                        }
                }
            }
            else
            {
                if (validStarts.find(sink->loopNext) != validStarts.end())
                {
                    // validStarts.erase(validStarts.find(sink->loopNext));
                    return sink->loopNext;
                }
            }
        }
        else if (isPatchEdg)
        {
            assert(sink->patchNext && "patch must exist");
            // std::cout<<"In patch edge \n";
            if (sink->loopNext)
            {
                // std::cout<<"Coming in loopNext \n";
                //Bifurcation exists INSERT....
                if (!sink->patchNext->twin())
                {
                        // std::cout<<"inserting "<<sink->patchNext->source->sID<<"--->"<<sink->patchNext->sink->sID<<std::endl;

                        validStarts.insert(sink->patchNext);
                }
                

                if (abs(sink->loopNext->wn) > abs(sink->patchNext->wn))
                {
                    if (validStarts.find(sink->loopNext) != validStarts.end())
                    {
                        // validStarts.erase(validStarts.find(sink->loopNext));
                        return sink->loopNext;
                    }
                    //                return sink->loopNext;
                }
                else
                {

                    if (!sink->patchNext->twin())
                    {
                        // std::cout<<"returning "<<sink->patchNext->source->sID<<"--->"<<sink->patchNext->sink->sID<<std::endl;
                        return sink->patchNext;
                    }
                }
            }
            else
            {
                validStarts.insert(sink->patchNext);
                return sink->patchNext;
            }
            
        }

        else if (isReversePatchEdg)
        {
            assert(sink->reversepatchNext && "patch must exist");
            if (sink->loopNext)
            {
                //Bifurcation exists INSERT....
                if (!sink->reversepatchNext->twin())
                {
                    validStarts.insert(sink->reversepatchNext);
                }
                

                if (abs(sink->loopNext->wn) > abs(sink->reversepatchNext->wn))
                {
                    if (validStarts.find(sink->loopNext) != validStarts.end())
                    {
                        // validStarts.erase(validStarts.find(sink->loopNext));
                        return sink->loopNext;
                    }
                    //                return sink->loopNext;
                }
                else
                {
                    if (!sink->reversepatchNext->twin())
                    {
                        return sink->reversepatchNext;
                    }
                }
            }
            else
            {
                validStarts.insert(sink->reversepatchNext);
                return sink->reversepatchNext;
            }
        }

        //Check validStarts to see if there exist some loop edge with nonzero winding number
        for (const auto& vStarts : validStarts)
        {
            if ( vStarts->wn!=0)
            {
                return vStarts;
            }
        }
        return nullptr;
    }

    const ClipperEdge* ClipperEdge::twin() const
    {
        if(isLoop())
        {
            if(sink->patchPrev)
            {
                if(sink->patchPrev->source==source)
                {// a boundary edge
                    return sink->patchPrev;
                }
                else if (sink->reversepatchPrev->source==source) //Added by Yash
                {
                    return sink->reversepatchPrev;
                }
            }
        }
        else if(isPatch())
        {
            if(sink->loopPrev)
            {
                if(sink->loopPrev->source==source)
                {// a boundary edge
                    return sink->loopPrev;
                }
            }
        }
        else if (isReversePatch()) //Added by Yash
        {
            if (sink->loopPrev)
            {
                if (sink->loopPrev->source == source)
                { // a boundary edge
                    return sink->loopPrev;
                }
            }
        }
        else
        {
            assert(false);
        }
        return nullptr;
    }

    std::string ClipperEdge::tag() const
    {
        assert(source.get() && "BAD SOURCE");
            assert(sink.get() && "BAD SINK");
        return std::to_string(source->sID)+"->"+std::to_string(sink->sID);
    }

    /******************************************************************************/
    LoopEdge::LoopEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
            const std::shared_ptr<LoopPathClipperNode>& sink_in,
            const int& wn_in) :
    ClipperEdge(source_in,sink_in,wn_in)
    {
        assert(source->loopNext==nullptr);
        source->loopNext=this;
        assert(sink->loopPrev==nullptr);
        sink->loopPrev=this;
    }

    /******************************************************************************/
    PatchEdge::PatchEdge(const std::shared_ptr<LoopPathClipperNode>& source_in,
                        const std::shared_ptr<LoopPathClipperNode>& sink_in,
                        const int& wn_in) :
    ClipperEdge(source_in,sink_in,wn_in)
    {
        assert(source->patchNext==nullptr);
        source->patchNext=this;
        assert(sink->patchPrev==nullptr);
        sink->patchPrev=this;
    }

    ReversePatchEdge::ReversePatchEdge (const std::shared_ptr<LoopPathClipperNode>& source_in,
                        const std::shared_ptr<LoopPathClipperNode>& sink_in,
                        const int& wn_in) :
                    ClipperEdge(source_in,sink_in,wn_in)
    {
        assert(source->reversepatchNext==nullptr);
        source->reversepatchNext=this;
        assert(sink->reversepatchPrev==nullptr);
        sink->reversepatchPrev=this;
    };

    /******************************************************************************/
    struct LoopPathClipper : private std::map<std::pair<size_t,size_t>,LoopEdge>
    /*                    */,private std::map<std::pair<size_t,size_t>,PatchEdge>
    /*                    */,private std::map<std::pair<size_t,size_t>,ReversePatchEdge>
    /*                    */,private std::map<size_t,const std::shared_ptr<LoopPathClipperNode>>
    /*                    */,private std::vector<std::vector<const ClipperEdge*>>
    {
        typedef Eigen::Matrix<double,2,1> VectorDim;
        typedef std::map<size_t,const std::shared_ptr<LoopPathClipperNode>> NodeCointainerType;

        //Member variable added by Yash
        bool patchCompletelyInside;

        NodeCointainerType& nodes()
        {
            return *this;
        }

        const NodeCointainerType& nodes() const
        {
            return *this;
        }
        
        const std::map<std::pair<size_t,size_t>,LoopEdge>& loopEdges() const
        {
            return *this;
        }

        std::map<std::pair<size_t,size_t>,LoopEdge>& loopEdges()
        {
            return *this;
        }
        
        const std::map<std::pair<size_t,size_t>,PatchEdge>& patchEdges() const
        {
            return *this;
        }
        
        std::map<std::pair<size_t,size_t>,PatchEdge>& patchEdges()
        {
            return *this;
        }

        const std::map<std::pair<size_t,size_t>,ReversePatchEdge>& reversepatchEdges() const
        {
            return *this;
        }
        
        std::map<std::pair<size_t,size_t>,ReversePatchEdge>& reversepatchEdges()
        {
            return *this;
        }

        const std::vector<std::vector<const ClipperEdge*>>& clippedPolygons() const
        {
            return *this;
        }

        
        std::vector<std::vector<const ClipperEdge*>>& clippedPolygons()
        {
            return *this;
        }
        
        
        
        LoopPathClipper(const std::vector<VectorDim>& loop,
                        const std::vector<VectorDim>& patch)
        {
            std::map<float,std::shared_ptr<LoopPathClipperNode>> loopIntersections;  // intersection points along loop. Key with float precision to round off errors
            std::map<float,std::shared_ptr<LoopPathClipperNode>> patchIntersections; // intersection points along patch. Key with float precision to round off errors
            // std::cout<<"Printing loop points \n";
            // for (const auto &lPoints : loop)
            // {
            //     std::cout<<lPoints.transpose()<<std::endl;
            // }

            // std::cout << "Printing patch points \n";
            // for (const auto &pPoints : patch)
            // {
            //     std::cout << pPoints.transpose() << std::endl;
            // }
            // std::cout<<"PRINTING DONE... \n";
            for (size_t i = 0; i < loop.size(); ++i)
            {
                const size_t i1(i + 1 < loop.size() ? i + 1 : 0);
                for (size_t j = 0; j < patch.size(); ++j)
                {
                    const size_t j1(j + 1 < patch.size() ? j + 1 : 0);
                    SegmentSegmentDistance<2> ssd(loop[i], loop[i1], patch[j], patch[j1]);
                    for (const auto &inter : ssd.intersectionSegment())
                    {
                        const double t(fmod(i + std::get<1>(inter), loop.size()));  // intersection parameter in loopIntersections
                        const double u(fmod(j + std::get<2>(inter), patch.size())); // intersection parameter in patchIntersections

                        //                    std::cout<<t<<" "<<u<<std::endl;
                        const auto loopIter(loopIntersections.find(t));
                        const auto patchIter(patchIntersections.find(u));

                        if (loopIter != loopIntersections.end() && patchIter != patchIntersections.end())
                        { // both u and t exist, so intersection was already found
                            assert(loopIter->second == patchIter->second && "INTERSECTION MUST BE SAME NODE");
                        }
                        else if (loopIter != loopIntersections.end() && patchIter == patchIntersections.end())
                        {   // u not found and t found. Intersection point exists on loop only
                            //                        assert(false && "Intersection point exists on patch only. FINISH HERE");
                            patchIntersections.emplace(u, loopIter->second);
                        }
                        else if (loopIter == loopIntersections.end() && patchIter != patchIntersections.end())
                        {   // u found and t not found. Intersection point exists on patch only
                            //                        assert(false && "Intersection point exists on loop only. FINISH HERE");
                            loopIntersections.emplace(t, patchIter->second);
                        }
                        else
                        { // neither u nor t exist. New intersection point
                            std::shared_ptr<LoopPathClipperNode> newNode(new LoopPathClipperNode(std::get<0>(inter)));
                            nodes().emplace(newNode->sID, newNode);
                            loopIntersections.emplace(t, newNode);
                            patchIntersections.emplace(u, newNode);
                        }
                    }
                }
            }
    //         for(size_t i=0;i<loop.size();++i)
    //         {
    //             const size_t i1(i+1<loop.size()? i+1 :0);
    //             for(size_t j=0;j<patch.size();++j)
    //             {
    //                 const size_t j1(j+1<patch.size()? j+1 :0);
    //                 SegmentSegmentDistance<2> ssd(loop[i],loop[i1],patch[j],patch[j1]);
    //                 for(const auto& inter : ssd.intersectionSegment())
    //                 {
    //                     const double t(fmod(i+std::get<1>(inter), loop.size())); // intersection parameter in loopIntersections
    //                     const double u(fmod(j+std::get<2>(inter),patch.size())); // intersection parameter in patchIntersections
                        
    // //                    std::cout<<t<<" "<<u<<std::endl;
    //                     const auto loopIter(  loopIntersections.find(t));
    //                     const auto patchIter(patchIntersections.find(u));

    //                     if(    loopIter!= loopIntersections.end()
    //                     && patchIter!=patchIntersections.end())
    //                     {// both u and t exist, so intersection was already found
    //                         assert(loopIter->second==patchIter->second && "INTERSECTION MUST BE SAME NODE");
                            
    //                     }
    //                     else if(    loopIter== loopIntersections.end()
    //                             && patchIter!=patchIntersections.end())
    //                     {// u not found and t found. Intersection point exists on loop only
    //                         std::ofstream ofsLoop("polyPoints.txt");
    //                         for (const auto &lpoints : loop)
    //                         {
    //                             //windingNumbers.push_back();
    //                             ofsLoop << std::setprecision(15) << lpoints.transpose() << "\n";
    //                         }
    //                         ofsLoop.close();
    //                         std::ofstream ofsPatch("patchPoints.txt");
    //                         for (const auto &ppoints : patch)
    //                         {
    //                             ofsPatch << std::setprecision(15) << ppoints.transpose() << "\n";
    //                         }

    //                         ofsPatch.close();
    //                         assert(false && "Intersection point exists on patch only. FINISH HERE");
    //                     }
    //                     else if(    loopIter!= loopIntersections.end()
    //                             && patchIter==patchIntersections.end())
    //                     {// u found and t not found. Intersection point exists on patch only
    //                         std::ofstream ofsLoop("polyPoints.txt");
    //                         for (const auto &lpoints : loop)
    //                         {
    //                             //windingNumbers.push_back();
    //                             ofsLoop << std::setprecision(15) << lpoints.transpose() << "\n";
    //                         }
    //                         ofsLoop.close();
    //                         std::ofstream ofsPatch("patchPoints.txt");
    //                         for (const auto &ppoints : patch)
    //                         {
    //                             ofsPatch << std::setprecision(15) << ppoints.transpose() << "\n";
    //                         }

    //                         ofsPatch.close();
    //                         assert(false && "Intersection point exists on loop only. FINISH HERE");
    //                     }
    //                     else
    //                     {// neither u nor t exist. New intersection point
    //                         std::shared_ptr<LoopPathClipperNode> newNode( new LoopPathClipperNode(std::get<0>(inter)));
    //                         nodes().emplace(newNode->sID,newNode);
    //                         loopIntersections.emplace(t,newNode);
    //                         patchIntersections.emplace(u,newNode);
    //                     }
    //                 }

    //             }
    //         }
            
            patchCompletelyInside=(!patchIntersections.size() && !loopIntersections.size()); //Location of this declaration is important

            for(size_t i=0;i<loop.size();++i)
            {// insert loop vertices as loopIntersections points, if not already existing
                if(loopIntersections.find(i+0.0)==loopIntersections.end())
                {
                    std::shared_ptr<LoopPathClipperNode> newNode( new LoopPathClipperNode(loop[i]));
                    nodes().emplace(newNode->sID,newNode);
                    loopIntersections.emplace(i+0.0,newNode);
                }
            }
            
            for(size_t j=0;j<patch.size();++j)
            {// insert patch vertices as patchintersection points, if not already existing
                if(patchIntersections.find(j+0.0)==patchIntersections.end())
                {
                    std::shared_ptr<LoopPathClipperNode> newNode( new LoopPathClipperNode(patch[j]));
                    nodes().emplace(newNode->sID,newNode);
                    patchIntersections.emplace(j+0.0,newNode);
                }
            }

            
            for(auto iter=loopIntersections.begin();iter!=loopIntersections.end();++iter)
            {
                auto iter1(iter);
                iter1++;
                if(iter1==loopIntersections.end())
                {
                    iter1=loopIntersections.begin();
                }
                const auto source(iter->second);
                const auto   sink(iter1->second);
                const VectorDim midPoint(0.5*(*source+*sink));
                loopEdges().emplace(std::piecewise_construct,
                                    std::forward_as_tuple(source->sID,sink->sID),
                                    std::forward_as_tuple(source,sink,Polygon2D::windingNumber(midPoint,patch))
                                    );
            }

            for(auto iter=patchIntersections.begin();iter!=patchIntersections.end();++iter)
            {
                auto iter1(iter);
                iter1++;
                if(iter1==patchIntersections.end())
                {
                    iter1=patchIntersections.begin();
                }
                const auto source(iter->second);
                const auto   sink(iter1->second);
                const VectorDim midPoint(0.5*(*source+*sink));
                patchEdges().emplace(std::piecewise_construct,
                                    std::forward_as_tuple(source->sID,sink->sID),
                                    std::forward_as_tuple(source,sink,Polygon2D::windingNumber(midPoint,loop))
                                    );
                reversepatchEdges().emplace(std::piecewise_construct,
                                    std::forward_as_tuple(sink->sID,source->sID),
                                    std::forward_as_tuple(sink,source,-1*Polygon2D::windingNumber(midPoint,loop))  //Added by Yash
                                    );
            }

        }

        void makePaths()
        {
            //        std::vector<std::vector<VectorDim>> clippedPolygons;
            // std::cout << "PatchCompletelyInside =" << patchCompletelyInside << std::endl;
            if (!patchCompletelyInside)
            {
                std::set<const ClipperEdge *> validStarts;
                for (const auto &edg : loopEdges())
                {
                    if (edg.second.wn != 0 || edg.second.twin())
                    {
                        validStarts.insert(&edg.second);
                    }
                }

            bool clipperboolIterator(true);
                while (clipperboolIterator)
                {
                    std::vector<const ClipperEdge *> path;
                    clipperboolIterator = false;

                    for (const auto &vstarts : validStarts)
                    {
                        if (vstarts->wn != 0 )
                        {
                            clipperboolIterator = true;
                            const ClipperEdge *edg(vstarts);
                            clipPolygon(edg, path, validStarts);
                            clipperboolIterator=true;
                            break;
                        }
                    }
        
                    if (path.size())
                    {
                        // std::cout<<"Creating Path created \n";
                        clippedPolygons().push_back(path);
                        // std::cout<<"Path created \n";
                        // for (const auto& cp : clippedPolygons())
                        // {
                        //     std::cout<<"Pathi \n";
                        //     for(const auto& cpEdge : cp)
                        //     {
                        //         std::cout<<cpEdge->source->sID<<"--->"<<cpEdge->sink->sID<<std::endl;
                        //     }
                        // }
                    }

                }
            }
            else
            {
                bool sameWNPatch;
                for (const auto &edg : patchEdges())
                {
                    sameWNPatch = (edg.second.wn == (patchEdges().begin())->second.wn);
                }

                bool sameWNLoop;
                for (const auto &edg : loopEdges())
                {
                    sameWNLoop = (edg.second.wn == (loopEdges().begin())->second.wn);
                }

                assert(sameWNPatch && "WN has to be same....Error"); //Check that the winding numbers are the same
                assert(sameWNLoop && "WN has to be same....Error"); //Check that the winding numbers are the same

    //Include a for loop for inserting the patch edges wn number of times
                if ((patchEdges().begin())->second.wn != 0)
                {
                        std::cout << "Printing the wn for the patches " << (patchEdges().begin())->second.wn << std::endl;
                    if ((patchEdges().begin())->second.wn > 0)
                    {
                        std::set<const ClipperEdge *> validStarts;
                        validStarts.insert(&((patchEdges().begin())->second));
                        const ClipperEdge *edg(*validStarts.begin());

                        std::vector<const ClipperEdge *> path;
                        clipPolygon(edg, path, validStarts);
                        if (path.size())
                        {
                            clippedPolygons().push_back(path);
                        }
                    }
                    else if ((reversepatchEdges().begin())->second.wn > 0)
                    {
                        std::cout << "Printing the wn for the reverse patches " << (reversepatchEdges().begin())->second.wn << std::endl;

                        std::set<const ClipperEdge *> validStarts;
                        validStarts.insert(&((reversepatchEdges().begin())->second));
                        const ClipperEdge *edg(*validStarts.begin());

                        std::vector<const ClipperEdge *> path;
                        clipPolygon(edg, path, validStarts);
                        if (path.size())
                        {
                            clippedPolygons().push_back(path);
                        }
                    }
                    else
                    {
                        assert(false && "Either the patch or the reverse patch edge has to be inserted");
                    }
                    
                }
                else if ((loopEdges().begin())->second.wn != 0)
                {
                    //Insert the loop
                    std::vector<const ClipperEdge *> path;
                    for (const auto& lEdges : loopEdges())
                    {
                        assert(lEdges.second.wn!=0);
                        path.push_back(&lEdges.second);
                    }
                    if (path.size())
                    {
                        clippedPolygons().push_back(path);
                    }

                }
                else
                {
                    if ((loopEdges().begin())->second.wn != 0 || (patchEdges().begin())->second.wn != 0)
                    {
                        assert(0 && "Either the loopEdges or the pathc edges have to be inserted ");
                    }
                }
                
            }
        }

        void clipPolygon(const ClipperEdge *edg,
                        std::vector<const ClipperEdge *> &path,
                        std::set<const ClipperEdge *> &validStarts) const
        {

            // std::cout << "edg= " << edg->tag() <<"is (Loop Patch or ReversePatch)"<<edg->isLoop()<<edg->isPatch()<<edg->isReversePatch()<< std::endl;
            // std::cout << "printing valid starts (source->sink) \n";
            //     if (validStarts.size() > 0)
            //     {
            //         for (const auto &validStarttemp : validStarts)
            //         {
            //             std::cout << validStarttemp->source->sID << "--->" << validStarttemp->sink->sID << "\n";
            //         }
            //     }
            assert(validStarts.find(edg) != validStarts.end());
            validStarts.erase(validStarts.find(edg));

            const auto twin(edg->twin());
            // need to determine if edg goes in path
            if (edg->wn != 0)
            { // add edg to path
                path.push_back(edg);
            }
            else if (edg->wn == 0)
            { // check if this is a boundary edge
                if (twin)
                {
                    path.push_back(twin);
                }
            }

            if (path.size())
            {
                if (path.back()->sink != path.front()->source)
                { // edg is same as beginning of path. Stop
                    const ClipperEdge* edg_next = edg->next(validStarts);
                    clipPolygon(edg_next, path, validStarts);
                    //                    clipPolygon(edg->next(),path,checkedEdges);
                }
                else
                {
                    //This will update valid start
                    const ClipperEdge* edg_next = edg->next(validStarts);
                    if(edg_next == path.front())
                    {
                        validStarts.erase(validStarts.find(edg_next));
                    }
                }
                
            }
            else
            {
                const ClipperEdge* edg_next = edg->next(validStarts);
                clipPolygon(edg_next,path, validStarts);
                //                clipPolygon(edg->next(),path,checkedEdges);
            }
        }

        void printRelevant()
        {
            std::cout<<"Print loopEdges \n";
            for (const auto& lEdges : loopEdges())
            {
                std::cout<<(*lEdges.second.source).transpose()<<"-->"<<(*lEdges.second.sink).transpose()<<std::endl;

            }

            std::cout<<"Print patcEdges \n";

            for (const auto &pEdges : patchEdges())
            {
                std::cout << (*pEdges.second.source).transpose() << "-->" << (*pEdges.second.sink).transpose() << std::endl;
            }
        }

         std::vector<std::vector<VectorDim>> getClippedPolygons()
         {
             std::vector<std::vector<VectorDim>> nodeSeq;
             for (size_t k = 0; k < clippedPolygons().size(); ++k)
             {
                 std::vector<VectorDim> nodeSeqK;
                 for (const auto &edg : this->clippedPolygons()[k])
                 {
                     nodeSeqK.emplace_back(*(edg->source));
                    //  std::cout << edg->tag() << std::endl;
                 }
                 assert(nodeSeqK.size() > 0 && "Clipped Container can't be empty");
                 nodeSeq.emplace_back(nodeSeqK);
                 // assert(this->clippedPolygons()[k].back()->sink == this->clippedPolygons()[k].front()->source && "ERROR, OPEN PATH");
             }
             return nodeSeq;
         }

        /**********************************************************************/
        template <class T>
        friend T &operator<<(T &os, const LoopPathClipper &lpc)
        {
            for (const auto &edge : lpc.loopEdges())
            {
                os << edge.second.source->transpose() << " " << edge.second.sink->transpose() << " " << edge.second.wn << "\n";
            }
            for (const auto &edge : lpc.patchEdges())
            {
                os << edge.second.source->transpose() << " " << edge.second.sink->transpose() << " " << edge.second.wn << "\n";
            }
            os << std::flush;
            return os;
            }
    };

}
#endif