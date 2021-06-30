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
                        const std::vector<VectorDim>& patch);

        void makePaths();

        void clipPolygon(const ClipperEdge *edg,
                        std::vector<const ClipperEdge *> &path,
                         std::set<const ClipperEdge *> &validStarts) const;

        void printRelevant();

        std::vector<std::vector<VectorDim>> getClippedPolygons();

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
