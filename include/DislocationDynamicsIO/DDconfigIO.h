/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2020 by Danny Perez <danny_perez@lanl.gov>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigIO_H_
#define model_DDconfigIO_H_

#include <type_traits>
#include <vector>
#include <map>
#include <chrono>
#include <iostream>     // std::cout
#include <fstream>      // std::ifstream
#include <cfloat>      // std::ifstream
#include <DDbaseIO.h>
#include <DislocationLoopIO.h>
#include <DislocationLoopNodeIO.h>
#include <DislocationLoopLinkIO.h>
#include <DislocationNodeIO.h>
#include <DislocationSegmentIO.h>
#include <EshelbyInclusionIO.h>
//#include <DislocationDynamicsModule.h>
//#include <DislocationNetwork.h>
//#include <PeriodicLoopIO.h>

#include <TerminalColors.h>
#include <set>


namespace model
{

    template <int dim>
    class DDconfigIO : public DDbaseIO
    /*              */,private std::vector<DislocationNodeIO<dim> >
    /*              */,private std::vector<DislocationLoopIO<dim> >
    /*              */,private std::vector<DislocationLoopNodeIO<dim> >
    /*              */,private std::vector<DislocationLoopLinkIO<dim> >
    /*              */,private std::vector<EshelbyInclusionIO<dim>>
    {
        
        
        std::map<size_t,const size_t> _loopNodeMap;
        std::map<size_t,const size_t> _nodeMap;
        std::map<size_t, const size_t> _loopMap;
        
        void make_maps();
        
        
    public:
        
        DDconfigIO(const std::string& folderName,const std::string& suffix="");
        
        
        const std::vector<DislocationNodeIO<dim> >& nodes() const;
        std::vector<DislocationNodeIO<dim> >& nodes();
        const std::vector<DislocationLoopNodeIO<dim> >& loopNodes() const;
        std::vector<DislocationLoopNodeIO<dim> >& loopNodes();
        const std::vector<DislocationLoopIO<dim> >& loops() const;
        std::vector<DislocationLoopIO<dim> >& loops();
        const std::vector<DislocationLoopLinkIO<dim> >& loopLinks() const;
        std::vector<DislocationLoopLinkIO<dim> >& loopLinks();
        const std::vector<EshelbyInclusionIO<dim> >& eshelbyInclusions() const;
        std::vector<EshelbyInclusionIO<dim> >& eshelbyInclusions();
        const std::map<size_t,const size_t>& nodeMap() const;
        std::map<size_t,const size_t>& nodeMap();
        const std::map<size_t,const size_t>& loopNodeMap() const;
        std::map<size_t,const size_t>& loopNodeMap();
        const std::map<size_t, const size_t>& loopMap() const;
        std::map<size_t, const size_t>& loopMap();
        std::map<std::pair<size_t,size_t>,DislocationSegmentIO<dim> > segments() const;
        void write(const size_t& runID,const bool& outputBinary);
        void writeTxt(const size_t& runID);
        void writeTxtStream(std::ostream &file);
        void writeBin(const size_t& runID);
        void read(const size_t& runID);
        void readBin(const size_t& runID);
        void readTxt(const size_t& runID);
        void readTxtStream(std::istream &infile);
        void bin2txt(const size_t& runID,const bool& writeSegments);
        std::map<size_t,std::vector<size_t>> loopNodeSequence() const;
        
        //        std::map<std::pair<size_t,size_t>,std::set<size_t> > segmentloopMap() const
        //        {
        //
        //           std::map<std::pair<size_t,size_t>,std::set<size_t> >temp;
        //
        //            for(const auto& link : loopLinks())
        //            {
        //
        //
        //                const auto loopIter=loopMap().find(link.loopID);
        //                assert(loopIter!=loopMap().end());
        //
        //                const auto loopSourceIter(loopNodeMap().find(link.sourceID));
        //                assert(loopSourceIter!=loopNodeMap().end());
        //
        //                const auto loopSinkIter(loopNodeMap().find(link.sinkID));
        //                assert(loopSinkIter!=loopNodeMap().end());
        //
        //
        //                if(link.hasNetworkLink)
        //                {
        //                    const size_t sourceID(std::min(loopSourceIter->second->networkNodeID,loopSinkIter->second->networkNodeID));
        //                    const size_t   sinkID(std::max(loopSourceIter->second->networkNodeID,loopSinkIter->second->networkNodeID));
        //                    const auto key=std::make_pair(sourceID,sinkID);
        //
        //                    temp[key].insert(link.loopID);
        //
        //                }
        //            }
        //
        //            return temp;
        //        }
        
        
        
};

}
#endif
