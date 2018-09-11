/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SweepPlane_H_
#define model_SweepPlane_H_

#include <cfloat>
#include <tuple>
#include <map>
#include <set>
#include <deque>
#include <algorithm>
#include <Eigen/Dense>

#include <model/Geometry/SegmentSegmentDistance.h>

namespace model
{
    
    /*! https://en.wikipedia.org/wiki/Bentley–Ottmann_algorithm
     */
    template <typename LineType,int dim>
    struct SweepPlane :
    /* base */ public std::map<double,std::pair<std::set<const LineType*>,std::set<const LineType*>>>,
    /* base */ public std::deque<std::pair<const LineType*,const LineType*>>
    {
        
        typedef std::set<const LineType*> EventSetType;
        typedef std::map<double,std::pair<EventSetType,EventSetType>> EventContainerType;
        typedef std::deque<std::pair<const LineType*,const LineType*>> IntersectionPairContainerType;
        /**********************************************************************/
        const EventContainerType& events() const
        {
            return *this;
        }
        
        /**********************************************************************/
        EventContainerType& events()
        {
            return *this;
        }
        
        const IntersectionPairContainerType& potentialIntersectionPairs() const
        {
            return *this;
        }
        
        const std::pair<const LineType*,const LineType*>& potentialIntersectionPair(const size_t& k) const
        {
            return potentialIntersectionPairs()[k];
        }
        
        /**********************************************************************/
        void addSegment(const double& xStart,const double& xEnd,const LineType& line)
        {
            
            const bool success1=events()[std::min(xStart,xEnd)].first.insert(&line).second;     // line starts at xStart
            const bool success2=events()[std::max(xStart,xEnd)].second.insert(&line).second;    // line ends at xEnd
            assert(success1==success2);
        }
        
        /**********************************************************************/
        void computeIntersectionPairs()
        {
            std::set<const LineType*> activeSegments;
            
            for(const auto& event : events())
            {
                
                for(const auto& startingSegment : event.second.first)
                {
                    for(const auto& activeSegment : activeSegments)
                    {
                        if(activeSegment!=startingSegment)
                        {
                            IntersectionPairContainerType::emplace_back(activeSegment,startingSegment);
                        }
                        
                    }
                    const bool success=activeSegments.insert(startingSegment).second;
                    assert(success);
                }
                
                for(const auto& endingSegment : event.second.second)
                {
                    const size_t nErased=activeSegments.erase(endingSegment);
                    assert(nErased==1);
                }
            }
            
            assert(activeSegments.size()==0); // no segments must be active when plane finishes sweep
            
        }
        
        
    };
    
}
#endif
