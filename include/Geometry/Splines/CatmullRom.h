/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_CatmullRom_H_
#define model_CatmullRom_H_

#include <set>
#include <Eigen/Dense>
#include <SplineSegment.h>
#include <LoopLink.h>

namespace model
{
    
    struct CatmullRom
    {        
        
        /**********************************************************************/
        template <typename T>
        static Eigen::Matrix<double,T::dim,1> loopTangent(const std::set<LoopLink<T>*>& segmentSet)
        {
            
            double cT=0.0;
            for(const auto& link : segmentSet)
            {
                const double cL=link->pLink->parametricChordLength();
                cT+=cL;
            }
            
            Eigen::Matrix<double,T::dim,1> temp=Eigen::Matrix<double,T::dim,1>::Zero();
            if(cT>0.0)
            {
                for(const auto& link : segmentSet)
                {
                    const double cL=link->pLink->parametricChordLength();
                    if(cL>0.0)
                    {
                        temp+=(link->sink()->get_P()-link->source()->get_P())/cL*(cT-cL)/cT;
                    }
                }
            }

            
            return temp;
        }
        
        /**********************************************************************/
        template <typename T>
        static std::map<size_t,std::pair<double,Eigen::Matrix<double,T::dim,1>>> loopTangentCoeffs(const std::set<LoopLink<T>*>& segmentSet)
        {/*!\returns a map<snID,pair>, where snID is the NetworkComponent ID of
          * a node, pair.second is the position of a neighbor node, and pair.first
          * is the scalar coefficient which multiplies the position of the 
          * neighbor node in the calculation of the loop tangent at the node snID.
          */
            
            double cT=0.0;
            for(const auto& link : segmentSet)
            {
                const double cL=link->pLink->parametricChordLength();
                cT+=cL;
            }
            
            //Eigen::Matrix<double,T::dim,1> temp=Eigen::Matrix<double,T::dim,1>::Zero();
            std::map<size_t,std::pair<double,Eigen::Matrix<double,T::dim,1>>> temp;
            
            for(const auto& link : segmentSet)
            {
                temp[link->source()->snID()]=std::make_pair(0.0,link->source()->get_P());
                temp[link->  sink()->snID()]=std::make_pair(0.0,link->  sink()->get_P());
            }
            
            if(cT>0.0)
            {
                for(const auto& link : segmentSet)
                {
                    const double cL=link->pLink->parametricChordLength();
                    if(cL>0.0)
                    {
                        temp[link->source()->snID()].first-=1.0/cL*(cT-cL)/cT;
                        temp[link->  sink()->snID()].first+=1.0/cL*(cT-cL)/cT;
                    }
                }
            }

            
            return temp;
        }
        
    };
    
}
#endif
