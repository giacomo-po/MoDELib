/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <deque>
#include <cfloat>
#include <experimental/filesystem>
#include <model/DislocationDynamics/IO/EVLio.h>

using namespace model;

namespace fs = std::experimental::filesystem;


int main(int argc, char** argv)
{
    const std::string path = "./evl";

    const std::string subS1="evl_";
    const std::string subS2=".bin";
    const std::string subS3=".txt";

    
    std::map<size_t,std::deque<double>> fileMap;
    
    for (const auto & p : fs::directory_iterator(path))
    {
        
        const size_t found1=p.path().string().find(subS1);
        const size_t found2=p.path().string().find(subS2);
        const size_t found3=p.path().string().find(subS3);
        
        if (   found1 != std::string::npos
            && (found2 != std::string::npos || found3 != std::string::npos))
        {
            const std::string stringID=(found2 != std::string::npos)? p.path().string().substr(found1+subS1.length(),found2-(found1+subS1.length()))
            /*                                                   */ : p.path().string().substr(found1+subS1.length(),found3-(found1+subS1.length()));
            const size_t runID=std::stol(stringID);
            
            if(fileMap.find(runID)==fileMap.end())
            {
                EVLio<3> evlio;
                if (EVLio<3>::isBinGood(runID))
                {
                    evlio.readBin(runID);
                }
                else if(EVLio<3>::isTxtGood(runID))
                {
                    evlio.readTxt(runID);
                }
                
                std::map<std::pair<size_t,size_t>,DislocationSegmentIO<3>> segmentMap=evlio.segments();
                
                
                std::map<size_t,const DislocationNodeIO<3>* const> nodeMap=evlio.nodeMap();
                
                double screwLength=0.0;
                double edgeLength=0.0;
                
                
                
                for(const auto& segmentPair : segmentMap)
                {
                    
                    const double bNorm(segmentPair.second.b.norm());
                    const double nNorm(segmentPair.second.n.norm());
                    
                    if(segmentPair.second.meshLocation==0
                       && bNorm>FLT_EPSILON
                       && nNorm>FLT_EPSILON
                       && fabs(segmentPair.second.b.dot(segmentPair.second.n))<FLT_EPSILON*bNorm*nNorm)
                    {
                        const size_t& sourceID(segmentPair.first.first);
                        const size_t&   sinkID(segmentPair.first.second);
                        
                        
                        const DislocationNodeIO<3>& source(*nodeMap.find(sourceID)->second);
                        const DislocationNodeIO<3>&   sink(*nodeMap.find(  sinkID)->second);
                        
                        const Eigen::Matrix<double,3,1> T(sink.P-source.P);
                        const double L(T.norm());

                        
                        const double cosTheta=T.dot(segmentPair.second.b)/L/bNorm;
                        const double cosTheta2=cosTheta*cosTheta;
                        const double sinTheta2=1.0-cosTheta2;
                        
                        screwLength+=cosTheta2*L;
                        edgeLength +=sinTheta2*L;
                        
                        
                    }
                    
                }
                
                fileMap[runID].push_back(screwLength);
                fileMap[runID].push_back(edgeLength);

            }
            
            
            
        }
    }
    
    std::ofstream lengthFile ("lengthFile.txt");
    if (lengthFile.is_open())
    {
        for(const auto& pair : fileMap)
        {
            lengthFile<<pair.first;
            for(const double& val : pair.second)
            {
                lengthFile<<std::setprecision(15)<<std::scientific<<" "<<val;
            }
            lengthFile<<"\n";
        }
        lengthFile.close();
    }

    return 0;
    
}


