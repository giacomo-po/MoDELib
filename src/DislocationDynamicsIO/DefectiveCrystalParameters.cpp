/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystalParameters_cpp_
#define model_DefectiveCrystalParameters_cpp_

#include <DefectiveCrystalParameters.h>

namespace model
{
    

        
                /**********************************************************************/
        void DefectiveCrystalParameters::manageRestart()
        {
            // Menage restart
            IDreader<1,200,double> vReader(traitsIO.fFolder+"/F");
            vReader.readLabelsFile(traitsIO.fFolder+"/F_labels.txt");
            if (vReader.isGood(0,true))
            {// F/F_0.txt exists
                vReader.read(0,true);
                if(runID<0)
                {// Restart from last available step
                    if(vReader.size())
                    {// at least a line is available
//                        vReader.readLabelsFile("F/F_labels.txt");
                        runID=vReader.rbegin()->first;
                        totalTime=vReader.last("time [b/cs]");
                        dt=vReader.last("dt [b/cs]");
                    }
                    else
                    {// file is empty, keep default initialization
                        runID=0;
                    }
                }
                else
                {// Restart from specific runID
//                    const auto iter=vReader.find(runID);
                    totalTime=vReader(runID,"time [b/cs]");
                    dt=vReader(runID,"dt [b/cs]");
                }
            }
            else
            {// F/F_0.txt is not there, keep default initialization
                std::cout<<"Unable to read F/F_0.txt"<<std::endl;
                runID=0;
            }
            
            std::cout<<"starting at time step "<<runID<<std::endl;
            std::cout<<"totalTime= "<<totalTime<<std::endl;
            std::cout<<"dt= "<<dt<<std::endl;
        }

        /**********************************************************************/
        std::set<int> DefectiveCrystalParameters::getSubCyclingSet(const std::vector<int> &inpVector)
        {
            std::set<int> temp;
            for (const auto &iv : inpVector)
            {
                temp.insert(iv);
            }
            return temp;
        }

        /**********************************************************************/
        DefectiveCrystalParameters::DefectiveCrystalParameters(const std::string& folderName) :
        /* init */ traitsIO(folderName)
        /* init */,simulationType(TextFileParser(traitsIO.ddFile).readScalar<int>("simulationType",true))
        /* init */,useDislocations(TextFileParser(traitsIO.ddFile).readScalar<int>("useDislocations",true))
        /* init */,useCracks(TextFileParser(traitsIO.ddFile).readScalar<int>("useCracks",true))
        /* init */,periodicImageSize(simulationType==DDtraitsIO::PERIODIC_IMAGES? TextFileParser(traitsIO.ddFile).readArray<int>("periodicImageSize",true) : std::vector<int>())
        /* init */,Nsteps(TextFileParser(traitsIO.ddFile).readScalar<size_t>("Nsteps",true))
        /* init */,timeIntegrationMethod(TextFileParser(traitsIO.ddFile).readScalar<int>("timeIntegrationMethod",true))
        /* init */,useSubCycling(TextFileParser(traitsIO.ddFile).readScalar<int>("useSubCycling",true))
        /* init */,subcyclingBins(getSubCyclingSet(TextFileParser(traitsIO.ddFile).readArray<int>("subcyclingBins",true)))
        /* init */,externalLoadControllerName(TextFileParser(traitsIO.ddFile).readString("externalLoadControllerName",true))
        /* init */,virtualSegmentDistance((simulationType==DDtraitsIO::FINITE_FEM || simulationType==DDtraitsIO::PERIODIC_FEM || simulationType==DDtraitsIO::PERIODIC_IMAGES)? TextFileParser(traitsIO.ddFile).readScalar<double>("virtualSegmentDistance",true) : 0.0)
        /* init */,use_stochasticForce(TextFileParser(traitsIO.ddFile).readScalar<int>("use_stochasticForce",true))
        /* init */,periodicFaceIDs(TextFileParser(traitsIO.polyFile).template readSet<int>("periodicFaceIDs",true))
        /* init */,runID(TextFileParser(traitsIO.ddFile).readScalar<long int>("startAtTimeStep",true))
        /* init */,totalTime(0.0)
        /* init */,dt(10.0)
        {
            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");

            manageRestart();

            
        }
        

        
        bool DefectiveCrystalParameters::isPeriodicSimulation() const
        {
            return simulationType==DDtraitsIO::PERIODIC_IMAGES || simulationType==DDtraitsIO::PERIODIC_FEM;
        }
        
  }
#endif
