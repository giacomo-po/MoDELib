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
            IDreader<'F',1,200,double> vReader;
            vReader.readLabelsFile("F/F_labels.txt");
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
        DefectiveCrystalParameters::DefectiveCrystalParameters(int& , char* []) :
        /* init */ simulationType(TextFileParser("inputFiles/DD.txt").readScalar<int>("simulationType",true))
        /* init */,useDislocations(TextFileParser("inputFiles/DD.txt").readScalar<int>("useDislocations",true))
        /* init */,useCracks(TextFileParser("inputFiles/DD.txt").readScalar<int>("useCracks",true))
        /* init */,periodicImages_x(simulationType==PERIODIC_IMAGES? TextFileParser("inputFiles/DD.txt").readScalar<int>("periodicImages_x",true) : 0)
        /* init */,periodicImages_y(simulationType==PERIODIC_IMAGES? TextFileParser("inputFiles/DD.txt").readScalar<int>("periodicImages_y",true) : 0)
        /* init */,periodicImages_z(simulationType==PERIODIC_IMAGES? TextFileParser("inputFiles/DD.txt").readScalar<int>("periodicImages_z",true) : 0)
        /* init */,Nsteps(TextFileParser("inputFiles/DD.txt").readScalar<size_t>("Nsteps",true))
        /* init */,timeIntegrationMethod(TextFileParser("inputFiles/DD.txt").readScalar<int>("timeIntegrationMethod",true))
        /* init */,externalLoadControllerName(TextFileParser("inputFiles/DD.txt").readString("externalLoadControllerName",true))
        /* init */,virtualSegmentDistance((simulationType==FINITE_FEM || simulationType==PERIODIC_FEM || simulationType==PERIODIC_IMAGES)? TextFileParser("inputFiles/DD.txt").readScalar<double>("virtualSegmentDistance",true) : 0.0)
        /* init */,runID(TextFileParser("inputFiles/DD.txt").readScalar<long int>("startAtTimeStep",true))
        /* init */,totalTime(0.0)
        /* init */,dt(10.0)
        {
            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");

            manageRestart();

            
        }
        

        
        bool DefectiveCrystalParameters::isPeriodicSimulation() const
        {
            return simulationType==PERIODIC_IMAGES || simulationType==PERIODIC_FEM;
        }
        
  }
#endif
