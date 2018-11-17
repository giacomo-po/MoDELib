/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystal_H_
#define model_DefectiveCrystal_H_

#include <model/DislocationDynamics/DislocationNetwork.h>


namespace model
{
    
    template <int _dim, short unsigned int _corder, typename InterpolationType>
    class DefectiveCrystal //: public GlidePlaneObserver<_dim>
    {
        
        
        long int runID;
        size_t Nsteps;

        
        SimplicialMesh<_dim> mesh;
        Polycrystal<_dim> poly;
        BVPsolver<_dim,2> bvpSolver;
        DislocationNetwork<_dim,_corder,InterpolationType> DN;
//        ExternalLoadControllerType extStressController;


    public:
        
        
        /**********************************************************************/
        DefectiveCrystal(int& argc, char* argv[]) :
        /* init */ runID(TextFileParser("inputFiles/DD.txt").readScalar<int>("startAtTimeStep",true))
        /* init */,Nsteps(TextFileParser("inputFiles/DD.txt").readScalar<size_t>("Nsteps",true))
        /* init */,mesh(TextFileParser("inputFiles/DD.txt").readScalar<int>("meshID",true))
        /* init */,poly("./inputFiles/polycrystal.txt",mesh)
        /* init */,bvpSolver(mesh)
        /* init */,DN(argc,argv,mesh,poly,bvpSolver,runID)
        {
            assert(mesh.simplices().size() && "MESH IS EMPTY.");
            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");

        }
        
        /**********************************************************************/
        void singleStep()
        {
            model::cout<<blueBoldColor<< "runID="<<runID<<" (of "<<Nsteps<<")"
//            /*                    */<< ", time="<<totalTime
//            /*                    */<< ": nodes="<<this->nodes().size()
//            /*                    */<< ", segments="<<this->links().size()
//            /*                    */<< ", loopSegments="<<this->loopLinks().size()
//            /*                    */<< ", loops="<<this->loops().size()
//            /*                    */<< ", components="<<this->components().size()
            /*                    */<< defaultColor<<std::endl;

            
            DN.singleStep(runID);
            ++runID;
        }
        
        /**********************************************************************/
        void runSteps()
        {/*! Runs Nsteps simulation steps
          */
            const auto t0= std::chrono::system_clock::now();
            while (runID<Nsteps)
            {
                model::cout<<std::endl; // leave a blank line
                singleStep();
            }
//            updateQuadraturePoints(); // necessary if quadrature data are necessary in main
            model::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<Nsteps<< " simulation steps completed in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" [sec]"<<defaultColor<<std::endl;
        }
        
        
    };
}
#endif
