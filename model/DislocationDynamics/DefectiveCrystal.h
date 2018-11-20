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

#include <iostream>
#include <vector>
#include <memory>

//#ifndef ExternalLoadControllerFile
//#define ExternalLoadControllerFile <model/DislocationDynamics/ExternalLoadControllers/DummyExternalLoadController.h>
//#endif
//#include ExternalLoadControllerFile

#include <model/DislocationDynamics/DefectiveCrystalParameters.h>
#include <model/DislocationDynamics/DislocationNetwork.h>
#include <model/DislocationDynamics/ExternalLoadControllers/UniformExternalLoadController.h>


namespace model
{
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    class DefectiveCrystal //: public GlidePlaneObserver<dim>
    {
        
        typedef DislocationNetwork<dim,corder,InterpolationType> DislocationNetworkType;
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef BVPsolver<dim,2> BVPsolverType;
        DefectiveCrystalParameters simulationParameters;
        
//        long int runID;
        


        const SimplicialMesh<dim> mesh;
        const std::vector<VectorDim,Eigen::aligned_allocator<VectorDim>> periodicShifts;
        const Polycrystal<dim> poly;
        DislocationNetworkType DN;
        const std::unique_ptr<BVPsolverType> bvpSolver;
        const std::unique_ptr<ExternalLoadControllerBase<dim>> externalLoadController;
        std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> straightSegmentsDeq;


        /**********************************************************************/
        static std::unique_ptr<ExternalLoadControllerBase<dim>> getExternalLoadController(const DefectiveCrystalParameters& params,
                                                                                          const DislocationNetworkType& dn,
                                                                                          const long int& rID)
        {
            // Set up periodic shifts
            if(params.externalLoadControllerName=="UniformExternalLoadController")
            {
                return std::unique_ptr<ExternalLoadControllerBase<dim>>(new UniformExternalLoadController<DislocationNetworkType>(dn,rID));
            }
            else if(params.externalLoadControllerName=="None")
            {
                return std::unique_ptr<ExternalLoadControllerBase<dim>>(nullptr);
            }
            else
            {
                model::cout<<"Unknown externalLoadController name "<<params.externalLoadControllerName<<". Use 'None' for no-controller. EXITING."<<std::endl;
                exit(EXIT_FAILURE);
            }
            
        }
        
        /**********************************************************************/
        static std::vector<VectorDim,Eigen::aligned_allocator<VectorDim>> getPeriodicShifts(const SimplicialMesh<dim>& m,
                                                                                            const DefectiveCrystalParameters& params)
        {
            // Set up periodic shifts
            std::vector<VectorDim,Eigen::aligned_allocator<VectorDim>> temp;
            const VectorDim meshDimensions(m.xMax()-m.xMin());
            model::cout<<"meshDimensions="<<meshDimensions.transpose()<<std::endl;
            for(int i=-params.periodicImages_x;i<=params.periodicImages_x;++i)
            {
                for(int j=-params.periodicImages_y;j<=params.periodicImages_y;++j)
                {
                    for(int k=-params.periodicImages_z;k<=params.periodicImages_z;++k)
                    {
                        const Eigen::Array<int,dim,1> cellID((Eigen::Array<int,dim,1>()<<i,j,k).finished());
                        temp.push_back((meshDimensions.array()*cellID.template cast<double>()).matrix());
                    }
                }
            }
            
            model::cout<<"periodic shift vectors:"<<std::endl;
            for(const auto& shift : temp)
            {
                model::cout<<shift.transpose()<<std::endl;
                
            }
            
            return temp;

        }
        
        /**********************************************************************/
        void updateStressStraightSegments()
        {
            
            const auto t0= std::chrono::system_clock::now();
            model::cout<<"		Collecting StressStraight objects: MUST CHANGE THIS BY SHITFING THE FIELD POINT RATHER THAN POPOLATING THE IMAGES"<<std::flush;
            
            straightSegmentsDeq.clear();
            //                std::deque<StressStraight<dim>,Eigen::aligned_allocator<StressStraight<dim>>> straightSegmentsDeq;
            size_t currentSize=0;
            //            if(computeDDinteractions)
            //            {
            for(const auto& link : DN.networkLinks())
            {
                link.second->addToStressStraight(straightSegmentsDeq);
            }
            
            currentSize=straightSegmentsDeq.size();
            
            for(const auto& shift : periodicShifts)
            {
                if(shift.squaredNorm()!=0.0)
                {
                    for (size_t c=0;c<currentSize;++c)
                    {
                        straightSegmentsDeq.emplace_back(straightSegmentsDeq[c].P0+shift,straightSegmentsDeq[c].P1+shift,straightSegmentsDeq[c].b);
                    }
                }
            }
            
            model::cout<< straightSegmentsDeq.size()<<" straight segments ("<<currentSize<<"+"<<straightSegmentsDeq.size()-currentSize<<" images)"<<std::flush;
            model::cout<<magentaColor<<std::setprecision(3)<<std::scientific<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]."<<defaultColor<<std::endl;
            
        }
        
        /**********************************************************************/
        void updateLoadControllers(const long int& runID)
        {/*! Updates bvpSolver using the stress and displacement fields of the
          *  current DD configuration.
          */
            const int quadraturePerTriangle=37;
            if(bvpSolver)
            {
                if (!(runID%bvpSolver->stepsBetweenBVPupdates))
                {// enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
                    model::cout<<"		Updating elastic bvp... "<<std::endl;
                    bvpSolver->template assembleAndSolve<DislocationNetworkType,quadraturePerTriangle>(DN);
                }
            }
            if (externalLoadController)
            {
                externalLoadController->update(DN,runID);
            }
        }
        
    public:
        
        
        /**********************************************************************/
        DefectiveCrystal(int& argc, char* argv[]) :
        /* init */ simulationParameters(argc,argv)
//        /* init */,runID(TextFileParser("inputFiles/DD.txt").readScalar<int>("startAtTimeStep",true))
        /* init */,mesh(TextFileParser("inputFiles/DD.txt").readScalar<int>("meshID",true))
        /* init */,periodicShifts(getPeriodicShifts(mesh,simulationParameters))
        /* init */,poly("./inputFiles/polycrystal.txt",mesh)
        /* init */,DN(argc,argv,simulationParameters,mesh,poly,bvpSolver,externalLoadController,periodicShifts,simulationParameters.runID)
        /* init */,bvpSolver(simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM? new BVPsolverType(mesh) : nullptr)
        /* init */,externalLoadController(getExternalLoadController(simulationParameters,DN,simulationParameters.runID))
        {
            assert(mesh.simplices().size() && "MESH IS EMPTY.");
            
            
            switch (simulationParameters.simulationType)
            {// Initilization based on type of simulation


                case DefectiveCrystalParameters::FINITE_NO_FEM:
                {
//                    externalLoadController->init(DN,runID);  // have to initialize it after mesh!
                    break;
                }
                    
                case DefectiveCrystalParameters::FINITE_FEM:
                {
                    bvpSolver->init(DN);
                    break;
                }
                   
                case DefectiveCrystalParameters::PERIODIC:
                {
                    assert(poly.grains().size()==1 && "ONLY SINGLE-CRYSTAL PERIODIC SIMULATIONS SUPPORTED.");
                    model::cout<<"Checking that mesh and lattice are commensurate"<<std::endl;
                    const VectorDim meshSize(mesh.xMax()-mesh.xMin());
                    for(int d=0;d<dim;++d)
                    {
                        VectorDim v(VectorDim::Zero());
                        v(d)=1*meshSize(d);
                        LatticeVector<dim> lv(v,poly.grains().begin()->second);
                    }
                    break;
                }
                    
                default:
                {
                    model::cout<<"simulationType MUST BE 0,1, or 2. EXITING."<<std::endl;
                    exit(EXIT_FAILURE);
                    break;
                }
            }
            

            
        }
        
        /**********************************************************************/
        void singleStep()
        {
            model::cout<<blueBoldColor<< "runID="<<simulationParameters.runID<<" (of "<<simulationParameters.Nsteps<<")"
            /*                    */<< ", time="<<simulationParameters.totalTime
            /*                    */<< ": nodes="<<DN.nodes().size()
            /*                    */<< ", segments="<<DN.links().size()
            /*                    */<< ", loopSegments="<<DN.loopLinks().size()
            /*                    */<< ", loops="<<DN.loops().size()
            /*                    */<< ", components="<<DN.components().size()
            /*                    */<< defaultColor<<std::endl;

            
            DN.updateGeometry(simulationParameters.dt);
            
            updateStressStraightSegments();

            updateLoadControllers(simulationParameters.runID);
            
            DN.assembleAndSolve(simulationParameters.runID,straightSegmentsDeq);
            simulationParameters.dt=DN.get_dt(); // TO DO: MAKE THIS std::min between DN and CrackSystem
            simulationParameters.totalTime+=simulationParameters.dt;
//            DN.updatePlasticDistortionRateFromVelocities();
            
            // output
            DN.io().output(simulationParameters.runID);

            // move
            DN.move(simulationParameters.dt);

            // menage discrete topological events
            DN.singleStepDiscreteEvents(simulationParameters.runID);
            
            ++simulationParameters.runID;
        }
        
        /**********************************************************************/
        void runSteps()
        {/*! Runs Nsteps simulation steps
          */
            const auto t0= std::chrono::system_clock::now();
            while (simulationParameters.runID<simulationParameters.Nsteps)
            {
                model::cout<<std::endl; // leave a blank line
                singleStep();
            }
//            updateQuadraturePoints(); // necessary if quadrature data are necessary in main
            model::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<simulationParameters.Nsteps<< " simulation steps completed in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" [sec]"<<defaultColor<<std::endl;
        }
        
        
    };
}
#endif
