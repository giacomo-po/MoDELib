/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2017 by Yinan Cui <cuiyinan@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DefectiveCrystal_cpp_
#define model_DefectiveCrystal_cpp_

#include <DefectiveCrystal.h>

namespace model
{
        
        
        /**********************************************************************/
        template <int _dim, short unsigned int corder>
        std::unique_ptr<ExternalLoadControllerBase<DefectiveCrystal<_dim,corder>::dim>> DefectiveCrystal<_dim,corder>::getExternalLoadController(const DefectiveCrystalParameters& params,
                                                                                          const DefectiveCrystalType& dc,
                                                                                          const long int& rID)
        {
                     
            if(params.simulationType==DefectiveCrystalParameters::FINITE_FEM)
            {
                return std::unique_ptr<ExternalLoadControllerBase<DefectiveCrystal<_dim,corder>::dim>>(nullptr);
            }
            else
            {
                if(params.externalLoadControllerName=="UniformExternalLoadController")
                {
                    return std::unique_ptr<ExternalLoadControllerBase<DefectiveCrystal<_dim,corder>::dim>>(new UniformExternalLoadController<DefectiveCrystalType>(dc,rID));
                }
//                else if(params.externalLoadControllerName=="None")
//                {
//                    return std::unique_ptr<ExternalLoadControllerBase<DefectiveCrystal<_dim,corder>::dim>>(nullptr);
//                }
                else
                {
                    std::cout<<"Unknown externalLoadController name "<<params.externalLoadControllerName<<"No controller applied."<<std::endl;
                    return std::unique_ptr<ExternalLoadControllerBase<DefectiveCrystal<_dim,corder>::dim>>(nullptr);
//
//                    exit(EXIT_FAILURE);
//                    return nullptr;
                }
            }
            
            
        }
        
        /**********************************************************************/
        template <int _dim, short unsigned int corder>
        std::vector<typename DefectiveCrystal<_dim,corder>::VectorDim> DefectiveCrystal<_dim,corder>::getPeriodicShifts(const SimplicialMesh<DefectiveCrystal<_dim,corder>::dim>& m,
                                                        const DefectiveCrystalParameters& params)
        {
            // Set up periodic shifts
            std::vector<VectorDim> temp;
            if(params.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES)
            {
                const VectorDim meshDimensions(m.xMax()-m.xMin());
                std::cout<<"meshDimensions="<<meshDimensions.transpose()<<std::endl;
                for(int i=-params.periodicImages_x;i<=params.periodicImages_x;++i)
                {
                    for(int j=-params.periodicImages_y;j<=params.periodicImages_y;++j)
                    {
                        for(int k=-params.periodicImages_z;k<=params.periodicImages_z;++k)
                        {
                            const Eigen::Array<int,DefectiveCrystal<_dim,corder>::dim,1> cellID((Eigen::Array<int,DefectiveCrystal<_dim,corder>::dim,1>()<<i,j,k).finished());
                            temp.push_back((meshDimensions.array()*cellID.template cast<double>()).matrix());
                        }
                    }
                }
            }
            else
            {
                temp.push_back(VectorDim::Zero());
            }
            
            std::cout<<"periodic shift vectors:"<<std::endl;
            for(const auto& shift : temp)
            {
                std::cout<<shift.transpose()<<std::endl;
                
            }
            
            return temp;
            
        }
        
        
        
        /**********************************************************************/
        template <int _dim, short unsigned int corder>
        void DefectiveCrystal<_dim,corder>::updateLoadControllers(const long int& runID, const bool& isClimbStep)
        {/*! Updates bvpSolver using the stress and displacement fields of the
          *  current DD configuration.
          */
            const int quadraturePerTriangle=37;
            if(bvpSolver)
            {
                if (!(runID%bvpSolver->stepsBetweenBVPupdates))
                {// enter the if statement if use_bvp!=0 and runID is a multiple of use_bvp
                    std::cout<<"		Updating elastic bvp... "<<std::endl;
                    bvpSolver->template assembleAndSolve<DislocationNetworkType,quadraturePerTriangle>(*DN, isClimbStep);
                }
            }
            if (externalLoadController)
            {
                externalLoadController->update(runID);
            }
        }
        
    
        
        /**********************************************************************/
        template <int _dim, short unsigned int corder>
        DefectiveCrystal<_dim,corder>::DefectiveCrystal(const std::string& folderName) :
        /* init */ simulationParameters(folderName)
        /* init */,periodicFaceIDs(TextFileParser(simulationParameters.traitsIO.polyFile).template readSet<int>("periodicFaceIDs",true))
        /* init */,mesh(simulationParameters.traitsIO.meshFile,
                        TextFileParser(simulationParameters.traitsIO.polyFile).readMatrix<double>("A",3,3,true),
                        TextFileParser(simulationParameters.traitsIO.polyFile).readMatrix<double>("x0",1,3,true).transpose(),periodicFaceIDs)
        /* init */,periodicShifts(getPeriodicShifts(mesh,simulationParameters))
        /* init */,poly(simulationParameters.traitsIO.polyFile,mesh)
        /* init */,DN(simulationParameters.useDislocations? new DislocationNetworkType(simulationParameters,mesh,poly,bvpSolver,externalLoadController,periodicShifts,simulationParameters.runID) : nullptr)
        /* init */,CS(simulationParameters.useCracks? new CrackSystemType() : nullptr)
        //        /* init */,DN(argc,argv,simulationParameters,mesh,poly,bvpSolver,externalLoadController,periodicShifts,simulationParameters.runID)
        /* init */,bvpSolver(simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM? new BVPsolverType(mesh,*DN) : nullptr)
        /* init */,externalLoadController(getExternalLoadController(simulationParameters,*this,simulationParameters.runID))
        {
            assert(mesh.simplices().size() && "MESH IS EMPTY.");
            

            //Commensurate check happens while creation of periodic glideplane factory
            
//             if(   simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_IMAGES
//                || simulationParameters.simulationType==DefectiveCrystalParameters::PERIODIC_FEM)
//             {
//                 assert(poly.grains().size()==1 && "ONLY SINGLE-CRYSTAL PERIODIC SIMULATIONS SUPPORTED.");
                
//                 for(const auto& rIter : mesh.regions())
//                 {
//                     for(const auto& pair : rIter.second->parallelFaces())
//                     {
//                         std::cout<<"Checking if parallel faces "<<pair.first<<"<->"<<pair.second<<" are commensurate"<<std::endl;
//                         const PlanarMeshFace<DefectiveCrystal<_dim,corder>::dim>& face1(*rIter.second->faces().at(pair.first));
//                         const PlanarMeshFace<DefectiveCrystal<_dim,corder>::dim>& face2(*rIter.second->faces().at(pair.second));
//                         const VectorDim cc(face1.center()-face2.center());
//                         const VectorDim ccc(cc.dot(face1.outNormal())*face1.outNormal());
                        
//                         const LatticeDirection<DefectiveCrystal<_dim,corder>::dim> ld(poly.grains().begin()->second.latticeDirection(face1.outNormal()));
//                         const double normRatio(ccc.norm()/ld.cartesian().norm());
//                         if(std::fabs(std::round(normRatio)-normRatio)>FLT_EPSILON)
//                         {
// //                            std::cout<<"Face outNormal="<<std::setprecision(15)<<std::scientific<<face1.outNormal().transpose()<<std::endl;
//                             std::cout<<"Mesh in direction "<< std::setprecision(15)<<std::scientific<<ld.cartesian().normalized().transpose()<<" is not commensurate for periodicity"<<std::endl;
//                             std::cout<<"Mesh size in that direction must be a multiple of "<< std::setprecision(15)<<std::scientific<<ld.cartesian().norm()<<std::endl;
//                             std::cout<<"Size detected="<< std::setprecision(15)<<std::scientific<<ccc.norm()<<std::endl;
//                             std::cout<<"Closest commensurate size="<< std::setprecision(15)<<std::scientific<<std::round(normRatio)*ld.cartesian().norm()<<std::endl;
//                             assert(false && "MESH NOT COMMENSURATE");
//                         }
// //                        LatticeVector<DefectiveCrystal<_dim,corder>::dim> lv(ccc,poly.grains().begin()->second);
//                     }
//                 }
//             }
            
          
        }

        template <int _dim, short unsigned int corder>
        double DefectiveCrystal<_dim,corder>::getMaxVelocity() const
        {
            double vmax = 0.0;

            for (const auto &nodeIter : DN->networkNodes())
            {
                    const double vNorm(nodeIter.second.lock()->get_V().norm());
                    if (vNorm > vmax)
                    {
                        vmax = vNorm;
                    }
            }
            return vmax;
        }
        
        /**********************************************************************/
        template <int _dim, short unsigned int corder>
        void DefectiveCrystal<_dim,corder>::singleGlideStep()
        {
            std::cout<<blueBoldColor<< "runID="<<simulationParameters.runID<<" (of "<<simulationParameters.Nsteps<<")"
            /*                    */<< ", time="<<simulationParameters.totalTime;
            if(DN)
            {
                std::cout<< ": networkNodes="<<DN->networkNodes().size()
                /*                    */<< ", networkSegments="<<DN->networkLinks().size()
                /*                    */<< ", loopNodes="<<DN->loopNodes().size()
                /*                    */<< ", loopSegments="<<DN->loopLinks().size()
                /*                    */<< ", loops="<<DN->loops().size();
//                /*                    */<< ", components="<<DN->components().size();
            }
            std::cout<< defaultColor<<std::endl;

            if(DN)
            {
                // for (const auto& netNode : DN->networkNodes())
                // {
                //     std::cout<<std::scientific<<std::setprecision(7)<<"networkNode sID"<<netNode.second.lock()->sID<<"-->"<<netNode.second.lock()->get_P().transpose()<<
                //     "-->"<<netNode.second.lock()->get_V().transpose()<<"-->"<<netNode.second.lock()->velocityReduction()<<std::endl;
                // }
                DislocationNode<dim,corder>::totalCappedNodes=0;
                DN->updateGeometry();
                // DN->io().output(simulationParameters.runID);

                updateLoadControllers(simulationParameters.runID, false);
                const double maxVelocity(getMaxVelocity());
                DN->assembleAndSolveGlide(simulationParameters.runID, maxVelocity);
                simulationParameters.dt=DN->timeIntegrator.getGlideTimeIncrement(*DN); // TO DO: MAKE THIS std::min between DN and CrackSystem
                // output
                DN->io().output(simulationParameters.runID);

                DislocationCrossSlip<DislocationNetworkType> crossSlip(*DN);
                // move
//                DN->dummyMove(simulationParameters.runID);
                DN->moveGlide(simulationParameters.dt);
                crossSlip.execute();
                // DN->io().output(simulationParameters.runID);


                // manage discrete topological events
                DN->singleGlideStepDiscreteEvents(simulationParameters.runID);
                // DN->io().output(simulationParameters.runID);
                if (true)
                {
                    std::cout<<redBoldColor<<"( "<<(DislocationNode<dim,corder>::totalCappedNodes)<<" total nodes capped "<<defaultColor<<std::endl;
                    std::cout<<redBoldColor<<", "<<(double(DislocationNode<dim,corder>::totalCappedNodes)/double(DN->networkNodes().size()))<<" fraction of nodes capped "
                    <<defaultColor<<" )"<<std::endl;
                }
            }
            simulationParameters.totalTime+=simulationParameters.dt;
            ++simulationParameters.runID;
        }
//
//        /**********************************************************************/
        template <int _dim, short unsigned int corder>
        void DefectiveCrystal<_dim,corder>::runGlideSteps()
        {/*! Runs a number of simulation time steps defined by simulationParameters.Nsteps
          */
            const auto t0= std::chrono::system_clock::now();
            while (simulationParameters.runID<simulationParameters.Nsteps)
            {
                std::cout<<std::endl; // leave a blank line
                singleGlideStep();
            }
            std::cout<<greenBoldColor<<std::setprecision(3)<<std::scientific<<simulationParameters.Nsteps<< " simulation steps completed in "<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" [sec]"<<defaultColor<<std::endl;
        }

        /**********************************************************************/
        template <int _dim, short unsigned int corder>
        typename DefectiveCrystal<_dim,corder>::MatrixDim DefectiveCrystal<_dim,corder>::plasticDistortion() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(MatrixDim::Zero());
            if(DN)
            {
                temp+=DN->plasticDistortion();
            }
            if(CS)
            {
                temp+=CS->plasticDistortion();
            }
            return temp;
        }

        /**********************************************************************/
        template <int _dim, short unsigned int corder>
        typename DefectiveCrystal<_dim,corder>::MatrixDim DefectiveCrystal<_dim,corder>::plasticDistortionRate() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(MatrixDim::Zero());
            if(DN)
            {
                temp+=DN->plasticDistortionRate();
            }
            if(CS)
            {
                temp+=CS->plasticDistortionRate();
            }
            return temp;
        }

        /**********************************************************************/
        template <int _dim, short unsigned int corder>
        typename DefectiveCrystal<_dim,corder>::MatrixDim DefectiveCrystal<_dim,corder>::plasticStrainRate() const
        {/*!\param[in] P position vector
          * \returns The stress field in the DefectiveCrystal at P
          * Note:
          */
            MatrixDim temp(plasticDistortionRate());
            return 0.5*(temp+temp.transpose());
        }
//
template class DefectiveCrystal <3,0>;
    
}
#endif
