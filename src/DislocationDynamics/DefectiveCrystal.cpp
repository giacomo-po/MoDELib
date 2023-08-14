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
                     
            if(params.simulationType==DDtraitsIO::FINITE_FEM)
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
                }
            }
            
            
        }
        
//        /**********************************************************************/
//        template <int _dim, short unsigned int corder>
//        std::vector<typename DefectiveCrystal<_dim,corder>::VectorDim> DefectiveCrystal<_dim,corder>::getPeriodicShifts(const SimplicialMesh<DefectiveCrystal<_dim,corder>::dim>& m,
//                                                        const DefectiveCrystalParameters& params)
//        {
//            // Set up periodic shifts
//            std::vector<VectorDim> temp;
//            temp.push_back(VectorDim::Zero());
//            if(params.simulationType==DDtraitsIO::PERIODIC_IMAGES)
//            {
//                const auto shiftVectors(m.periodicBasis());
//                std::cout<<"Box periodicity vectors ("<<shiftVectors.size()<<"):"<<std::endl;
//                for(const auto& shift : shiftVectors)
//                {
//                    std::cout<<shift.transpose()<<std::endl;
//                }
//                
//                if(shiftVectors.size()!=params.periodicImageSize.size())
//                {
//                    std::cout<<"shiftVectors.size()="<<shiftVectors.size()<<std::endl;
//                    std::cout<<"periodicImageSize.size()="<<params.periodicImageSize.size()<<std::endl;
//                    throw std::runtime_error("shiftVectors.size() must equal periodicImageSize.size()");
//                }
//                
//                
//                for(size_t k=0;k<shiftVectors.size();++k)
//                {
//                    const auto shiftVector(shiftVectors[k]);
//                    const int imageSize(std::abs(params.periodicImageSize[k]));
//                    std::vector<VectorDim> newTemp;
//                    for(const auto& v:temp)
//                    {// grab existing shift vectors
//                        for(int i=-imageSize;i<=imageSize;++i)
//                        {// add current shift times corresponding image size
//                            newTemp.push_back(v+i*shiftVector);
//                        }
//                    }
//                    temp.swap(newTemp);
//                }
//            }
//
//
//
//            
//            std::cout<<"Image shift vectors ("<<temp.size()<<"):"<<std::endl;
//            for(const auto& shift : temp)
//            {
//                std::cout<<shift.transpose()<<std::endl;
//            }
//            
//            return temp;
//            
//        }
        
        
        
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
                    std::cout<<"Updating bvpSolver ... "<<std::endl;
                    bvpSolver->template assembleAndSolve<DislocationNetworkType,quadraturePerTriangle>(*DN, isClimbStep);
                }
            }
            if (externalLoadController)
            {
                std::cout<<"Updating externalLoadController... "<<std::endl;
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
//        /* init */,periodicShifts(getPeriodicShifts(mesh,simulationParameters))
        /* init */,periodicShifts(mesh.periodicShifts(simulationParameters.periodicImageSize))
        /* init */,poly(simulationParameters.traitsIO.polyFile,mesh)
        /* init */,DN(simulationParameters.useDislocations? new DislocationNetworkType(simulationParameters,mesh,poly,bvpSolver,externalLoadController,periodicShifts,simulationParameters.runID) : nullptr)
        /* init */,CS(simulationParameters.useCracks? new CrackSystemType() : nullptr)
        //        /* init */,DN(argc,argv,simulationParameters,mesh,poly,bvpSolver,externalLoadController,periodicShifts,simulationParameters.runID)
        /* init */,bvpSolver(simulationParameters.simulationType==DDtraitsIO::FINITE_FEM? new BVPsolverType(mesh,*DN) : nullptr)
        /* init */,externalLoadController(getExternalLoadController(simulationParameters,*this,simulationParameters.runID))
        {
            
            if(!mesh.simplices().size())
            {
                throw std::runtime_error("Mesh is empty");
            }
          
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
            }
            std::cout<< defaultColor<<std::endl;

            if(DN)
            {
                DislocationNode<dim,corder>::totalCappedNodes=0;
                DN->updateGeometry();
                updateLoadControllers(simulationParameters.runID, false);
                const double maxVelocity(getMaxVelocity());
                DN->assembleGlide(simulationParameters.runID, maxVelocity);
                DN->storeSingleGlideStepDiscreteEvents(simulationParameters.runID);
                DN->solveGlide(simulationParameters.runID);
                simulationParameters.dt=DN->timeIntegrator.getGlideTimeIncrement(*DN); // TO DO: MAKE THIS std::min between DN and CrackSystem
                DN->updateRates();
                DN->io().output(simulationParameters.runID);
                DN->moveGlide(simulationParameters.dt);
                DN->executeSingleGlideStepDiscreteEvents(simulationParameters.runID);
                if (DN->capMaxVelocity)
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
            if(DN)
            {
                DN->updateGeometry();
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

template <int _dim, short unsigned int corder>
typename DefectiveCrystal<_dim,corder>::MatrixDim DefectiveCrystal<_dim,corder>::plasticStrain() const
{/*!\param[in] P position vector
  * \returns The stress field in the DefectiveCrystal at P
  * Note:
  */
    MatrixDim temp(plasticDistortion());
    return 0.5*(temp+temp.transpose());
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
