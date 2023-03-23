/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONNETWORKIO_H_
#define model_DISLOCATIONNETWORKIO_H_

#include <chrono>
#include <string>
#include <IDreader.h>


//#include <UniqueOutputFile.h>
//#include <SequentialOutputFile.h>
#include <TerminalColors.h>

#include <PeriodicGlidePlane.h>
#include <LatticeModule.h>
#include <LatticeModule.h>
//#include <BoundaryDisplacementPoint.h>
//#include <DislocationNodeIO.h>
//#include <DDtimeIntegrator.h>
//#include <DislocationNodeContraction.h>
//#include <GrainBoundaryTransmission.h>
#include <DislocationLinkingNumber.h>
//#include <EVLio.h>
#include <DDconfigIO.h>
#include <DDauxIO.h>

#include <SphericalInclusion.h>
#include <PolyhedronInclusion.h>
//#include <DisplacementPoint.h>
#include <FEMnodeEvaluation.h>



namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename DislocationNetworkType>
    struct DislocationNetworkIO
    {
        
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;
        
        const DislocationNetworkType& DN;
        std::ofstream f_file;
        std::ofstream F_labels;

        /**********************************************************************/
        DislocationNetworkIO(const DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        /* init */,f_file(DN.simulationParameters.traitsIO.fFile,std::ios_base::app)
        /* init */,F_labels(DN.simulationParameters.traitsIO.flabFile,std::ios_base::app)
        {
            if (!f_file.is_open())
              {
                  throw std::runtime_error("Cannot open file "+DN.simulationParameters.traitsIO.fFile);
              }
            if (!F_labels.is_open())
              {
                  throw std::runtime_error("Cannot open file "+DN.simulationParameters.traitsIO.flabFile);
              }

        }

        /**********************************************************************/
        void output(const size_t& runID)
        {
            if (!(runID%DN.outputFrequency))
            {
                const auto t0=std::chrono::system_clock::now();
#ifdef _MODEL_DD_MPI_
                if(ModelMPIbase::mpiRank()==0)
                {
                    outputFiles(runID);
                }
#else
                outputFiles(runID);
#endif
                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            
        }
        
        /**********************************************************************/
        DDconfigIO<dim> configIO() const
        {
            DDconfigIO<dim> temp(DN.simulationParameters.traitsIO.evlFolder);
            for(const auto& loop : DN.loops())
            {
                temp.loops().emplace_back(*loop.second.lock());
            }
            
            // Store LoopNodes
            for(const auto& node : DN.loopNodes())
            {
                temp.loopNodes().emplace_back(*node.second.lock());
            }
            
            // Store LoopLinks
            for(const auto& link : DN.loopLinks())
            {
                temp.loopLinks().emplace_back(link.second);
            }
            
            // Store NetworkNodes
            for(const auto& node : DN.networkNodes())
            {
                temp.nodes().emplace_back(*node.second.lock());
            }
            
            for(const auto& node : DN.polyhedronInclusionNodes())
            {
                temp.polyhedronInclusionNodes().emplace_back(node.second);
            }
            
            // Store Eshelby Inclusions
            for(const auto& ei : DN.eshelbyInclusions())
            {
                
                auto* sphericalDerived = dynamic_cast<SphericalInclusion<dim>*>(ei.second.get());
                if (sphericalDerived)
                {
                    temp.sphericalInclusions().emplace_back(*sphericalDerived);
                }
                
                auto* polyhedronDerived = dynamic_cast<PolyhedronInclusion<dim>*>(ei.second.get());
                if (polyhedronDerived)
                {
                    temp.polyhedronInclusions().emplace_back(*polyhedronDerived);
                    for(const auto& face : polyhedronDerived->faces)
                    {
                        for(size_t k=0;k<face.second.size();++k)
                        {
                            const size_t k1(k<face.second.size()-1? k+1 : 0);
                            temp.polyhedronInclusionEdges().emplace_back(polyhedronDerived->sID,face.first,face.second[k],face.second[k1]);
                        }
                    }
                }
            }
            

            
            
            return temp;
        }
        
        /**********************************************************************/
        DDauxIO<dim> auxIO() const
        {
            DDauxIO<dim> temp(DN.simulationParameters.traitsIO.auxFolder);
            
            if(DN.outputMeshDisplacement)
            {
                std::vector<FEMnodeEvaluation<typename DislocationNetworkType::ElementType,dim,1>> fieldPoints; // the container of field points
                fieldPoints.reserve(DN.mesh.template observer<0>().size());
                for (const auto& sIter : DN.mesh.template observer<0>())
                {
                    if(sIter.second->isBoundarySimplex())
                    {
                        fieldPoints.emplace_back(sIter.second->xID(0),sIter.second->P0);
                    }
                }
                temp.meshNodes().reserve(fieldPoints.size());
                
                DN.displacement(fieldPoints);
                
                for(auto& node : fieldPoints)
                {// add FEM solution and output
                    if(DN.simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM)
                    {
                        const size_t femID=DN.bvpSolver->finiteElement().mesh2femIDmap().at(node.pointID)->gID;
                        node+=DN.bvpSolver->displacement().dofs(femID);
                    }
                    temp.meshNodes().emplace_back(node);
                }
            }
            
            if (DN.outputQuadraturePoints)
            {
                for (const auto& link : DN.networkLinks())
                {
                    for(const auto& qPoint : link.second.lock()->quadraturePoints())
                    {
                        temp.quadraturePoints().push_back(qPoint);
                    }
                }
            }
            
            if(DN.outputGlidePlanes)
            {
                for(const auto& pair : DN.glidePlaneFactory.glidePlanes())
                {
                    const auto glidePlane(pair.second.lock());
                    if(glidePlane)
                    {
                        temp.glidePlanes().emplace_back(glidePlane->key);
//                        for(const auto& seg : glidePlane->meshIntersections)
//                        {
//                            glidePlanesBoundaries().emplace_back(glidePlane->sID,*seg);
//                        }
                    }
                }
                
//                temp.setGlidePlaneBoundaries(DN.glidePlaneFactory);
            }
            
            return temp;
        }
        
        /* outputTXT **********************************************************/
        void outputFiles(const size_t& runID)
        {/*! Outputs DislocationNetwork data to the following files (x is the runID):
          */
            configIO().write(runID,DN.outputBinary);
            auxIO().write(runID,DN.outputBinary);

 
            
            if(DN.outputSegmentPairDistances)
            {
                const std::string outFileName(DN.simulationParameters.traitsIO.simulationFolder+"H/H_"+std::to_string(runID));
                std::ofstream h_file(outFileName);
                if (h_file.is_open())
                {
                    const auto t0=std::chrono::system_clock::now();
                    std::cout<<"writing to "<<outFileName<<std::flush;
                    
                    for(auto linkIter1=DN.networkLinks().begin();linkIter1!=DN.networkLinks().end();++linkIter1)
                    {
                        for(auto linkIter2=linkIter1;linkIter2!=DN.networkLinks().end();++linkIter2)
                        {
                            if(   !linkIter1->second.lock()->isBoundarySegment()
                               && !linkIter2->second.lock()->isBoundarySegment()
                               && !linkIter1->second.lock()->hasZeroBurgers()
                               && !linkIter2->second.lock()->hasZeroBurgers())
                            {
                                SegmentSegmentDistance<dim> ssi(linkIter1->second.lock()->source->get_P(),
                                                                linkIter1->second.lock()->sink->get_P(),
                                                                linkIter2->second.lock()->source->get_P(),
                                                                linkIter2->second.lock()->sink->get_P());
                                
                                h_file<<sqrt(ssi.D1)<<" "<<sqrt(ssi.D2)<<" "<<ssi.dMin<<"\n";
                            }
                        }
                    }
                    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                }
                else
                {
                    throw std::runtime_error("Cannot open "+outFileName);
                }

                
            }
            
//            if (DN.bvpSolver && DN.outputFEMsolution )
//            {
//                if(!(runID%DN.bvpSolver->stepsBetweenBVPupdates))
//                {// Output displacement and stress on external mesh faces
////                    const auto t0=std::chrono::system_clock::now();
////                    model::SequentialOutputFile<'U',1>::set_count(runID); // Vertices_file;
////                    model::SequentialOutputFile<'U',1>::set_increment(DN.outputFrequency); // Vertices_file;
////                    model::SequentialOutputFile<'U',true> u_file;
////                    std::cout<<"		writing to U/U_"<<u_file.sID<<".txt"<<std::flush;
////                    u_file<<DN.bvpSolver->displacement().onBoundary();
////                    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
////
////                    const auto t1=std::chrono::system_clock::now();
////                    model::SequentialOutputFile<'S',1>::set_count(runID); // Vertices_file;
////                    model::SequentialOutputFile<'S',1>::set_increment(DN.outputFrequency); // Vertices_file;
////                    model::SequentialOutputFile<'S',true> s_file;
////                    std::cout<<"		writing to S/S_"<<s_file.sID<<".txt"<<std::flush;
////                    s_file<<DN.bvpSolver->stress().onBoundary();
////                    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;
//                }
//            }
            

            
            if (DN.outputLinkingNumbers)
            {
                DislocationLinkingNumber<DislocationNetworkType> LN(DN);
                const std::string outFileName(DN.simulationParameters.traitsIO.simulationFolder+"Z/Z_"+std::to_string(runID));
                std::ofstream z_file(outFileName);
                if (z_file.is_open())
                {
                    z_file<<LN;
                }
                else
                {
                    throw std::runtime_error("Cannot open "+outFileName);
                }
            }
            
            // Output to F file
            std::cout<<"Writing F/F_0.txt"<<std::flush;
            
            f_file<< runID<<" "<<std::setprecision(15)<<std::scientific<<DN.simulationParameters.totalTime<<" "<<DN.simulationParameters.dt<<" ";
            if(runID==0)
            {
                F_labels<<"runID\n";
                F_labels<<"time [b/cs]\n";
                F_labels<<"dt [b/cs]"<<std::endl;
            }
            
            
            const Eigen::Matrix<double,dim,dim>& pD(DN.plasticDistortion());
            f_file<<pD.row(0)<<" "<<pD.row(1)<<" "<<pD.row(2)<<" "<<pD.trace()<<" "<<pD.norm()<<" ";
            if(runID==0)
            {
                F_labels<<"betaP_11\n";
                F_labels<<"betaP_12\n";
                F_labels<<"betaP_13\n";
                F_labels<<"betaP_21\n";
                F_labels<<"betaP_22\n";
                F_labels<<"betaP_23\n";
                F_labels<<"betaP_31\n";
                F_labels<<"betaP_32\n";
                F_labels<<"betaP_33\n";
                F_labels<<"tr(betaP)\n";
                F_labels<<"norm(betaP)"<<std::endl;

            }
            
            if(DN.outputPlasticDistortionRate)
            {
                const Eigen::Matrix<double,dim,dim>& pDR(DN.plasticDistortionRate());
                f_file<<pDR.row(0)<<" "<<pDR.row(1)<<" "<<pDR.row(2)<<" "<<pDR.trace()<<" "<<pDR.norm()<<" ";
                if(runID==0)
                {
                    F_labels<<"dotBetaP_11 [cs/b]\n";
                    F_labels<<"dotBetaP_12 [cs/b]\n";
                    F_labels<<"dotBetaP_13 [cs/b]\n";
                    F_labels<<"dotBetaP_21 [cs/b]\n";
                    F_labels<<"dotBetaP_22 [cs/b]\n";
                    F_labels<<"dotBetaP_23 [cs/b]\n";
                    F_labels<<"dotBetaP_31 [cs/b]\n";
                    F_labels<<"dotBetaP_32 [cs/b]\n";
                    F_labels<<"dotBetaP_33 [cs/b]\n";
                    F_labels<<"tr(dotBetaP) [cs/b]\n";
                    F_labels<<"norm(dotBetaP) [cs/b]"<<std::endl;
                }
            }
            
            if(DN.outputDislocationLength)
            {
                const std::tuple<double,double,double,double> length=DN.networkLength();
                f_file<<std::get<0>(length)<<" "<<std::get<1>(length)<<" "<<std::get<2>(length)<<" "<<std::get<3>(length)<<" ";
                if(runID==0)
                {
                    F_labels<<"glissile length [b]\n";
                    F_labels<<"sessile length [b]\n";
                    F_labels<<"boundary length [b]\n";
                    F_labels<<"grain boundary length [b]"<<std::endl;
                }
            }
            
            if(DN.computeElasticEnergyPerLength)
            {
                double eE(0.0);
                for(const auto& linkIter : DN.networkLinks())
                {// Collect LoopLinks by loop IDs
                    const auto link(linkIter.second.lock());
                    for(const auto& qPoint : link->quadraturePoints())
                    {
                        eE+=qPoint.elasticEnergyPerLength*qPoint.dL;
                    }
                }
                f_file<<eE<<" ";
                if(runID==0)
                {
                    F_labels<<"elastic energy [mu b^3]\n";
                }
            }
            
            if (DN.externalLoadController)
            {
                DN.externalLoadController->output(runID,f_file,F_labels);
            }
            
            if(DN.bvpSolver)
            {
                DN.bvpSolver->loadController().output(DN,runID,f_file,F_labels);
            }
            
#ifdef userOutputFile
#include userOutputFile
#endif
            
            f_file<<std::endl;
            if(runID==0)
            {
                F_labels<<std::flush;
            }
            
        }
        
    };
    
}
#endif

