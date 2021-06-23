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


#include <UniqueOutputFile.h>
#include <SequentialOutputFile.h>
#include <TerminalColors.h>

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

#include <EshelbyInclusion.h>
//#include <DisplacementPoint.h>
#include <FEMnodeEvaluation.h>



namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename DislocationNetworkType>
    struct DislocationNetworkIO
    {
        
//        enum {dim=DislocationNetworkType::dim};
        static constexpr int dim=TypeTraits<DislocationNetworkType>::dim;

        //    public:
//        typedef typename DislocationNetworkType::VectorDim VectorDim;
//        typedef typename DislocationNetworkType::MatrixDim MatrixDim;
//        typedef typename DislocationNetworkType::NodeType NodeType;
//        typedef typename DislocationNetworkType::BvpSolverType::FiniteElementType FiniteElementType;
//        typedef typename FiniteElementType::ElementType ElementType;
//        typedef typename DislocationNetworkType::BvpSolverType::TrialFunctionType TrialFunctionType;
//        typedef LatticeVector<dim> LatticeVectorType;
//        typedef typename DislocationNetworkType::LoopType LoopType;
//        typedef typename DislocationNetworkType::LinkType LinkType;
        //        typedef typename DislocationNetworkType::StressField StressField;
//        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
        
//        enum {NdofXnode=NodeType::NdofXnode};
        
        
        const DislocationNetworkType& DN;
        const std::string suffix;
        
        /**********************************************************************/
        DislocationNetworkIO(const DislocationNetworkType& DN_in,
                             const std::string& suffix_in="") :
        /* init */ DN(DN_in),
        /* init */ suffix(suffix_in)
        {
            
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
            return DDconfigIO<dim>(DN,suffix);
        }
        
        /**********************************************************************/
        DDauxIO<dim> auxIO() const
        {
            return DDauxIO<dim>(DN,suffix);
        }
        
        /* outputTXT **********************************************************/
        void outputFiles(const size_t& runID) const
        {/*! Outputs DislocationNetwork data to the following files (x is the runID):
          */
            configIO().write(runID,DN.outputBinary);
            auxIO().write(runID,DN.outputBinary);

            if(DN.outputElasticEnergy)
            {
                //                this->template computeNeighborField<ElasticEnergy>();
                
                assert(0 && "RE-IMPLEMENT THIS FOR STRAIGHT SEGMENTS");
                //
                //                
                //                if(outputElasticEnergy)
                //                {
                //                    //                typedef typename DislocationParticleType::ElasticEnergy ElasticEnergy;
                //                }
                
            }
            
            //            typedef BoundaryDisplacementPoint<DislocationNetworkType> FieldPointType;
            //            typedef typename FieldPointType::DisplacementField DisplacementField;
            
//            if(DN.outputMeshDisplacement)
//            {
//                
//                const auto t0=std::chrono::system_clock::now();
//                model::SequentialOutputFile<'D',1>::set_count(runID); // Vertices_file;
//                model::SequentialOutputFile<'D',1>::set_increment(DN.outputFrequency); // Vertices_file;
//                model::SequentialOutputFile<'D',true> d_file;
//                std::cout<<"        writing to D/D_"<<d_file.sID<<std::flush;
//                
//                std::vector<FEMnodeEvaluation<ElementType,dim,1>> fieldPoints; // the container of field points
//                fieldPoints.reserve(DN.mesh.template observer<0>().size());
//                for (const auto& sIter : DN.mesh.template observer<0>())
//                {
//                    if(sIter.second->isBoundarySimplex())
//                    {
//                        fieldPoints.emplace_back(sIter.second->xID(0),sIter.second->P0);
//                    }
//                }
//                
//                DN.displacement(fieldPoints);
//                
//                for(auto& node : fieldPoints)
//                {// add FEM solution and output
//                    if(DN.simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM)
//                    {
//                        const size_t femID=DN.bvpSolver->finiteElement().mesh2femIDmap().at(node.pointID)->gID;
//                        node+=DN.bvpSolver->displacement().dofs(femID);
//                    }
//                    d_file<<node.pointID<<" "<<node.transpose()<<"\n";
//                }
//                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
//            }
            
            if(DN.outputSegmentPairDistances)
            {
                const auto t0=std::chrono::system_clock::now();
                model::SequentialOutputFile<'H',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'H',1>::set_increment(DN.outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'H',true> h_file;
                std::cout<<"		writing to H/H_"<<h_file.sID<<".txt"<<std::flush;
                
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
            
            if (DN.bvpSolver && DN.outputFEMsolution )
            {
                if(!(runID%DN.bvpSolver->stepsBetweenBVPupdates))
                {// Output displacement and stress on external mesh faces
                    const auto t0=std::chrono::system_clock::now();
                    model::SequentialOutputFile<'U',1>::set_count(runID); // Vertices_file;
                    model::SequentialOutputFile<'U',1>::set_increment(DN.outputFrequency); // Vertices_file;
                    model::SequentialOutputFile<'U',true> u_file;
                    std::cout<<"		writing to U/U_"<<u_file.sID<<".txt"<<std::flush;
                    u_file<<DN.bvpSolver->displacement().onBoundary();
                    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                    
                    const auto t1=std::chrono::system_clock::now();
                    model::SequentialOutputFile<'S',1>::set_count(runID); // Vertices_file;
                    model::SequentialOutputFile<'S',1>::set_increment(DN.outputFrequency); // Vertices_file;
                    model::SequentialOutputFile<'S',true> s_file;
                    std::cout<<"		writing to S/S_"<<s_file.sID<<".txt"<<std::flush;
                    s_file<<DN.bvpSolver->stress().onBoundary();
                    std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;
                }
            }
            

            
            if (DN.outputLinkingNumbers)
            {
                DislocationLinkingNumber<DislocationNetworkType> LN(DN);
                model::SequentialOutputFile<'Z',1>::set_count(runID);                   // Vertices_file;
                model::SequentialOutputFile<'Z',1>::set_increment(DN.outputFrequency);  // Vertices_file;
                model::SequentialOutputFile<'Z',true> z_file;
                z_file<<LN;
                
            }
            
            // Output to F file
            UniqueOutputFile<'F'> f_file;
            std::cout<<"Writing F/F_0.txt"<<std::flush;
            
//            std::ofstream F_labels ("F/F_labels.txt", std::ios::out | std::ios::app);
            std::ofstream F_labels;
            if(runID==0)
            {
                F_labels.open("F/F_labels.txt");
            }
            
            f_file<< runID<<" "<<std::setprecision(15)<<std::scientific<<DN.simulationParameters.totalTime<<" "<<DN.simulationParameters.dt<<" ";
            //            int labelCol=0;
            if(runID==0)
            {
                F_labels<<"runID\n";
                F_labels<<"time [b/cs]\n";
                F_labels<<"dt [b/cs]\n";
                //                labelCol+=3;
            }
            
            
            //            if(DN.outputPlasticDistortion)
            //            {
            const Eigen::Matrix<double,dim,dim>& pD(DN.plasticDistortion());
            f_file<<pD.row(0)<<" "<<pD.row(1)<<" "<<pD.row(2)<<" ";
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
                //                    labelCol+=9;
            }
            //            }
            
            if(DN.outputPlasticDistortionRate)
            {
                const Eigen::Matrix<double,dim,dim>& pDR(DN.plasticDistortionRate());
                f_file<<pDR.row(0)<<" "<<pDR.row(1)<<" "<<pDR.row(2)<<" ";
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
                    //                    labelCol+=9;
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
                    F_labels<<"grain boundary length [b]\n";
                    //                    labelCol+=4;
                }
            }
            
            if (DN.externalLoadController)
            {
                DN.externalLoadController->output(runID,f_file,F_labels);
            }
            
            if(DN.bvpSolver)
            {
                DN.bvpSolver->loadController().output(DN,runID,f_file,F_labels);
                
                //                f_file<<std::setprecision(15)<<std::scientific<<DN.bvpSolver->loadController().output(DN,runID,f_file,F_labels);
                //                if(runID==0)
                //                {
                //                    assert(0 && "FINISH HERE, pass F_labels to loadController.output()");
                //                }
            }
            
#ifdef userOutputFile
#include userOutputFile
#endif
            
            f_file<<std::endl;
            
        }
        
    };
    
}
#endif

