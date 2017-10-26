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
#include <model/IO/IDreader.h>

#include <model/IO/UniqueOutputFile.h>
#include <model/IO/SequentialOutputFile.h>
#include <model/IO/SequentialBinFile.h>

#include <model/Utilities/TerminalColors.h>
#include <model/DislocationDynamics/GlidePlanes/GlidePlaneObserver.h>
#include <model/MPI/MPIcout.h>
#include <model/LatticeMath/LatticeMath.h>
#include <model/LatticeMath/LatticeMath.h>
#include <model/DislocationDynamics/BVP/BoundaryDisplacementPoint.h>
#include <model/DislocationDynamics/IO/DislocationNodeIO.h>


namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename DislocationNetworkType>
    struct DislocationNetworkIO
    {
        
        enum {dim=DislocationNetworkType::dim};
        
        //    public:
        typedef typename DislocationNetworkType::VectorDimD VectorDimD;
        typedef typename DislocationNetworkType::NodeType NodeType;
//        typedef typename DislocationNetworkType::NetworkLinkContainerType NetworkLinkContainerType;
//        typedef typename DislocationNetworkType::NetworkNodeContainerType NetworkNodeContainerType;
        typedef typename DislocationNetworkType::GlidePlaneObserverType GlidePlaneObserverType;
        typedef typename DislocationNetworkType::SpatialCellObserverType SpatialCellObserverType;
        typedef typename SpatialCellObserverType::CellMapType CellMapType;
        typedef typename DislocationNetworkType::BvpSolverType::FiniteElementType FiniteElementType;
        typedef typename DislocationNetworkType::BvpSolverType::TrialFunctionType TrialFunctionType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef typename DislocationNetworkType::LoopType LoopType;
        
        enum {NdofXnode=NodeType::NdofXnode};
        
        static int  outputFrequency;
        static bool outputBinary;
        static bool outputGlidePlanes;
        static bool outputSpatialCells;
        static bool outputPKforce;
        static bool outputElasticEnergy;
        static bool outputMeshDisplacement;
        static bool outputFEMsolution;
        static bool outputDislocationLength;
        static bool outputPlasticDistortion;
        static bool outputPlasticDistortionRate;
        static bool outputQuadratureParticles;
        static unsigned int _userOutputColumn;
        
        
        
        /* readVertices *******************************************************/
        static void readVertices(DislocationNetworkType& DN, const size_t& fileID)
        {/*! Reads file V/V_0.txt and creates DislocationNodes
          */
//            typedef VertexReader<'V',7,double> VertexReaderType;
            typedef IDreader<'V',1,10,double> VertexReaderType;
            VertexReaderType  vReader;	// sID,Px,Py,Pz,Tx,Ty,Tz,snID,meshLocation,grainID
            if (vReader.isGood(fileID,false)) // bin file exists
            {
                vReader.read(fileID,false);
                
            }
            else
            {
                if (vReader.isGood(fileID,true)) // txt file exists
                {
                    vReader.read(fileID,true);
                    
                }
                else
                {
                    assert(0 && "UNABLE TO READ VERTEX FILE V/V_x (x is the requested file ID).");
                    
                }
            }
            
            size_t kk(1);
//            for (VertexReaderType::iterator vIter=vReader.begin();vIter!=vReader.end();++vIter)
            for (const auto& vIter : vReader)
            {
                const size_t nodeIDinFile(vIter.first);
                NodeType::set_count(nodeIDinFile);
                model::cout << "\r \r" << "Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<vReader.size()<<")"<<std::flush;
                
                const Eigen::Map<const Eigen::Matrix<double,1,9>> row(vIter.second.data());
                const VectorDimD P(row.template segment<NdofXnode>(0));
                const VectorDimD V(row.template segment<NdofXnode>(NdofXnode));
                const double velocityReductionCoeff(row(6));
//                const double snID(row(7));
//                const bool onMeshBoundary(8);
                
                //                const int grainID(row(9));
                
                const size_t nodeID=DN.insertDanglingNode(P,V,velocityReductionCoeff).first->first;
                
                //LatticeVectorType L(VectorDimD());
                //const size_t nodeID(DN.insertVertex(L).first->first);
//                const size_t nodeID(DN.insertVertex(P,grainID).first->first);
                assert(nodeID==nodeIDinFile);
                kk++;
            }
            model::cout<<std::endl;
        }
        
        /* readEdges **********************************************************/
        static void readEdges(DislocationNetworkType& DN, const size_t& fileID)
        {/*! Reads file E/E_0.txt and creates DislocationSegments
          */
            
            // Reading loops
//            typedef VertexReader<'L',11,double> VertexReaderType;
            typedef IDreader<'L',1,10,double> VertexReaderType;
            VertexReaderType  vReader;	// sID,Px,Py,Pz,Tx,Ty,Tz,snID,meshLocation,grainID
            if (vReader.isGood(fileID,false)) // bin file exists
            {
                vReader.read(fileID,false);
                
            }
            else
            {
                if (vReader.isGood(fileID,true)) // txt file exists
                {
                    vReader.read(fileID,true);
                    
                }
                else
                {
                    assert(0 && "UNABLE TO READ VERTEX FILE L/L_x (x is the requested file ID).");
                    
                }
            }
            
            // Reading LoopLinks
            typedef IDreader<'E',3,0,double> LoopLinkReaderType;
            LoopLinkReaderType llreader;
            if (llreader.isGood(fileID,false)) // bin file exists
            {
                llreader.read(fileID,false);
                
            }
            else
            {
                if (llreader.isGood(fileID,true)) // txt file exists
                {
                    llreader.read(fileID,true);
                    
                }
                else
                {
                    assert(0 && "UNABLE TO READ VERTEX FILE E/E_x (x is the requested file ID).");
                    
                }
            }
            
            
            std::map<size_t,std::map<size_t,size_t>> loopMap;
            
            for(const auto& looplink : llreader)
            {
                loopMap[looplink.first[0]].emplace(looplink.first[1],looplink.first[2]);
            }


            
            assert(loopMap.size()==vReader.size());
            
            size_t loopLumber=1;
            for(const auto& loop : vReader)
            {
                Eigen::Map<const Eigen::Matrix<double,1,10>> row(loop.second.data());

                const size_t loopID=loop.first;
                const size_t grainID=row(9);

                const auto loopFound=loopMap.find(loopID);
                assert(loopFound!=loopMap.end());
                

                
                std::vector<size_t> nodeIDs;
                nodeIDs.push_back(loopFound->second.begin()->first);
                for(size_t k=0;k<loopFound->second.size();++k)
                {
                    const auto nodeFound=loopFound->second.find(*nodeIDs.rbegin());
                    if(k<loopFound->second.size()-1)
                    {
                        nodeIDs.push_back(nodeFound->second);
                    }
                    else
                    {
                        assert(nodeFound->second==nodeIDs[0]);
                    }
                }
                

                
                LoopType::set_count(loopID);
//                const LatticeVector<dim> B=DN.poly.grain(grainID).latticeVector(row.template segment<dim>(0*dim).transpose());
//                const LatticePlaneBase N(DN.poly.grain(grainID).reciprocalLatticeDirection(row.template segment<dim>(1*dim).transpose())); // BETTER TO CONSTRUCT WITH PRIMITIVE VECTORS ON THE PLANE
//                const LatticeVector<dim> P=DN.poly.grain(grainID).latticeVector(row.template segment<dim>(2*dim).transpose());
                const VectorDimD B=row.template segment<dim>(0*dim).transpose();
                const VectorDimD N=row.template segment<dim>(1*dim).transpose(); // BETTER TO CONSTRUCT WITH PRIMITIVE VECTORS ON THE PLANE
                const VectorDimD P=row.template segment<dim>(2*dim).transpose();
                

                
                model::cout<<"Creating Dislocation Loop "<<loopID<<" ("<<loopLumber<<" of "<<vReader.size()<<")"<<std::endl;
                const size_t newLoopID=DN.insertLoop(nodeIDs,B,N,P,grainID)->sID;
                assert(loopID==newLoopID);
                loopLumber++;
            }
            
            DN.clearDanglingNodes();
            
            
//            typedef EdgeReader  <'E',11,double>	EdgeReaderType;
//            EdgeReaderType    eReader;	// sourceID,sinkID,Bx,By,Bz,Nx,Ny,Nz
//            if (eReader.isGood(fileID,false)) // bin file exists
//            {
//                eReader.read(fileID,false);
//                
//            }
//            else
//            {
//                if (eReader.isGood(fileID,true)) // txt file exists
//                {
//                    eReader.read(fileID,true);
//                    
//                }
//                else
//                {
//                    assert(0 && "UNABLE TO READ EDGE FILE E/E_x (x is the requested file ID).");
//                }
//            }
//            
//            unsigned int kk(1);
//            for (EdgeReaderType::iterator eIter=eReader.begin();eIter!=eReader.end();++eIter)
//            {
//                VectorDimD B(eIter->second.template segment<dim>(0  ).transpose()); // Burgers vector
//                VectorDimD N(eIter->second.template segment<dim>(dim).transpose()); // Glide plane normal
//                const size_t sourceID(eIter->first.first );
//                const size_t   sinkID(eIter->first.second);
//                model::cout << "\r \r" << "Creating DislocationSegment "<<sourceID<<"->"<<sinkID<<" ("<<kk<<" of "<<eReader.size()<<")              "<<std::flush;
//                
//                const auto isSource=DN.node(sourceID);
//                assert(isSource.first && "Source node does not exist");
//                const auto isSink=DN.node(sinkID);
//                assert(isSink.first && "Sink node does not exist");
//                
//                assert(isSource.second->grain.grainID==isSink.second->grain.grainID && "Source and Sink are in different Grains");
//                
////                const bool success=DN.connect(sourceID,sinkID,isSource.second->grain.latticeVector(B));
////               assert(success && "UNABLE TO CREATE CURRENT DISLOCATION SEGMENT.");
//                kk++;
//            }
            model::cout<<std::endl;
        }
        
        /* outputTXT **********************************************************/
        static void output(const DislocationNetworkType& DN, const unsigned int& runID)
        {/*! Outputs DislocationNetwork data to the following files (x is the runID):
          * ./E/E_x.txt (DislocationSegment(s) are always outputted)
          * ./V/V_x.txt (DislocationNode(s) are always outputted)
          * ./C/C_x.txt (DislocationCell(s) only if outputSpatialCells==true)
          * ./G/G_x.txt (GlidePlane(s) only if outputGlidePlanes==true)
          * ./P/P_x.txt (PK forces only if outputPKforce==true)
          * ./W/W_x.txt (elastic energy only if outputElasticEnergy==true)
          * ./D/D_x.txt (mesh displacement only if outputMeshDisplacement==true)
          */
            model::cout<<"		Writing to "<<std::flush;
            
            if (outputBinary)
            {
                assert(0 && "FINISH BIN OUTPUT");

                typedef DislocationNodeIO<dim> BinVertexType;
                SequentialBinFile<'V',BinVertexType>::set_count(runID);
                SequentialBinFile<'V',BinVertexType>::set_increment(outputFrequency);
                SequentialBinFile<'V',BinVertexType> binVertexFile;
                //                for (typename NetworkNodeContainerType::const_iterator nodeIter=DN.nodeBegin();nodeIter!=DN.nodeEnd();++nodeIter)
                for (const auto& node : DN.nodes())
                {
                    binVertexFile.write(BinVertexType(*node.second));
                }
                model::cout<<" V/V_"<<binVertexFile.sID<<".bin"<<std::flush;

            }
            else
            {
                SequentialOutputFile<'V',1>::set_count(runID); // vertexFile;
                SequentialOutputFile<'V',1>::set_increment(outputFrequency); // vertexFile;
                SequentialOutputFile<'V',1> vertexFile;
                //vertexFile << *(const NetworkNodeContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
                for (const auto& node : DN.nodes())
                {
                    vertexFile << *node.second<<"\n";
                }
                model::cout<<", V/V_"<<vertexFile.sID<<".txt"<<std::flush;
                
                SequentialOutputFile<'L',1>::set_count(runID); // vertexFile;
                SequentialOutputFile<'L',1>::set_increment(outputFrequency); // vertexFile;
                SequentialOutputFile<'L',1> loopFile;
                //vertexFile << *(const NetworkNodeContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
                for (const auto& loop : DN.loops())
                {
                    loopFile << *loop.second<<"\n";
                }
                model::cout<<", L/L_"<<loopFile.sID<<".txt"<<std::flush;
                
                SequentialOutputFile<'E',1>::set_count(runID); // vertexFile;
                SequentialOutputFile<'E',1>::set_increment(outputFrequency); // vertexFile;
                SequentialOutputFile<'E',1> loopLinkFile;
                //vertexFile << *(const NetworkNodeContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
                for (const auto& loopLink : DN.loopLinks())
                {
                    loopLinkFile << loopLink.second <<"\n";
                }
                model::cout<<", E/E_"<<loopLinkFile.sID<<".txt"<<std::flush;
                
                SequentialOutputFile<'K',1>::set_count(runID); // linkFile;
                SequentialOutputFile<'K',1>::set_increment(outputFrequency); // linkFile;
                SequentialOutputFile<'K',1> linkFile;
                //linkFile << *(const NetworkLinkContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
                for (const auto& linkIter : DN.networkLinks())
                {
                    linkFile<< *linkIter.second<<"\n";
                }
                model::cout<<" K/K_"<<linkFile.sID<<".txt"<<std::flush;

            }
            
            
            if(outputSpatialCells)
            {
                //! 3- Outputs the nearest neighbor Cell structures to file C_*.txt where * is the current simulation step
                SequentialOutputFile<'C',1>::set_count(runID); // Cell_file;
                SequentialOutputFile<'C',1>::set_increment(outputFrequency); // Cell_file;
                SequentialOutputFile<'C',1> Cell_file;
                //              SpatialCellObserverType SPC;
                int cID(0);
                //for (typename CellMapType::const_iterator cellIter=SpatialCellObserverType::cellBegin();cellIter!=SpatialCellObserverType::cellEnd();++cellIter)
                for (const auto& cell : SpatialCellObserverType::cells())
                {
                    Cell_file<<cID<<"\t"<<cell.second->cellID.transpose()<<"\t"<<SpatialCellObserverType::cellSize()
#ifndef _MODEL_ENABLE_CELL_VERTEX_ALPHA_TENSORS_
                    /*     */<<"\t"<<std::get<0>(*cell.second).row(0)
                    /*     */<<"\t"<<std::get<0>(*cell.second).row(1)
                    /*     */<<"\t"<<std::get<0>(*cell.second).row(2)
#endif
                    <<"\n";
                    ++cID;
                }
                model::cout<<", C/C_"<<Cell_file.sID<<std::flush;
            }
            
            if(outputGlidePlanes)
            {
                //! 4- Outputs the glide planes
                SequentialOutputFile<'G',1>::set_count(runID); // GlidePlanes_file;
                SequentialOutputFile<'G',1>::set_increment(outputFrequency); // GlidePlanes_file;
                SequentialOutputFile<'G',1> glide_file;
                glide_file << *dynamic_cast<const GlidePlaneObserverType*>(&DN);
                model::cout<<", G/G_"<<glide_file.sID<<std::flush;
            }
            
            if(outputPKforce)
            {
                SequentialOutputFile<'P',1>::set_count(runID); // Edges_file;
                SequentialOutputFile<'P',1>::set_increment(outputFrequency); // Edges_file;
                SequentialOutputFile<'P',1> p_file;
                for (const auto& linkIter : DN.links())
                {
                    const int qOrder(linkIter.second->rgauss.cols());
                    for (int q=0;q<qOrder;++q)
                    {
                        p_file << linkIter.second->source->sID<<" "
                        /*  */ << linkIter.second->sink->sID<<" "
                        /*  */ <<q<<" "
                        /*  */ << linkIter.second->rgauss.col(q).transpose()<<" "
                        /*  */ <<linkIter.second->pkGauss.col(q).transpose()<<"\n";
                    }
                }
                model::cout<<", P/P_"<<p_file.sID<<std::flush;
            }
            
            if(outputElasticEnergy)
            {
                typedef typename DislocationNetworkType::DislocationParticleType::ElasticEnergy ElasticEnergy;
                SequentialOutputFile<'W',1>::set_count(runID);
                SequentialOutputFile<'W',1>::set_increment(outputFrequency);
                SequentialOutputFile<'W',1> w_file; //energy_file
                int ll=0;
                for (const auto& linkIter : DN.links())
                {
                    const int qOrder(linkIter.second->rgauss.cols());
                    for (size_t q=0;q<linkIter.second->quadratureParticleContainer.size();++q)
                    {
                        w_file << ll*qOrder+q<<" "<< linkIter.second->rgauss.col(q).transpose()<<" "<< linkIter.second->quadratureParticleContainer[q]->template field<ElasticEnergy>()<<"\n";
                    }
                    ll++;
                }
                model::cout<<", W/W_"<<w_file.sID<<std::flush;
            }
            
            typedef BoundaryDisplacementPoint<DislocationNetworkType> FieldPointType;
            typedef typename FieldPointType::DisplacementField DisplacementField;

            if(outputMeshDisplacement)
            {
                if(DN.use_bvp)
                {
                    //                    if (!(DN.runningID()%DN.use_bvp))
                    //                    {
                    const auto t0=std::chrono::system_clock::now();
                    model::SequentialOutputFile<'D',1>::set_count(runID); // Vertices_file;
                    model::SequentialOutputFile<'D',1>::set_increment(outputFrequency); // Vertices_file;
                    model::SequentialOutputFile<'D',true> d_file;
                    model::cout<<"		writing to D/D_"<<d_file.sID<<std::flush;
                    
                    std::deque<FieldPointType> fieldPoints; // the container of field points
                    for (const auto& sIter : DN.mesh.template observer<0>())
                    {
                        if(sIter.second->isBoundarySimplex())
                        {
                            fieldPoints.emplace_back(*(sIter.second));
                        }
                    }
                    DN.template computeField<FieldPointType,DisplacementField>(fieldPoints);
                    
                    
                    for(auto node : fieldPoints)
                    {
                        Eigen::Matrix<double,dim,1> nodeDisp = node.template field<DisplacementField>();
                        
                        // Sum solid angle jump
                        if (DN.use_virtualSegments) // solid-angle jump of virtual segments
                        {
                            for(const auto& segment : DN.links())
                            {
                                segment.second->addToSolidAngleJump(node.P,node.S,nodeDisp);
                            }
                        }
                        
                        // Sum FEM solution
                        const size_t femID=DN.bvpSolver.finiteElement().mesh2femIDmap().at(node.gID)->gID;
                        nodeDisp+=DN.bvpSolver.displacement().dofs(femID);
                        
                        // output
                        d_file<<node.gID<<" "<<nodeDisp.transpose()<<"\n";
                    }
                    model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                }
                
            }

            
            if (DN.use_bvp && outputFEMsolution && !(DN.runningID()%DN.use_bvp))
            {
                /**************************************************************************/
                // Output displacement and stress on external mesh faces
                const auto t0=std::chrono::system_clock::now();
                model::SequentialOutputFile<'U',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'U',1>::set_increment(outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'U',true> u_file;
                model::cout<<"		writing to U/U_"<<u_file.sID<<".txt"<<std::flush;
                u_file<<DN.bvpSolver.displacement().onBoundary();
                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                
                const auto t1=std::chrono::system_clock::now();
                model::SequentialOutputFile<'S',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'S',1>::set_increment(outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'S',true> s_file;
                model::cout<<"		writing to S/S_"<<s_file.sID<<".txt"<<std::flush;
                s_file<<DN.bvpSolver.stress().onBoundary();
                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            
            if (outputQuadratureParticles)
            {
                model::SequentialOutputFile<'Q',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'Q',1>::set_increment(outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'Q',true> q_file;
                for (const auto& particle : DN.particles())
                {
                    q_file<<particle<<"\n";
                }
                
                
            }
            
            // Output to F file
            UniqueOutputFile<'F'> f_file;
            model::cout<<" F/F_0.txt"<<std::flush;
            f_file<< DN.runningID()<<" "<<std::setprecision(15)<<std::scientific<<DN.get_totalTime()<<" "<<DN.get_dt()<<" ";
            
            if(outputPlasticDistortion)
            {
                const Eigen::Matrix<double,dim,dim>& pD(DN.plasticDistortion());
                f_file<<pD.row(0)<<" "<<pD.row(1)<<" "<<pD.row(2)<<" ";
            }
            
            if(outputPlasticDistortionRate)
            {
                const Eigen::Matrix<double,dim,dim>& pDR(DN.plasticDistortionRate());
                f_file<<pDR.row(0)<<" "<<pDR.row(1)<<" "<<pDR.row(2)<<" ";
            }
            
            if(outputDislocationLength)
            {
                const std::tuple<double,double,double> length=DN.networkLength();
                //                const auto length=DN.networkLength();
                //                f_file<<length.first<<" "<<length.second<<" ";
                f_file<<std::get<0>(length)<<" "<<std::get<1>(length)<<" "<<std::get<2>(length)<<" ";
                
            }
            
            if(DN.use_bvp)
            {
                f_file<<std::setprecision(15)<<std::scientific<<DN.bvpSolver.loadController().output(DN);
            }
            
#ifdef userOutputFile
#include userOutputFile
#endif
            
            f_file<<std::endl;
            
        }
        
    };
    
    // Declare static data
    template <typename DislocationNetworkType>
    int DislocationNetworkIO<DislocationNetworkType>::outputFrequency=1;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputBinary=false;
    
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputGlidePlanes=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputSpatialCells=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputPKforce=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputElasticEnergy=true;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputMeshDisplacement=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputFEMsolution=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputDislocationLength=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortion=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputPlasticDistortionRate=false;
    
    template <typename DislocationNetworkType>
    bool DislocationNetworkIO<DislocationNetworkType>::outputQuadratureParticles=false;
    
    template <typename DislocationNetworkType>
    unsigned int DislocationNetworkIO<DislocationNetworkType>::_userOutputColumn=3;
    
    /**************************************************************************/
    /**************************************************************************/
} // namespace model
#endif

