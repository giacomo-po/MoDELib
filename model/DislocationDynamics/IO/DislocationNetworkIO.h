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
#include <model/DislocationDynamics/DDtimeIntegrator.h>
#include <model/DislocationDynamics/DislocationNodeContraction.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundaryTransmission.h>
#include <model/DislocationDynamics/IO/DislocationLinkingNumber.h>


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
        typedef typename DislocationNetworkType::GlidePlaneObserverType GlidePlaneObserverType;
        typedef typename DislocationNetworkType::SpatialCellObserverType SpatialCellObserverType;
        typedef typename SpatialCellObserverType::CellMapType CellMapType;
        typedef typename DislocationNetworkType::BvpSolverType::FiniteElementType FiniteElementType;
        typedef typename DislocationNetworkType::BvpSolverType::TrialFunctionType TrialFunctionType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef typename DislocationNetworkType::LoopType LoopType;
        typedef typename DislocationNetworkType::LinkType LinkType;
        typedef typename DislocationNetworkType::StressField StressField;
        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
        
        enum {NdofXnode=NodeType::NdofXnode};
        
        
        DislocationNetworkType& DN;
        
        /**********************************************************************/
        DislocationNetworkIO(DislocationNetworkType& DN_in) :
        /* init */ DN(DN_in)
        {
            
        }

        
        /**********************************************************************/
        void read(const std::string& inputDirectoryName_in, std::string inputFileName)
        { // TO DO: move this to DislocationNetworkIO.h
            
            std::ostringstream fullName;
            fullName<<inputDirectoryName_in<<inputFileName;
            
            model::cout<<greenBoldColor<<"Reading "<<fullName.str()<<"..."<<defaultColor<<std::endl;
            
            
            // Create a file-reader object
            EigenDataReader EDR;
            
            // IO
            EDR.readScalarInFile(fullName.str(),"outputFrequency",DN.outputFrequency);
            EDR.readScalarInFile(fullName.str(),"outputBinary",DN.outputBinary);
            EDR.readScalarInFile(fullName.str(),"outputGlidePlanes",DN.outputGlidePlanes);
            EDR.readScalarInFile(fullName.str(),"outputSpatialCells",DN.outputSpatialCells);
            EDR.readScalarInFile(fullName.str(),"outputPKforce",DN.outputPKforce);
            EDR.readScalarInFile(fullName.str(),"outputMeshDisplacement",DN.outputMeshDisplacement);
            EDR.readScalarInFile(fullName.str(),"outputElasticEnergy",DN.outputElasticEnergy);
            EDR.readScalarInFile(fullName.str(),"outputLinkingNumbers",DN.outputLinkingNumbers);
            EDR.readScalarInFile(fullName.str(),"outputLoopLength",DN.outputLoopLength);
            
            
            EDR.readScalarInFile(fullName.str(),"outputPlasticDistortion",DN.outputPlasticDistortion);
            if(DN.outputPlasticDistortion)
            {
                DN._userOutputColumn+=9;
            }
            EDR.readScalarInFile(fullName.str(),"outputPlasticDistortionRate",DN.outputPlasticDistortionRate);
            if(DN.outputPlasticDistortionRate)
            {
                DN._userOutputColumn+=9;
            }
            
            EDR.readScalarInFile(fullName.str(),"outputDislocationLength",DN.outputDislocationLength);
            if(DN.outputDislocationLength)
            {
                DN._userOutputColumn+=3;
            }
            
            EDR.readScalarInFile(fullName.str(),"outputQuadratureParticles",DN.outputQuadratureParticles);
            
            // Parametrization exponent
            EDR.readScalarInFile(fullName.str(),"parametrizationExponent",LinkType::alpha);
            assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
            assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
            
            // Temperature. Make sure you initialize before calling Material<Isotropic>::select()
            EDR.readScalarInFile(fullName.str(),"temperature",Material<Isotropic>::T); // temperature
            
            // Material and crystal orientation
            unsigned int materialZ;
            EDR.readScalarInFile("./polyCrystalInput.txt","material",materialZ); // material by atomic number Z
            Material<Isotropic>::select(materialZ);
            
            // quadPerLength
            EDR.readScalarInFile(fullName.str(),"quadPerLength",LinkType::quadPerLength); // quadPerLength
            
            // core size
            EDR.readScalarInFile(fullName.str(),"coreSize",StressField::a); // core-width
            assert((StressField::a)>0.0 && "coreSize MUST BE > 0.");
            StressField::a2=StressField::a*StressField::a;
            

            // multipole expansion
            double cellSize(0.0);
            EDR.readScalarInFile(fullName.str(),"dislocationCellSize",cellSize);
            SpatialCellObserverType::setCellSize(cellSize);
            EDR.readScalarInFile(fullName.str(),"use_DisplacementMultipole",DislocationDisplacement<dim>::use_multipole);
            EDR.readScalarInFile(fullName.str(),"use_StressMultipole",DislocationStress<dim>::use_multipole);
            EDR.readScalarInFile(fullName.str(),"use_EnergyMultipole",DislocationEnergy<dim>::use_multipole);
            
           
            //dt=0.0;
            EDR.readScalarInFile(fullName.str(),"timeIntegrationMethod",DN.timeIntegrationMethod);
            EDR.readScalarInFile(fullName.str(),"dxMax",DDtimeIntegrator<0>::dxMax);
            assert(DDtimeIntegrator<0>::dxMax>0.0);
            //            EDR.readScalarInFile(fullName.str(),"shearWaveSpeedFraction",shearWaveSpeedFraction);
            //            assert(shearWaveSpeedFraction>=0.0);
            EDR.readScalarInFile(fullName.str(),"use_velocityFilter",NodeType::use_velocityFilter);
            EDR.readScalarInFile(fullName.str(),"velocityReductionFactor",NodeType::velocityReductionFactor);
            assert(NodeType::velocityReductionFactor>0.0 && NodeType::velocityReductionFactor<=1.0);

            EDR.readScalarInFile(fullName.str(),"computeDDinteractions",DN.computeDDinteractions);

            
            
            // Eternal Stress
           // EDR.readMatrixInFile(fullName.str(),"externalStress",DN.externalStress);
            EDR.readScalarInFile("./loadInput.txt","use_externalStress",DN.use_externalStress);
            if (DN.use_externalStress)
            {
                 DN._userOutputColumn+=18;  //put here in order for right bvp restart
            } 
            // Use Changing external stress field induced by straight dislocations. 
            EDR.readScalarInFile("./loadInput.txt","use_externaldislocationstressfield",DN.use_externaldislocationstressfield); 
            if (DN.use_externaldislocationstressfield)
            {
				DN.ssdeq.clear();
				typedef IDreader<'B',1,10,double> IDreaderType;
				IDreaderType vReader;
				if (vReader.isGood(0,true))
				{ 
					vReader.read(0,true);
					 for (const auto& vIter : vReader)
                    {
						Eigen::Map<const Eigen::Matrix<double,1,9>> row(vIter.second.data());
                        VectorDimD P0(row.template segment<dim>(0));// P0 position
                        VectorDimD P1(row.template segment<dim>(dim)); // P1 position
                        VectorDimD B(row.template segment<dim>(dim*2));  // Burgers vector 
                        DN.ssdeq.emplace_back(StressStraight<dim>(P0,P1,B));                      
					}
				}
				else
				{
					model::cout<<"could not read runID from B/B_0.txt"<<std::endl;
				}					
			}
                        
            // Restart
            EDR.readScalarInFile(fullName.str(),"startAtTimeStep",DN.runID);
            //            VertexReader<'F',201,double> vReader;
            IDreader<'F',1,200,double> vReader;
            Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
            
            
            if (vReader.isGood(0,true))
            {
                
                vReader.read(0,true);
                
                if(DN.runID<0)
                {
                    if(vReader.size())
                    {
                        DN.runID=vReader.rbegin()->first;
                        temp=Eigen::Map<Eigen::Matrix<double,1,200>>(vReader.rbegin()->second.data());
                    }
                    else
                    {
                        DN.runID=0;
                    }
                }
                else
                {
                    const auto iter=vReader.find(DN.runID);
                    if(iter!=vReader.end())
                    {// runID has been found
                        temp=Eigen::Map<Eigen::Matrix<double,1,200>>(iter->second.data());
                    }
                    else
                    {
                        assert(0 && "runID NOT FOUND IN F/F_0.txt");
                    }
                }
                
                DN.dt=temp(1);
                
            }
            else
            {
                model::cout<<"could not read runID from F/F_0.txt"<<std::endl;
                DN.runID=0;
                DN.dt=10.0;
            }
            
            model::cout<<"dt="<<DN.dt<<std::endl;
            
            size_t curCol=0;
            DN.totalTime=temp(curCol);
            curCol+=2;
            
            if (DN.outputPlasticDistortion)
            {
                std::cout<<"reading PD"<<std::endl;
                
                for(int r=0;r<3;++r)
                {
                    for(int c=0;c<3;++c)
                    {
                        DN._plasticDistortion(r,c)=temp(curCol);
                        curCol+=1;
                    }
                }
            }
            
            model::cout<<"starting at time step "<<DN.runID<<std::endl;
            model::cout<<"totalTime= "<<DN.totalTime<<std::endl;
            model::cout<<"plasticDistortion=\n "<<DN._plasticDistortion<<std::endl;
            
            // time-stepping
            
            //            EDR.readScalarInFile(fullName.str(),"useImplicitTimeIntegration",useImplicitTimeIntegration);
            EDR.readScalarInFile(fullName.str(),"ddSolverType",DN.ddSolverType);
            
            
            EDR.readScalarInFile(fullName.str(),"Nsteps",DN.Nsteps);
            assert(DN.Nsteps>=0 && "Nsteps MUST BE >= 0");
            
            //            EDR.readScalarInFile(fullName.str(),"timeWindow",timeWindow);
            //            assert(timeWindow>=0.0 && "timeWindow MUST BE >= 0");
            
            // Check balance
            //            EDR.readScalarInFile(fullName.str(),"check_balance",check_balance);
            
            // JUNCTION FORMATION
            EDR.readScalarInFile(fullName.str(),"maxJunctionIterations",DN.maxJunctionIterations);


            
            //            EDR.readScalarInFile(fullName.str(),"collisionTol",DislocationJunctionFormation<DislocationNetworkType>::collisionTol);
            
            // VERTEX REDISTRIBUTION
            EDR.readScalarInFile(fullName.str(),"use_redistribution",DislocationNetworkRemesh<DislocationNetworkType>::use_redistribution);
            EDR.readScalarInFile(fullName.str(),"Lmin",DislocationNetworkRemesh<DislocationNetworkType>::Lmin);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=0.0);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=2.0*DDtimeIntegrator<0>::dxMax && "YOU MUST CHOOSE Lmin>2*dxMax.");
            EDR.readScalarInFile(fullName.str(),"Lmax",DislocationNetworkRemesh<DislocationNetworkType>::Lmax);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmax>DislocationNetworkRemesh<DislocationNetworkType>::Lmin);
            
            // Cross-Slip
            EDR.readScalarInFile(fullName.str(),"use_crossSlip",DN.use_crossSlip);
            if(DN.use_crossSlip)
            {
                EDR.readScalarInFile(fullName.str(),"crossSlipDeg",DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg);
                assert(DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg>=0.0 && DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg <= 90.0 && "YOU MUST CHOOSE 0.0<= crossSlipDeg <= 90.0");
                //                EDR.readScalarInFile(fullName.str(),"crossSlipLength",DislocationCrossSlip<DislocationNetworkType>::crossSlipLength);
                //                assert(DislocationCrossSlip<DislocationNetworkType>::crossSlipLength>=DislocationNetworkRemesh<DislocationNetworkType>::Lmin && "YOU MUST CHOOSE crossSlipLength>=Lmin.");
                EDR.readScalarInFile(fullName.str(),"verboseCrossSlip",DislocationCrossSlip<DislocationNetworkType>::verboseCrossSlip);

            }
            
            // Mesh and BVP
            EDR.readScalarInFile(fullName.str(),"use_boundary",DN.use_boundary);
            if (DN.use_boundary)
            {
                //                EDR.readScalarInFile(fullName.str(),"use_meshRegions",use_meshRegions);
                
                int meshID(0);
                EDR.readScalarInFile(fullName.str(),"meshID",meshID);
                DN.mesh.readMesh(meshID);
                assert(DN.mesh.simplices().size() && "MESH IS EMPTY.");
                
                // Initialize Polycrystal
                DN.poly.init(DN,"./polyCrystalInput.txt");
                
                EDR.readScalarInFile(fullName.str(),"use_virtualSegments",DN.use_virtualSegments);
                if(DN.use_virtualSegments)
                {
                    EDR.readScalarInFile(fullName.str(),"virtualSegmentDistance",LinkType::virtualSegmentDistance);
                }
                
                EDR.readScalarInFile(fullName.str(),"use_bvp",DN.use_bvp);
                if(DN.use_bvp)
                {
                    EDR.readScalarInFile(fullName.str(),"use_directSolver_FEM",DN.bvpSolver.use_directSolver);
                    EDR.readScalarInFile(fullName.str(),"solverTolerance",DN.bvpSolver.tolerance);
                    DN.bvpSolver.init(DN);
                }
            }
            else{ // no boundary is used, DislocationNetwork is in inifinite medium
                DN.use_bvp=0;	// never comupute boundary correction
            }
            
            if (DN.use_externalStress)
            {
                 DN.extStressController.init(DN);  // have to initialize it after mesh!
            } 
            
            // Verbose levels
            if(DN.maxJunctionIterations>0)
            {
                EDR.readScalarInFile(fullName.str(),"verboseJunctions",DislocationJunctionFormation<DislocationNetworkType>::verboseJunctions);
            }
            
            EDR.readScalarInFile(fullName.str(),"verboseNodeContraction",DislocationNodeContraction<DislocationNetworkType>::verboseNodeContraction);
            EDR.readScalarInFile(fullName.str(),"verboseDislocationNode",NodeType::verboseDislocationNode);

            // GrainBoundary model
            EDR.readScalarInFile(fullName.str(),"grainBoundaryTransmissionModel",GrainBoundaryTransmission<DislocationNetworkType>::grainBoundaryTransmissionModel);

            
            
            // Grain Boundary flags
            //            EDR.readScalarInFile(fullName.str(),"use_GBdissociation",GrainBoundaryDissociation<DislocationNetworkType>::use_GBdissociation);
            //            EDR.readScalarInFile(fullName.str(),"use_GBtransmission",GrainBoundaryTransmission<DislocationNetworkType>::use_GBtransmission);
            //            EDR.readScalarInFile(fullName.str(),"use_GBdislocations",GrainBoundary<dim>::use_GBdislocations);
            
            // Read Vertex and Edge information
            readVertices(DN.runID); // this requires mesh to be up-to-date
            readEdges(DN.runID);    // this requires mesh to be up-to-date
            
            //#ifdef DislocationNucleationFile
            //            EDR.readScalarInFile(fullName.str(),"nucleationFreq",nucleationFreq);
            //#endif
            
#ifdef _MODEL_MPI_
            // Avoid that a processor starts writing before other are reading
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            
            // Initializing configuration
            //            DN.move(0.0);	// initial configuration
        }
        
        /* readVertices *******************************************************/
        void readVertices(const size_t& fileID)
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
        void readEdges(const size_t& fileID)
        {/*! Reads file E/E_0.txt and creates DislocationSegments
          */
            
            // Reading loops
            //            typedef VertexReader<'L',11,double> VertexReaderType;
            typedef IDreader<'L',1,13,double> VertexReaderType;
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
            {// for each line of the L files
                Eigen::Map<const Eigen::Matrix<double,1,10>> row(loop.second.data());
                
                const size_t loopID=loop.first;
                const size_t grainID=row(9);
                
                const auto loopFound=loopMap.find(loopID); // there must be an entry with key loopID in loopMap
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
                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            
        }
        
        /* outputTXT **********************************************************/
        void outputFiles(const size_t& runID) const
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
            
            if (DN.outputBinary)
            {
                assert(0 && "FINISH BIN OUTPUT");
                
                typedef DislocationNodeIO<dim> BinVertexType;
                SequentialBinFile<'V',BinVertexType>::set_count(runID);
                SequentialBinFile<'V',BinVertexType>::set_increment(DN.outputFrequency);
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
                SequentialOutputFile<'V',1>::set_increment(DN.outputFrequency); // vertexFile;
                SequentialOutputFile<'V',1> vertexFile;
                //vertexFile << *(const NetworkNodeContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
                for (const auto& node : DN.nodes())
                {
                    vertexFile << *node.second<<"\n";
                }
                model::cout<<", V/V_"<<vertexFile.sID<<".txt"<<std::flush;
                
                SequentialOutputFile<'L',1>::set_count(runID); // vertexFile;
                SequentialOutputFile<'L',1>::set_increment(DN.outputFrequency); // vertexFile;
                SequentialOutputFile<'L',1> loopFile;
                //vertexFile << *(const NetworkNodeContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
                for (const auto& loop : DN.loops())
                {
                    loopFile << *loop.second<<"\n";
                }
                model::cout<<", L/L_"<<loopFile.sID<<".txt"<<std::flush;
                
                SequentialOutputFile<'E',1>::set_count(runID); // vertexFile;
                SequentialOutputFile<'E',1>::set_increment(DN.outputFrequency); // vertexFile;
                SequentialOutputFile<'E',1> loopLinkFile;
                //vertexFile << *(const NetworkNodeContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
                for (const auto& loopLink : DN.loopLinks())
                {
                    loopLinkFile << loopLink.second <<"\n";
                }
                model::cout<<", E/E_"<<loopLinkFile.sID<<".txt"<<std::flush;
                
                SequentialOutputFile<'K',1>::set_count(runID); // linkFile;
                SequentialOutputFile<'K',1>::set_increment(DN.outputFrequency); // linkFile;
                SequentialOutputFile<'K',1> linkFile;
                //linkFile << *(const NetworkLinkContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
                for (const auto& linkIter : DN.networkLinks())
                {
                    linkFile<< *linkIter.second<<"\n";
                }
                model::cout<<" K/K_"<<linkFile.sID<<".txt"<<std::flush;
                
            }
            
            
            if(DN.outputSpatialCells)
            {
                //! 3- Outputs the nearest neighbor Cell structures to file C_*.txt where * is the current simulation step
                SequentialOutputFile<'C',1>::set_count(runID); // Cell_file;
                SequentialOutputFile<'C',1>::set_increment(DN.outputFrequency); // Cell_file;
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
            
            if(DN.outputGlidePlanes)
            {
                //! 4- Outputs the glide planes
                SequentialOutputFile<'G',1>::set_count(runID); // GlidePlanes_file;
                SequentialOutputFile<'G',1>::set_increment(DN.outputFrequency); // GlidePlanes_file;
                SequentialOutputFile<'G',1> glide_file;
                glide_file << *dynamic_cast<const GlidePlaneObserverType*>(&DN);
                model::cout<<", G/G_"<<glide_file.sID<<std::flush;
            }
            
            if(DN.outputPKforce)
            {
                SequentialOutputFile<'P',1>::set_count(runID); // Edges_file;
                SequentialOutputFile<'P',1>::set_increment(DN.outputFrequency); // Edges_file;
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
            
            if(DN.outputElasticEnergy)
            {
                typedef typename DislocationNetworkType::DislocationParticleType::ElasticEnergy ElasticEnergy;
                SequentialOutputFile<'W',1>::set_count(runID);
                SequentialOutputFile<'W',1>::set_increment(DN.outputFrequency);
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
            
            if(DN.outputMeshDisplacement)
            {
                if(DN.use_bvp)
                {
                    //                    if (!(DN.runningID()%DN.use_bvp))
                    //                    {
                    const auto t0=std::chrono::system_clock::now();
                    model::SequentialOutputFile<'D',1>::set_count(runID); // Vertices_file;
                    model::SequentialOutputFile<'D',1>::set_increment(DN.outputFrequency); // Vertices_file;
                    model::SequentialOutputFile<'D',true> d_file;
                    model::cout<<"		writing to D/D_"<<d_file.sID<<std::flush;
                    
                    std::deque<FieldPointType,Eigen::aligned_allocator<FieldPointType>> fieldPoints; // the container of field points
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
            
            
            if (DN.use_bvp && DN.outputFEMsolution && !(DN.runningID()%DN.use_bvp))
            {
                /**************************************************************************/
                // Output displacement and stress on external mesh faces
                const auto t0=std::chrono::system_clock::now();
                model::SequentialOutputFile<'U',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'U',1>::set_increment(DN.outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'U',true> u_file;
                model::cout<<"		writing to U/U_"<<u_file.sID<<".txt"<<std::flush;
                u_file<<DN.bvpSolver.displacement().onBoundary();
                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                
                const auto t1=std::chrono::system_clock::now();
                model::SequentialOutputFile<'S',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'S',1>::set_increment(DN.outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'S',true> s_file;
                model::cout<<"		writing to S/S_"<<s_file.sID<<".txt"<<std::flush;
                s_file<<DN.bvpSolver.stress().onBoundary();
                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            
            if (DN.outputQuadratureParticles)
            {
                model::SequentialOutputFile<'Q',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'Q',1>::set_increment(DN.outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'Q',true> q_file;
                for (const auto& particle : DN.particles())
                {
                    q_file<<particle<<"\n";
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
            model::cout<<" F/F_0.txt"<<std::flush;
            f_file<< DN.runningID()<<" "<<std::setprecision(15)<<std::scientific<<DN.get_totalTime()<<" "<<DN.get_dt()<<" ";
            
            if(DN.outputPlasticDistortion)
            {
                const Eigen::Matrix<double,dim,dim>& pD(DN.plasticDistortion());
                f_file<<pD.row(0)<<" "<<pD.row(1)<<" "<<pD.row(2)<<" ";
            }
            
            if(DN.outputPlasticDistortionRate)
            {
                const Eigen::Matrix<double,dim,dim>& pDR(DN.plasticDistortionRate());
                f_file<<pDR.row(0)<<" "<<pDR.row(1)<<" "<<pDR.row(2)<<" ";
            }
            
            if(DN.outputDislocationLength)
            {
                const std::tuple<double,double,double> length=DN.networkLength();
                //                const auto length=DN.networkLength();
                //                f_file<<length.first<<" "<<length.second<<" ";
                f_file<<std::get<0>(length)<<" "<<std::get<1>(length)<<" "<<std::get<2>(length)<<" ";
                
            }
            if (DN.use_externalStress)
            {
               f_file<<DN.extStressController.output();
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
  /*  template <typename DislocationNetworkType>
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
    unsigned int DislocationNetworkIO<DislocationNetworkType>::_userOutputColumn=3;*/
    
}
#endif

