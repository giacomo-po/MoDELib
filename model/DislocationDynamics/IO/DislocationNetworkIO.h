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
#include <model/DislocationDynamics/IO/EVLio.h>
#include <model/DislocationDynamics/ElasticFields/EshelbyInclusion.h>



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
        typedef typename DislocationNetworkType::MatrixDimD MatrixDimD;
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
        const std::string suffix;
        
        /**********************************************************************/
        DislocationNetworkIO(DislocationNetworkType& DN_in,
                             const std::string& suffix_in="") :
        /* init */ DN(DN_in),
        suffix(suffix_in)
        {
            
        }
        
        
        /**********************************************************************/
        void read(const std::string& inputDirectoryName_in, std::string inputFileName) __attribute__ ((deprecated))
        { //
            
//            std::ostringstream fullName;
//            fullName<<inputDirectoryName_in<<inputFileName;
            
//            model::cout<<greenBoldColor<<"Reading "<<fullName.str()<<"..."<<defaultColor<<std::endl;
            
            
            // Create a file-reader object
//            EigenDataReader EDR;
            
            // IO
//            EDR.readScalarInFile(fullName.str(),"outputFrequency",DN.outputFrequency);
//            EDR.readScalarInFile(fullName.str(),"outputBinary",DN.outputBinary);
//            EDR.readScalarInFile(fullName.str(),"outputGlidePlanes",DN.outputGlidePlanes);
//            EDR.readScalarInFile(fullName.str(),"outputSpatialCells",DN.outputSpatialCells);
//            EDR.readScalarInFile(fullName.str(),"outputPKforce",DN.outputPKforce);
//            EDR.readScalarInFile(fullName.str(),"outputMeshDisplacement",DN.outputMeshDisplacement);
//            EDR.readScalarInFile(fullName.str(),"outputElasticEnergy",DN.outputElasticEnergy);
//            EDR.readScalarInFile(fullName.str(),"outputLinkingNumbers",DN.outputLinkingNumbers);
//            EDR.readScalarInFile(fullName.str(),"outputLoopLength",DN.outputLoopLength);
            
            
//            EDR.readScalarInFile(fullName.str(),"outputPlasticDistortion",DN.outputPlasticDistortion);
//            if(DN.outputPlasticDistortion)
//            {
//                DN._userOutputColumn+=9;
//            }
//            EDR.readScalarInFile(fullName.str(),"outputPlasticDistortionRate",DN.outputPlasticDistortionRate);
//            if(DN.outputPlasticDistortionRate)
//            {
//                DN._userOutputColumn+=9;
//            }
//            
//            EDR.readScalarInFile(fullName.str(),"outputDislocationLength",DN.outputDislocationLength);
//            if(DN.outputDislocationLength)
//            {
//                DN._userOutputColumn+=3;
//            }
            
//            EDR.readScalarInFile(fullName.str(),"outputQuadratureParticles",DN.outputQuadratureParticles);
            
            // Parametrization exponent
//            EDR.readScalarInFile(fullName.str(),"parametrizationExponent",LinkType::alpha);
//            assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
//            assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
            
//            // Temperature. Make sure you initialize before calling Material<Isotropic>::select()
//            EDR.readScalarInFile(fullName.str(),"temperature",Material<Isotropic>::T); // temperature
            
//            // Material and crystal orientation
//            unsigned int materialZ;
//            EDR.readScalarInFile("./polyCrystalInput.txt","material",materialZ); // material by atomic number Z
//            Material<Isotropic>::select(materialZ);
            
            // quadPerLength
//            EDR.readScalarInFile(fullName.str(),"quadPerLength",LinkType::quadPerLength); // quadPerLength
            
            // core size
//            EDR.readScalarInFile(fullName.str(),"coreSize",StressField::a); // core-width
//            assert((StressField::a)>0.0 && "coreSize MUST BE > 0.");
//            StressField::a2=StressField::a*StressField::a;
            
            
            // multipole expansion
//            double cellSize(0.0);
//            EDR.readScalarInFile(fullName.str(),"dislocationCellSize",cellSize);
//            SpatialCellObserverType::setCellSize(cellSize);
//            EDR.readScalarInFile(fullName.str(),"use_DisplacementMultipole",DislocationDisplacement<dim>::use_multipole);
//            EDR.readScalarInFile(fullName.str(),"use_StressMultipole",DislocationStress<dim>::use_multipole);
//            EDR.readScalarInFile(fullName.str(),"use_EnergyMultipole",DislocationEnergy<dim>::use_multipole);
            
            
            //dt=0.0;
//            EDR.readScalarInFile(fullName.str(),"dxMax",DDtimeIntegrator<0>::dxMax);
//            assert(DDtimeIntegrator<0>::dxMax>0.0);
            //            EDR.readScalarInFile(fullName.str(),"shearWaveSpeedFraction",shearWaveSpeedFraction);
            //            assert(shearWaveSpeedFraction>=0.0);
            
//            EDR.readScalarInFile(fullName.str(),"use_stochasticForce",DN.use_stochasticForce);
//            int stochasticForceSeed=-1;
//            EDR.readScalarInFile(fullName.str(),"stochasticForceSeed",stochasticForceSeed);
//            if(stochasticForceSeed<0)
//            {
//                StochasticForceGenerator::init(std::chrono::system_clock::now().time_since_epoch().count());
//            }
//            else
//            {
//                StochasticForceGenerator::init(stochasticForceSeed);
//            }
            
            
            
            // Use Changing external stress field induced by straight dislocations.
//            EDR.readScalarInFile("./loadInput.txt","use_externaldislocationstressfield",DN.use_externaldislocationstressfield);
//            if (DN.use_externaldislocationstressfield)
//            {
//                DN.ssdeq.clear();
//                typedef IDreader<'B',1,10,double> IDreaderType;
//                IDreaderType vReader;
//                if (vReader.isGood(0,true))
//                {
//                    vReader.read(0,true);
//                    for (const auto& vIter : vReader)
//                    {
//                        Eigen::Map<const Eigen::Matrix<double,1,9>> row(vIter.second.data());
//                        VectorDimD P0(row.template segment<dim>(0));// P0 position
//                        VectorDimD P1(row.template segment<dim>(dim)); // P1 position
//                        VectorDimD B(row.template segment<dim>(dim*2));  // Burgers vector
//                        DN.ssdeq.emplace_back(StressStraight<dim>(P0,P1,B));
//                    }
//                }
//                else
//                {
//                    model::cout<<"could not read runID from B/B_0.txt"<<std::endl;
//                }
//            }
            
            // Restart
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
                        DN._plasticDistortionFromVelocities(r,c)=temp(curCol);
                        curCol+=1;
                    }
                }
            }
            

            model::cout<<"starting at time step "<<DN.runID<<std::endl;
            model::cout<<"totalTime= "<<DN.totalTime<<std::endl;
            model::cout<<"plasticDistortionFromVelocities=\n "<<DN._plasticDistortionFromVelocities<<std::endl;
            
            // time-stepping
            
//            EDR.readScalarInFile(fullName.str(),"outputDislocationStiffnessAndForce",DislocationNetworkComponentType::outputKF);
            
            
            
            //            EDR.readScalarInFile(fullName.str(),"timeWindow",timeWindow);
            //            assert(timeWindow>=0.0 && "timeWindow MUST BE >= 0");
            
            // Check balance
            //            EDR.readScalarInFile(fullName.str(),"check_balance",check_balance);
            
            
            
            
            //            EDR.readScalarInFile(fullName.str(),"collisionTol",DislocationJunctionFormation<DislocationNetworkType>::collisionTol);
            
            
            // Cross-Slip
//            EDR.readScalarInFile(fullName.str(),"crossSlipModel",DN.crossSlipModel);
//            if(DN.crossSlipModel)
//            {
//                DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg=TextFileParser("inputFiles/DD.txt").readScalar<double>("crossSlipDeg",true);
////                EDR.readScalarInFile(fullName.str(),"crossSlipDeg",DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg);
//                assert(DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg>=0.0 && DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg <= 90.0 && "YOU MUST CHOOSE 0.0<= crossSlipDeg <= 90.0");
//                //                EDR.readScalarInFile(fullName.str(),"crossSlipLength",DislocationCrossSlip<DislocationNetworkType>::crossSlipLength);
//                //                assert(DislocationCrossSlip<DislocationNetworkType>::crossSlipLength>=DislocationNetworkRemesh<DislocationNetworkType>::Lmin && "YOU MUST CHOOSE crossSlipLength>=Lmin.");
////                EDR.readScalarInFile(fullName.str(),"verboseCrossSlip",DislocationCrossSlip<DislocationNetworkType>::verboseCrossSlip);
//                DislocationCrossSlip<DislocationNetworkType>::verboseCrossSlip=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseCrossSlip",true);
//            }
            
            // Mesh and BVP
            if (DN.use_boundary)
            {
//                EDR.readScalarInFile(fullName.str(),"surfaceAttractionDistance",DN.surfaceAttractionDistance);
                
//                EDR.readScalarInFile(fullName.str(),"dislocationImages_x",DN.dislocationImages_x);
//                EDR.readScalarInFile(fullName.str(),"dislocationImages_y",DN.dislocationImages_y);
//                EDR.readScalarInFile(fullName.str(),"dislocationImages_z",DN.dislocationImages_z);
                
                //                EDR.readScalarInFile(fullName.str(),"use_meshRegions",use_meshRegions);
                
//                int meshID(0);
//                EDR.readScalarInFile(fullName.str(),"meshID",meshID);
//                DN.mesh.readMesh(meshID);
//                assert(DN.mesh.simplices().size() && "MESH IS EMPTY.");
                
                // Initialize Polycrystal
//                DN.poly.init(DN,"./polyCrystalInput.txt");
                
//                EDR.readScalarInFile(fullName.str(),"use_virtualSegments",DN.use_virtualSegments);
//                if(DN.use_virtualSegments)
//                {
//                    EDR.readScalarInFile(fullName.str(),"virtualSegmentDistance",LinkType::virtualSegmentDistance);
//                }
                
                if(DN.use_bvp)
                {
                    DN.bvpSolver.use_directSolver=TextFileParser("inputFiles/DD.txt").readScalar<int>("use_directSolver_FEM",true);
                    DN.bvpSolver.tolerance=TextFileParser("inputFiles/DD.txt").readScalar<double>("solverTolerance",true);
//                    
//                    EDR.readScalarInFile(fullName.str(),"use_directSolver_FEM",DN.bvpSolver.use_directSolver);
//                    EDR.readScalarInFile(fullName.str(),"solverTolerance",DN.bvpSolver.tolerance);
                    DN.bvpSolver.init(DN);
                }
            }
            else{ // no boundary is used, DislocationNetwork is in inifinite medium
                DN.use_bvp=0;	// never comupute boundary correction
            }
            
            DN.extStressController.init(DN);  // have to initialize it after mesh!


            
            // VERTEX REDISTRIBUTION
            DislocationNetworkRemesh<DislocationNetworkType>::remeshFrequency=TextFileParser("inputFiles/DD.txt").readScalar<int>("remeshFrequency",true);
//            EDR.readScalarInFile(fullName.str(),"remeshFrequency",DislocationNetworkRemesh<DislocationNetworkType>::remeshFrequency);
            double Lmin=TextFileParser("inputFiles/DD.txt").readScalar<double>("Lmin",true);
//            EDR.readScalarInFile(fullName.str(),"Lmin",Lmin);
            double Lmax=TextFileParser("inputFiles/DD.txt").readScalar<double>("Lmax",true);
            double nodeRemoveAngleDeg=TextFileParser("inputFiles/DD.txt").readScalar<double>("nodeRemoveAngleDeg",true);

            //            EDR.readScalarInFile(fullName.str(),"Lmax",Lmax);
            if(DN.use_boundary)
            {
                const double minMeshSize=std::min(DN.mesh.xMax(0)-DN.mesh.xMin(0),std::min(DN.mesh.xMax(1)-DN.mesh.xMin(1),DN.mesh.xMax(2)-DN.mesh.xMin(2)));
                assert(Lmax<1.0 && "IF USING A BOUNDARY Lmax MUST BE RELATIVE TO BOX SIZE (Lmax<1)");
                assert(Lmin<=Lmax);
                DislocationNetworkRemesh<DislocationNetworkType>::Lmax=Lmax*minMeshSize;
                DislocationNetworkRemesh<DislocationNetworkType>::Lmin=Lmin*minMeshSize;
                DislocationNetworkRemesh<DislocationNetworkType>::cosRemove=cos(nodeRemoveAngleDeg*M_PI/180.0);

            }
            else
            {
                DislocationNetworkRemesh<DislocationNetworkType>::Lmax=Lmax;
                DislocationNetworkRemesh<DislocationNetworkType>::Lmin=Lmin;
            }
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmax>3.0*DislocationNetworkRemesh<DislocationNetworkType>::Lmin);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=0.0);
            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=2.0*DDtimeIntegrator<0>::dxMax && "YOU MUST CHOOSE Lmin>2*dxMax.");
            
            
            
            // Verbose levels
            if(DN.maxJunctionIterations>0)
            {
//                EDR.readScalarInFile(fullName.str(),"verboseJunctions",DislocationJunctionFormation<DislocationNetworkType>::verboseJunctions);
                DislocationJunctionFormation<DislocationNetworkType>::verboseJunctions=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseJunctions",true);
            }
            
//            EDR.readScalarInFile(fullName.str(),"verboseNodeContraction",DislocationNodeContraction<DislocationNetworkType>::verboseNodeContraction);
//            EDR.readScalarInFile(fullName.str(),"verboseDislocationNode",NodeType::verboseDislocationNode);
            DislocationNodeContraction<DislocationNetworkType>::verboseNodeContraction=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseNodeContraction",true);
            
            // GrainBoundary model
//            EDR.readScalarInFile(fullName.str(),"grainBoundaryTransmissionModel",GrainBoundaryTransmission<DislocationNetworkType>::grainBoundaryTransmissionModel);
            GrainBoundaryTransmission<DislocationNetworkType>::grainBoundaryTransmissionModel=TextFileParser("inputFiles/DD.txt").readScalar<int>("grainBoundaryTransmissionModel",true);

//            EDR.readScalarInFile(fullName.str(),"outputSegmentPairDistances",DN.outputSegmentPairDistances);
            
            
            
            // Grain Boundary flags
            //            EDR.readScalarInFile(fullName.str(),"use_GBdissociation",GrainBoundaryDissociation<DislocationNetworkType>::use_GBdissociation);
            //            EDR.readScalarInFile(fullName.str(),"use_GBtransmission",GrainBoundaryTransmission<DislocationNetworkType>::use_GBtransmission);
            //            EDR.readScalarInFile(fullName.str(),"use_GBdislocations",GrainBoundary<dim>::use_GBdislocations);
            
            // Read Vertex and Edge information
            EVLio<dim> evl;
            if(EVLio<dim>::isBinGood(DN.runID,suffix))
            {
                evl.readBin(DN.runID,suffix);
            }
            else
            {
                if(EVLio<dim>::isTxtGood(DN.runID,suffix))
                {
                    evl.readTxt(DN.runID,suffix);
                }
                else
                {
                    std::cout<<"COULD NOT FIND INPUT FILEs evl/evl_"<<DN.runID<<".bin or evl/evl_"<<DN.runID<<".txt"<<std::endl;
                    assert(0 && "COULD NOT FIND INPUT FILEs.");
                }
            }
            createVertices(evl);
            createEdges(evl);
            DN.updatePlasticDistortionFromAreas();
            
            IDreader<'E',1,14,double> inclusionsReader;
            inclusionsReader.read(0,true);
            
            for(const auto& pair : inclusionsReader)
            {
                
                const size_t& inclusionID(pair.first);
                Eigen::Map<const Eigen::Matrix<double,1,14>> row(pair.second.data());
                
                const VectorDimD C(row.template segment<dim>(0));
                const double a(row(dim+0));
                MatrixDimD eT(MatrixDimD::Zero());
                const int typeID(row(13));
                int k=dim+1;
                for(int i=0;i<dim;++i)
                {
                    for(int j=0;j<dim;++j)
                    {
                        eT(i,j)=row(k);
                        k++;
                    }
                }
                
                
                
                EshelbyInclusion<dim>::set_count(inclusionID);
                DN.eshelbyInclusions().emplace(std::piecewise_construct,
                                               std::make_tuple(inclusionID),
                                        std::make_tuple(C,a,eT,DN.poly.nu,DN.poly.mu,typeID) );
            }
            
            //            readVertices(DN.runID); // this requires mesh to be up-to-date
            //            readEdges(DN.runID);    // this requires mesh to be up-to-date
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
        void createVertices(const EVLio<dim>& evl)
        {/*!Creates DislocationNode(s) based on the data read by the EVLio<dim>
          * object.
          */
            size_t kk(1);
            for (const auto& node : evl.nodes())
            {
                const size_t nodeIDinFile(node.sID);
                model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<")"<<std::endl;
                NodeType::set_count(nodeIDinFile);
                const size_t nodeID=DN.insertDanglingNode(node.P,node.V,node.velocityReduction).first->first;
                assert(nodeID==nodeIDinFile);
                kk++;
            }
        }
        
        /* readEdges **********************************************************/
        void createEdges(const EVLio<dim>& evl)
        {/*!
          */
            
            
            
            std::map<size_t,std::map<size_t,size_t>> loopMap;
            for(const auto& looplink : evl.links())
            {
                loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
            }
            
            
            
            assert(loopMap.size()==evl.loops().size());
            
            size_t loopLumber=1;
            for(const auto& loop : evl.loops())
            {// for each loop in the EVLio<dim> object
                
                const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
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
                
                
                
                LoopType::set_count(loop.sID);
                //                const LatticeVector<dim> B=DN.poly.grain(grainID).latticeVector(row.template segment<dim>(0*dim).transpose());
                //                const LatticePlaneBase N(DN.poly.grain(grainID).reciprocalLatticeDirection(row.template segment<dim>(1*dim).transpose())); // BETTER TO CONSTRUCT WITH PRIMITIVE VECTORS ON THE PLANE
                //                const LatticeVector<dim> P=DN.poly.grain(grainID).latticeVector(row.template segment<dim>(2*dim).transpose());
                //                const VectorDimD B=row.template segment<dim>(0*dim).transpose();
                //                const VectorDimD N=row.template segment<dim>(1*dim).transpose(); // BETTER TO CONSTRUCT WITH PRIMITIVE VECTORS ON THE PLANE
                //                const VectorDimD P=row.template segment<dim>(2*dim).transpose();
                
                
                
                model::cout<<"Creating Dislocation Loop "<<loop.sID<<" ("<<loopLumber<<" of "<<evl.loops().size()<<")"<<std::endl;
                const size_t newLoopID=DN.insertLoop(nodeIDs,loop.B,loop.N,loop.P,loop.grainID)->sID;
                assert(loop.sID==newLoopID);
                loopLumber++;
            }
            
            DN.clearDanglingNodes();
            
            
            
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
            
            //            EVLio<dim> evl;
            if (DN.outputBinary)
            {
                EVLio<dim>::writeBin(DN,suffix);
            }
            else
            {
                EVLio<dim>::writeTxt(DN,suffix);
            }
            
            //            if (DN.outputBinary)
            //            {
            //                assert(0 && "FINISH BIN OUTPUT");
            //
            //                typedef DislocationNodeIO<dim> BinVertexType;
            //                SequentialBinFile<'V',BinVertexType>::set_count(runID);
            //                SequentialBinFile<'V',BinVertexType>::set_increment(DN.outputFrequency);
            //                SequentialBinFile<'V',BinVertexType> binVertexFile;
            //                //                for (typename NetworkNodeContainerType::const_iterator nodeIter=DN.nodeBegin();nodeIter!=DN.nodeEnd();++nodeIter)
            //                for (const auto& node : DN.nodes())
            //                {
            //                    binVertexFile.write(BinVertexType(*node.second));
            //                }
            //                model::cout<<" V/V_"<<binVertexFile.sID<<".bin"<<std::flush;
            //
            //            }
            //            else
            //            {
            //                SequentialOutputFile<'V',1>::set_count(runID); // vertexFile;
            //                SequentialOutputFile<'V',1>::set_increment(DN.outputFrequency); // vertexFile;
            //                SequentialOutputFile<'V',1> vertexFile;
            //                //vertexFile << *(const NetworkNodeContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
            //                for (const auto& node : DN.nodes())
            //                {
            //                    vertexFile << *node.second<<"\n";
            //                }
            //                model::cout<<", V/V_"<<vertexFile.sID<<".txt"<<std::flush;
            //
            //                SequentialOutputFile<'L',1>::set_count(runID); // vertexFile;
            //                SequentialOutputFile<'L',1>::set_increment(DN.outputFrequency); // vertexFile;
            //                SequentialOutputFile<'L',1> loopFile;
            //                //vertexFile << *(const NetworkNodeContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
            //                for (const auto& loop : DN.loops())
            //                {
            //                    loopFile << *loop.second<<"\n";
            //                }
            //                model::cout<<", L/L_"<<loopFile.sID<<".txt"<<std::flush;
            //
            //                SequentialOutputFile<'E',1>::set_count(runID); // vertexFile;
            //                SequentialOutputFile<'E',1>::set_increment(DN.outputFrequency); // vertexFile;
            //                SequentialOutputFile<'E',1> loopLinkFile;
            //                //vertexFile << *(const NetworkNodeContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
            //                for (const auto& loopLink : DN.loopLinks())
            //                {
            //                    loopLinkFile << loopLink.second <<"\n";
            //                }
            //                model::cout<<", E/E_"<<loopLinkFile.sID<<".txt"<<std::flush;
            //
            //                SequentialOutputFile<'K',1>::set_count(runID); // linkFile;
            //                SequentialOutputFile<'K',1>::set_increment(DN.outputFrequency); // linkFile;
            //                SequentialOutputFile<'K',1> linkFile;
            //                //linkFile << *(const NetworkLinkContainerType*)(&DN); // intel compiler doesn't accept this, so use following loop
            //                for (const auto& linkIter : DN.networkLinks())
            //                {
            //                    linkFile<< *linkIter.second<<"\n";
            //                }
            //                model::cout<<" K/K_"<<linkFile.sID<<".txt"<<std::flush;
            //
            //            }
            
            
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
            
            if(DN.outputSegmentPairDistances)
            {
                const auto t0=std::chrono::system_clock::now();
                model::SequentialOutputFile<'H',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'H',1>::set_increment(DN.outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'H',true> h_file;
                model::cout<<"		writing to H/H_"<<h_file.sID<<".txt"<<std::flush;
                
                for(auto linkIter1=DN.links().begin();linkIter1!=DN.links().end();++linkIter1)
                {
                    for(auto linkIter2=linkIter1;linkIter2!=DN.links().end();++linkIter2)
                    {
                        if(   !linkIter1->second->isBoundarySegment()
                           && !linkIter2->second->isBoundarySegment()
                           && !linkIter1->second->hasZeroBurgers()
                           && !linkIter2->second->hasZeroBurgers())
                        {
                            SegmentSegmentDistance<dim> ssi(linkIter1->second->source->get_P(),
                                                            linkIter1->second->sink->get_P(),
                                                            linkIter2->second->source->get_P(),
                                                            linkIter2->second->sink->get_P());
                            
                            h_file<<sqrt(ssi.D1)<<" "<<sqrt(ssi.D2)<<" "<<ssi.dMin<<"\n";
                        }
                    }
                }
                
                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
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
            
            std::ofstream F_labels ("F/F_labels.txt", std::ios::out | std::ios::app);
            
            f_file<< DN.runningID()<<" "<<std::setprecision(15)<<std::scientific<<DN.get_totalTime()<<" "<<DN.get_dt()<<" ";
            int labelCol=0;
            if(DN.runningID()==0)
            {
                F_labels<<labelCol+0<<"    step #\n";
                F_labels<<labelCol+1<<"    time [b/cs]\n";
                F_labels<<labelCol+2<<"    dt [b/cs]\n";
                labelCol+=3;
            }
            
            
            if(DN.outputPlasticDistortion)
            {
                const Eigen::Matrix<double,dim,dim>& pD(DN.plasticDistortion());
                f_file<<pD.row(0)<<" "<<pD.row(1)<<" "<<pD.row(2)<<" ";
                if(DN.runningID()==0)
                {
                    F_labels<<labelCol+0<<"    betaP_11\n";
                    F_labels<<labelCol+1<<"    betaP_12\n";
                    F_labels<<labelCol+2<<"    betaP_13\n";
                    F_labels<<labelCol+3<<"    betaP_21\n";
                    F_labels<<labelCol+4<<"    betaP_22\n";
                    F_labels<<labelCol+5<<"    betaP_23\n";
                    F_labels<<labelCol+6<<"    betaP_31\n";
                    F_labels<<labelCol+7<<"    betaP_32\n";
                    F_labels<<labelCol+8<<"    betaP_33\n";
                    labelCol+=9;
                }
            }
            
            if(DN.outputPlasticDistortionRate)
            {
                const Eigen::Matrix<double,dim,dim>& pDR(DN.plasticDistortionRate());
                f_file<<pDR.row(0)<<" "<<pDR.row(1)<<" "<<pDR.row(2)<<" ";
                if(DN.runningID()==0)
                {
                    F_labels<<labelCol+0<<"    dotBetaP_11 [cs/b]\n";
                    F_labels<<labelCol+1<<"    dotBetaP_12 [cs/b]\n";
                    F_labels<<labelCol+2<<"    dotBetaP_13 [cs/b]\n";
                    F_labels<<labelCol+3<<"    dotBetaP_21 [cs/b]\n";
                    F_labels<<labelCol+4<<"    dotBetaP_22 [cs/b]\n";
                    F_labels<<labelCol+5<<"    dotBetaP_23 [cs/b]\n";
                    F_labels<<labelCol+6<<"    dotBetaP_31 [cs/b]\n";
                    F_labels<<labelCol+7<<"    dotBetaP_32 [cs/b]\n";
                    F_labels<<labelCol+8<<"    dotBetaP_33 [cs/b]\n";
                    labelCol+=9;
                }
            }
            
            if(DN.outputDislocationLength)
            {
                const std::tuple<double,double,double,double> length=DN.networkLength();
                f_file<<std::get<0>(length)<<" "<<std::get<1>(length)<<" "<<std::get<2>(length)<<" "<<std::get<3>(length)<<" ";
                if(DN.runningID()==0)
                {
                    F_labels<<labelCol+0<<"    glissile length [b]\n";
                    F_labels<<labelCol+1<<"    sessile length [b]\n";
                    F_labels<<labelCol+2<<"    boundary length [b]\n";
                    F_labels<<labelCol+3<<"    grain boundary length [b]\n";
                    labelCol+=4;
                }
            }
            
            if (DN.use_externalStress)
            {
                DN.extStressController.output(DN.runningID(),f_file,F_labels,labelCol);
            }
            
            if(DN.use_bvp)
            {
                f_file<<std::setprecision(15)<<std::scientific<<DN.bvpSolver.loadController().output(DN);
                if(DN.runningID()==0)
                {
                    assert(0 && "FINISH HERE, pass F_labels to loadController.output()");
                }
            }
            
#ifdef userOutputFile
#include userOutputFile
#endif
            
            f_file<<std::endl;
            
        }
        
    };
    
}
#endif

