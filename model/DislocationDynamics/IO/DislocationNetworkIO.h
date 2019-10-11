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
#include <SequentialBinFile.h>
#include <TerminalColors.h>
#include <MPIcout.h>
#include <LatticeMath.h>
#include <LatticeMath.h>
//#include <BoundaryDisplacementPoint.h>
#include <DislocationNodeIO.h>
#include <DDtimeIntegrator.h>
#include <DislocationNodeContraction.h>
#include <GrainBoundaryTransmission.h>
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
        
        enum {dim=DislocationNetworkType::dim};
        
        //    public:
//        typedef typename DislocationNetworkType::VectorDim VectorDim;
//        typedef typename DislocationNetworkType::MatrixDim MatrixDim;
        typedef typename DislocationNetworkType::NodeType NodeType;
        typedef typename DislocationNetworkType::BvpSolverType::FiniteElementType FiniteElementType;
        typedef typename FiniteElementType::ElementType ElementType;
        typedef typename DislocationNetworkType::BvpSolverType::TrialFunctionType TrialFunctionType;
        typedef LatticeVector<dim> LatticeVectorType;
        typedef typename DislocationNetworkType::LoopType LoopType;
        typedef typename DislocationNetworkType::LinkType LinkType;
        //        typedef typename DislocationNetworkType::StressField StressField;
        typedef DislocationNetworkComponent<NodeType,LinkType> DislocationNetworkComponentType;
        
        enum {NdofXnode=NodeType::NdofXnode};
        
        
        const DislocationNetworkType& DN;
        const std::string suffix;
        
        /**********************************************************************/
        DislocationNetworkIO(const DislocationNetworkType& DN_in,
                             const std::string& suffix_in="") :
        /* init */ DN(DN_in),
        /* init */ suffix(suffix_in)
        {
            
        }
        
        
//        /**********************************************************************/
//        void read(const std::string& inputDirectoryName_in, std::string inputFileName,
//                  long int& runID) __attribute__ ((deprecated))
//        { //
//            
//            //            std::ostringstream fullName;
//            //            fullName<<inputDirectoryName_in<<inputFileName;
//            
//            //            model::cout<<greenBoldColor<<"Reading "<<fullName.str()<<"..."<<defaultColor<<std::endl;
//            
//            
//            // Create a file-reader object
//            //            EigenDataReader EDR;
//            
//            // IO
//            //            EDR.readScalarInFile(fullName.str(),"outputFrequency",DN.outputFrequency);
//            //            EDR.readScalarInFile(fullName.str(),"outputBinary",DN.outputBinary);
//            //            EDR.readScalarInFile(fullName.str(),"outputGlidePlanes",DN.outputGlidePlanes);
//            //            EDR.readScalarInFile(fullName.str(),"outputSpatialCells",DN.outputSpatialCells);
//            //            EDR.readScalarInFile(fullName.str(),"outputPKforce",DN.outputPKforce);
//            //            EDR.readScalarInFile(fullName.str(),"outputMeshDisplacement",DN.outputMeshDisplacement);
//            //            EDR.readScalarInFile(fullName.str(),"outputElasticEnergy",DN.outputElasticEnergy);
//            //            EDR.readScalarInFile(fullName.str(),"outputLinkingNumbers",DN.outputLinkingNumbers);
//            //            EDR.readScalarInFile(fullName.str(),"outputLoopLength",DN.outputLoopLength);
//            
//            
//            //            EDR.readScalarInFile(fullName.str(),"outputPlasticDistortion",DN.outputPlasticDistortion);
//            //            if(DN.outputPlasticDistortion)
//            //            {
//            //                DN._userOutputColumn+=9;
//            //            }
//            //            EDR.readScalarInFile(fullName.str(),"outputPlasticDistortionRate",DN.outputPlasticDistortionRate);
//            //            if(DN.outputPlasticDistortionRate)
//            //            {
//            //                DN._userOutputColumn+=9;
//            //            }
//            //            
//            //            EDR.readScalarInFile(fullName.str(),"outputDislocationLength",DN.outputDislocationLength);
//            //            if(DN.outputDislocationLength)
//            //            {
//            //                DN._userOutputColumn+=3;
//            //            }
//            
//            //            EDR.readScalarInFile(fullName.str(),"outputQuadraturePoints",DN.outputQuadraturePoints);
//            
//            // Parametrization exponent
//            //            EDR.readScalarInFile(fullName.str(),"parametrizationExponent",LinkType::alpha);
//            //            assert((LinkType::alpha)>=0.0 && "parametrizationExponent MUST BE >= 0.0");
//            //            assert((LinkType::alpha)<=1.0 && "parametrizationExponent MUST BE <= 1.0");
//            
//            //            // Temperature. Make sure you initialize before calling Material<Isotropic>::select()
//            //            EDR.readScalarInFile(fullName.str(),"temperature",Material<Isotropic>::T); // temperature
//            
//            //            // Material and crystal orientation
//            //            unsigned int materialZ;
//            //            EDR.readScalarInFile("./polyCrystalInput.txt","material",materialZ); // material by atomic number Z
//            //            Material<Isotropic>::select(materialZ);
//            
//            // quadPerLength
//            //            EDR.readScalarInFile(fullName.str(),"quadPerLength",LinkType::quadPerLength); // quadPerLength
//            
//            // core size
//            //            EDR.readScalarInFile(fullName.str(),"coreSize",StressField::a); // core-width
//            //            assert((StressField::a)>0.0 && "coreSize MUST BE > 0.");
//            //            StressField::a2=StressField::a*StressField::a;
//            
//            
//            // multipole expansion
//            //            double cellSize(0.0);
//            //            EDR.readScalarInFile(fullName.str(),"dislocationCellSize",cellSize);
//            //            SpatialCellObserverType::setCellSize(cellSize);
//            //            EDR.readScalarInFile(fullName.str(),"use_DisplacementMultipole",DislocationDisplacement<dim>::use_multipole);
//            //            EDR.readScalarInFile(fullName.str(),"use_StressMultipole",DislocationStress<dim>::use_multipole);
//            //            EDR.readScalarInFile(fullName.str(),"use_EnergyMultipole",DislocationEnergy<dim>::use_multipole);
//            
//            
//            //dt=0.0;
//            //            EDR.readScalarInFile(fullName.str(),"dxMax",DDtimeIntegrator<0>::dxMax);
//            //            assert(DDtimeIntegrator<0>::dxMax>0.0);
//            //            EDR.readScalarInFile(fullName.str(),"shearWaveSpeedFraction",shearWaveSpeedFraction);
//            //            assert(shearWaveSpeedFraction>=0.0);
//            
//            //            EDR.readScalarInFile(fullName.str(),"use_stochasticForce",DN.use_stochasticForce);
//            //            int stochasticForceSeed=-1;
//            //            EDR.readScalarInFile(fullName.str(),"stochasticForceSeed",stochasticForceSeed);
//            //            if(stochasticForceSeed<0)
//            //            {
//            //                StochasticForceGenerator::init(std::chrono::system_clock::now().time_since_epoch().count());
//            //            }
//            //            else
//            //            {
//            //                StochasticForceGenerator::init(stochasticForceSeed);
//            //            }
//            
//            
//            
//            // Use Changing external stress field induced by straight dislocations.
//            //            EDR.readScalarInFile("./loadInput.txt","use_externaldislocationstressfield",DN.use_externaldislocationstressfield);
//            //            if (DN.use_externaldislocationstressfield)
//            //            {
//            //                DN.ssdeq.clear();
//            //                typedef IDreader<'B',1,10,double> IDreaderType;
//            //                IDreaderType vReader;
//            //                if (vReader.isGood(0,true))
//            //                {
//            //                    vReader.read(0,true);
//            //                    for (const auto& vIter : vReader)
//            //                    {
//            //                        Eigen::Map<const Eigen::Matrix<double,1,9>> row(vIter.second.data());
//            //                        VectorDim P0(row.template segment<dim>(0));// P0 position
//            //                        VectorDim P1(row.template segment<dim>(dim)); // P1 position
//            //                        VectorDim B(row.template segment<dim>(dim*2));  // Burgers vector
//            //                        DN.ssdeq.emplace_back(StressStraight<dim>(P0,P1,B));
//            //                    }
//            //                }
//            //                else
//            //                {
//            //                    model::cout<<"could not read runID from B/B_0.txt"<<std::endl;
//            //                }
//            //            }
//            
//            //            // Restart
//            //            IDreader<'F',1,200,double> vReader;
//            //            vReader.readLabelsFile("F/F_labels.txt");
//            //            Eigen::Matrix<double,1,200> temp(Eigen::Matrix<double,1,200>::Zero());
//            //            
//            //            
//            //            if (vReader.isGood(0,true))
//            //            {
//            //                
//            //                vReader.read(0,true);
//            //                
//            //                if(runID<0)
//            //                {
//            //                    if(vReader.size())
//            //                    {
//            //                        runID=vReader.rbegin()->first;
//            //                        temp=Eigen::Map<Eigen::Matrix<double,1,200>>(vReader.rbegin()->second.data());
//            //                    }
//            //                    else
//            //                    {
//            //                        runID=0;
//            //                    }
//            //                }
//            //                else
//            //                {
//            //                    const auto iter=vReader.find(runID);
//            //                    if(iter!=vReader.end())
//            //                    {// runID has been found
//            //                        temp=Eigen::Map<Eigen::Matrix<double,1,200>>(iter->second.data());
//            //                    }
//            //                    else
//            //                    {
//            //                        assert(0 && "runID NOT FOUND IN F/F_0.txt");
//            //                    }
//            //                }
//            //                
//            //                DN.dt=temp(1);
//            //                
//            //            }
//            //            else
//            //            {
//            //                model::cout<<"could not read runID from F/F_0.txt"<<std::endl;
//            //                runID=0;
//            //                DN.dt=10.0;
//            //            }
//            //            
//            //            model::cout<<"dt="<<DN.dt<<std::endl;
//            //            
//            //            size_t curCol=0;
//            //            DN.totalTime=temp(curCol);
//            //            curCol+=2;
//            //            
//            //            if (DN.outputPlasticDistortion)
//            //            {
//            //                std::cout<<"reading PD"<<std::endl;
//            //                
//            //                for(int r=0;r<3;++r)
//            //                {
//            //                    for(int c=0;c<3;++c)
//            //                    {
//            //                        DN._plasticDistortionFromVelocities(r,c)=temp(curCol);
//            //                        curCol+=1;
//            //                    }
//            //                }
//            //            }
//            //            
//            //
//            //            model::cout<<"starting at time step "<<runID<<std::endl;
//            //            model::cout<<"totalTime= "<<DN.totalTime<<std::endl;
//            //            model::cout<<"plasticDistortionFromVelocities=\n "<<DN._plasticDistortionFromVelocities<<std::endl;
//            
//            // time-stepping
//            
//            //            EDR.readScalarInFile(fullName.str(),"outputDislocationStiffnessAndForce",DislocationNetworkComponentType::outputKF);
//            
//            
//            
//            //            EDR.readScalarInFile(fullName.str(),"timeWindow",timeWindow);
//            //            assert(timeWindow>=0.0 && "timeWindow MUST BE >= 0");
//            
//            // Check balance
//            //            EDR.readScalarInFile(fullName.str(),"check_balance",check_balance);
//            
//            
//            
//            
//            //            EDR.readScalarInFile(fullName.str(),"collisionTol",DislocationJunctionFormation<DislocationNetworkType>::collisionTol);
//            
//            
//            // Cross-Slip
//            //            EDR.readScalarInFile(fullName.str(),"crossSlipModel",DN.crossSlipModel);
//            //            if(DN.crossSlipModel)
//            //            {
//            //                DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg=TextFileParser("inputFiles/DD.txt").readScalar<double>("crossSlipDeg",true);
//            ////                EDR.readScalarInFile(fullName.str(),"crossSlipDeg",DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg);
//            //                assert(DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg>=0.0 && DislocationCrossSlip<DislocationNetworkType>::crossSlipDeg <= 90.0 && "YOU MUST CHOOSE 0.0<= crossSlipDeg <= 90.0");
//            //                //                EDR.readScalarInFile(fullName.str(),"crossSlipLength",DislocationCrossSlip<DislocationNetworkType>::crossSlipLength);
//            //                //                assert(DislocationCrossSlip<DislocationNetworkType>::crossSlipLength>=DislocationNetworkRemesh<DislocationNetworkType>::Lmin && "YOU MUST CHOOSE crossSlipLength>=Lmin.");
//            ////                EDR.readScalarInFile(fullName.str(),"verboseCrossSlip",DislocationCrossSlip<DislocationNetworkType>::verboseCrossSlip);
//            //                DislocationCrossSlip<DislocationNetworkType>::verboseCrossSlip=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseCrossSlip",true);
//            //            }
//            
//            // Mesh and BVP
//            //            if (DN.use_boundary)
//            //            {
//            //                EDR.readScalarInFile(fullName.str(),"surfaceAttractionDistance",DN.surfaceAttractionDistance);
//            
//            //                EDR.readScalarInFile(fullName.str(),"dislocationImages_x",DN.dislocationImages_x);
//            //                EDR.readScalarInFile(fullName.str(),"dislocationImages_y",DN.dislocationImages_y);
//            //                EDR.readScalarInFile(fullName.str(),"dislocationImages_z",DN.dislocationImages_z);
//            
//            //                EDR.readScalarInFile(fullName.str(),"use_meshRegions",use_meshRegions);
//            
//            //                int meshID(0);
//            //                EDR.readScalarInFile(fullName.str(),"meshID",meshID);
//            //                DN.mesh.readMesh(meshID);
//            //                assert(DN.mesh.simplices().size() && "MESH IS EMPTY.");
//            
//            // Initialize Polycrystal
//            //                DN.poly.init(DN,"./polyCrystalInput.txt");
//            
//            //                EDR.readScalarInFile(fullName.str(),"useVirtualExternalLoops",DN.useVirtualExternalLoops);
//            //                if(DN.useVirtualExternalLoops)
//            //                {
//            //                    EDR.readScalarInFile(fullName.str(),"virtualSegmentDistance",LinkType::virtualSegmentDistance);
//            //                }
//            
//            //                if(DN.use_bvp)
//            //                {
//            ////                    DN.bvpSolver->use_directSolver=TextFileParser("inputFiles/DD.txt").readScalar<int>("use_directSolver_FEM",true);
//            //                    DN.bvpSolver->tolerance=TextFileParser("inputFiles/DD.txt").readScalar<double>("solverTolerance",true);
//            ////                    
//            ////                    EDR.readScalarInFile(fullName.str(),"use_directSolver_FEM",DN.bvpSolver->use_directSolver);
//            ////                    EDR.readScalarInFile(fullName.str(),"solverTolerance",DN.bvpSolver->tolerance);
//            //                    DN.bvpSolver->init(DN);
//            //                }
//            //            }
//            //            else{ // no boundary is used, DislocationNetwork is in inifinite medium
//            //                DN.use_bvp=0;    // never comupute boundary correction
//            //            }
//            
//            //            DN.externalLoadController.init(DN,runID);  // have to initialize it after mesh!
//            
//            
//            
//            // VERTEX REDISTRIBUTION
////            DislocationNetworkRemesh<DislocationNetworkType>::remeshFrequency=TextFileParser("inputFiles/DD.txt").readScalar<int>("remeshFrequency",true);
//            //            EDR.readScalarInFile(fullName.str(),"remeshFrequency",DislocationNetworkRemesh<DislocationNetworkType>::remeshFrequency);
////            double Lmin=TextFileParser("inputFiles/DD.txt").readScalar<double>("Lmin",true);
//            //            EDR.readScalarInFile(fullName.str(),"Lmin",Lmin);
////            double Lmax=TextFileParser("inputFiles/DD.txt").readScalar<double>("Lmax",true);
////            double nodeRemoveAngleDeg=TextFileParser("inputFiles/DD.txt").readScalar<double>("nodeRemoveAngleDeg",true);
//            
//            //            EDR.readScalarInFile(fullName.str(),"Lmax",Lmax);
//            //            if(DN.use_boundary)
//            //            {
////            const double minMeshSize=std::min(DN.mesh.xMax(0)-DN.mesh.xMin(0),std::min(DN.mesh.xMax(1)-DN.mesh.xMin(1),DN.mesh.xMax(2)-DN.mesh.xMin(2)));
////            assert(Lmax<1.0 && "IF USING A BOUNDARY Lmax MUST BE RELATIVE TO BOX SIZE (Lmax<1)");
////            assert(Lmin<=Lmax);
////            DislocationNetworkRemesh<DislocationNetworkType>::Lmax=Lmax*minMeshSize;
////            DislocationNetworkRemesh<DislocationNetworkType>::Lmin=Lmin*minMeshSize;
////            DislocationNetworkRemesh<DislocationNetworkType>::cosRemove=cos(nodeRemoveAngleDeg*M_PI/180.0);
////            //
////            //            }
////            //            else
////            //            {
////            //                DislocationNetworkRemesh<DislocationNetworkType>::Lmax=Lmax;
////            //                DislocationNetworkRemesh<DislocationNetworkType>::Lmin=Lmin;
////            //            }
////            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmax>3.0*DislocationNetworkRemesh<DislocationNetworkType>::Lmin);
////            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=0.0);
////            assert(DislocationNetworkRemesh<DislocationNetworkType>::Lmin>=2.0*DDtimeIntegrator<0>::dxMax && "YOU MUST CHOOSE Lmin>2*dxMax.");
//            
//            
//            
//            // Verbose levels
////            if(DN.maxJunctionIterations>0)
////            {
////                //                EDR.readScalarInFile(fullName.str(),"verboseJunctions",DislocationJunctionFormation<DislocationNetworkType>::verboseJunctions);
////                DislocationJunctionFormation<DislocationNetworkType>::verboseJunctions=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseJunctions",true);
////            }
//            
//            //            EDR.readScalarInFile(fullName.str(),"verboseNodeContraction",DislocationNodeContraction<DislocationNetworkType>::verboseNodeContraction);
//            //            EDR.readScalarInFile(fullName.str(),"verboseDislocationNode",NodeType::verboseDislocationNode);
////            DislocationNodeContraction<DislocationNetworkType>::verboseNodeContraction=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseNodeContraction",true);
//            //NodeType::verboseDislocationNode=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseDislocationNode",true);
//            //            LinkType::verboseDislocationSegment=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseDislocationSegment",true);
//            
//            // GrainBoundary model
//            //            EDR.readScalarInFile(fullName.str(),"grainBoundaryTransmissionModel",GrainBoundaryTransmission<DislocationNetworkType>::grainBoundaryTransmissionModel);
////            GrainBoundaryTransmission<DislocationNetworkType>::grainBoundaryTransmissionModel=TextFileParser("inputFiles/DD.txt").readScalar<int>("grainBoundaryTransmissionModel",true);
//            
//            //            EDR.readScalarInFile(fullName.str(),"outputSegmentPairDistances",DN.outputSegmentPairDistances);
//            
//            
//            
//            // Grain Boundary flags
//            //            EDR.readScalarInFile(fullName.str(),"use_GBdissociation",GrainBoundaryDissociation<DislocationNetworkType>::use_GBdissociation);
//            //            EDR.readScalarInFile(fullName.str(),"use_GBtransmission",GrainBoundaryTransmission<DislocationNetworkType>::use_GBtransmission);
//            //            EDR.readScalarInFile(fullName.str(),"use_GBdislocations",GrainBoundary<dim>::use_GBdislocations);
//            
//            // Read Vertex and Edge information
//            EVLio<dim> evl;
//            if(EVLio<dim>::isBinGood(runID,suffix))
//            {
//                evl.readBin(runID,suffix);
//            }
//            else
//            {
//                if(EVLio<dim>::isTxtGood(runID,suffix))
//                {
//                    evl.readTxt(runID,suffix);
//                }
//                else
//                {
//                    std::cout<<"COULD NOT FIND INPUT FILEs evl/evl_"<<runID<<".bin or evl/evl_"<<runID<<".txt"<<std::endl;
//                    assert(0 && "COULD NOT FIND INPUT FILEs.");
//                }
//            }
////            createVertices(evl);
////            createEdges(evl);
////            DN.updatePlasticDistortionFromAreas(DN.simulationParameters.dt);
//            
////            IDreader<'E',1,14,double> inclusionsReader;
////            inclusionsReader.read(0,true);
////
////            const std::vector<double> inclusionsMobilityReduction(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsMobilityReduction",true));
////            for(const auto& pair : inclusionsReader)
////            {
////
////                const size_t& inclusionID(pair.first);
////                Eigen::Map<const Eigen::Matrix<double,1,14>> row(pair.second.data());
////
////                const VectorDim C(row.template segment<dim>(0));
////                const double a(row(dim+0));
////                MatrixDimD eT(MatrixDimD::Zero());
////                const int typeID(row(13));
////                int k=dim+1;
////                for(int i=0;i<dim;++i)
////                {
////                    for(int j=0;j<dim;++j)
////                    {
////                        eT(i,j)=row(k);
////                        k++;
////                    }
////                }
////
////
////
////                EshelbyInclusion<dim>::set_count(inclusionID);
////                DN.eshelbyInclusions().emplace(std::piecewise_construct,
////                                               std::make_tuple(inclusionID),
////                                               std::make_tuple(C,a,eT,DN.poly.nu,DN.poly.mu,inclusionsMobilityReduction[typeID],typeID) );
////            }
//            
//            //            readVertices(runID); // this requires mesh to be up-to-date
//            //            readEdges(runID);    // this requires mesh to be up-to-date
//            //#ifdef DislocationNucleationFile
//            //            EDR.readScalarInFile(fullName.str(),"nucleationFreq",nucleationFreq);
//            //#endif
//            
////#ifdef _MODEL_MPI_
////            // Avoid that a processor starts writing before other are reading
////            MPI_Barrier(MPI_COMM_WORLD);
////#endif
//            
//            // Initializing configuration
//            //            DN.move(0.0);    // initial configuration
//        }
        
//        /* readVertices *******************************************************/
//        void createVertices(const EVLio<dim>& evl)
//        {/*!Creates DislocationNode(s) based on the data read by the EVLio<dim>
//          * object.
//          */
//            size_t kk(1);
//            for (const auto& node : evl.nodes())
//            {
//                const size_t nodeIDinFile(node.sID);
//                NodeType::set_count(nodeIDinFile);
//                if(node.sID==node.masterID)
//                {// a regular node is created
//                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<")"<<std::endl;
//                    const size_t nodeID=DN.insertDanglingNode(node.P,node.V,node.velocityReduction).first->first;
//                    assert(nodeID==nodeIDinFile);
//                }
//                else
//                {
//                    model::cout<<"Creating DislocationNode "<<nodeIDinFile<<" ("<<kk<<" of "<<evl.nodes().size()<<"), virtual of "<<node.masterID<<std::endl;
//                    const auto isNode(DN.node(node.masterID));
//                    assert(isNode.first);
//                    isNode.second->resetVirtualBoundaryNode();
//                }
//                kk++;
//            }
//        }
//
//        /* readEdges **********************************************************/
//        void createEdges(const EVLio<dim>& evl)
//        {/*!
//          */
//
//
//
//            std::map<size_t,std::map<size_t,size_t>> loopMap;
//            for(const auto& looplink : evl.links())
//            {
//                loopMap[looplink.loopID].emplace(looplink.sourceID,looplink.sinkID);
//            }
//
//
//
//            assert(loopMap.size()==evl.loops().size());
//
//            size_t loopLumber=1;
//            for(const auto& loop : evl.loops())
//            {// for each loop in the EVLio<dim> object
//
//                const auto loopFound=loopMap.find(loop.sID); // there must be an entry with key loopID in loopMap
//                assert(loopFound!=loopMap.end());
//                std::vector<size_t> nodeIDs;
//                nodeIDs.push_back(loopFound->second.begin()->first);
//                for(size_t k=0;k<loopFound->second.size();++k)
//                {
//                    const auto nodeFound=loopFound->second.find(*nodeIDs.rbegin());
//                    if(k<loopFound->second.size()-1)
//                    {
//                        nodeIDs.push_back(nodeFound->second);
//                    }
//                    else
//                    {
//                        assert(nodeFound->second==nodeIDs[0]);
//                    }
//                }
//
//                LoopType::set_count(loop.sID);
//
//                if(loop.N.squaredNorm()>FLT_EPSILON)
//                {
//                    model::cout<<"Creating Dislocation Loop "<<loop.sID<<" ("<<loopLumber<<" of "<<evl.loops().size()<<")"<<std::endl;
//                    const size_t newLoopID=DN.insertLoop(nodeIDs,loop.B,loop.N,loop.P,loop.grainID)->sID;
//                    assert(loop.sID==newLoopID);
//                }
//                else
//                {
//                    model::cout<<"Creating Dislocation Loop "<<loop.sID<<" ("<<loopLumber<<" of "<<evl.loops().size()<<") virtual"<<std::endl;
//                    std::vector<std::shared_ptr<NodeType>> sharedNodes;
//                    for(const size_t nodeID : nodeIDs)
//                    {// collect shared_ptrs to nodes
//                        const auto isNode(DN.node(nodeID));
//                        assert(isNode.first);
//                        if(isNode.second->masterNode)
//                        {// a virtual node
//                            sharedNodes.push_back(isNode.second->masterNode->virtualBoundaryNode());
//                        }
//                        else
//                        {
//                            const auto isSharedNode(DN.danglingNode(nodeID));
//                            if(!isSharedNode.first)
//                            {
//                                model::cout<<"node "<<nodeID<<" not found"<<std::endl;
//                                assert(false && "node shared pointer not found");
//                            }
//                            sharedNodes.push_back(isSharedNode.second);
//                        }
//                    }
//                    const size_t newLoopID=DN.insertLoop(sharedNodes,loop.B,loop.grainID)->sID;
//                    assert(loop.sID==newLoopID);
//                }
//                loopLumber++;
//            }
//
//            DN.clearDanglingNodes();
//
//
//
//            model::cout<<std::endl;
//        }
        
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
            model::cout<<"		Writing to "<<std::flush;
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
                
                assert(0 && "FINISH BINARY OUTPUT OF QUADRATURE POINTS");
                //                typedef typename DislocationNetworkType::DislocationParticleType::ElasticEnergy ElasticEnergy;
                //                SequentialOutputFile<'W',1>::set_count(runID);
                //                SequentialOutputFile<'W',1>::set_increment(DN.outputFrequency);
                //                SequentialOutputFile<'W',1> w_file; //energy_file
                //                int ll=0;
                //                for (const auto& linkIter : DN.links())
                //                {
                //                    const int qOrder(linkIter.second->rgauss.cols());
                //                    for (size_t q=0;q<linkIter.second->quadratureParticleContainer.size();++q)
                //                    {
                //                        w_file << ll*qOrder+q<<" "<< linkIter.second->rgauss.col(q).transpose()<<" "<< linkIter.second->quadratureParticleContainer[q]->template field<ElasticEnergy>()<<"\n";
                //                    }
                //                    ll++;
                //                }
                //                model::cout<<", W/W_"<<w_file.sID<<std::flush;
            }
            
            //            typedef BoundaryDisplacementPoint<DislocationNetworkType> FieldPointType;
            //            typedef typename FieldPointType::DisplacementField DisplacementField;
            
            if(DN.outputMeshDisplacement)
            {
                
                const auto t0=std::chrono::system_clock::now();
                model::SequentialOutputFile<'D',1>::set_count(runID); // Vertices_file;
                model::SequentialOutputFile<'D',1>::set_increment(DN.outputFrequency); // Vertices_file;
                model::SequentialOutputFile<'D',true> d_file;
                model::cout<<"		writing to D/D_"<<d_file.sID<<std::flush;
                
                std::vector<FEMnodeEvaluation<ElementType,dim,1>> fieldPoints; // the container of field points
                fieldPoints.reserve(DN.mesh.template observer<0>().size());
                for (const auto& sIter : DN.mesh.template observer<0>())
                {
                    if(sIter.second->isBoundarySimplex())
                    {
                        fieldPoints.emplace_back(sIter.second->xID(0),sIter.second->P0);
                    }
                }
                
                DN.displacement(fieldPoints);
                
                for(auto& node : fieldPoints)
                {// add FEM solution and output
                    if(DN.simulationParameters.simulationType==DefectiveCrystalParameters::FINITE_FEM)
                    {
                        const size_t femID=DN.bvpSolver->finiteElement().mesh2femIDmap().at(node.pointID)->gID;
                        node+=DN.bvpSolver->displacement().dofs(femID);
                    }
                    d_file<<node.pointID<<" "<<node.transpose()<<"\n";
                }
                model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
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
            
            if (DN.bvpSolver && DN.outputFEMsolution )
            {
                if(!(runID%DN.bvpSolver->stepsBetweenBVPupdates))
                {// Output displacement and stress on external mesh faces
                    const auto t0=std::chrono::system_clock::now();
                    model::SequentialOutputFile<'U',1>::set_count(runID); // Vertices_file;
                    model::SequentialOutputFile<'U',1>::set_increment(DN.outputFrequency); // Vertices_file;
                    model::SequentialOutputFile<'U',true> u_file;
                    model::cout<<"		writing to U/U_"<<u_file.sID<<".txt"<<std::flush;
                    u_file<<DN.bvpSolver->displacement().onBoundary();
                    model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
                    
                    const auto t1=std::chrono::system_clock::now();
                    model::SequentialOutputFile<'S',1>::set_count(runID); // Vertices_file;
                    model::SequentialOutputFile<'S',1>::set_increment(DN.outputFrequency); // Vertices_file;
                    model::SequentialOutputFile<'S',true> s_file;
                    model::cout<<"		writing to S/S_"<<s_file.sID<<".txt"<<std::flush;
                    s_file<<DN.bvpSolver->stress().onBoundary();
                    model::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t1)).count()<<" sec]"<<defaultColor<<std::endl;
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

