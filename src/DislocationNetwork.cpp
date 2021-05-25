/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationNetwork_cpp_
#define model_DislocationNetwork_cpp_


#include <DislocationNetwork.h>

namespace model
{

    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    DislocationNetwork<dim,corder,InterpolationType>::DislocationNetwork(int& argc, char* argv[],
                       const DefectiveCrystalParameters& _simulationParameters,
                       const SimplicialMesh<dim>& _mesh,
                       const Polycrystal<dim>& _poly,
                       const std::unique_ptr<BVPsolver<dim,2>>& _bvpSolver,
                       const std::unique_ptr<ExternalLoadControllerBase<dim>>& _externalLoadController,
                       const std::vector<VectorDim>& _periodicShifts,
                       long int& runID) :
//    /* init */ LoopNetworkType(_simulationParameters.isPeriodicSimulation()? std::shared_ptr<NetworkComponentType>(new NetworkComponentType()) : nullptr)
    /* init */ simulationParameters(_simulationParameters)
    /* init */,mesh(_mesh)
    /* init */,poly(_poly)
    /* init */,glidePlaneFactory(poly)
//    /* init */,periodicDislocationLoopFactory(simulationParameters.isPeriodicSimulation()? new PeriodicDislocationLoopFactory<DislocationNetworkType>(poly,glidePlaneFactory) : nullptr)
    /* init */,bvpSolver(_bvpSolver)
    /* init */,externalLoadController(_externalLoadController)
    /* init */,periodicShifts(_periodicShifts)
//    /* init */,networkRemesher(*this)
//    /* init */,junctionsMaker(*this)
//    /* init */,nodeContractor(*this)
//    /* init */,gbTransmission(*this)
    //        /* init */,timeIntegrationMethod(TextFileParser("inputFiles/DD.txt").readScalar<int>("timeIntegrationMethod",true))
    ///* init */,maxJunctionIterations(TextFileParser("inputFiles/DD.txt").readScalar<int>("maxJunctionIterations",true))
    //        /* init */,runID(TextFileParser("inputFiles/DD.txt").readScalar<int>("startAtTimeStep",true)),
    //        /* init */,totalTime(0.0),
    //        /* init */ dt(0.0),
    //        /* init */ vMax(0.0),
    //        /* init */ Nsteps(TextFileParser("inputFiles/DD.txt").readScalar<size_t>("Nsteps",true)),
    //        /* init */,_plasticDistortionFromVelocities(MatrixDim::Zero())
    /* init */,_plasticDistortionFromAreas(std::make_pair(0.0,MatrixDim::Zero()))
    /* init */,_plasticDistortionRateFromVelocities(MatrixDim::Zero())
    /* init */,_plasticDistortionRateFromAreas(MatrixDim::Zero())
    /* init */,ddSolverType(TextFileParser("inputFiles/DD.txt").readScalar<int>("ddSolverType",true))
    /* init */,computeDDinteractions(TextFileParser("inputFiles/DD.txt").readScalar<int>("computeDDinteractions",true))
    /* init */,crossSlipModel(TextFileParser("inputFiles/DD.txt").readScalar<int>("crossSlipModel",true))
    /* init */,outputFrequency(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFrequency",true))
    /* init */,outputBinary(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputBinary",true))
    /* init */,outputGlidePlanes(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputGlidePlanes",true))
    /* init */,outputElasticEnergy(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputElasticEnergy",true))
    /* init */,outputMeshDisplacement(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputMeshDisplacement",true))
    /* init */,outputFEMsolution(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputFEMsolution",true))
    /* init */,outputDislocationLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputDislocationLength",true))
    //        /* init */,outputPlasticDistortion(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortion",true))
    /* init */,outputPlasticDistortionRate(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPlasticDistortionRate",true))
    /* init */,outputQuadraturePoints(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputQuadraturePoints",true))
    /* init */,outputLinkingNumbers(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLinkingNumbers",true))
    /* init */,outputLoopLength(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputLoopLength",true))
    /* init */,outputSegmentPairDistances(TextFileParser("inputFiles/DD.txt").readScalar<int>("outputSegmentPairDistances",true))
    /* init */,outputPeriodicConfiguration(simulationParameters.isPeriodicSimulation()? TextFileParser("inputFiles/DD.txt").readScalar<int>("outputPeriodicConfiguration",true) : false)
    //        /* init */ _userOutputColumn(3)
    /* init */,use_stochasticForce(TextFileParser("inputFiles/DD.txt").readScalar<int>("use_stochasticForce",true))
    /* init */,surfaceAttractionDistance(TextFileParser("inputFiles/DD.txt").readScalar<double>("surfaceAttractionDistance",true))
    //        /* init */,computePlasticDistortionRateFromVelocities(TextFileParser("inputFiles/DD.txt").readScalar<int>("computePlasticDistortionRateFromVelocities",true))
    /* init */,folderSuffix("")
    {
        
        // Some sanity checks
        //            assert(Nsteps>=0 && "Nsteps MUST BE >= 0");
        
        // Initialize static variables
        LoopNetworkType::verboseLevel=TextFileParser("inputFiles/DD.txt").readScalar<int>("verboseLoopNetwork",true);
        LoopType::initFromFile("inputFiles/DD.txt");
        LoopNodeType::initFromFile("inputFiles/DD.txt");
        LoopLinkType::initFromFile("inputFiles/DD.txt");
        NetworkLinkType::initFromFile("inputFiles/DD.txt");
        NetworkNodeType::initFromFile("inputFiles/DD.txt");
//        PeriodicDislocationBase::initFromFile("inputFiles/DD.txt");
//        DislocationNetworkComponentType::initFromFile("inputFiles/DD.txt");
        DislocationStressBase<dim>::initFromFile("inputFiles/DD.txt");
//        DDtimeIntegrator<0>::initFromFile("inputFiles/DD.txt");
//        DislocationCrossSlip<DislocationNetworkType>::initFromFile("inputFiles/DD.txt");
        int stochasticForceSeed=TextFileParser("inputFiles/DD.txt").readScalar<int>("stochasticForceSeed",true);
        if(stochasticForceSeed<0)
        {
            StochasticForceGenerator::init(std::chrono::system_clock::now().time_since_epoch().count());
        }
        else
        {
            StochasticForceGenerator::init(stochasticForceSeed);
        }
        
        if(argc>1)
        {
            folderSuffix=argv[1];
            //                std::cout<<"folderSuffix="<<folderSuffix<<std::endl;
        }

//        DDconfigIO<dim> evl(folderSuffix);
//        evl.read(runID);
//        setConfiguration(evl);
        createEshelbyInclusions();
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationNetwork<dim,corder,InterpolationType>::EshelbyInclusionContainerType& DislocationNetwork<dim,corder,InterpolationType>::eshelbyInclusions() const
    {
        return *this;
    }
    
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNetwork<dim,corder,InterpolationType>::EshelbyInclusionContainerType& DislocationNetwork<dim,corder,InterpolationType>::eshelbyInclusions()
    {
        return *this;
    }

    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    void DislocationNetwork<dim,corder,InterpolationType>::createEshelbyInclusions()
    {
        for(const auto& grain : poly.grains())
        {
            EshelbyInclusion<dim>::addSlipSystems(grain.second.slipSystems());
        }
        
        
        IDreader<'E',1,14,double> inclusionsReader;
        inclusionsReader.read(0,true);
        
        const std::vector<double> inclusionsMobilityReduction(TextFileParser("./inputFiles/initialMicrostructure.txt").readArray<double>("inclusionsMobilityReduction",true));
        for(const auto& pair : inclusionsReader)
        {
            
            const size_t& inclusionID(pair.first);
            Eigen::Map<const Eigen::Matrix<double,1,14>> row(pair.second.data());
            
            const VectorDim C(row.template segment<dim>(0));
            const double a(row(dim+0));
            MatrixDim eT(MatrixDim::Zero());
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
            eshelbyInclusions().emplace(std::piecewise_construct,
                                        std::make_tuple(inclusionID),
                                        std::make_tuple(C,a,eT,poly.nu,poly.mu,inclusionsMobilityReduction[typeID],typeID) );
        }
    }

    

    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNetwork<dim,corder,InterpolationType>::MatrixDim DislocationNetwork<dim,corder,InterpolationType>::plasticDistortion() const
    {
        MatrixDim temp(MatrixDim::Zero());
        assert(0 && "FINISH HERE");
//        for(const auto& loop : this->loops())
//        {
//            temp+= loop.second->plasticDistortion();
//        }
        return temp;
    }
    
    //        /**********************************************************************/
    //        const MatrixDimD& plasticDistortion() const
    //        {
    //            return computePlasticDistortionRateFromVelocities? _plasticDistortionFromVelocities : _plasticDistortionFromAreas;
    //        }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    const typename DislocationNetwork<dim,corder,InterpolationType>::MatrixDim& DislocationNetwork<dim,corder,InterpolationType>::plasticDistortionRate() const
    {
        return  _plasticDistortionRateFromAreas;
    }
    
    /**********************************************************************/
    template <int dim, short unsigned int corder, typename InterpolationType>
    typename DislocationNetwork<dim,corder,InterpolationType>::MatrixDim DislocationNetwork<dim,corder,InterpolationType>::plasticStrainRate() const
    {/*!\returns the plastic strain rate tensor generated during the last time step.
      */
        //const MatrixDimD temp(plasticDistortionRate());
        return (plasticDistortionRate()+plasticDistortionRate().transpose())*0.5;
    }
    
}
#endif
