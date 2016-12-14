/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Tamer Crosby     <tcrosby@ucla.edu>.
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>,
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_TRANSMITSEGMENT_H
#define model_TRANSMITSEGMENT_H

#include <math.h> // isfinite
#include <float.h>
#include <list>
#include <stdlib.h> // rand()
#include <ctime>
#include <time.h>
#include <iterator>
#include <vector>
#include <set>
#include <map>
#include <memory> // unique_ptr
#include <tuple>
#include <model/DislocationDynamics/DislocationSharedObjects.h>
#include <model/DislocationDynamics/Polycrystals/GrainBoundary.h>
#include <model/DislocationDynamics/Polycrystals/Grain.h>
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/LatticeMath/LatticeMath.h>


namespace model
{
    
    template <typename DislocationSegmentType>
    struct TransmitSegment
    {
        
        static constexpr int dim=DislocationSegmentType::dim;
        typedef typename DislocationSegmentType::VectorDim VectorDim;
        typedef LatticeVector<DislocationSegmentType::dim> LatticeVectorType;
        typedef typename DislocationSegmentType::NodeType NodeType;
        typedef std::map<double,std::pair<SlipSystem,const GrainBoundary<dim>* const>> MapType;
        
//        DislocationSharedObjects<dim> shared;
        const size_t sourceID;
        const size_t sinkID;
        bool isValidTransmission;
        VectorDim transmitSourceP;
        VectorDim transmitSinkP;
        const Grain<dim>* transmitGrain;
        std::unique_ptr<LatticeVectorType> transmitBurgers;
        std::unique_ptr<LatticeVectorType> transmitMidpoint;

//        std::unique_ptr<LatticeVectorType> residualBurgers;
        
        /**********************************************************************/
        TransmitSegment(const DislocationSegmentType& ds) :
        /* init list */ sourceID(ds.source->sID),
        /* init list */ sinkID(ds.sink->sID),
//        /* init list */ grain(ds.grain),
        /* init list */ isValidTransmission(false),
        /* init list */ transmitSourceP(ds.source->get_P()),
        /* init list */ transmitSinkP(ds.sink->get_P())
//        /* init list */ transmitMidpoint(ds.source->get_L())
//        /* init list */ transmitBurgers(ds.flow),
//        /* init list */ residualBurgers(transmitBurgers)
        {
            
            const VectorDim chord(ds.sink->get_P()-ds.source->get_P());
 //           const double chorNorm=chord.norm();
            
            MapType slipSystemMap;

            // Nucleation into other grains
            for(const auto& gb : ds.grainBoundarySet)
            {

                const size_t otherGrainID=(gb->grainBndID.first==ds.grain.grainID)? gb->grainBndID.second : gb->grainBndID.first;
                for(const auto& slipSystem : gb->grain(otherGrainID).slipSystems())
                {
                    //                const std::pair<LatticeVectorType,LatticeVectorType>& primitiveVectors(gb->latticePlane(ds.grain.grainID).n.primitiveVectors);
                    if(fabs(slipSystem.n.cartesian().dot(chord))<FLT_EPSILON) // slip system contains chord
                    {
                        const VectorDim pkForce=(ds.midPointStress()*slipSystem.s.cartesian()).cross(-chord);
                        if(pkForce.dot(gb->latticePlane(otherGrainID).n.cartesian())<0.0) // PK force points inside the grain
                        {
                            slipSystemMap.emplace(ds.Burgers.dot(slipSystem.s.cartesian())+slipSystem.s.cartesian().squaredNorm(),std::make_pair(slipSystem,gb));
                        }
                    }
                }

                
            }
            
            if(slipSystemMap.size())
            {
//                if(slipSystemMap.begin()->first<0.0)
                    if(true)
                {
                    const GrainBoundary<dim>* const gb=slipSystemMap.begin()->second.second;
                    const size_t otherGrainID=(gb->grainBndID.first==ds.grain.grainID)? gb->grainBndID.second : gb->grainBndID.first;
                    isValidTransmission=true;
                    transmitGrain=&(gb->grain(otherGrainID));
                    transmitBurgers.reset(new LatticeVectorType(slipSystemMap.begin()->second.first.s));
                    transmitSourceP=gb->latticePlane(otherGrainID).snapToLattice(transmitSourceP).cartesian();
                    transmitSinkP=gb->latticePlane(otherGrainID).snapToLattice(transmitSinkP).cartesian();
                    
                    LatticeVectorType P1(transmitGrain->snapToLattice(ds.source->get_P()));
//                    std::cout<<"Here 0.25"<<std::endl;
//                    std::cout<<"ds.source->get_P()="<<std::setprecision(15)<<std::scientific<<ds.source->get_P().transpose()<<std::endl;
//                    std::cout<<"P1="<<std::setprecision(15)<<std::scientific<<P1.cartesian().transpose()<<std::endl;

//                    LatticeVectorType  P2(ds.source->get_P(),transmitGrain->covBasis(),transmitGrain->contraBasis());
//                    std::cout<<"Here 0.5"<<std::endl;

                    LatticePlane snapPlane(P1,(slipSystemMap.begin()->second).first.n);
                    
                    

                    VectorDim midD=0.5*(transmitSourceP+transmitSinkP)-chord.cross(slipSystemMap.begin()->second.first.n.cartesian()).normalized()*10.0;
                    auto temp=DislocationSharedObjects<dim>::mesh.searchWithGuess(midD,ds.source->includingSimplex());

                    if(temp.second->region->regionID!=otherGrainID)
                    {
                        midD=0.5*(transmitSourceP+transmitSinkP)+chord.cross(slipSystemMap.begin()->second.first.n.cartesian()).normalized()*10.0;
                    }
                    
                    temp=DislocationSharedObjects<dim>::mesh.searchWithGuess(midD,ds.source->includingSimplex());
                    assert(temp.first);
                    
                    transmitMidpoint.reset(new LatticeVectorType(snapPlane.snapToLattice(midD)));
                
                    transmitSourceP=transmitGrain->snapToLattice(transmitSourceP).cartesian();
                    transmitSinkP=transmitGrain->snapToLattice(transmitSinkP).cartesian();

                }

                
            }
            
        }
        
    };
    
} // namespace model
#endif



