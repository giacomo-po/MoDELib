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
        
        DislocationSharedObjects<dim> shared;
        const size_t sourceID;
        const size_t sinkID;
        const Grain<dim>& grain;
        bool isValidTransmission;
//        VectorDim transmitSourceP;
//        VectorDim transmitSinkP;
        LatticeVectorType transmitMidpoint;
        LatticeVectorType transmitBurgers;
        LatticeVectorType residualBurgers;
        
        /**********************************************************************/
        TransmitSegment(const DislocationSegmentType& ds) :
        /* init list */ sourceID(ds.source->sID),
        /* init list */ sinkID(ds.sink->sID),
        /* init list */ grain(ds.grain),
        /* init list */ isValidTransmission(false),
//        /* init list */ transmitSourceP(ds.source->get_P()),
//        /* init list */ transmitSinkP(ds.sink->get_P()),
        /* init list */ transmitMidpoint(ds.source->get_L()),
        /* init list */ transmitBurgers(ds.flow),
        /* init list */ residualBurgers(transmitBurgers)
        {
            
            const VectorDim chord(ds.sink->get_P()-ds.source->get_P());
 //           const double chorNorm=chord.norm();
            
            std::map<double,std::pair<LatticeVectorType,const GrainBoundary<dim>* const>> primitiveMap;

            for(const auto& gb : ds.grainBoundarySet)
            {
                const std::pair<LatticeVectorType,LatticeVectorType>& primitiveVectors(gb->latticePlane(ds.grain.grainID).n.primitiveVectors);
                
                
//                primitiveMap.emplace(ds.Burgers.dot(primitiveVectors.first.cartesian()),std::make_pair(primitiveVectors.first,gb));
//                primitiveMap.emplace(ds.Burgers.dot(-primitiveVectors.first.cartesian()),std::make_pair(primitiveVectors.first*-1,gb));
//                primitiveMap.emplace(ds.Burgers.dot(primitiveVectors.second.cartesian()),std::make_pair(primitiveVectors.second,gb));
//                primitiveMap.emplace(ds.Burgers.dot(-primitiveVectors.second.cartesian()),std::make_pair(primitiveVectors.second*-1,gb));
                
                primitiveMap.emplace(ds.Burgers.dot(primitiveVectors.first.cartesian())+primitiveVectors.first.cartesian().squaredNorm(),std::make_pair(primitiveVectors.first,gb));
                primitiveMap.emplace(ds.Burgers.dot(-primitiveVectors.first.cartesian())+primitiveVectors.first.cartesian().squaredNorm(),std::make_pair(primitiveVectors.first*-1,gb));
                primitiveMap.emplace(ds.Burgers.dot(primitiveVectors.second.cartesian())+primitiveVectors.second.cartesian().squaredNorm(),std::make_pair(primitiveVectors.second,gb));
                primitiveMap.emplace(ds.Burgers.dot(-primitiveVectors.second.cartesian())+primitiveVectors.second.cartesian().squaredNorm(),std::make_pair(primitiveVectors.second*-1,gb));

                
            }
            
            transmitBurgers=primitiveMap.begin()->second.first;
            residualBurgers=transmitBurgers+ds.flow;
            if(primitiveMap.begin()->first<0.0 && residualBurgers.squaredNorm()>0)
            {
                isValidTransmission=true;
                VectorDim dir=chord.cross(primitiveMap.begin()->second.second->latticePlane(ds.grain.grainID).n.cartesian()).normalized();
                VectorDim pkForce=(ds.midPointStress()*transmitBurgers.cartesian()).cross(-chord);
                if(dir.dot(pkForce)<0.0)
                {
                    dir*=-1.0;
                }
                transmitMidpoint=primitiveMap.begin()->second.second->latticePlane(ds.grain.grainID).snapToLattice(0.5*(ds.source->get_P()+ds.sink->get_P())+dir*10.0);
            }
            
        }
        
    };
    
} // namespace model
#endif



