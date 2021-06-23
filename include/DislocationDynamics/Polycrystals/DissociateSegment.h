/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 * Copyright (C) 2016 by Nathaniel Burbery <nbur049@aucklanduni.ac.nz>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DissociateSegment_H
#define model_DissociateSegment_H

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
#include <DislocationSharedObjects.h>
#include <GrainBoundary.h>
#include <Grain.h>
#include <Material.h>
#include <LatticeModule.h>


namespace model
{
    
    template <typename DislocationSegmentType>
    struct DissociateSegment
    {
        
        static constexpr int dim=DislocationSegmentType::dim;
        static constexpr double dissociateDistance=5.0;

        typedef typename DislocationSegmentType::VectorDim VectorDim;
        typedef LatticeVector<DislocationSegmentType::dim> LatticeVectorType;
        
        const size_t sourceID;
        const size_t sinkID;
        const Grain<dim>& grain;
        bool isValidDissociation;
        LatticeVectorType dissociateMidpoint;
        LatticeVectorType dissociateBurgers;
        LatticeVectorType residualBurgers;
        
        /**********************************************************************/
        DissociateSegment(const DislocationSegmentType& ds) :
        /* init list */ sourceID(ds.source->sID),
        /* init list */ sinkID(ds.sink->sID),
        /* init list */ grain(ds.grain),
        /* init list */ isValidDissociation(false),
        /* init list */ dissociateMidpoint(ds.source->get_L()),
        /* init list */ dissociateBurgers(ds.flow),
        /* init list */ residualBurgers(dissociateBurgers)
        {
            
            assert(ds.grainBoundarySet.size() && "grainBoundarySet has zero size.");
            
            const VectorDim chord(ds.sink->get_P()-ds.source->get_P());
 //           const double chorNorm=chord.norm();
            
            std::map<double,std::pair<LatticeVectorType,const GrainBoundary<dim>* const>> primitiveMap;

            
            for(const auto& gb : ds.grainBoundarySet)
            {
                const std::pair<LatticeVectorType,LatticeVectorType>& primitiveVectors(gb->latticePlane(ds.grain.grainID).n.primitiveVectors);
                
                primitiveMap.emplace(ds.Burgers.dot(primitiveVectors.first.cartesian())+primitiveVectors.first.cartesian().squaredNorm(),std::make_pair(primitiveVectors.first,gb));
                primitiveMap.emplace(ds.Burgers.dot(-primitiveVectors.first.cartesian())+primitiveVectors.first.cartesian().squaredNorm(),std::make_pair(primitiveVectors.first*-1,gb));
                primitiveMap.emplace(ds.Burgers.dot(primitiveVectors.second.cartesian())+primitiveVectors.second.cartesian().squaredNorm(),std::make_pair(primitiveVectors.second,gb));
                primitiveMap.emplace(ds.Burgers.dot(-primitiveVectors.second.cartesian())+primitiveVectors.second.cartesian().squaredNorm(),std::make_pair(primitiveVectors.second*-1,gb));

            }
            
            
            dissociateBurgers=primitiveMap.begin()->second.first;
            residualBurgers=dissociateBurgers+ds.flow;
            
            
            if(primitiveMap.begin()->first<0.0 && residualBurgers.squaredNorm()>0)
            {
                isValidDissociation=true;
                VectorDim dir=chord.cross(primitiveMap.begin()->second.second->latticePlane(ds.grain.grainID).n.cartesian()).normalized();
                const VectorDim pkForce=(ds.midPointStress()*dissociateBurgers.cartesian()).cross(-chord);
                if(dir.dot(pkForce)<0.0)
                {
                    dir*=-1.0;
                }
                dissociateMidpoint=primitiveMap.begin()->second.second->latticePlane(ds.grain.grainID).snapToLattice(0.5*(ds.source->get_P()+ds.sink->get_P())+dir*dissociateDistance);
            }
            
        }
        
    };
    
} // namespace model
#endif



