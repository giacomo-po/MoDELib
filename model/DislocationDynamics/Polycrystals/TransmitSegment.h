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
#include <model/DislocationDynamics/Materials/Material.h>
#include <model/LatticeMath/LatticeMath.h>


namespace model 
{
    
    template <typename DislocationSegmentType>
	struct TransmitSegment
    {
		
    		typedef typename DislocationSegmentType::VectorDim VectorDim;
        typedef LatticeVector<DislocationSegmentType::dim> LatticeVectorType;
        typedef typename DislocationSegmentType::NodeType NodeType;
        
        DislocationSharedObjects<3> shared;
        NodeType& source;
        NodeType& sink;
        int transmitRID;
        const Grain<3>& grain;
        std::pair<bool,GrainBoundary<3>> grainBoundary;
        const std::vector<LatticePlaneBase> planeNormals;
        const std::vector<SlipSystem> slipSystems; 
        const  LatticeVectorType Burgers;
        LatticeVectorType midPointO;
        LatticeVectorType midPoint;
        const  VectorDim pkForce;
		    const  Eigen::Matrix<double, 3, 3> midPointStress;
		    
        //LatticePlane gbPlane;
        
				/* Constructor *******************************************************/		
				TransmitSegment(DislocationSegmentType& ds) :
				/*init list   */ source(*ds.source),
				/*init list   */ sink(*ds.sink),
				/*init list   */ transmitRID(sink.grainBoundary_rID2),
				/*init list   */ grain(shared.poly.grain(transmitRID)),
				/*init list   */ grainBoundary(transmitRID>0&&source.grainBoundary_rID2>0 ? std::make_pair(true,shared.poly.grainBoundary(source.includingSimplex()->region->regionID,transmitRID)) :  std::make_pair(false,shared.poly.grainBoundary(source.sID,sink.sID))),
				/*init list   */ planeNormals(grainBoundary.second.grain(transmitRID).planeNormals()),
				/*init list   */ slipSystems(grainBoundary.second.grain(transmitRID).slipSystems()),
				/*init list   */ Burgers(ds.flow),  //This is not needed, but is useful for debug check! 
				/*init list   */ midPointO(grainBoundary.second.grain(source.grain.grainID).snapToLattice(ds.get_r(0.5))),
				/*init list   */ midPoint(grainBoundary.second.grain(transmitRID).snapToLattice(ds.get_r(0.5))),
				/*init list   */ pkForce(ds.integralPK()),
				/*init list   */ midPointStress(ds.midPointStress())
        {
            assert(transmitRID!=source.grain.grainID);
        }	
        
    };
    
} // namespace model
#endif



