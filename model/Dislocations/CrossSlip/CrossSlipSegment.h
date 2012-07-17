/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Can Erel         <canerel55@gmail.com>,
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_CROSSSLIPSEGMENT_H
#define model_CROSSSLIPSEGMENT_H

#include <math.h>
#include <float.h>

namespace model {


	template <typename DislocationSegmentType>
	struct CrossSlipSegment {
		
		typedef typename DislocationSegmentType::VectorDim VectorDim;

		/*const*/ size_t sourceID;
		/*const*/ size_t   sinkID;
		/*const*/ VectorDim chord;
		/*const*/ VectorDim midPoint;
		/*const*/ VectorDim Burgers;
		/*const*/ bool sourceOnMeshBoundary;
		/*const*/ bool sinkOnMeshBoundary;
		/*const*/ VectorDim pkForce;
		/*const*/ VectorDim normalPrimary;
        /*const*/ bool isSessile;
		/*const*/ VectorDim normalConjugate;
		/*const*/ double rssPrimary;
		/*const*/ double rssConjugate;
		/*const*/ double crossSlipFactor;
		/*const*/ bool isCrossSlipSegment;
		
		/* Constructor *******************************************************/		
		CrossSlipSegment(const DislocationSegmentType& ds, const double& sinThetaCrossSlipCr,const double& crossSlipLength) : 
		/*init list   */ sourceID(ds.source->sID),
		/*init list   */   sinkID(ds.sink->sID),
		/*init list   */ chord(ds.chord()),
		/*init list   */ midPoint(ds.get_r(0.5)),
		/*init list   */ Burgers(ds.Burgers),
		/*init list   */ sourceOnMeshBoundary(ds.source->nodeMeshLocation == onMeshBoundary),
		/*init list   */   sinkOnMeshBoundary(ds.sink  ->nodeMeshLocation == onMeshBoundary),
		/*init list   */ pkForce(ds.integralPK()),
		/*init list   */ normalPrimary(ds.glidePlaneNormal),
        /*init list   */ isSessile(std::fabs(Burgers.dot(normalPrimary))>FLT_EPSILON),
		/*init list   */ normalConjugate(isSessile? VectorDim::Zero() : ds.conjugatePlaneNormal()),
		/*init list   */ rssPrimary  ((pkForce-pkForce.dot( normalPrimary)* normalPrimary  ).norm()),
		/*init list   */ rssConjugate((pkForce-pkForce.dot(normalConjugate)*normalConjugate).norm()),
		/*init list   */ crossSlipFactor(1.5),
		/*init list   */ isCrossSlipSegment(chord.norm()>1.1*crossSlipLength
		/*                               */ && chord.normalized().cross(Burgers.normalized()).norm()<=sinThetaCrossSlipCr
		/*                               */ && !sourceOnMeshBoundary && !sinkOnMeshBoundary // not on the boundary
		/*                               */ && crossSlipFactor*rssPrimary<rssConjugate
        /*                               */ && !isSessile) {
		
		}
		
		
	};
		
	//////////////////////////////////////////////////////////////s
} // namespace model
#endif

