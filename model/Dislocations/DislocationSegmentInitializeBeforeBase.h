/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po       <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


#ifndef model_DISLOCATIONSEGMENTINITIALIZEBEFOREBASE_H
#define model_DISLOCATIONSEGMENTINITIALIZEBEFOREBASE_H

#include <assert.h>
#include <set>
#include <Eigen/Core>
#include <Eigen/Geometry> // cross
#include <model/Network/Operations/EdgeExpansion.h>
#include <model/Dislocations/DislocationSharedObjects.h>

namespace model {
	
	
	template <short unsigned int dim, typename MaterialType>
	struct DislocationSegmentInitializeBeforeBase{
		
		DislocationSharedObjects<dim,MaterialType> shared;
		
		
		typedef Eigen::Matrix<double,dim,1> VectorDim;

		/********************************************/
		VectorDim find_planeNormal(const VectorDim& chord, const VectorDim& Burgers/*,const VectorDim& T1,const VectorDim& T2*/){
			enum {Nslips=MaterialType::Nslips};
			std::set<SlipSystem<dim,Nslips> > allowedSlipSystems;
			shared.material.find_slipSystem(chord,Burgers,/*T1,T2,*/allowedSlipSystems);
			return allowedSlipSystems.begin()->normal; // DON't LIKE THIS
		}

		/********************************************/
		VectorDim get_sessileNormal(const VectorDim& chord, const VectorDim& Burgers){
			assert(chord.norm()>FLT_EPSILON && "CHORD TOO SMALL");
			assert(Burgers.norm()>FLT_EPSILON && "Burgers TOO SMALL");
			VectorDim temp(chord.normalized().cross(Burgers));
			double tempNorm(temp.norm());
			if (tempNorm>FLT_EPSILON){
				temp.normalize();
			}
			else{
				temp.setZero();
			}
			return temp;
		}
		
		
		//! The glide plane normal
		const VectorDim glidePlaneNormal;
		const VectorDim sessilePlaneNormal;
		
		
		
		// Constructor with chord and Burgers
		DislocationSegmentInitializeBeforeBase(const VectorDim& chord, const VectorDim& Burgers) : glidePlaneNormal(find_planeNormal(chord,Burgers).normalized()),
		/*                                                                                      */ sessilePlaneNormal(get_sessileNormal(chord,Burgers)){}
		
		// Constructor with plane normal
		DislocationSegmentInitializeBeforeBase(const VectorDim& normal_in, const VectorDim& chord, const VectorDim& Burgers) : glidePlaneNormal(normal_in.normalized()),
		/*                                                                                                                  */ sessilePlaneNormal(get_sessileNormal(chord,Burgers)){}
			
	};
	
	//////////////////////////////////////////////////////////////s
} // namespace model
#endif

