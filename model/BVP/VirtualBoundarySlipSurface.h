/* This file is part of model, the Mechanics of Material Defects Library.
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_VIRTUALBOUNDARYSLIPSURFACE_H_
#define  model_VIRTUALBOUNDARYSLIPSURFACE_H_


#include <Eigen/Dense>

namespace model {
		
		template <short unsigned int dim>
		class VirtualBoundarySlipSurface {

		typedef Eigen::Matrix<double,dim,dim> MatrixDim;			
		typedef Eigen::Matrix<double,dim,1>   VectorDim;


	VectorDim getOutNormal(){

	}


	public:

		const static double L;
		const VectorDim sourceP;
		const VectorDim sinkP;
		const VectorDim sourceN;
		const VectorDim sinkN;
		const VectorDim Burgers;

		const VectorDim sourcePatL;
		const VectorDim   sinkPatL;

	/* VirtualBoundarySlipSurface ********************************/
//	VirtualBoundarySlipSurface(const VectorDim& source_in, const VectorDim& sink_in,
//	/*                      */ const VectorDim& sourceN_in, const VectorDim& sinkN_in,
//	/*                      */ const VectorDim& B_in) : 
//   /*                                               */ sourceP(source_in),
//   /*                                               */   sinkP(  sink_in),
//   /*                                               */ sourceN(sourceN_in.normalized()),
//   /*                                               */   sinkN(  sinkN_in.normalized()),
//	/*                                               */   sourcePatL(sourceP+L*sourceN),
//	/*                                               */     sinkPatL(  sinkP+L*sinkN){
//		std::cout<<"Creating VirtualBoundarySlipSurface "<<std::endl;
//	}

	/* VirtualBoundarySlipSurface ********************************/
	template <typename DislocationSegmentType>
	VirtualBoundarySlipSurface(const DislocationSegmentType& ds) : 
   /*                                               */ sourceP(ds.source->get_P()),
   /*                                               */   sinkP(ds.sink->get_P()),
   /*                                               */ sourceN(getOutNormal(ds.chord().cross(ds.glidePlaneNormal))),
   /*                                               */   sinkN(dsourceN),
	/*                                               */   sourcePatL(sourceP+L*sourceN),
	/*                                               */     sinkPatL(  sinkP+L*sinkN){
		std::cout<<"Creating VirtualBoundarySlipSurface "<<std::endl;
	}

	/* stress ****************************************************/
	MatrixDim stress(const VectorDim& Rfield){

	}

	/* displacement **********************************************/
	VectorDim displacement(const VectorDim& Rfield, const VectorDim& S) const{
	VectorDim temp(VectorDim::Zero());
	if (isPiercing(Rfield,S)){
		temp+=...;	
	}
	return temp;
	}

	/* isPiercing **********************************************/
	bool isPiercing(const VectorDim& Rfield, const VectorDim& S) const {

	}

		};
		
// init static data
template <short unsigned int dim>
const double VirtualBoundarySlipSurface<dim>::L=3000.0;

} // namespace model
#endif
