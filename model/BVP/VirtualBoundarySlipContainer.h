/* This file is part of mmdl, the Mechanics of Material Defects Library.
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 *
 * mmdl is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  mmdl_VIRTUALBOUNDARYSLIPCONTAINER_H_
#define  mmdl_VIRTUALBOUNDARYSLIPCONTAINER_H_

#include <boost/ptr_container/ptr_vector.hpp>
#include <Eigen/Dense>

#include <mmdl/BVP/VirtualBoundarySlipSurface.h>

namespace mmdl {
		
		template <short unsigned int dim>
		class VirtualBoundarySlipSurface : public boost::ptr_vector<VirtualBoundarySlipSurface<dim> >{

		typedef boost::ptr_vector<VirtualBoundarySlipSurface<dim> BaseContainerType;

		typedef Eigen::Matrix<double,dim,dim> MatrixDim;			
		typedef Eigen::Matrix<double,dim,1>   VectorDim;



		public:

		/* stress ****************************************************/
	MatrixDim stress(const VectorDim& Rfield) const {
	MatrixDim temp(MatrixDim::Zero());
	for(typename BaseContainerType::const_iterator sIter=this->begin();sIter!=this->end();++sIter){
		temp+=sIter->stress(Rfield);	
	}
	return temp;
	}

	/* displacement **********************************************/
	VectorDim displacement(const VectorDim& Rfield, const VectorDim& S) const {
	VectorDim temp(VectorDim::Zero());
	for(typename BaseContainerType::const_iterator sIter=this->begin();sIter!=this->end();++sIter){
		temp+=sIter->displacement(Rfield,S);	
	}
	return temp;
	}		



		};
		


} // namespace mmdl
#endif
