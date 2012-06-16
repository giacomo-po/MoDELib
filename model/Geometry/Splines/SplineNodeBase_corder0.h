/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */


typedef NetworkNode<Derived> NetworkNodeType;

typedef typename NetworkNodeType::LinkType LinkType;

typedef boost::tuple<Derived*,LinkType*,short int> NeighborType;
typedef std::map<size_t,NeighborType> NeighborContainerType;
typedef typename NeighborContainerType::iterator NeighborIteratorType;
typedef typename NeighborContainerType::const_iterator constNeighborIteratorType;

typedef Eigen::Matrix<double, dim, 1> VectorDim;

private:

	VectorDim P;

public:
EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	const VectorDim & get_P() const{
		return P;
	}

	///////////////////////////////
	bool is_removable() const {
		return true;
	}
