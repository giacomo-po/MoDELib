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

//	///////////////////////////////
//	bool is_removable() const {
//		return true;
//	}



/* isCommonNeighborAt *************************************************/
std::pair<bool,size_t> isNeighborAt(const VectorDim& P0) const {
    std::pair<bool,size_t> temp(false,0);
    for (typename Derived::NeighborContainerType::const_iterator nIiter=this->neighborhood().begin();nIiter!=this->neighborhood().end();++nIiter){ // loop over neighborhood of source
        if (boost::tuples::get<0>(nIiter->second)->sID!=this->sID){ // neighbor is neither source nor sink
            if((boost::tuples::get<0>(nIiter->second)->get_P()-P0).norm()<FLT_EPSILON){ // a neighbor of I exists at P0
                //const size_t p0ID(D); // the sID of the neighbor at P0
                temp=std::pair<bool,size_t>(true,boost::tuples::get<0>(nIiter->second)->sID);
            }
        }
    }
    return temp;
}