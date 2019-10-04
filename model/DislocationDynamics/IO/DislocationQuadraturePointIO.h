/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DislocationQuadraturePointIO_H_
#define model_DislocationQuadraturePointIO_H_

#include <Eigen/Dense>


namespace model
{
    template<int dim>
    struct DislocationQuadraturePointIO
    {

        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        
        size_t sourceID;
        size_t sinkID;
        int qID;
        VectorDim r;                          // position
        double j;                             // jacobian dl/du
        VectorDim rl;                         // unit tangent dr/dl
        
        MatrixDim stress;
        VectorDim pkForce;
        VectorDim glideVelocity;
        
        /**********************************************************************/
        DislocationQuadraturePointIO() :
        /* init */ sourceID(0)
        /* init */,sinkID(0)
        /* init */,qID(0)
        /* init */,r(VectorDim::Zero())
        /* init */,j(0.0)
        /* init */,rl(VectorDim::Zero())
        /* init */,stress(MatrixDim::Zero())
        /* init */,pkForce(VectorDim::Zero())
        /* init */,glideVelocity(VectorDim::Zero())
        {
        }
        
        /**********************************************************************/
        template<typename DislocationQuadraturePointType>
        DislocationQuadraturePointIO(const DislocationQuadraturePointType& qPoint) :
        /* init */ sourceID(qPoint.sourceID)
        /* init */,sinkID(qPoint.sinkID)
        /* init */,qID(qPoint.qID)
        /* init */,r(qPoint.r)
        /* init */,j(qPoint.j)
        /* init */,rl(qPoint.rl)
        /* init */,stress(qPoint.stress)
        /* init */,pkForce(qPoint.pkForce)
        /* init */,glideVelocity(qPoint.glideVelocity)
        {
            
        }
        
        /**********************************************************************/
        DislocationQuadraturePointIO(std::stringstream& ss) :
        /* init */ sourceID(0)
        /* init */,sinkID(0)
        /* init */,qID(0)
        /* init */,r(VectorDim::Zero())
        /* init */,j(0.0)
        /* init */,rl(VectorDim::Zero())
        /* init */,stress(MatrixDim::Zero())
        /* init */,pkForce(VectorDim::Zero())
        /* init */,glideVelocity(VectorDim::Zero())
        {
            
            ss>>sourceID;
            ss>>sinkID;
            ss>>qID;
            for(int d=0;d<dim;++d)
            {
                ss>>r(d);
            }
            ss>>j;
            for(int i=0;i<dim;++i)
            {
                for(int j=0;j<dim;++j)
                {
                    ss>>stress(i,j);
                }
            }
            for(int d=0;d<dim;++d)
            {
                ss>>pkForce(d);
            }
            for(int d=0;d<dim;++d)
            {
                ss>>glideVelocity(d);
            }
        }
        
        /**********************************************************************/
        template <class T>
        friend T& operator << (T& os, const DislocationQuadraturePointIO<dim>& ds)
        {
            os  << ds.sourceID<<"\t"
            /**/<< ds.sinkID<<"\t"
            /**/<< ds.qID<<"\t"
            /**/<< ds.r.transpose()<<"\t"
            /**/<< ds.j<<"\t"
            /**/<< ds.rl<<"\t";
            for(int d=0;d<dim;++d)
            {
                os  << ds.stress.row(d)<<" ";
            }
            os  << ds.pkForce.transpose()<<"\t"
            /**/<< ds.glideVelocity.transpose();
            return os;
        }
        
    };
    
}
#endif
