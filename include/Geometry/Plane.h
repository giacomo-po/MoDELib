/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_Plane_H_
#define model_Plane_H_

#include <cfloat>
#include <tuple>
//#include <map>
#include <Eigen/Dense>
#include <iostream>


namespace model
{
    
    template <int dim>
    struct Plane
    {
        
        typedef Eigen::Matrix<double,dim,1> VectorDim;
        typedef Eigen::Matrix<double,dim-1,1> VectorLowerDim;
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        const VectorDim P;
        const VectorDim unitNormal;
        const MatrixDim L2G;

        /**********************************************************************/
        Plane(const VectorDim& p,const VectorDim& n);
        
        /**********************************************************************/
        bool contains(const VectorDim& P0) const;

        /**********************************************************************/
        bool isAbove(const VectorDim& P0) const;

        /**********************************************************************/
        bool isBelow(const VectorDim& P0) const;

        /**********************************************************************/
        VectorDim snapToPlane(const VectorDim& P0) const;
        
        /**********************************************************************/
        double distanceTo(const VectorDim& P0) const;
        
        /**********************************************************************/
        VectorLowerDim localPosition(const VectorDim& point) const;
        
        /**********************************************************************/
        VectorDim globalPosition(const VectorLowerDim& point) const;
        
        /**********************************************************************/
        static MatrixDim getL2G(VectorDim z);

    };
    
}
#endif
