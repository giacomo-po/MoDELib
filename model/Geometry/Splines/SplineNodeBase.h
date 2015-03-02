/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SPLINENODEBASE_H_
#define model_SPLINENODEBASE_H_

#include <Eigen/Dense>
//#include <model/Math/RoundEigen.h>
//#include <model/DislocationDynamics/Materials/CrystalOrientation.h>

namespace model
{
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename Derived,	short unsigned int dim, short unsigned int corder>
    class SplineNodeBase
    {
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename Derived, short unsigned int dim>
    class SplineNodeBase<Derived, dim,0>
    {
    public:
        typedef Eigen::Matrix<double, dim, 1>   VectorDim;
        typedef Eigen::Matrix<double, dim, dim> MatrixDim;
        
    private:
        VectorDim P;
        
//        VectorDim roundP(const VectorDim& P_in)
//        {
//            return CrystalOrientation<dim>::lattMat()*RoundEigen<double,dim>::round(CrystalOrientation<dim>::invLattMat()*P_in);
//        }
        
    public:
        
        /*************************************************/
        SplineNodeBase(const VectorDim& P_in) :
        /* init list */ P(P_in)
//        /* init list */ P(roundP(P_in))
        {

        }

        /*************************************************/
        void set_P(const VectorDim& P_in)
        {
            P=P_in;
//            P= roundP(P_in);
        }
        
        /*************************************************/
        const VectorDim& get_P() const
        {
            return P;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename Derived, short unsigned int dim>
    class SplineNodeBase<Derived, dim,1> : public SplineNodeBase<Derived, dim,0>
    {
        
    public:
        typedef SplineNodeBase<Derived, dim,0> Base;
        typedef typename Base::VectorDim VectorDim;
        typedef typename Base::MatrixDim MatrixDim;
        
    private:
        VectorDim T;
        
    public:
        
        MatrixDim prjM;	//! the projection matrix. THIS SHOULD BE PRIVATE
        
        SplineNodeBase(const VectorDim& P_in,const VectorDim& T_in) :
        /* init list */ Base(P_in),
        /* init list */ T(T_in),
        /* init list */ prjM(MatrixDim::Identity())
        {
            
        }
        
        /*************************************************/
        const VectorDim& get_T() const
        {
            return T;
        }
        
        /*************************************************/
        void set_T(const VectorDim& T_in)
        {
            T=T_in;
        }
        
    };
    
    /**************************************************************************/
    /**************************************************************************/
    template <typename Derived, short unsigned int dim>
    class SplineNodeBase<Derived, dim,2> : public SplineNodeBase<Derived, dim,1>
    {
        
        typedef Eigen::Matrix<double, dim, 1>	VectorDim;
        typedef Eigen::Matrix<double, dim, dim> MatrixDim;
        
        VectorDim K;
        
    public:
        
        const VectorDim& get_K() const
        {
            return K;
        }
        
    };
    
}
#endif

