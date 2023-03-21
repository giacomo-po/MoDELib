/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2018 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2018 by Yinan Cui <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_SphericalInclusionIO_H_
#define model_SphericalInclusionIO_H_

#ifndef _MODEL_NON_SINGULAR_DD_
#define _MODEL_NON_SINGULAR_DD_ 1
#endif

#include <Eigen/Dense>
//#include <Material.h>
#include <StaticID.h>

//#include <DislocationStress.h>
// http://solidmechanics.org/text/Chapter5_4/Chapter5_4.htm
#include <SlipSystem.h>

namespace model
{
    
    
    template <int dim>
    struct SphericalInclusionIO
    {
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;
        typedef Eigen::Matrix<double,dim,1>   VectorDim;
        
        size_t inclusionID;
        VectorDim C;
        double a; // size (radius)
        MatrixDim eT;
        double mobilityReduction;
        int phaseID;

        /**********************************************************************/
        template<typename SphericalInclusionType>
        SphericalInclusionIO(const SphericalInclusionType& ei) :
        /* init */ inclusionID(ei.sID)
        /* init */,C(ei.C)
        /* init */,a(ei.a)
        /* init */,eT(ei.eT)
        /* init */,mobilityReduction(ei.mobilityReduction)
        /* init */,phaseID(ei.phaseID)
        {
            
        }
        
        SphericalInclusionIO(const size_t& inclusionID_in,
                           const VectorDim& C_in,
                           const double& a_in, // size (radius)
                           const MatrixDim& eT_in,
                           const double& mobilityReduction_in,
                           const int& phaseID_in) :
        /* init */ inclusionID(inclusionID_in)
        /* init */,C(C_in)
        /* init */,a(a_in)
        /* init */,eT(eT_in)
        /* init */,mobilityReduction(mobilityReduction_in)
        /* init */,phaseID(phaseID_in)
        {
            
        }
        
        SphericalInclusionIO() :
        /* init */ inclusionID(0)
        /* init */,C(VectorDim::Zero())
        /* init */,a(0)
        /* init */,eT(MatrixDim::Zero())
        /* init */,mobilityReduction(1.0)
        /* init */,phaseID(0)
        {
            
        }
        
        SphericalInclusionIO(std::stringstream& ss) :
        /* init */ inclusionID(0)
        /* init */,C(VectorDim::Zero())
        /* init */,a(0)
        /* init */,eT(MatrixDim::Zero())
        /* init */,mobilityReduction(1.0)
        /* init */,phaseID(0)
        {
            ss>>inclusionID;
            for(int d=0;d<dim;++d)
            {
                ss>>C(d);
            }
            ss>>a;
            for(int i=0;i<dim;++i)
            {
                for(int j=0;j<dim;++j)
                {
                    double temp;
                    ss>>temp;
                    eT(i,j)=temp;
                }
            }
            ss>>mobilityReduction;
            ss>>phaseID;
        }
        
        template <class T>
        friend T& operator << (T& os, const SphericalInclusionIO<dim>& ei)
        {
            os  << ei.inclusionID<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ei.C.transpose()<<"\t"
            /**/<< std::setprecision(15)<<std::scientific<<ei.a<<"\t";
            for(int d=0;d<dim;++d)
            {
                os  <<std::setprecision(15)<<std::scientific<<ei.eT.row(d)<<"\t";
            }
            os  <<std::setprecision(15)<<std::scientific<<ei.mobilityReduction<<"\t"
            /**/<<ei.phaseID;
            return os;
        }
        
        
    };
    
    
}
#endif
