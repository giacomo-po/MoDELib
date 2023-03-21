/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PolyhedronInclusionIO_H_
#define model_PolyhedronInclusionIO_H_

#include <tuple>

namespace model
{
    
    template<short unsigned int dim>
    struct PolyhedronInclusionIO
    {
        
        typedef Eigen::Matrix<double,dim,dim> MatrixDim;

        
        size_t inclusionID;
        MatrixDim eT;
        double mobilityReduction;
        int phaseID;


        template<typename PolyhedronInclusionType>
        PolyhedronInclusionIO(const PolyhedronInclusionType& ei) :
        /* init */ inclusionID(ei.sID)
        /* init */,eT(ei.eT)
        /* init */,mobilityReduction(ei.mobilityReduction)
        /* init */,phaseID(ei.phaseID)
        {
            
        }
        
        PolyhedronInclusionIO(const size_t& inclusionID_in,
                           const MatrixDim& eT_in,
                           const double& mobilityReduction_in,
                           const int& phaseID_in) :
        /* init */ inclusionID(inclusionID_in)
        /* init */,eT(eT_in)
        /* init */,mobilityReduction(mobilityReduction_in)
        /* init */,phaseID(phaseID_in)
        {
            
        }
        
        PolyhedronInclusionIO() :
        /* init */ inclusionID(0)
        /* init */,eT(MatrixDim::Zero())
        /* init */,mobilityReduction(1.0)
        /* init */,phaseID(0)
        {
            
        }
        
        PolyhedronInclusionIO(std::stringstream& ss) :
        /* init */ inclusionID(0)
        /* init */,eT(MatrixDim::Zero())
        /* init */,mobilityReduction(1.0)
        /* init */,phaseID(0)
        {
            ss>>inclusionID;
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
        friend T& operator << (T& os, const PolyhedronInclusionIO<dim>& ei)
        {
            os  << ei.inclusionID<<"\t";
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

