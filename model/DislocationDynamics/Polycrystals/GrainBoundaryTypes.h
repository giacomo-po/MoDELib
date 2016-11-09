/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GrainBoundaryTypes_H_
#define model_GrainBoundaryTypes_H_




namespace model
{
    
    
    
    template <int dim>
    struct GrainBoundaryType :
    /* base */ private std::map<int,const Grain<dim>* const>,
    /* base */ private std::map<int,LatticePlane>
    {
        
        typedef Eigen::Matrix<double,2*dim+1,1> BGkeyType;
        typedef Eigen::Matrix<  double,dim,1> VectorDimD;
//        typedef Eigen::Matrix<double,dim,dim> MatrixDimD;

        const VectorDimD t;
        const VectorDimD n;
        const double theta;
        const VectorDimD v1;
        const VectorDimD v2;
        

        GrainBoundaryType(const BGkeyType& key,
                          const VectorDimD& v1_in,
                          const VectorDimD& v2_in) :
        /* init list */ t(key.segment<dim>(0)),
        /* init list */ n(key.segment<dim>(dim)),
        /* init list */ theta(key(2*dim)),
        /* init list */ v1(v1_in),
        /* init list */ v2(v2_in)
        {
        
        }
        
        
    };
    
} // end namespace
#endif

