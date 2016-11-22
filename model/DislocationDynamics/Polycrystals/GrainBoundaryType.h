/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GrainBoundaryType_H_
#define model_GrainBoundaryType_H_

#include <string>
#include <math.h>       /* fabs */
#include <cfloat>       /* FLT_EPSILON */
#include <deque>
#include <model/LatticeMath/LatticeDirection.h>
#include <model/Math/CompileTimeMath/PermutationWithoutRepetition.h>


namespace model
{
    
    
    
    template <int dim>
    class GrainBoundaryType
    {
        
        
        static std::deque<Eigen::Matrix<double,dim,1>> reducePermutations(const Eigen::Matrix<double,dim,Eigen::Dynamic>& m)
        {
            std::deque<Eigen::Matrix<double,dim,1>> temp;
            
            for (int c=0;c<m.cols();++c)
            {
                bool isUnique=true;
                for(const auto& v : temp)
                {
                    isUnique*=( (m.col(c)-v).norm()>FLT_EPSILON && (m.col(c)+v).norm()>FLT_EPSILON);
                    
                    if(!isUnique)
                    {
                        break;
                    }
                }
                
                if(isUnique)
                {
                    temp.push_back(m.col(c));
                }
            }
            return temp;
        }
        
    public:
        
        typedef Eigen::Matrix<double,dim,1> VectorDimD;
        
        const std::string name;
        const double energyDensity;
        const double dislocationSpacing;
        const VectorDimD Burgers;
        
        const std::deque<Eigen::Matrix<double,dim,1>> axisPermutations;
        const std::deque<Eigen::Matrix<double,dim,1>> n1Permutations;
        const std::deque<Eigen::Matrix<double,dim,1>> n2Permutations;
        
        
        
        /**********************************************************************/
        GrainBoundaryType(const std::string& name_in,
                          const VectorDimD& axis,
                          const VectorDimD& n1,
                          const VectorDimD& n2,
                          const double& energy_in,
                          const double& spacing_in,
                          const VectorDimD& b_in) :
        /* init list */ name(name_in),
        /* init list */ energyDensity(energy_in),
        /* init list */ dislocationSpacing(spacing_in),
        /* init list */ Burgers(b_in),
        /* init list */ axisPermutations(reducePermutations(PermutationWithoutRepetition<dim>::permuteWithPlusMinusSign(axis))),
        /* init list */ n1Permutations(reducePermutations(PermutationWithoutRepetition<dim>::permuteWithPlusMinusSign(n1))),
        /* init list */ n2Permutations(reducePermutations(PermutationWithoutRepetition<dim>::permuteWithPlusMinusSign(n2)))
        {
            assert(energyDensity>=0.0);
            assert(dislocationSpacing>0.0);
            
            std::cout<<"AxisPermutations:"<<std::endl;
            for(const auto& v : axisPermutations)
            {
                std::cout<<v.transpose()<<std::endl;
            }
            
            std::cout<<"n1Permutations:"<<std::endl;
            for(const auto& v : n1Permutations)
            {
                std::cout<<v.transpose()<<std::endl;
            }
            
            std::cout<<"n2Permutations:"<<std::endl;
            for(const auto& v : n2Permutations)
            {
                std::cout<<v.transpose()<<std::endl;
            }

        }
        
        
    };
    
} // end namespace
#endif

