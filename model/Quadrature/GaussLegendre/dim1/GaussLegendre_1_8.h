 /* This file is part of MODEL, the Mechanics Of Defect Evolution Library. 
 * 
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>. 
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>. 
 */

#ifndef model_GAUSSLEGENDRE_1_8_H_ 
#define model_GAUSSLEGENDRE_1_8_H_ 

namespace model{

/*** This file is automatically generated by generateGaussLegendre.cpp ***/
   template<>
   struct GaussLegendre<1,8>{
       static Eigen::Matrix<double,2,8> abcsissasAndWeights(){
           Eigen::Matrix<double,8,2> aw;
           aw<<1.985507175123208e-02, 5.061426814518828e-02, 
               1.016667612931867e-01, 1.111905172266873e-01, 
               2.372337950418354e-01, 1.568533229389436e-01, 
               4.082826787521751e-01, 1.813418916891810e-01, 
               5.917173212478249e-01, 1.813418916891810e-01, 
               7.627662049581645e-01, 1.568533229389438e-01, 
               8.983332387068137e-01, 1.111905172266870e-01, 
               9.801449282487682e-01, 5.061426814518853e-02; 
               
       return aw.transpose();
       } 
   }; 
/*************************************************/
} 
#endif 
