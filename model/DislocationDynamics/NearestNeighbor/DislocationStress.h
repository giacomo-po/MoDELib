/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_DislocationStress_h_
#define _model_DislocationStress_h_

#include <model/ParticleInteraction/FieldBase.h>
#include <model/ParticleInteraction/FieldPoint.h>


namespace model
{
	
	/**************************************************************************/
	/**************************************************************************/
	template<short unsigned int _dim>
	class DislocationStress
    /* inheritance */ : public FieldBase<double,_dim,_dim>
    {
        
//        typedef typename DislocationParticleType::VectorDimD VectorDimD;
//        typedef typename DislocationParticleType::MatrixDim MatrixDim;
        
//        enum{dim=DislocationParticleType::dim};
//        enum{dataSize=dim*dim};
        
        
    public:
        

        
        typedef DislocationStress<_dim> DislocationStressType;
        typedef FieldBase<double,_dim,_dim> FieldBaseType;
        typedef typename FieldBaseType::MatrixType MatrixType;
        
        static  double a2;
		static const Eigen::Matrix<double,_dim,_dim> I;

        
        template <typename DislocationParticleType>
        static MatrixType compute(const DislocationParticleType& source,const DislocationParticleType& field)
        {/*!@param[in] source the DislocationParticle that is source of stress 
          * @param[in] field  the DislocationParticle on which stress is computed
          *\returns the stress field produced by source on field
          */
            
            Eigen::Matrix<double,_dim,1> R(field.P-source.P);
			double RaSquared (R.squaredNorm() + a2);
			return   (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*source.T*(source.B.cross(R)).transpose()
                               + 	R*(source.T.cross(source.B)).transpose()
                               +   0.5* R.cross(source.B).dot(source.T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
                               )/std::pow(sqrt(RaSquared),3)*source.quadWeight;
            
            
//            const double c2 (Material<Isotropic>::C2*source.quadWeight/std::pow(sqrt(RaSquared),3));
//            const double c1(Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared));
//            
//            const Eigen::Matrix<double,_dim,1> BcrossR(source.B.cross(R));
//            
//            const double BcrossRdotT(BcrossR.dot(source.T));
//            const double c3(-BcrossRdotT*(1.0+3.0*a2/RaSquared));
//            const double c4(-BcrossRdotT*3.0/RaSquared);
//            
//            
//            MatrixType temp(c2*c3*I);
//            temp. template selfadjointView<Eigen::Upper>().rankUpdate(source.T,BcrossR,c2*c1);
//            temp. template selfadjointView<Eigen::Upper>().rankUpdate(R,source.T.cross(source.B),c2);
//            temp. template selfadjointView<Eigen::Upper>().rankUpdate(R,R,c2*c4);
//            
//            return temp.template selfadjointView<Eigen::Upper>();
            
        }
        
        template <typename DislocationParticleType, typename OtherParticleType>
        static MatrixType compute(const DislocationParticleType& source, const OtherParticleType& field)
        {/*!@param[in] source the DislocationParticle that is source of stress
          * @param[in] field  the DislocationParticle on which stress is computed
          *\returns the stress field produced by source on field
          */
            
            Eigen::Matrix<double,_dim,1> R(field.P-source.P);
			double RaSquared (R.squaredNorm() + a2);
			return   (Material<Isotropic>::C1*(1.0+1.5*a2/RaSquared)*source.T*(source.B.cross(R)).transpose()
                      + 	R*(source.T.cross(source.B)).transpose()
                      +   0.5* R.cross(source.B).dot(source.T) * (I*(1.0+3.0*a2/RaSquared) + 3.0/RaSquared*R*R.transpose())
                      )/std::pow(sqrt(RaSquared),3)*source.quadWeight;
            
        }
        
        
				
	};
    
    // Static data members
	template<short unsigned int _dim>
    double DislocationStress<_dim>::a2=1.0;  // square of core size a
    
	template<short unsigned int _dim>
	const Eigen::Matrix<double,_dim,_dim> DislocationStress<_dim>::I=Eigen::Matrix<double,_dim,_dim>::Identity();  // square of core size a

    
    /**************************************************************************/
    /**************************************************************************/
}	// close namespace
#endif
