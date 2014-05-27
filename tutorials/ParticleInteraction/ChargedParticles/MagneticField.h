/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _MagneticField_h
#define _MagneticField_h

#include <model/ParticleInteraction/FieldBase.h>


namespace model {
    
    template <short unsigned int dim>
    struct MagneticField : public FieldBase<double,dim,1>
    {
        
        typedef MagneticField<dim> MagneticFieldType;
        typedef FieldBase<double,dim,1> FieldBaseType;
        typedef typename FieldBaseType::MatrixType MatrixType;
        
        template <typename ChargeParticleType>
        static MatrixType compute(const ChargeParticleType& source,const ChargeParticleType& field)
        {
            
            MatrixType temp(MatrixType::Zero());
            MatrixType R(field.P()-source.P());
            const double r2( R.squaredNorm() ); // squared distance between particles
            if(r2!=0.0)
            {
                temp = source.q/r2*source.V().cross(R.normalized()); // the force on particle 2
            }
            
            return temp;
        }
        
    };
        
}
#endif
