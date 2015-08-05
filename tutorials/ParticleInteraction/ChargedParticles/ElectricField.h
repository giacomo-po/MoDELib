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

#ifndef _ElectricField_h
#define _ElectricField_h

#include <model/ParticleInteraction/FieldBase.h>


namespace model {
    
    template <short unsigned int dim>
    struct ElectricField : public FieldBase<double,dim,1>
    {
        
        typedef ElectricField<dim> ElectricFieldType;
        typedef FieldBase<double,dim,1> FieldBaseType;
        typedef typename FieldBaseType::MatrixType MatrixType;
        
        static bool use_multipole;

        /**********************************************************************/
        static const MatrixType& get(const MatrixType& temp)
        {
            return temp;
        }
        
        /**********************************************************************/
        template <typename ChargeParticleType, typename CellContainerType>
        static MatrixType multipole(const ChargeParticleType& field,const CellContainerType& farCells)
        {/*!@param[in] field  the FieldPoint at which stress is computed
          * @param[in] farCells container of SpatialCell(s) that are not neighbors of field
          *\returns not miplemente yet
          *
          */
            assert(0 && "multipole not implemented");
            return MatrixType::Zero(); //
        }
        
        template <typename ChargeParticleType>
        static MatrixType compute(const ChargeParticleType& source,const ChargeParticleType& field)
        {
            
            MatrixType temp(MatrixType::Zero());
            MatrixType R(field.P-source.P);
            const double r2( R.squaredNorm() ); // squared distance between particles
            if(r2!=0.0)
            {
                temp = source.q/r2*R.normalized(); // the force on particle 2
            }

            return temp;
        }
        
    };
    
    
    template <short unsigned int dim>
    bool ElectricField<dim>::use_multipole=false;

}
#endif

