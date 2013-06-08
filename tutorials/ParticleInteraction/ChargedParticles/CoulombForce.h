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

#ifndef _CoulombForce_h
#define _CoulombForce_h

#include <iostream>
#include <model/ParticleInteraction/InteractionBase.h>
//#include <tutorials/ParticleInteraction/ChargedParticles/ChargedParticle.h>



namespace model {
    
    // Class predeclaration
//    class ChargedParticle;
    
    template <typename ChargedParticleType>
    struct CoulombForce :
    /* inheritance   */ public InteractionBase<double,3>
    { // change this 3 to ChargedParticle::dim
        
        typedef typename ChargedParticleType::PositionType PositionType;
//        typedef Eigen::Matrix<double,3,1> PositionType;
//        typedef Eigen::Matrix<double,3,1> ForceType;
//        typedef Eigen::Matrix<double,3,1> ResultType;

        typedef typename ChargedParticleType::ForceType ForceType;
        typedef ForceType ResultType;
        
        
#ifdef _MODEL_MPI_
        /**********************************************************************/
        CoulombForce(ChargedParticleType& cp1,const ChargedParticleType& cp2)
        {/*! @param[in] cp1 A const reference to ChargedParticle 1
          *  @param[in] cp2 A const reference to ChargedParticle 2
          *
          *  This constructor works with ParticleSystem<true,ChargedParticle>,
          *  that is the paralle version of ParticleSystem. The result of the 
          *  interaction between cp1 and cp2 is stored in this->resultVector,
          *  therefore cp1 and cp2 can be const.
          */
            
            PositionType R(cp1.P()-cp2.P());
            double r2( R.squaredNorm() ); // squared distance between particles
            if(r2!=0.0)
            {
                ForceType f(cp1.q*cp2.q/r2*R.normalized()); // the force on particle 2
                this->resultVector[cp1.mpiID*3+0]+=f(0);
                this->resultVector[cp1.mpiID*3+1]+=f(1);
                this->resultVector[cp1.mpiID*3+2]+=f(2);
                cp1.force+=f;
            }
            
        }
        
        /**********************************************************************/
        static ResultType get(const ChargedParticleType& cp1)
        {
            return ResultType((ResultType()<<InteractionBase<double,3>::resultVector[cp1.mpiID*3+0],
                               /*         */ InteractionBase<double,3>::resultVector[cp1.mpiID*3+1],
                               /*         */ InteractionBase<double,3>::resultVector[cp1.mpiID*3+2]).finished());
        }
        
#else
        /**********************************************************************/
        CoulombForce(ChargedParticleType& cp1, const ChargedParticleType& cp2)
        {/*! @param[in] cp1 A reference to ChargedParticle 1
          *  @param[in] cp2 A reference to ChargedParticle 2
          *
          *  This constructor works with ParticleSystem<false,ChargedParticle>,
          *  that is the serial version of ParticleSystem. The result of the
          *  interaction between cp1 and cp2 is stored in each particle.
          */
            
            PositionType R(cp1.P()-cp2.P());
            double r2( R.squaredNorm() ); // squared distance between particles
            if(r2!=0.0)
            {
                ForceType f(cp1.q*cp2.q/r2*R.normalized()); // the force on particle 1
                cp1.force+=f;
            }
            
        }
        
        /**********************************************************************/
        static const ResultType& get(const ChargedParticleType& cp1)
        {
            return cp1.force;
        }
#endif

        
        /*****************************************/
        static void reset(ChargedParticleType& cp)
        {/*! Set force of cp to zero
          */
            cp.force.setZero();
        }
        
        /*****************************************/
        static int dataPerParticle()
        {
            return 3;
        }
        
//        /*****************************************/
//        static ResultType getResult(const ChargedParticle& p){
//            const int mpiID(p.mpiID);
//            
//            return (ResultType()<<resultVector[mpiID*3+0],resultVector[mpiID*3+1],resultVector[mpiID*3+2]).finished();
//            
//            //resultVector[mpiID*3+1]
//            //resultVector[mpiID*3+2]
//        }
        
    };
        
}
#endif

