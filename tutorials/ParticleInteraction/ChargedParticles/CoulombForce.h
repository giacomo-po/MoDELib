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
#include <tutorials/ParticleInteraction/ChargedParticles/ChargedParticle.h>

namespace model {
    
    
    struct CoulombForce : public pil::InteractionBase<double,3> { // change this 3 to ChargedParticle::dim
        
        typedef ChargedParticle::PositionType PositionType;
        typedef ChargedParticle::ForceType ForceType;
        typedef ForceType ResultType;
        
        //   static std::vector<double> resultVector;
        //typedef MPI_DOUBLE ResultType;
        //enum{ResultType=MPI_DOUBLE};
        
        /*****************************************/
        CoulombForce(ChargedParticle& cp1,ChargedParticle& cp2)
        {/*! Constructor with two ChargedParticle(s). Compu
          */
            //std::cout<<"Computing CoulombForce"<<std::endl;
            
            PositionType R(cp1.P()-cp2.P());
            double r2( R.squaredNorm() ); // squared distance between particles
            if(r2!=0.0)
            {
                ForceType f(cp1.q()*cp2.q()/r2*R.normalized()); // the force on particle 2
                //cp2.force+=f;
                //cp1.force-=f;
                this->resultVector[cp1.rID*3+0]+=f(0);
                this->resultVector[cp1.rID*3+1]+=f(1);
                this->resultVector[cp1.rID*3+2]+=f(2);
            }
            
        }
        
        /*****************************************/
        static void reset(ChargedParticle& cp)
        {/*! Set force of cp to zero
          */
            cp.force.setZero();
        }
        
        /*****************************************/
        static int dataPerParticle()
        {
            return 3;
        }
        
        /*****************************************/
        static ResultType getResult(const ChargedParticle& p){
            const int rid(p.rID);
            
            return (ResultType()<<resultVector[rid*3+0],resultVector[rid*3+1],resultVector[rid*3+2]).finished();
            
            //resultVector[rid*3+1]
            //resultVector[rid*3+2]
        }
        
    };
    
    //    // declare statica data
    //    std::vector<double> CoulombForce::resultVector;
    
}
#endif

