/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_EVLio_H_
#define model_EVLio_H_

#include <model/DislocationDynamics/IO/DislocationNodeIO.h>


namespace model
{
    
    
    
    struct EVLio : public StaticID<EVLio>
    {
        

        /**********************************************************************/
        template<typename DislocationNodeType>
        void write(const DislocationNodeType& dn)
        {
         
            std::string filename("bin/EVL"+std::to_string(dn.runningID())+".bin");
            std::ofstream of(filename.c_str(), std::ios::out  | std::ios::binary);
            if(of.is_open())
            {
                const size_t nV(dn.nodes().size());
                of.write(nV);
                for(const auto& node : dn.nodes())
                {
                    of.write(DislocationNodeIO<DislocationNodeType::dim>(*node.second));
                }
            }
            else
            {
                model::cout<<"CANNOT OPEN "<<filename<<std::endl;
            }
        }
        

        
	};
	
}
#endif

