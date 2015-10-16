/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGenerator_H_
#define model_MicrostructureGenerator_H_


#include <assert.h>
#include <model/Utilities/EigenDataReader.h>
#include <model/Mesh/SimplicialMesh.h> // defines mode::cout

namespace model
{
    
    class MicrostructureGenerator
    {
    
        SimplicialMesh<3> mesh;
        
    public:
        
        MicrostructureGenerator()
        {
            int meshID(0);
            EigenDataReader EDR;
            bool use_boundary=false;
            EDR.readScalarInFile("./DDinput.txt","use_boundary",use_boundary);
            
            
            
            double percentSessile=0.0;
            EDR.readScalarInFile("./DDinput.txt","percentSessile",percentSessile);
            
            if (use_boundary)
            {
                double targetDensity=0.0;
                EDR.readScalarInFile("./DDinput.txt","targetDensity",targetDensity);
                
                EDR.readScalarInFile("./DDinput.txt","meshID",meshID);
                mesh.readMesh(meshID);
                
                
                
                
            }
            else
            {
                assert(0 && "MICROSTRUCTURE GENERATION IN INFINITE SPACE NOT SUPPORTED YET.");
            }
            
        }
        
    
    };

}
#endif
