/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_RGBmap_H_
#define model_RGBmap_H_


// taken from:
// http://stackoverflow.com/questions/7706339/grayscale-to-red-green-blue-matlab-jet-color-scale


namespace model {
		
    
    struct RGBcolor{
         double r;
         double g;
         double b;
    };
    
	/*************************************************************/
	/* RGBmap **************************************************/
	/*************************************************************/
	struct RGBmap  {
        
        /*! Use:
         c = GetColour(v,-1.0,1.0);
         
         */
        static RGBcolor getColor(double v,double vmin,double vmax){
            RGBcolor c = {1.0,1.0,1.0}; // white
            double dv;
            
            if (v < vmin)
                v = vmin;
            if (v > vmax)
                v = vmax;
            dv = vmax - vmin;
            
            if (v < (vmin + 0.25 * dv)) {
                c.r = 0;
                c.g = 4 * (v - vmin) / dv;
            } else if (v < (vmin + 0.5 * dv)) {
                c.r = 0;
                c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
            } else if (v < (vmin + 0.75 * dv)) {
                c.r = 4 * (v - vmin - 0.5 * dv) / dv;
                c.b = 0;
            } else {
                c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
                c.b = 0;
            }
            
            return(c);
        
        }


	};
	
	
} // namespace model
#endif
