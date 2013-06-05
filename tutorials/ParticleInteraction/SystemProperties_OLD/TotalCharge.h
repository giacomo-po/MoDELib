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

#ifndef _TotalCharge_h
#define _TotalCharge_h


struct TotalCharge{
    
    typedef double PropertyType;
    
    
    PropertyType totalCharge;
        
    /*****************************************/
    TotalCharge() : totalCharge(0.0)
    {/*! Constructor initializes totalCharge to 0.0
      */
    }
    
    /*****************************************/
    operator const PropertyType& () const
    {/*! Cast operator returns totalCharge
      */
        return totalCharge;
    }

};


#endif

