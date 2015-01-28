/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_CellPropertiesBase_h_
#define _model_CellPropertiesBase_h_

#include <tuple> // c++11


namespace model
{
	
    /**************************************************************************/
    /**************************************************************************/
    template <typename...T>
    struct CellPropertiesBase
    {
        typedef std::tuple<T...> SpatialCellProperties;
    };

}
#endif

