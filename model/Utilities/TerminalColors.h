/* This file is part of PIL, the Particle Interaction Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * PIL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _TerminalColors_h
#define _TerminalColors_h

#include <iostream>
#include <iomanip>
#include <string>

namespace model {
    std::string defaultColor    = "\033[0m";	  // the default color for the console
    std::string redBoldColor    = "\033[1;31m";   // a bold red color
    std::string greenBoldColor  = "\033[1;32m";   // a bold green color
    std::string blueBoldColor   = "\033[1;34m";   // a bold blue color
    std::string magentaColor    = "\033[0;35m";   // a magenta color
}


#endif

