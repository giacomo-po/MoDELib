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

#ifndef _TerminalColors_h
#define _TerminalColors_h

#include <iostream>
#include <iomanip>
#include <string>

namespace model
{
    static std::string defaultColor    = "\033[0m";	  // the default color for the console

    static std::string redColor        = "\033[0;31m";
    static std::string redBoldColor    = "\033[1;31m";

    static std::string greenColor      = "\033[0;32m";   // a green color
    static std::string greenBoldColor  = "\033[1;32m";   // a bold green color
    
    static std::string yellowColor       = "\033[0;33m";   // a yellow color
    static std::string yellowBoldColor   = "\033[1;33m";   // a bold yellow color

    static std::string blueColor       = "\033[0;34m";   // a blue color
    static std::string blueBoldColor   = "\033[1;34m";   // a bold blue color
    
    static std::string magentaColor    = "\033[0;35m";   // a magenta color
    static std::string magentaBoldColor= "\033[1;35m";   // a bold magenta color

    static std::string cyanColor    = "\033[0;36m";   // a cyan color
    static std::string cyanBoldColor= "\033[1;36m";   // a bold cyan color

    
}


#endif

