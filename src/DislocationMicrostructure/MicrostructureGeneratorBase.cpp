/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_MicrostructureGeneratorBase_cpp_
#define model_MicrostructureGeneratorBase_cpp_

#include <iostream>
#include <string>

#include <TextFileParser.h>
#include <TerminalColors.h>


#include <MicrostructureGeneratorBase.h>

namespace model
{

    PeriodicPlanePatch<3>* PolyPoint::periodicPlanePatch() const
    {
        return nullptr;
    }

    MicrostructureGeneratorBase::MicrostructureGeneratorBase(const std::string& fileName):
    /* init */ parser(fileName)
    /* init */,tag(parser.readString("tag",true))
//    /* init */,microstructureType(parser.readString("microstructureType",true))
    {
        std::cout<<greenColor<<"Generating microstructure from "<<fileName<<defaultColor<<std::endl;
    }


}
#endif
