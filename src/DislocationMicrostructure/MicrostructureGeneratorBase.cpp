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

    std::shared_ptr<PeriodicPlanePatch<3>> PolyPoint::periodicPlanePatch() const
    {        
        return nullptr;
    }

    MicrostructureGeneratorBase::MicrostructureGeneratorBase(const std::string& fileName):
    /* init */ microstructureFileName(fileName)
    /* init */,parser(microstructureFileName)
    /* init */,type(parser.readString("type",false))
    /* init */,style(parser.readString("style",false))
    /* init */,tag(parser.readString("tag",false))
    {
        std::cout<<greenColor<<tag<<" -> "<<type<<" "<<style<<defaultColor<<" ("<<fileName<<")"<<std::endl;
    }


}
#endif
