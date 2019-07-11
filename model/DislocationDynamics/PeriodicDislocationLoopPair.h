/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2019 by Yash Pachaury <ypachaur@purdue.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_PeriodicDislocationLoopPair_H_
#define model_PeriodicDislocationLoopPair_H_

#include <string>
#include <TextFileParser.h>
#include <IDreader.h>
#include <MPIcout.h>


namespace model
{
    template<typename LinkSequenceType>
    struct PeriodicDislocationLoopPair : public std::set<const LinkSequenceType*>
    {
        
        ~PeriodicDislocationLoopPair()
        {
            assert(this->empty());

        }
        
        void insertLoop(const LinkSequenceType* linkSequence)
        {
            this->insert(linkSequence);
            assert(this->size()<=2);
        }
        
        void eraseLoop(const LinkSequenceType* linkSequence)
        {
            this->erase(linkSequence);
        }
        
    };
}
#endif
