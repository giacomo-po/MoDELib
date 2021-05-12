/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigVtk_h_
#define model_DDconfigVtk_h_


#include <DDconfigIO.h>
#include <DDconfigVtkBase.h>
#include <SimplicialMesh.h>
#include <DislocationNodeActor.h>
#include <DislocationSegmentActor.h>
#include <DislocationLoopActor.h>

namespace model
{
    struct DDconfigVtk : public DDconfigVtkBase
    /*                */,public DDconfigIO<3>
    {
        
        DislocationNodeActor nodes;
        DislocationSegmentActor segments;
        DislocationLoopActor loops;
        
        /**********************************************************************/
        DDconfigVtk(vtkRenderer* const ren,const SimplicialMesh<3>& mesh) :
        /* init */ nodes(ren)
        /* init */,segments(ren)
        /* init */,loops(ren)
        {
            
        }
        
        /**********************************************************************/
        void updateConfiguration(const size_t& frameID)
        {
            if(DDconfigIO<3>::isBinGood(frameID))
            {
                this->readBin(frameID);
            }
            else
            {
                this->readTxt(frameID);
            }
            
            nodes.updateConfiguration(*this);
            segments.updateConfiguration(*this,nodes.nodePolyData);
            loops.updateConfiguration(*this,nodes.nodePolyData);
        }
        
        /**********************************************************************/
        void modify()
        {
            nodes.modify();
            segments.modify();
            loops.modify();
        }
        
    };
    
}
#endif
