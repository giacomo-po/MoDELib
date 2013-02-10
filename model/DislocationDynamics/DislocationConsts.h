/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DISLOCATIONCONSTS_H_
#define model_DISLOCATIONCONSTS_H_

namespace model {
	
//enum {freeNode=0, fixedNode=1, moveLine=2, boundaryNode=3};
//enum {freeNode=0, fixedNode=1, moveLine=2};	
enum {outsideMesh=0, insideMesh=1, onMeshBoundary=2};
enum {noBoundary=0, softBoundary=1, hardBoundary=2};
	
}

#endif