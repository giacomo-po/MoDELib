/* This file is part of finite element solution of BVP attached with model "the Mechanics of Material Defects Library".
 *
 * Copyright (C) 2011 by Mamdouh Mohamed <mamdouh.s.mohamed@gmail.com>, 
 * Copyright (C) 2011 by Giacomo Po <giacomopo@gmail.com>.
 * 
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef bvpfe_face_H_
#define bvpfe_face_H_

#include "model/BVP/Triangle.h"
#include <map>
#include <boost/ptr_container/ptr_map.hpp>
//#include <memory>

namespace bvpfe{

	class Face  {
	  
		 typedef boost::ptr_map<size_t,Triangle>    triContainerType;

		public:
		  
		  
			std::map<size_t,bvpfe::Triangle*> triContainer;
			//triContainerType triContainer;

			//-------------------------------------------------------
			//void insertTri (std::auto_ptr<Triangle> triPtr)
			void insertTri (Triangle* triPtr)
			{
				assert(triContainer.insert(std::make_pair(triPtr->sID,triPtr)).second);
				
				//assert(triContainer.insert(triPtr->sID,triPtr).second);
				
				//assert(triContainer.insert(triPtr->sID,triPtr).second);
		
				
			}
				
			//-------------------------------------------------------
			size_t faceSize()
			{
				return triContainer.size();
			}		

	};
		
}  //  namespace bvpfe
#endif
