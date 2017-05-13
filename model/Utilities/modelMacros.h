/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2013 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_macros_h_
#define _model_macros_h_



// see http://eigen.tuxfamily.org/dox-devel/Macros_8h_source.html


#ifdef _MODEL_NO_DEBUG_ // RELEASE MODE

#else // DEBUG MODE
#include <assert.h>
#include <iostream>
#include <cstdlib>   // abort
//namespace model
//{
//    inline void assert_fail(const char *condition, const char *function, const char *file, int line)
//    {
//        std::cerr << "assertion failed: " << condition << " in function " << function << " at " << file << ":" << line << std::endl;
//        abort();
//    }
//}

#endif


/* model_removeAssert *********************************************************/
#ifndef model_removeAssert
#ifdef _MODEL_NO_DEBUG_
// in release mode, model_removeAssert(x,y) removes x from the code
#define model_removeAssert(x)
#else
// in debug mode, model_removeAssert(x,y)=assert(x)
// TO DO: CHANGE THIS DEFINITION BECAUSE assert(x) CAN BE REMOVED BY NDEBUG!!! MUST USE mode::assert_fail
#define model_removeAssert(x) assert(x)
#endif
#endif

/* model_execAssert ***********************************************************/
#ifndef model_execAssert
#ifdef _MODEL_NO_DEBUG_
// in release mode, model_execAssert(x,y,z) executes x
#define model_execAssert(x,y,z) x
#else
// in debug mode, model_execAssert(x,y) executes xy and asserts (xy && z),
// where z is typically a string
// TO DO: CHANGE THIS DEFINITION BECAUSE assert(x) CAN BE REMOVED BY NDEBUG!!! MUST USE mode::assert_fail
#define model_execAssert(x,y,z) assert(x y && z)
//#define model_execAssert(x,y) model::assert_fail(x && y)
#endif
#endif

/* model_execAssert ***********************************************************/
//#ifndef model_checkInput
//#ifdef _MODEL_NO_BAD_INPUT_CHECK_
//// if function inputs are not checked, model_checkInput(x) removes x from the code
//#define model_checkInput(x)
//#else
//// if function inputs are checked, model_checkInput(x)=assert(x)
//// TO DO: CHANGE THIS DEFINITION BECAUSE assert(x) CAN BE REMOVED BY NDEBUG!!! MUST USE mode::assert_fail
//#define model_checkInput(x) assert(x)
//#endif
//#endif

#endif // _model_macros_h_
