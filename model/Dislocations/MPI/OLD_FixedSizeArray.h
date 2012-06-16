#ifndef mmdl_FIXEDSIZEARRAY_H_
#define mmdl_FIXEDSIZEARRAY_H_

#include <iostream>
#include <assert.h>

namespace mmdl {
	
	/****************************************************/
	/*FixedSizeArray<T,N> general case ******************/
	/****************************************************/
	template <typename T,unsigned int N>
	class FixedSizeArray : private FixedSizeArray<T,N-1> {
		
		T data;
		
		public:
				
		/************************************************/
		const T& operator()(const unsigned int& k) const {
			/*! Read-only acces to the k-th element in the array.
			 */
			assert(k<N && "INDEX MUST BE IN [0,N-1]");
			return (k==N-1)? data : FixedSizeArray<T,N-1>::operator()(k);
		}
		
		/************************************************/		
		T& operator()(const unsigned int& k) {
			/*! Read-write acces to the k-th element in the array.
			 */
			assert(k<N && "INDEX MUST BE IN [0,N-1]");
			return (k==N-1)? data : FixedSizeArray<T,N-1>::operator()(k);
		}
		
//		/************************************************/
//		template<unsigned int k>
//		const T& operator()() const {
//			/*! Read-only acces to the k-th element in the array.
//			 */
//			assert(k<N && "INDEX MUST BE IN [0,N-1]");
//			return (k==N-1)? data : FixedSizeArray<T,N-1>::operator()(k);
//		}
//		
//		/************************************************/		
//		T& operator()(const unsigned int& k) {
//			/*! Read-write acces to the k-th element in the array.
//			 */
//			assert(k<N && "INDEX MUST BE IN [0,N-1]");
//			return (k==N-1)? data : FixedSizeArray<T,N-1>::operator()(k);
//		}
		
		/************************************************/		
		unsigned int size() const {
			/*! The container size, i.d. N.
			 */
			return N;
		}
		
	};
	
	
	/****************************************************/
	/*FixedSizeArray<T,1> template specialization *******/
	/****************************************************/
	template <typename T>
	class FixedSizeArray<T,1> {
		
		T data;
			
	public:
		
//		const unsigned int size;
		
//		FixedSizeArray() : size(1){}
		
		/************************************************/
		const T& operator()(const unsigned int& k) const {
			assert(k==0 && "INDEX MUST BE IN [0,N-1]");
			return data;
		}
		
		/************************************************/
		T& operator()(const unsigned int& k) {
			assert(k==0 && "INDEX MUST BE IN [0,N-1]");
			return data;
		}
		
		/************************************************/		
		unsigned int size() const {
			/*! The container size, i.d. 1.
			 */
			return 1;
		}
		
	};
	
	
	/****************************************************/
	/*FixedSizeArray<T,0> template specialization *******/
	/****************************************************/
	template <typename T>
	class FixedSizeArray<T,0> {
		// Avoid infinite recursion if 
		// FixedSizeArray<T,0> is found by compiler
	};
	
}
#endif