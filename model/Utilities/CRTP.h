/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef  model_CRTP_H_
#define model_CRTP_H_

namespace model {
		
		/*! \brief Implementation of the Curiously Recurring Template Pattern (CRTP).
		 * CRTP achieves static polyphormism  by offerings access to the derived class within 
		 * the base class.
		 */
		template <typename Derived>
		class CRTP
    {
			
		protected:
			
			//! A reference to the derived object
			Derived&   derived() { return *static_cast<Derived*>(this); }
			
			//! A pointer to the derived object
			Derived* p_derived() { return  static_cast<Derived*>(this); }
			
			//! A  reference to the const derived object
			const Derived&   derived() const { return *static_cast<const Derived*>(this); }
			
            //! A pointer to the const derived object
            const Derived* p_derived() const { return  static_cast<const Derived*>(this); }
			
		public:
			//! The default constructor
			CRTP(){
			/*! Sample use: \code
			 * #include <iostream>
			 * #include <CRTP.h>
			 *
			 * template <typename Derived>
			 * class BaseClass : public CRTP<Derived>{
			 *	public:
			 *	void print() const {std::cout<<"Base class of DerivedClass"<< this->p_derived()->id()<<std::endl;}
			 * };
		 
			 * class DerivedClass1 : public BaseClass<DerivedClass1>{
			 *	public:
			 *	int id() const {return 1;}
			 * };
		         *
			 * class DerivedClass2 : public BaseClass<DerivedClass2>{
			 *	public:
			 *	int id() const {return 2;}
			 * };
			 *
			 * int main(){
			 * 	DerivedClass1 dc1;
			 * 	DerivedClass2 dc2;
			 * 	dc1.print();
			 * 	dc2.print();
			 * }	
			* \endcode
			*/	
			}
				
		};
		
		//////////////////////////////////////////////////////////////
} // namespace model
#endif
