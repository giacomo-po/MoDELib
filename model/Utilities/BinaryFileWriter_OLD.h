/* This file is part of mmdl, the Mechanics of Material Defects Library.
 *
 * Copyright (C) 2011 by Benjamin Ramirez <ramirezbrf@gmail.com>.
 *
 * mmdl is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef mmdl_BINARYFILEWRITER_H_
#define mmdl_BINARYFILEWRITER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>

namespace mmdl {

template <typename scalar>
class BinaryFileWriter  {
	
/*
 
 #include <iostream>
 #include <mmdl/Utilities/BinaryFileReader.h>
 
 
 class SameClass{
 
 public:
 int a;
 
 void print() const{
 std::cout<<a<<std::endl;
 }
 
 };
 
 
 int main(){
 
 
 SameClass fnum[3];
 //int i;
 fnum[0].a=0;
 fnum[1].a=1;
 fnum[2].a=2;
 
 std::ofstream out("numbers", std::ios::out | std::ios::binary);
 if(!out) {
 std::cout << "Cannot open file.";
 return 1;
 }
 
 //std
 
 out.write((char *) &fnum, (sizeof fnum));
 
 out.close();
 
 
 mmdl::BinaryFileReader<SameClass> f("numbers");
 
 for (int k=0;k<f.size();++k){
 f[k].print();
 }
 //std::cout<<f.size()<<std::endl;
 
 
 //double fnum[4] = {9.5, -3.4, 1.0, 2.1};
 
 //std::cout<<sizeof fnum<<std::endl;
 return 0;
 }
 
 
 */
	
	
	
};

} // namespace mmdl
#endif
