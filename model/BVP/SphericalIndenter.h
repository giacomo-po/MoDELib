#ifndef bvpfe_SphericalIndenter_H_
#define bvpfe_SphericalIndenter_H_

#include<cmath>

//namespace bvpfe{

	class SphericalIndenter  {
	  
		double R;
	  
		public:
		  
		  SphericalIndenter(double R_in) : R(R_in) {}
		  
		  double localDepth (double depth, double r){
		    return (depth - R + sqrt(pow(R,2)-pow(r,2)) );
		  }
		  
		  double contactRadius (double depth){
		    return sqrt(pow(R,2)-pow(R-depth,2));
		  }

	};
		
//}  //  namespace bvpfe
#endif