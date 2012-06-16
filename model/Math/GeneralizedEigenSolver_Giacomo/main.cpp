#include <iostream>

//#include <float.h>

//#include <set>
//

#include <model/Math/GeneralizedEigenSolver_Giacomo/GivensHT.h>


template<int N>
class GivensQZiteration{

public:
GivensQZiteration(Eigen::Matrix<double,N,N>& H, Eigen::Matrix<double,N,N>& T){

	int n=N-1; // zero-based
	int m=n-1; // zero-based
	
	double HT00(H(0,0)/T(0,0));
	double a10=((H(m,m)/T(m,m)-HT00)*(H(n,n),T(n,n)-HT00) - H(m,n)/T(n,n)*H(n,m)/T(m,m) + H(n,m)/T(m,m)*T(m,n)/T(n,n)*HT00)*H(0,0)/H(1,0) + H(0,1)/T(1,1) - HT00*T(0,1)/T(2,2);
	double a20=H(1,1)/T(1,1)-HT00 - H(1,0)/T(0,0)*T(0,1)/T(1,1) - H(m,m)/T(m,m)-HT00 - H(n,n)/T(n,n)-HT00 + H(n,m)/T(m,m)*H(m,n)/T(n,n);
	double a30=H(2,1)/T(1,1);

	std::cout<<"a10="<<a10<<std::endl;
		std::cout<<"a20="<<a20<<std::endl;
			std::cout<<"a30="<<a30<<std::endl;

			Eigen::Matrix<double,2,2> G(model::Givens::Q<double>(a20,a30));
			H.template block<2,N>(1,0)=(G*H.template block<2,N>(1,0)).eval();	// use .eval() to avoid aliasing
			T.template block<2,N>(1,0)=(G*T.template block<2,N>(1,0)).eval();	// use .eval() to avoid aliasing
	//		Q.template block<N,2>(0,I)=(Q.template block<N,2>(0,I)*G.transpose()).eval();   // use .eval() to avoid aliasing
			Eigen::Matrix<double,2,1> temp(G*(Eigen::Matrix<double,2,1>()<<a20,a30).finished());	
			G=model::Givens::Q<double>(a10,temp(0));
			H.template block<2,N>(0,0)=(G*H.template block<2,N>(0,0)).eval();	// use .eval() to avoid aliasing
			T.template block<2,N>(0,0)=(G*T.template block<2,N>(0,0)).eval();	// use .eval() to avoid aliasing

//std::cout<<temp;
//		G=Givens::Q<double>(-T(I+1,I+1),T(I+1,I)).transpose();
	//		T.template block<I+2,2>(0,I)=(T.template block<I+2,2>(0,I)*G).eval();	// use .eval() to avoid aliasing
	//		H.template block<N,2>(0,I)=(H.template block<N,2>(0,I)*G).eval();	    // use .eval() to avoid aliasing
	//		Z.template block<2,N>(I,0)=(G.transpose()*Z.template block<2,N>(I,0)).eval();	    // use .eval() to avoid aliasing



}

};


int main(){
	
	const int N(5);
	
	Eigen::Matrix<double,N,N> A=Eigen::Matrix<double,N,N>::Random();
	Eigen::Matrix<double,N,N> B=Eigen::Matrix<double,N,N>::Random();
//	A=random
		
	std::cout<<A<<std::endl<<std::endl;
	std::cout<<B<<std::endl<<std::endl;
	
	
	model::GivensHT<N> Ght(A,B);	
	
	std::cout<<"Decomposition is succesful? "<<Ght.success<<std::endl;
	
	
	
	Eigen::Matrix<double,N,N> H=Ght.matrixH();
	Eigen::Matrix<double,N,N> T=Ght.matrixT();
	
	
//	int N=N;

std::cout<<"H="<<std::endl<<H<<std::endl;
		std::cout<<"T="<<std::endl<<T<<std::endl;

	GivensQZiteration<N>(H,T);

	std::cout<<"H="<<std::endl<<H<<std::endl;
	std::cout<<"T="<<std::endl<<T<<std::endl;
	
	return 0;
}




