/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <EshelbyInclusion.h>
#include <Eigen/Dense>

using namespace model;

int main (int argc, char * const argv[])
{
    
    const std::string inputFile= argc > 1 ? std::string(argv[1]) : "inputFile.txt";
    
    const Eigen::Matrix<double,3,1> C(TextFileParser(inputFile).readMatrix<double>("C",1,3,true).transpose());
    const double a(TextFileParser(inputFile).readScalar<double>("a",true));
    const Eigen::Matrix<double,3,3> eT(TextFileParser(inputFile).readMatrix<double>("eT",3,3,true));
    const double nu(TextFileParser(inputFile).readScalar<double>("nu",true));
    const double mu(TextFileParser(inputFile).readScalar<double>("mu",true));
    const Eigen::MatrixXd X0(TextFileParser(inputFile).readMatrixRows<double>("x",1,true));
    const Eigen::MatrixXd X1(TextFileParser(inputFile).readMatrixRows<double>("y",1,true));
    const Eigen::MatrixXd X2(TextFileParser(inputFile).readMatrixRows<double>("z",1,true));

    EshelbyInclusion<3> ei(C,a,eT,nu,mu,1,1);

    std::ofstream outFile("inclusion.txt");

    for(int k=0;k<X2.cols();++k)
    {
        for(int i=0;i<X0.cols();++i)
        {
            for(int j=0;j<X1.cols();++j)
            {
                const Eigen::Matrix<double,3,1> x(X0(i),X1(j),X2(k));
                const auto stress(ei.stress(x));
                outFile<< x.transpose()<<" "<< stress.row(0)<<" "<<stress.row(1)<<" "<<stress.row(2)<<"\n";
//                std::cout<<x.transpose()<<std::endl;
            }
        }
    }
    
    std::ofstream vS_file("velocityS.txt");
//    std::ofstream vE_file("velocityE.txt");

    
    return 0;
}


