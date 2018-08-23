#include <iostream>
#include <model/IO/TextFileParser.h>

using namespace model;

int main()
{
    
    
    // Create parser
    TextFileParser parser("in.txt");
    
    // Read a bool
    bool readVerbose=parser.readScalar<int>("readVerbose",true);

    
    // Read a string
    std::string a=parser.readString("a",readVerbose);
    
    // Read an int
    int b=parser.readScalar<int>("b",readVerbose);
    
    // Read an double
    double c=parser.readScalar<double>("c",readVerbose);
    
    // Read an vector of floats
    std::vector<float> d=parser.readArray<float>("d",readVerbose);
    
    // Read an 3x2 matrix
    Eigen::MatrixXd e32=parser.readMatrix<double>("e",3,2,readVerbose);
    
    // Read an 2x3 matrix
    Eigen::MatrixXd e23=parser.readMatrix<double>("e",2,3,readVerbose);
    
    return 0;
}
