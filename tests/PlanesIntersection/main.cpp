
#include <fstream>

#include <TextFileParser.h>
#include <PlanesIntersection.h>



using namespace model;


template<int dim>
void snapToPlanes(const Eigen::Matrix<double,dim,Eigen::Dynamic>& N,
                                         const Eigen::Matrix<double,dim,Eigen::Dynamic>& P,
                                         const Eigen::Matrix<double,dim,1>& x,
                                         const std::string& outputFileName)
{
    
    PlanesIntersection<dim> pInt(N,P);
    const auto snapped(pInt.snap(x));
    
    std::ofstream outFile(outputFileName);
    outFile<<snapped.first<<" "<<snapped.second.transpose()<<std::endl;
    
}


int main(int argc, char * argv[])
{
    
    const std::string inputFileName(argc > 1 ? argv[1] : "inputFile.txt");
    const std::string outputFileName(argc > 2 ? argv[2] : "outputFile.txt");

    const int dim(TextFileParser(inputFileName).readScalar<int>("dim",true));
    const int numPlanes(TextFileParser(inputFileName).readScalar<int>("numPlanes",true));

    
    const Eigen::MatrixXd N(TextFileParser(inputFileName).readMatrix<double>("N",dim,numPlanes,true));
    const Eigen::MatrixXd P(TextFileParser(inputFileName).readMatrix<double>("P",dim,numPlanes,true));
    const Eigen::MatrixXd x(TextFileParser(inputFileName).readMatrix<double>("x",dim,1,true));


    
    if(   P.rows()==dim && P.cols()==numPlanes
       && x.rows()==dim && x.cols()==1)
    {
        
        
        switch (dim)
        {
            case 1:
            {
                snapToPlanes<1>(N,P,x,outputFileName);
                break;
            }
                
            case 2:
            {
                snapToPlanes<2>(N,P,x,outputFileName);
                break;
            }
                
            case 3:
            {
                snapToPlanes<3>(N,P,x,outputFileName);
                break;
            }
                
            default:
                throw std::runtime_error("Test in dim>3 not implemented");
                break;
        }
        
    }
    else
    {
        throw std::runtime_error("Input data size mismatch");
    }
    



    return 0;
}

