#include <iostream>
#include <string>
#include <sstream>      // std::istringstream

//#include <istringstream>



int main()
{
//    std::string str = "Hello, world";
//    std::istringstream in(str);
//    std::string word1, word2;
//    
//    in >> word1;
//    in.seekg(0,in.end); // rewind
//    int length = in.tellg();
////    in.seekg(0); // rewind
//    in.seekg(length - 1, in.beg);
//
//    in >> word2;
//    
//    std::cout << "word1 = " << word1 << '\n'
//    << "word2 = " << word2 << '\n';
    
    std::string str;
//    std::cin.getline(&str,2);
//    std::cout<<str<<std::endl;
    
    while(!std::cin.eof())
    {
        
        std::cin>>str;
//        std::cout
//        std::cin.seekg( 0);

//        std::cin.seekg( std::ios_base::seekdir::end);
//        std::cin.clear();
//        std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        //                    std::cout <<"   invalid input, try again: "<<std::flush;
        //std::cin>>frameID;
        // code'
        std::cout<<str<<" "<<std::cin.fail()<<std::endl;
    }
    return 0;
}
