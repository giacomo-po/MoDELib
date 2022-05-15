#include <iostream>
#include <filesystem>
#include <iomanip>

int main()
{

    std::filesystem::path p1("/usr/local");
    std::cout<<  (p1 / "include/bar.txt").string() <<std::endl;
    std::cout<<  (p1 / "include/bar.txt").parent_path().string() <<std::endl;
    std::cout<<  (p1 / "include/bar.txt").filename().string() <<std::endl;
    std::cout<<  (p1 / "include/bar.txt").stem().string() <<std::endl;
    std::cout<<  (p1 / "include/bar.txt").extension().string() <<std::endl;

    return 0;
}
