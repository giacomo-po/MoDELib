/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <SutherlandHodgman.h>

using namespace model;

template<typename PointType>
void writePolygon(const std::vector<PointType>& poly, const std::string& filename)
{
    std::ofstream outputFile(filename);
    for(const auto& pt : poly)
    {
        outputFile<<std::setprecision(15)<<std::scientific<<pt(0)<<" "<<pt(1)<<"\n";
    }
}

int main(int argc, char** argv)
{
    // subject polygon
    std::vector<Eigen::Vector2d> subjectPolygon{
        (Eigen::Vector2d()<<50,150).finished(),
        (Eigen::Vector2d()<<200,50).finished(),
        (Eigen::Vector2d()<<350,150).finished(),
        (Eigen::Vector2d()<<350,300).finished(),
        (Eigen::Vector2d()<<250,300).finished(),
        (Eigen::Vector2d()<<200,250).finished(),
        (Eigen::Vector2d()<<150,350).finished(),
        (Eigen::Vector2d()<<100,250).finished(),
        (Eigen::Vector2d()<<100,200).finished()}; // must be counter-clockwise
    
    // clipping polygon
    std::vector<Eigen::Vector2d> clipPolygon{
        (Eigen::Vector2d()<<100,100).finished(),
        (Eigen::Vector2d()<<300,100).finished(),
        (Eigen::Vector2d()<<300,300).finished(),
        (Eigen::Vector2d()<<100,300).finished()}; // must be counter-clockwise
    
    // The new clipped polygon
    std::vector<Eigen::Vector2d> newPolygon(SutherlandHodgman::clip(subjectPolygon,clipPolygon));
    
    // print clipped polygon points
    std::cout << "Clipped polygon points:" << std::endl;
    for(const auto& pt : newPolygon)
    {
        std::cout << "(" << pt(0) << ", " << pt(1) << ")" << std::endl;
    }
    
    writePolygon(subjectPolygon,"subject.txt");
    writePolygon(clipPolygon,"clipper.txt");
    writePolygon(newPolygon,"result.txt");

    
    return 0;
}
