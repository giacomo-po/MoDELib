
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

#include <vector>
#include <Eigen/Dense>
#include "gpc.c"
#include <stdio.h>



std::vector<Eigen::Matrix<double,2,1>> readPolyVector(const std::string& fileName)
{
    std::vector<Eigen::Matrix<double,2,1>> temp;
    std::ifstream file ( fileName.c_str() , std::ifstream::in );
    if(file.is_open())
    {
        
        std::string line;
        double x,y;
        
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            
            ss >> x >>y;
            temp.push_back((Eigen::Matrix<double,2,1>()<<x,y).finished());
            //std::cout<<x<<","<<y<<std::endl;
            
            //poly[0] <<ClipperLib::DoublePoint(x,y);
        }
        
        std::cout<<"Polygon:"<<std::endl;
        for(const auto& point : temp)
        {
            std::cout<<point.transpose()<<std::endl;
        }
    }
    else
    {
        std::cout<<"CANNOT READ "+fileName<<std::endl;
    }
    
    
    
    return temp;
}


struct GPCpolygon : public gpc_polygon
{
    
    //    gpc_polygon poly;
    
    GPCpolygon(const std::vector<Eigen::Matrix<double,2,1>>& v)
    {// constructor specialized for no holes
        this->num_contours=1; // no holes
        this->contour= (gpc_vertex_list*) malloc(this->num_contours*sizeof(gpc_vertex_list));
        this->contour[0].num_vertices=v.size();
        this->contour[0].vertex=(gpc_vertex*) malloc(this->contour[0].num_vertices*sizeof(gpc_vertex));
        for(int k=0;k<v.size();++k)
        {
            this->contour[0].vertex[k].x=v[k](0);
            this->contour[0].vertex[k].y=v[k](1);
        }
        
    }
    
    ~GPCpolygon()
    {
        for(int k=0;k<this->num_contours;++k)
        {
            free(this->contour[k].vertex);
        }
        free(this->contour);
    }
    
    std::vector<std::vector<Eigen::Matrix<double,2,1>>> clip(GPCpolygon& clip_p)
    {
        
        gpc_polygon temp_p;
        
        gpc_polygon_clip(GPC_INT,
                         this,
                         &clip_p,
                         &temp_p);
        
        std::vector<std::vector<Eigen::Matrix<double,2,1>>> temp;
        std::cout<<"num_contours="<<temp_p.num_contours<<std::endl;
        for(int k=0;k<temp_p.num_contours;++k)
        {
            temp.push_back(std::vector<Eigen::Matrix<double,2,1>>());
            temp.back().reserve(temp_p.contour[k].num_vertices);
            for(int v=0;v<temp_p.contour[k].num_vertices;++v)
            {
                temp.back().push_back((Eigen::Matrix<double,2,1>()<<temp_p.contour[k].vertex[v].x,temp_p.contour[k].vertex[v].y).finished());
            }
        }
        //        if(temp_p.num_contours>0)
        //        {
        //
        //        }
        return temp;
    }
    
    std::vector<Eigen::Matrix<double,2,1>> unionPoly(GPCpolygon& union_p)
    {
        
        gpc_polygon temp_p;
        
        gpc_polygon_clip(GPC_UNION,
                         this,
                         &union_p,
                         &temp_p);
        
        std::vector<Eigen::Matrix<double,2,1>> temp;
        std::cout<<"num_contours="<<temp_p.num_contours<<std::endl;
        for(int k=0;k<temp_p.num_contours;++k)
        {
            // temp.push_back(std::vector<Eigen::Matrix<double,2,1>>());
            temp.reserve(temp_p.contour[k].num_vertices);
            for(int v=0;v<temp_p.contour[k].num_vertices;++v)
            {
                temp.push_back((Eigen::Matrix<double,2,1>()<<temp_p.contour[k].vertex[v].x,temp_p.contour[k].vertex[v].y).finished());
            }
        }
        //        if(temp_p.num_contours>0)
        //        {
        //
        //        }
        return temp;
    }
    
    std::vector<Eigen::Matrix<double, 2, 1>> XORPoly(GPCpolygon &union_p)
    {
        
        gpc_polygon temp_p;
        
        gpc_polygon_clip(GPC_XOR,
                         this,
                         &union_p,
                         &temp_p);
        
        std::vector<Eigen::Matrix<double, 2, 1>> temp;
        std::cout << "num_contours=" << temp_p.num_contours << std::endl;
        for (int k = 0; k < temp_p.num_contours; ++k)
        {
            // temp.push_back(std::vector<Eigen::Matrix<double,2,1>>());
            temp.reserve(temp_p.contour[k].num_vertices);
            for (int v = 0; v < temp_p.contour[k].num_vertices; ++v)
            {
                temp.push_back((Eigen::Matrix<double, 2, 1>() << temp_p.contour[k].vertex[v].x, temp_p.contour[k].vertex[v].y).finished());
            }
        }
        //        if(temp_p.num_contours>0)
        //        {
        //
        //        }
        return temp;
    }
};



int main(int argc, char * argv[])
{
    
    //    gpc_polygon subject,clip,result;
    //
    //    FILE* file0=fopen("poly0.txt","r");
    //    gpc_read_polygon(file0,0,&subject);
    //
    //    FILE* file1=fopen("poly1.txt","r");
    //    gpc_read_polygon(file1,0,&clip);
    //
    //
    //    gpc_polygon_clip(GPC_INT,
    //                          &subject,
    //                          &clip,
    //                          &result);
    //
    //    FILE* file2=fopen("result.txt","w");
    //    gpc_write_polygon(file2,0,&result);
    
    
    
    GPCpolygon subject(readPolyVector("poly0.txt"));
    GPCpolygon clip(readPolyVector("poly1.txt"));
    std::cout<<"Clipping\n";
    auto result=subject.clip(clip);
    std::cout<<"Clipped size is "<<result.size()<<std::endl;
    // auto result=clip.clip(subject);
    
    if (true /*result.size()==1*/)
    {
        std::ofstream ofs("result.txt");
        ofs << result.size() << "\n";
        for (const auto &v : result)
        {
            ofs << v.size() << "\n";
            for (const auto &pt : v)
            {
                ofs << pt.transpose() << "\n";
            }
        }
    }
    else
    {
        std::vector<Eigen::Matrix<double,2,1>> resulttemp(result[0]);
        std::ofstream ofs("result.txt");
        std::cout << "XORing\n";
        
        for (size_t i=1;i<result.size();i++)
        {
            GPCpolygon fullPoly(resulttemp);
            GPCpolygon temp(result[i]);
            resulttemp=fullPoly.XORPoly(temp);
        }
        ofs << 1 << "\n";
        ofs << resulttemp.size() << "\n";
        
        for (const auto &pt : resulttemp)
        {
            ofs << pt.transpose() << "\n";
        }
        
        std::ofstream ofsOrig("resultOrig.txt");
        ofsOrig << result.size() << "\n";
        for (const auto &v : result)
        {
            ofsOrig << v.size() << "\n";
            for (const auto &pt : v)
            {
                ofsOrig << pt.transpose() << "\n";
            }
        }
    }
    
    
    
    
    
    
    return 0;
}

//void readPoly(const std::string& fileName,ClipperLib::Paths &poly)
//{
//
//    std::ifstream file ( fileName.c_str() , std::ifstream::in );
//    if(file.is_open())
//    {
//        std::string line;
//        double x,y;
//
//        while (std::getline(file, line))
//        {
//            std::stringstream ss(line);
//
//            ss >> x >>y;
//            std::cout<<x<<","<<y<<std::endl;
//
//            poly[0] <<ClipperLib::DoublePoint(x,y);
//        }
//
//        std::cout<<"Polygon:"<<std::endl;
//        for(const auto& point : poly[0])
//        {
//            std::cout<<point.X<<","<<point.Y<<std::endl;
//        }
//    }
//    else
//    {
//        std::cout<<"CANNOT READ "+fileName<<std::endl;
//    }
//
//}
//
//bool SaveToFile(const std::string& filename, ClipperLib::Paths &ppg, double scale = 1.0, unsigned decimal_places = 0)
//{
//    std::ofstream ofs(filename);
//    if (!ofs) return false;
//
//    if (decimal_places > 8) decimal_places = 8;
//    ofs << std::setprecision(decimal_places) << std::fixed;
////    ofs<<"hi"<<std::endl;
////    ClipperLib::Path pg;
//    for (size_t i = 0; i < ppg.size(); ++i)
//    {
//        for (size_t j = 0; j < ppg[i].size(); ++j)
//            ofs << ppg[i][j].X / scale << ", " << ppg[i][j].Y / scale << "," << std::endl;
//        ofs << std::endl;
//    }
//    ofs.close();
//    return true;
//}
//
//int main(int argc, char * argv[])
//{
//    ClipperLib::Paths subj(2), clip(1), solution;
//
//    readPoly("poly0.txt",subj);
//    readPoly("poly1.txt",clip);
//    ClipperLib::Clipper c;
//    c.AddPaths(subj, ClipperLib::ptSubject, true);
//    c.AddPaths(clip, ClipperLib::ptClip, true);
//    c.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
//    SaveToFile("intersection.txt",solution,1,7);
//    return 0;
//}

