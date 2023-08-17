/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDauxIO_cpp_
#define model_DDauxIO_cpp_

#include <DDauxIO.h>

namespace model
{
    
        template<int dim>
        DDauxIO<dim>::DDauxIO(const std::string& folderName,const std::string& suffix) :
        /* init */ DDbaseIO(folderName,"ddAux",suffix)
        {
            
        }

template<int dim>
void DDauxIO<dim>::clear()
{
    meshNodes().clear();
    periodicGlidePlanePatches().clear();
    quadraturePoints().clear();
}

        template<int dim>
        void DDauxIO<dim>::addPeriodicGlidePlane(const PeriodicGlidePlane<dim>& pgp)
        {
            for(const auto& patch : pgp.patches())
            {
                periodicGlidePlanePatches().emplace_back(*patch.second.lock());
            }
        }

        template<int dim>
        const std::vector<MeshNodeIO<dim>>& DDauxIO<dim>::meshNodes() const
        {
            return *this;
        }
        
        template<int dim>
        std::vector<MeshNodeIO<dim>>& DDauxIO<dim>::meshNodes()
        {
            return *this;
        }
        
//        template<int dim>
//        const std::vector<GlidePlaneIO<dim>>& DDauxIO<dim>::glidePlanes() const
//        {
//            return *this;
//        }
//
//        template<int dim>
//        std::vector<GlidePlaneIO<dim>>& DDauxIO<dim>::glidePlanes()
//        {
//            return *this;
//        }
        
        template<int dim>
        const std::vector<PeriodicPlanePatchIO<dim>>& DDauxIO<dim>::periodicGlidePlanePatches() const
        {
            return *this;
        }
        
        template<int dim>
        std::vector<PeriodicPlanePatchIO<dim>>& DDauxIO<dim>::periodicGlidePlanePatches()
        {
            return *this;
        }
        
        template<int dim>
        const std::vector<DislocationQuadraturePointIO<dim>>& DDauxIO<dim>::quadraturePoints() const
        {
            return *this;
        }
        
        template<int dim>
        std::vector<DislocationQuadraturePointIO<dim>>& DDauxIO<dim>::quadraturePoints()
        {
            return *this;
        }
        
        template<int dim>
        void DDauxIO<dim>::write(const size_t& runID,const bool& outputBinary)
        {
            if(outputBinary)
            {
                writeBin(runID);
            }
            else
            {
                writeTxt(runID);
            }
        }
        
        template<int dim>
        void DDauxIO<dim>::writeTxt(const long int& runID)
        {
            const auto t0=std::chrono::system_clock::now();
            const std::string filename(this->getTxtFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                std::cout<<"Writing "<<filename<<std::flush;
                // Write header
                file<<meshNodes().size()<<"\n";
                file<<quadraturePoints().size()<<"\n";
//                file<<glidePlanes().size()<<"\n";
                file<<periodicGlidePlanePatches().size()<<"\n";
//                file<<periodicLoopNodes().size()<<"\n";
//                file<<periodicLoopLinks().size()<<"\n";

                // Write body
                for(const auto& mPoint : meshNodes())
                {
                    file<<mPoint<<"\n";
                }
                for(const auto& qPoint : quadraturePoints())
                {
                    file<<qPoint<<"\n";
                }
//                for(const auto& gpBnd : glidePlanes())
//                {
//                    file<<gpBnd<<"\n";
//                }
                for(const auto& patch : periodicGlidePlanePatches())
                {
                    file<<patch<<"\n";
                }

                file.close();
                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            else
            {
                std::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
        }
        
        template<int dim>
        void DDauxIO<dim>::writeBin(const size_t& runID)
        {
            const auto t0=std::chrono::system_clock::now();
            const std::string filename(this->getBinFilename(runID));
            std::ofstream file(filename.c_str(), std::ios::out  | std::ios::binary);
            if(file.is_open())
            {
                std::cout<<"Writing "<<filename<<std::flush;
                // Write header
                binWrite(file,meshNodes().size());
                binWrite(file,quadraturePoints().size());
//                binWrite(file,glidePlanes().size());
                binWrite(file,periodicGlidePlanePatches().size());

                // Write body
                for(const auto& mPoint : meshNodes())
                {
                    binWrite(file,mPoint);
                }
                for(const auto& qPoint : quadraturePoints())
                {
                    binWrite(file,qPoint);
                }
//                for(const auto& gpb : glidePlanes())
//                {
//                    binWrite(file,gpb);
//                }
                for(const auto& patch : periodicGlidePlanePatches())
                {
                    binWrite(file,patch);
                }
                
                file.close();
                std::cout<<magentaColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
            }
            else
            {
                std::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
        }
        
        template<int dim>
        void DDauxIO<dim>::readBin(const size_t& runID)
        {
            clear();
            const std::string filename(this->getBinFilename(runID));
            std::ifstream infile (filename.c_str(), std::ios::in|std::ios::binary);
            if(infile.is_open())
            {
                const auto t0=std::chrono::system_clock::now();
                std::cout<<"reading "<<filename<<std::flush;

                // Read header
                size_t sizeMP;
                infile.read (reinterpret_cast<char*>(&sizeMP), 1*sizeof(sizeMP));
                size_t sizeQP;
                infile.read (reinterpret_cast<char*>(&sizeQP), 1*sizeof(sizeQP));
//                size_t sizeGP;
//                infile.read (reinterpret_cast<char*>(&sizeGP), 1*sizeof(sizeGP));
                size_t sizePPP;
                infile.read (reinterpret_cast<char*>(&sizePPP), 1*sizeof(sizePPP));

                // Read body
                meshNodes().resize(sizeMP);
                infile.read (reinterpret_cast<char*>(meshNodes().data()),meshNodes().size()*sizeof(MeshNodeIO<dim>));
                quadraturePoints().resize(sizeQP);
                infile.read (reinterpret_cast<char*>(quadraturePoints().data()),quadraturePoints().size()*sizeof(DislocationQuadraturePointIO<dim>));
//                glidePlanes().resize(sizeGP);
//                infile.read (reinterpret_cast<char*>(glidePlanes().data()),glidePlanes().size()*sizeof(GlidePlaneIO<dim>));
                periodicGlidePlanePatches().resize(sizePPP);
                infile.read (reinterpret_cast<char*>(periodicGlidePlanePatches().data()),periodicGlidePlanePatches().size()*sizeof(PeriodicPlanePatchIO<dim>));

                infile.close();
                printLog(t0);
            }
            else
            {
                std::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
            
        }
        
        template<int dim>
        void DDauxIO<dim>::readTxt(const size_t& runID)
        {
            clear();
            const std::string filename(this->getTxtFilename(runID));
            std::ifstream infile (filename.c_str(), std::ios::in);
            if(infile.is_open())
            {

                const auto t0=std::chrono::system_clock::now();
                std::cout<<"reading "<<filename<<std::flush;
                std::string line;
                std::stringstream ss;

                size_t sizeMP;
                std::getline(infile, line);
                ss<<line;
                ss >> sizeMP;
                ss.str("");
                ss.clear();
                
                size_t sizeQP;
                std::getline(infile, line);
                ss<<line;
                ss >> sizeQP;
                ss.str("");
                ss.clear();
                
//                size_t sizeGP;
//                std::getline(infile, line);
//                ss<<line;
//                ss >> sizeGP;
//                ss.str("");
//                ss.clear();

                size_t sizePPP;
                std::getline(infile, line);
                ss<<line;
                ss >> sizePPP;
                ss.str("");
                ss.clear();

//                meshNodes().clear();
                meshNodes().reserve(sizeMP);
                for(size_t k=0;k<sizeMP;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    meshNodes().emplace_back(ss);
                    ss.str("");
                    ss.clear();
                }
                
//                quadraturePoints().clear();
                quadraturePoints().reserve(sizeQP);
                for(size_t k=0;k<sizeQP;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    quadraturePoints().emplace_back(ss);
                    ss.str("");
                    ss.clear();
                }
                
//                glidePlanes().clear();
//                glidePlanes().reserve(sizeGP);
//                for(size_t k=0;k<sizeGP;++k)
//                {
//                    std::getline(infile, line);
//                    ss<<line;
//                    glidePlanes().emplace_back(ss);
//                    ss.str("");
//                    ss.clear();
//                }
                
//                periodicGlidePlanePatches().clear();
                periodicGlidePlanePatches().reserve(sizePPP);
                for(size_t k=0;k<sizePPP;++k)
                {
                    std::getline(infile, line);
                    ss<<line;
                    periodicGlidePlanePatches().emplace_back(ss);
                    ss.str("");
                    ss.clear();
                }
                
                infile.close();
                printLog(t0);
            }
            else
            {
                std::cout<<"CANNOT OPEN "<<filename<<std::endl;
                assert(false && "CANNOT OPEN FILE.");
            }
            
        }
        
        template<int dim>
        void DDauxIO<dim>::printLog(const std::chrono::time_point<std::chrono::system_clock>& t0) const
        {
            std::cout<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<std::endl;
            std::cout<<"  "<<meshNodes().size()<<" meshNodes "<<std::endl;
//            std::cout<<"  "<<glidePlanes().size()<<" glidePlanes "<<std::endl;
            std::cout<<"  "<<periodicGlidePlanePatches().size()<<" periodicPlanePatches"<<std::endl;
            std::cout<<"  "<<quadraturePoints().size()<<" quadraturePoints"<<std::endl;
        }

        template<int dim>
        void DDauxIO<dim>::read(const size_t& runID)
        {
            if(isBinGood(runID))
            {
                readBin(runID);
            }
            else
            {
                if(isTxtGood(runID))
                {
                    readTxt(runID);
                }
                else
                {
                    std::cout<<"COULD NOT FIND INPUT FILEs evl/evl_"<<runID<<".bin or evl/evl_"<<runID<<".txt"<<std::endl;
                    assert(0 && "COULD NOT FIND INPUT FILEs.");
                }
            }
        }
    
template struct DDauxIO<3>;

}
#endif

