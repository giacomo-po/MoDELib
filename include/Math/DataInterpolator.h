/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DataInterpolator_H_
#define model_DataInterpolator_H_

#include <assert.h>
#include <string>
#include <stdio.h>
#include <map>

namespace model{
    
    class DataInterpolator
    {
        
        typedef std::map<double,double> mapType;
        
        /**********************************************************************/
        mapType data2map(const std::string& filename)
        {
            FILE * pFile=fopen(filename.c_str(), "r");
            assert(pFile!=NULL && "DataInterpolator CANNOT READ FILE.");
            mapType temp;
            while (true)
            {
                double x(0.0);
                double y(0.0);
                const int nRead(fscanf (pFile, "%lf%lf",  &x, &y));
                if(nRead==2) // read x and y
                {
                    assert(temp.insert(std::make_pair(x,y)).second && "DATA BREAKPOINTS MUST BE UNIQUE.");
                }
                else // stop reading
                {
                    break;
                }
            }
            
            fclose (pFile);
            return temp;
        }
        
    public:
        
        const mapType dataMap;
        const double xMin;
        const double xMax;
        
        /**********************************************************************/
        DataInterpolator(const std::string& filename) :
        /* init list*/ dataMap(data2map(filename)),
        /* init list*/ xMin(dataMap.begin()->first),
        /* init list*/ xMax(dataMap.rbegin()->first)
        {/*! Constructure reads data abd initializes dataMap.
          */
            //std::cout<<"range is ["<<xMin<<","<<xMax<<"]"<<std::endl;
        }
        
        
        /**********************************************************************/
        double operator()(const double& x) const
        {
            assert(x>=xMin && "OUT OF RANGE.");
            assert(x<=xMax && "OUT OF RANGE.");
            double y(0.0);
            const mapType::const_iterator iterLower(dataMap.lower_bound(x));
            const mapType::const_iterator iterUpper(dataMap.upper_bound(x));
            if(iterLower==iterUpper) // x is in between breakpoints. 
            {
                // Because of assertion x<=xMax, iter1 cannot be dataMap.end()
                // Moreover, iterUpper can never be dataMap.begin()
                // Terefore it is safe to use element before iterLower
                mapType::const_iterator iter1(iterLower);
                --iter1;
                y=iter1->second+(iterLower->second-iter1->second)/(iterLower->first-iter1->first)*(x-iter1->first);
            }
            else  // x is on breakpoint. iterUpper may be dataMap.end().
            {
                y=iterLower->second;
            }
            return y;
        }
        
        
    };
    
    
    
} // close namespace model
#endif