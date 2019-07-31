/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2012 by Giacomo Po <gpo@ucla.edu>
 * Copyright (C) 2012 by Shao-Ching Huang <sch@ucla.edu>
 * Copyright (C) 2012 by Tajendra Singh <tvsingh@ucla.edu>
 * Copyright (C) 2012 by Tamer Crosby <tcrosby@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef _model_LPTpartitioner_h_
#define _model_LPTpartitioner_h_

#include <map>
#include <vector>
#include <deque>
#include <algorithm>    // std::min_element, std::max_element
#include <assert.h>
#include <iostream>
#include <iomanip>

namespace model {
    
    
    /*!\brief Template implementation of the Longest Processing Time (LPT) 
     * algorithm for task partiotioning.
     */
    template <typename T>
    class LPTpartitioner :
    /* inheritance    */ private std::multimap<double,T*>,
    /* inheritance    */ private std::vector<std::deque<T*> >,
    /* inheritance    */ private std::vector<double>
    {
        
    public:

        typedef std::multimap<double,T*> BaseMapType;
        typedef std::deque<T*> BaseDequeType;
        typedef std::vector<BaseDequeType> BaseDequeVectorType;
        typedef std::vector<double> BaseWeightVectorType;
        
        
        static bool partitionIsValid;

        
        void clear()
        {
            BaseMapType::clear();
            BaseDequeVectorType::clear();
            BaseWeightVectorType::clear();
            partitionIsValid=false;
        }
        
        /**********************************************************************/
        void insert(const double& w,T* pT)
        {
            assert(w>=0 && "WEIGHT MUST BE POSITIVE");
            BaseMapType::insert(std::make_pair(w,pT));
            partitionIsValid=false;
        }
        
        /**********************************************************************/
        void partition(const int& nProc)
        {
            // Resize
            BaseDequeVectorType::resize(nProc);
            BaseWeightVectorType::resize(nProc,0.0);
            
            // Traverse BaseMapType in descending order (highest to lowest weight)
            for (typename BaseMapType::const_reverse_iterator rIter =BaseMapType::rbegin();
                 /*                                        */ rIter!=BaseMapType::rend();
                 /*                                      */ ++rIter)
            {
                // Compute the index of minimum weight in BaseWeightVectorType
                const size_t argMin(min_element(BaseWeightVectorType::begin(), BaseWeightVectorType::end())-BaseWeightVectorType::begin());
                // Add the current weigth to the bin of minimum weigth
                BaseWeightVectorType::operator[](argMin)+=rIter->first;
                // Add the current pointer to the bin of minimum weight
                BaseDequeVectorType::operator[](argMin).push_back(rIter->second);
            }
            
            partitionIsValid=true;

        }
        
        /**********************************************************************/
        double totalWeight() const
        {
            double temp(0.0);
            for (unsigned int k=0;k<BaseWeightVectorType::size();++k)
            {
                temp+=weight(k);
            }
            return temp;
        }


        /**********************************************************************/
        size_t bins() const
        {
            return BaseDequeVectorType::size();
        }
        
        /**********************************************************************/
        const BaseDequeType& bin(const int& k) const
        {
            return BaseDequeVectorType::operator[](k);
        }
        
        /**********************************************************************/
        const double& weight(const int& k) const
        {
            return BaseWeightVectorType::operator[](k);
        }
        
        /**********************************************************************/
        void show() const
        {
 
            std::cout<<"Processor loads [%]: ";
            const double wT(totalWeight());
            for (unsigned int k=0;k<BaseWeightVectorType::size();++k)
            {
                double wr(weight(k)/wT);

                std::cout<<std::fixed << std::setprecision(2)<<wr*100.0<<" ";

//                std::cout<<"processor "<<std::setw(3)<<k<<std::setw(5)<<" ("<<wr*100<<"%): ";
//                int nC(wr*20*BaseWeightVectorType::size());
//                for (int c=0;c<nC;++c)
//                {
//                    std::cout<<"*";
//                }
//                std::cout<<"\n";
            }
            std::cout<<std::endl;

        }
        
        
    };
    
    // declare static data
    template <typename T>
    bool LPTpartitioner<T>::partitionIsValid=false;
    
    
} // end namespace
#endif