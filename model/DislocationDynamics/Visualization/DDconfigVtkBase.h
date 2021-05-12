/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DDconfigVtkBase_h_
#define model_DDconfigVtkBase_h_



namespace model
{
    struct DDconfigVtkBase
    {
        
        enum ColorScheme {colorBurgers=0,colorSessile=1,colorNormal=2,colorEdgeScrew=3,colorComponent=4};
        
        static float alpha;
        static float tubeRadius;
        static bool scaleRadiusByBurgers;
        static ColorScheme clr;
        static size_t Np;      // No. of vertices per line
        static bool showBoundarySegments;
        static bool blackGrainBoundarySegments;
        static bool showVelocities;
        static bool showNodeIDs;
        static float velocityFactor;
        static bool showZeroBuergers;
        static bool showSingleNode;
        static size_t singleNodeID;
        static bool showNodes;
        static bool showSlippedArea;
        static float slippedAreaOpacity;
        static unsigned char nodeClr[4][3];
        
    };
    
    float DDconfigVtkBase::tubeRadius=5.0;
    bool DDconfigVtkBase::scaleRadiusByBurgers=false;
    float DDconfigVtkBase::alpha=0.5;
    DDconfigVtkBase::ColorScheme DDconfigVtkBase::clr=DDconfigVtkBase::colorBurgers;
    size_t DDconfigVtkBase::Np=2;
    bool DDconfigVtkBase::showBoundarySegments=false;
    bool DDconfigVtkBase::blackGrainBoundarySegments=false;
    bool DDconfigVtkBase::showVelocities=false;
    bool DDconfigVtkBase::showNodeIDs=false;
    float DDconfigVtkBase::velocityFactor=100.0;
    bool DDconfigVtkBase::showZeroBuergers=false;
    bool DDconfigVtkBase::showSingleNode=false;
    size_t DDconfigVtkBase::singleNodeID=0;
    bool DDconfigVtkBase::showNodes=false;
    bool DDconfigVtkBase::showSlippedArea=false;
    float DDconfigVtkBase::slippedAreaOpacity=0.1;
    unsigned char DDconfigVtkBase::nodeClr[4][3]={{100,100,100},{0,255,255},{255,0,255},{1,1,1}};
}
#endif

