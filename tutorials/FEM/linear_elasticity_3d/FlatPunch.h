
struct FlatPunch
{
    //! z-coordinate of face at which FlatPunch is applied
    const double zFace;

    //! Height of sample in z direction
    const double Lz;

    //! Height of sample in z direction
    const double geometricTol;
    
    //! side length of the FlatPunch
    static double punchSize;
    
    //! Initial value of the punch displacement
    static double displacement;

    
    /**************************************/
    template <typename FiniteElementType>
    FlatPunch(const FiniteElementType& fe) :
    /* init list */ zFace(fe.xMax()(2)), // zFace is the one at x(2)= max of finite element node coordinates
    /* init list */ Lz(fe.xMax()(2)-fe.xMin()(2)),
    /* init list */ geometricTol(0.001)
    {
    }
    
    /**************************************/
    template <typename NodeType>
    bool operator()(const NodeType& node) const
    {/*!@param[in] node a finite element node
      *\returns true if the node coordinate is under the FlatPunch
      */
        return (   std::fabs(node.p0(0))<(0.5*punchSize+geometricTol)
        /*   */ && std::fabs(node.p0(1))<(0.5*punchSize+geometricTol)
        /*   */ && std::fabs(node.p0(2)-zFace)<geometricTol);
        
    }
    
    /**************************************/
    template <typename NodeType>
    double at(const NodeType&) const
    {
        
        return displacement;
    }
    
    
};

// Static data
double FlatPunch::punchSize=400.0;

double FlatPunch::displacement=-1.0;
