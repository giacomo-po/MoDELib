

MatrixDim userStress(const VectorDim& x) const
{
    
    // Center of point defect in simulation box
    const VectorDim P0((VectorDim()<<0.25*(this->network().mesh.xMax()(0)-this->network().mesh.xMin()(0)),0.0,0.0).finished());
    
    
    
    const VectorDim meshSize(this->network().mesh.xMax()-this->network().mesh.xMin());
    
    double tau_13=0.0;
    
    // Consider stress of images
    for(int i=-this->network().dislocationImages_x;i<=this->network().dislocationImages_x;++i)
    {
        for(int j=-this->network().dislocationImages_y;j<=this->network().dislocationImages_y;++j)
        {
            for(int k=-this->network().dislocationImages_z;k<=this->network().dislocationImages_z;++k)
            {
                const Eigen::Matrix<int,3,1> cellID((Eigen::Matrix<int,3,1>()<<i,j,k).finished());
                const VectorDim P=P0+(meshSize.array()*cellID.cast<double>().array()).matrix(); // image of P0 in cell (i,j,k)
                tau_13+=1.0/(x-P).norm(); // dummy stress using x and P
            }
        }
    }
    
    
    
    
    
    MatrixDim temp(MatrixDim::Zero());
    temp(0,2)=tau_13;
    temp(2,0)=temp(0,2);
    
    
    
    return temp;
    
}
