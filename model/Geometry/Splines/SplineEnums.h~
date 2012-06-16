/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2011 by Giacomo Po <gpo@ucla.edu>.
 * Copyright (C) 2011 by Benjamin Ramirez<ramirezbrf@gmail.com>.
 *
 * model is distributed without any warranty under the 
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

/////////////////////////////////////////////////
// Common Public Enumerators
public:
enum  {Ncoeff= 2*(corder+1)};
enum  {pOrder= 2*corder+1};
enum  {Ndof  = dim*Ncoeff};
enum  {eigenSize=pOrder*pOrder};

typedef Eigen::Matrix<double, 1, Ncoeff> RowNcoeff;
typedef Eigen::Matrix<double, 1, Ncoeff-1> RowNcoeffu;
typedef Eigen::Matrix<double, 1, Ncoeff-2> RowNcoeffuu;
typedef Eigen::Matrix<double, Ncoeff, Ncoeff> MatrixNcoeff;
typedef Eigen::Matrix<double, Ncoeff, dim> MatrixNcoeffDim;
typedef Eigen::Matrix<double, Ndof,1> VectorNdof;
typedef Eigen::Matrix<double,dim,dim> MatrixDim;
typedef Eigen::Matrix<double,dim,1>   VectorDim;
typedef Eigen::Matrix<double, dim, Ndof> MatrixDimNdof;

typedef Eigen::Matrix<double,Ndof,Ndof>	MatrixNdof;

