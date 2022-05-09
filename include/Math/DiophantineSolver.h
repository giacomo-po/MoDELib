/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>
 *
 * MODEL is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_DiophantineSolver_h_
#define model_DiophantineSolver_h_

#include <Eigen/Dense>
//#include <LatticeGCD.h>
#include <IntegerMath.h>
#include<climits>
namespace model
{
    class DiophantineSolver
    {
        public:
        // Finds such x and y, that a  x + b  y = gcd(a, b)
        //https://github.com/ADJA/algos/blob/master/NumberTheory/DiophantineEquation.cpp
        static long int extended_gcd(long int a, long int b, long int &x, long int &y)
        {
            if (b == 0)
            {
                x = 1;
                y = 0;
                return a;
            }
            long int x1, y1;
            long int g = extended_gcd(b, a % b, x1, y1);
            x = y1;
            y = x1 - (a / b) * y1;
            return g;
        }

        // Solves equation a  x + b  y = c, writes answer to x and y
        static void solveDiophantine2vars(long int a, long int b, long int c, long int &x, long int &y)
        {

            long int g = extended_gcd(a, b, x, y);

            if (c % g != 0)
            {
                puts("Impossible");
                exit(0);
            }

            c /= g;
            // long int temp;
            // if (a<b)
            // {
            //     temp=0;
            // }
            // else
            // {
            //     temp=LLONG_MAX;
            // }
            // x = x * c + temp*b/g;
            // y = y * c - temp*a/g;

            x = x * c ;
            y = y * c ;
        }
        // http://mathafou.free.fr/exe_en/exedioph3.html
        // https://math.stackexchange.com/questions/514105/how-do-i-solve-a-linear-diophantine-equation-with-three-unknowns
        static void solveDiophantine3vars(const Eigen::Matrix<long int, 3, 1> &alphas, const long int t1, Eigen::Matrix<long int, 3, 1> &sols)
        {
            // const auto gcd(LatticeGCD<dim>::gcd(key.reciprocalDirectionComponents().transpose()*N));
//            const long int gcd(LatticeGCD<3>::gcd(alphas));
            const long int gcd(IntegerMath<long int>::gcd(alphas));


            if (t1 % gcd != 0)
            {
                assert(0 && "Integral solution for shifts does not exist");
            }

            const long int a(alphas(0) / gcd);
            const long int b(alphas(1) / gcd);
            const long int c(alphas(2) / gcd);
            const long int t(t1 / gcd);

//            const long int gcdab(LatticeGCD<3>::gcd(a, b));
            const long int gcdab(IntegerMath<long int>::gcd(a, b));

            const long int ap(a / gcdab);
            const long int bp(b / gcdab);
            long int u0, v0;
            solveDiophantine2vars(ap, bp, c, u0, v0);
            long int z0, t0;
            solveDiophantine2vars(c, gcdab, t, z0, t0);
            long int x0, y0;
            solveDiophantine2vars(ap, bp, t0, x0, y0);

            //k and m are arbitrary variables
            const long int k = 0;
            const long int m = 0;

            const long int x(x0 + bp * k - u0 * m);
            const long int y(y0 - ap * k + v0 * m);
            const long int z(z0 + gcdab * m);

            sols << x, y, z;
        }
    };
}

#endif
