/*
 * -------------------------------------------------------------------
 *
 *                        Wavefront-related code.
 *
 *        Copyright (c) 2015, Juan Manuel Solano-Altamirano
 *                                 <jmsolanoalt@gmail.com>
 *
 * -------------------------------------------------------------------
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ---------------------------------------------------------------------
 *
 */


#ifndef _LEGENDRE_POLYNOMIAL_H_
#define _LEGENDRE_POLYNOMIAL_H_

namespace legendrepolynomial {

/** Returns the normalized Associated Legendre polynomial Pnm evaluated at x.  */
double NormalPnm(int n,int m,double x);
/** Returns the normalized Legendre polynomial Pn evaluated at x.  */
inline double NormalPn(int n,double x) {return NormalPnm(n,0,x);}
/** Returns dPn/dt, evaluated at t=x  */
double FirstDerivativeNormalPn(int n,double x);
/** Auxiliar function (used for saving computing time)  */
double FirstDerivativeNormalPnFromPnAndPnplus1(int n,\
      double x,double pn,double pnp1);
/** Auxiliar function (used for saving computing time)  */
double FirstDerivativeNormalPnFromPnAndPnminus1(int n,\
      double x,double pn,double pnm1);
}

#endif /* defined(_LEGENDRE_POLYNOMIAL_H_) */
