/*  
 * This file is part of the project Multilevel Evolution.
 *
 * This program uses the GNU Scientific Library (GSL), which is licensed 
 * under the GNU GPL v3 or later. Therefore, this file and any binary 
 * linked with GSL must be distributed under the terms of the GPL v3+.
 *
 * Copyright (C) 2025 Alessandra F. Lütz
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef BAC_H
#define BAC_H
/********************************************
*     GSL-compatible RHS function           *
*********************************************/
int bac_rhs(double t, const double y[], double f[], void *params);

/********************************************
*     GSL-compatible jacobian               *
*********************************************/
int bac_jac(double t, const double y[], double *dfdy,double *dfdt, void *params);
/********************************************************************
*              Euler integration for migration terms                *
*********************************************************************/
void eulerBacMig(int idh,int nviz,int *nlist, double dt,double **y);
/********************************************
*     Buil a GSL  ODE system struct         *
*********************************************/
void bac_make_system(gsl_odeiv2_system *sys, SysParams *P);
/*************************************************************
 *          Bacteria Dynamics                                *
 *************************************************************/
void bacDynamics(gsl_odeiv2_driver *driver,Event *event);
#endif
