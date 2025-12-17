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


#ifndef MEASURES_H
#define MEASURES_H
/**********************************************
*   stores the frequency of beneficial       *
*   bacteria in a single host                *
**********************************************/
void densB1Xt(double tnow,double ti,int idh);
/**********************************************
*   stores the frequency of beneficial       *
*   bacteria in a single host                *
**********************************************/
void averInvestmentXt(int ntnow,double tnow,int nti);
/**************************************
*  snapshot of host network           *
*  The colors indicate the acumulated *
*  investment of each host            *
***************************************/
void save_config(int t, int tf);
/**********************************************
 * calculates investment distribution among   *
 * hosts                                      *
 **********************************************/
void calcInvDist(double *hist_inv, double binsize);
/**********************************************
*   this routine united all time measures     *
*   done inside the time loop                 *
**********************************************/
void measures(TimeMeasures meas);
#endif
