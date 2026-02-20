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


#ifndef INIT_H
#define INIT_H
/***********************************************************
*     allocate memory for global arrays and structs        *
*     initialize system parameters                         *
************************************************************/
void allocateMemory(Event *event,TimeMeasures *meas);
/******************************************
 *  Open Global Files                     *
 ******************************************/
void openFiles(TimeMeasures *meas);
/*****************************************************
*   Populates Host and Microbial layer               *
*   Bacterial layer: types of bacteria               *
*   start with fixed fractions equal to 1/TYPES      *
******************************************************/
void initialStateFixedFrac(void);
/*****************************************************
*   Populates Host and Microbial layer               *
*   Bacterial layer: types of bacteria               *
*   are uniformly distributed                        *
*****************************************************/
void initialStateUniD(void);
/**********************************************************
*       Populates Host and Microbial layer                *
*       Bacterial layer: frequencies of each type j       *
*       is the probability of normal                      *
*       distribution with x=(investment[j]-mean)/stdinv   *
***********************************************************/
void initialStateNormD(void);
/*****************************************************
*   Populates Host layer with a single host          *
*   Bacteria yypes are uniformly distributed         *
*****************************************************/
void initialStateSingleH(TimeMeasures *meas);
/*************************************************
* Set Investments for a System with only helpers *
* (different types of helpers,  0<investment<1   *
* (paper version)                                *
**************************************************/
void setInvestmentsPaper(void);
/*************************************************
* Set Investments (0<=investment<=1)             *
**************************************************/
void setInvestments(void);
/****************************************************************************
*                     Build System                                          *
*****************************************************************************/
void setSystem(Event *event,TimeMeasures *meas);
#endif
