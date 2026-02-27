#ifndef INIT_H
#define INIT_H
/***********************************************************
*     allocate memory for global arrays and structs        *
*     initialize system parameters                         *
************************************************************/
void allocateMemory(Event *event,TimeMeasures *meas);
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
 *                     Set Initial Conditions                               *
 ***************************************************************************/
void setCI(TimeMeasures *meas);
/****************************************************************************
*                     Build System                                          *
*****************************************************************************/
void setSystem(Event *event,TimeMeasures *meas);
#endif
