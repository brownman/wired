/******************************************************************************
 * Sukhchander Khanna [ khannas@rpi.edu ]
 * Computer Graphics Fall 04
 * Final Project Component One -- 'WIRED' -- My Roller Coaster Creation
 *
 * File: RollerCoaster.h
 * Desc: Defines functions to build the roller coaster
 *		 Defines variables to be exported into the data loading component
 *****************************************************************************/

#ifndef ROLLER_COASTER_H
#define ROLLER_COASTER_H

#include "Point.h"

void ResizeScene(int,int);

void ChangeRollerCoasterParameters(float,float,float,int,int);
void InitializeRollerCoaster(int,int,int,float,float,float);
void DrawRollerCoaster();

void InitColumns();
void DrawColumns();
void DrawPlatform();
void DrawTrees();

extern point *pPointControl;
extern int nbPointControl;

extern int startSegment;
extern int brakeSegment;
extern float avgSegmentLength;
extern float twistFactor;

extern point *pColumnCoord;
extern int *pColumnCoordinate;
extern int nbColumnCoord;

extern point *pTree;
extern int nbTree;

extern float metalLength;
extern float metalAngle;
extern point metalPosition;

#endif