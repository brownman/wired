/******************************************************************************
 * Sukhchander Khanna [ khannas@rpi.edu ]
 * Computer Graphics Fall 04
 * Final Project Component One -- 'WIRED' -- My Roller Coaster Creation
 *
 * File: Point.h
 * Desc: Defines functions that are associated with matrices and vectors
 *		 Defines 2D coordinate --> struct coord2d
 *		 Defines 3D coordinate --> struct point
 *****************************************************************************/

#ifndef POINT_H
#define POINT_H

// some math libraries don't define this
#ifndef M_PI
#define M_PI 3.14159265359f
#endif

typedef struct {
	float x,y;
} coord2d;

typedef struct {
	float x,y,z;
} point;

void Multiply(point*, float, point*);
void AdditiveMultiply(point*, float, point*);

void MatrixMultiply(float*, float*, float*);
void VectorMultiply(point*, float*, point*);

float Magnitude(point*);
float Length(point*);
void NormalizePoint(point*);

void Add(point*, point*, point*);
void Subtract(point*, point*, point*);

float ScalarProduct(point*, point*);
void VectorProduct(point*, point*, point*);

void RotatePoint(point*, point*, float, point*);

int Intersection(point*, point*, point*, point*, point*);

#endif