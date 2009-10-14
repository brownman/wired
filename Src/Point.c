/******************************************************************************
 * Sukhchander Khanna [ khannas@rpi.edu ]
 * Computer Graphics Fall 04
 * Final Project Component One -- 'WIRED' -- My Roller Coaster Creation
 *
 * File: Point.c
 * Desc: Implements functions that are associated with matrices and vectors
 *****************************************************************************/


#include <math.h>
#include "Point.h"


/******************************************************************************
 * multiply by scaling p by t and storing in r
 *****************************************************************************/
void Multiply(point *r, float t, point *p)
{
	r->x = t * p->x;
	r->y = t * p->y;
	r->z = t * p->z;
}


/******************************************************************************
 * multiply by scaling p by t and adding to the previous point r
 *****************************************************************************/
void AdditiveMultiply(point *r, float t, point *p)
{
	r->x += t * p->x;
	r->y += t * p->y;
	r->z += t * p->z;
}


/******************************************************************************
 * multiply by scaling a by b and storing in r
 *****************************************************************************/
void MatrixMultiply(float *r, float* a, float *b)
{
	int i,j;
	for (i=0; i<16; i++) r[i]=0.0f;		// initialize
	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)				// scale values in a by b
		{
			r[i+4*j] =  a[i]	*	b[4*j]		+ 
						a[i+4]	*	b[4*j+1]	+ 
						a[i+8]	*	b[4*j+2]	+ 
						a[i+12] *	b[4*j+3];
		}
	}
}


/******************************************************************************
 * multiply by normalizing m by w and storing in r
 *****************************************************************************/
void VectorMultiply(point *r, float *m, point *v)
{
	float w;
	// magnitude of vector
	w = m[3]*v->x + m[7]*v->y + m[11]*v->z + m[15];
	r->x = (m[0]*v->x + m[4]*v->y + m[8]*v->z + m[12]) / w;
	r->y = (m[1]*v->x + m[5]*v->y + m[9]*v->z + m[13]) / w;
	r->z = (m[2]*v->x + m[6]*v->y + m[10]*v->z + m[14]) / w;
}



/******************************************************************************
 * return the length of point
 *****************************************************************************/
float Length(point *p)
{
	return p->x*p->x + p->y*p->y + p->z*p->z;
}


/******************************************************************************
 * return the magnitude of a point
 *****************************************************************************/
float Magnitude(point *p)
{
	return (float)sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
}


/******************************************************************************
 * normalize point by the absolute magnitude
 *****************************************************************************/
void NormalizePoint(point *p)
{
	float n = Magnitude(p);
	if(n == 0.0f) return;
	p->x = p->x / n;
	p->y = p->y / n;
	p->z = p->z / n;
}


/******************************************************************************
 * subtract b from a and store in r
 *****************************************************************************/
void Subtract(point *r, point *a, point *b)
{
	r->x = a->x - b->x;
	r->y = a->y - b->y;
	r->z = a->z - b->z;
}


/******************************************************************************
 * add b to a and store in r
 *****************************************************************************/
void Add(point *r, point *a, point *b)
{
	r->x = a->x + b->x;
	r->y = a->y + b->y;
	r->z = a->z + b->z;
}


/******************************************************************************
 * scalar product / scale a by b
 *****************************************************************************/
float ScalarProduct(point *a, point *b)
{
	return a->x * b->x + a->y * b->y + a->z * b->z;
}


/******************************************************************************
 * vector product / cross product ( r = a x b )
 *****************************************************************************/
void VectorProduct(point *r, point *a, point *b)
{
	r->x = a->y * b->z - a->z * b->y;
	r->y = a->z * b->x - a->x * b->z;
	r->z = a->x * b->y - a->y * b->x;
}


/******************************************************************************
 * rotate point v around p scale by t and store in r
 *****************************************************************************/
void RotatePoint(point *r, point *v, float t, point *p)
{
	r->x = ( v->x * v->x + cos(t) * (1.0f - (v->x * v->x)) ) * p->x +
	       ( v->x * v->y * (1.0f - cos(t)) - v->z * sin(t) ) * p->y +
	       ( v->z * v->x * (1.0f - cos(t)) + v->y * sin(t) ) * p->z;
	r->y = ( v->x * v->y * (1.0f - cos(t)) + v->z * sin(t) ) * p->x +
	       ( v->y * v->y + cos(t) * (1.0f - (v->y * v->y)) ) * p->y +
	       ( v->y * v->z * (1.0f - cos(t)) - v->x * sin(t) ) * p->z;
	r->z = ( v->z * v->x * (1.0f - cos(t)) - v->y * sin(t) ) * p->x +
	       ( v->y * v->z * (1.0f - cos(t)) + v->x * sin(t) ) * p->y +
	       ( v->z * v->z + cos(t) * (1.0f - (v->z * v->z)) ) * p->z;
}


/******************************************************************************
 * calcuclates if the line a in direction v intersects plane p with normal n
 * return 0 if the line is parallel with the plane
 *****************************************************************************/
int Intersection(point *r, point *a, point *v, point *p, point *n)
{
	point b;
	float t;
	Subtract(&b,p,a);
	t = ScalarProduct(v,n);
	if(t == 0.0f) return 0;
	t = ScalarProduct(&b,n) / t;
	r->x = a->x + t * v->x;
	r->y = a->y + t * v->y;
	r->z = a->z + t * v->z;
	return 1;
}