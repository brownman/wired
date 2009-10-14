/******************************************************************************
 * Sukhchander Khanna [ khannas@rpi.edu ]
 * Computer Graphics Fall 04
 * Final Project Component One -- 'WIRED' -- My Roller Coaster Creation
 *		I'm a huge roller coaster fan.
 *		One winter break I tried designing a roller coaster in SolidWorks.
 *		I created a very simple roller coaster which had the correct geometry.
 *		I wanted to take that idea further by programming it in OpenGL.
 *		This project is the culmination of that idea.
 *		If I were to design a roller coaster for Bolliger & Mallibard, the
 *		company that has created rides like Batman and Hulk, I would design
 *		the following ride. The ride would be called WIRED. After riding it
 *		you will be wired. The roller coaster reaches new heights, new speeds,
 *		and a new experience, that most coasters don't achieve today.
 *		I'm sure a coaster similar to this one will be made.
 *
 *      "SK Productions Presents 'WIRED'"
 *      SK = my name
 *      WIRED = the new roller coaster ride that should be coming to an 
 *				amusement park near you in the future!
 *
 * File: RollerCoaster.c
 * Desc: all components of the scene are created here
 *		 the sky is drawn
 *		 the grass ground is created
 *		 the closed bezier curve is created (major coaster rail)
 *		 the bonds on the closed bezier curve are created
 *		 the outer bezier curves are created (outer rails joined by bonds)
 *		 the platform is drawn
 *		 the columns/support structures are drawn
 *		 the trees are drawn
 *		 the panorama is started
 *		 the ride begins upon completing of the panorama
 *****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include <windows.h>
#include <mmsystem.h>

#include "Point.h"
#include "ImageLoader.h"
#include "RollerCoaster.h"
#include "RollerCoasterLoader.h"


/******************************************************************************
 * if on a *nix system use this function to get the current system time
 * this returns the time in msecs
 *****************************************************************************/
#ifndef _WIN32
#include <sys/time.h>
unsigned long timeGetTime()
{
	struct timeval tv;
	gettimeofday(&tv,NULL);
	return tv.tv_sec*1000 + tv.tv_usec/1000;
}
#endif

// window info
int wndWidth, wndHeight;

// scene panning info
GLfloat roll_angle, pitch_angle;

// comments to the right
float aperture = 90.0f;				// field of view angle in degree
float near_clip = 0.01f;			// near clip plane coord
float far_clip = 90.0f;				// distance clip plane coord

int startSegment = -10;             // index of segment where the train starts
int brakeSegment = -30;             // index of segment where the train starts braking
float avgSegmentLength = 0.15f;		// length of segments that the track is made of
float twistFactor = 7.0f;           // amount of bank

float eye_sep = 0.5f;				// distance between eyes
float focallength = 20.0f;			// depth length to the object
int doIntro = 1;					// Show coaster before starting the ride

int nbPointControl;					// Bezier surface control points
int nbLine;							// number of total segments
int nbDimensions = 6;				// dimensions of the cylinder of the column, this will make a smooth cylinder
int nbBonds;						// number of bonds between columns

point *pPointControl;				// table of control points that are in the file
point *pControl;					// table of control points of the Bezier curve
point *pLine;						// table of the segment ends
point *pPos;						// table of the vector positions of the curve
point *pTrajectory;					// table of the vector trajectories of the curve
point *pCylinder;					// table of the coordinates of the verticies of the cylinders
point *pCylindern;					// table of the normal coordinates to the verticies of the cylinders
float *pCurve;						// table of the curves projected on the xy plan
point *pLine1, *pLine2;				// table of points for the outer lines
point *pRail1, *pRail2;				// table of points for the outer rails
point *pRail1n, *pRail2n;			// table of the normals of the outer coaster rails
point *pBonds;						// table of points where bonds are to be placed
point *pnBonds;						// table of number of bonds
point *pBondsn;						// table of the normals to the 
point *pTangent;					// tangents of the trajectory of each point
float *pNormal;						// lengths of the normal segments

// light definitions
float lightCylinder = 0.02f;
float lightBond = 0.01f;
float lengthBond = 0.1f;

float maxZ = 0.0f;					// maximum z coordinate of the curve
float maxDist = 0.0f;				// maximum distance of the curve from the origin

// textures for the coaster platform and the grass field
glBmpImage grass;
glBmpImage metal;
glBmpImage tree;


/******************************************************************************
 * load the necessary textures -- grass for ground, metal for platform,
 *								  trees for trees
 * this function relies on the ImageLoader library
 * this function also reduces the manual mapping and maintaining a list of the
 * textures loaded into the application
 *****************************************************************************/
GLvoid LoadTextures()
{
	// Grass
	glBmpInit(&grass);
	if (!glBmpLoadImage(&grass,"Data\\Grass.bmp"))
	{
		fprintf(stderr,"Error Loading 'Grass'\n");
		exit(5);
	}
	glBmpSetFilter(&grass,GL_LINEAR_MIPMAP_NEAREST,GL_LINEAR);
	glBmpSetTextureWrap(&grass,GL_REPEAT,GL_REPEAT);
	glBmpGenTextureMipMap(&grass);

	// Metal Platform
	glBmpInit(&metal);
	if (!glBmpLoadImage(&metal,"Data\\Metal.bmp"))
	{
		fprintf(stderr,"Error Loading 'Metal'\n");
		exit(5);
	}
	glBmpSetFilter(&metal,GL_LINEAR_MIPMAP_NEAREST,GL_LINEAR);
	glBmpSetTextureWrap(&metal,GL_REPEAT,GL_REPEAT);
	glBmpGenTextureMipMap(&metal);

	// Tree
	glBmpInit(&tree);
	if (!glBmpLoadImage(&tree,"Data\\Tree.tga"))
	{
		fprintf(stderr,"Error Loading 'Tree'\n");
		exit(5);
	}
	glBmpSetFilter(&tree,GL_LINEAR_MIPMAP_NEAREST,GL_LINEAR);
	glBmpSetTextureWrap(&tree,GL_CLAMP,GL_CLAMP);
	glBmpGenTextureMipMap(&tree);
}


/******************************************************************************
 * initialize the scene
 * basic steps
 *		set the bg color to that of sky color
 *		enable shading
 *		enable culling to show material properties on the coaster rails
 *****************************************************************************/
void InitScene()
{
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glClearDepth(1.0);

	// sky color
	glClearColor(0.3f, 0.3f, 0.7f, 0.0f);
	
	// enable some options
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_BLEND);
	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);

	glDisable(GL_LIGHTING);
	glDisable(GL_DITHER);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glPolygonMode(GL_FRONT,GL_FILL);
	glFrontFace(GL_CCW);

	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


/******************************************************************************
 * resize the scene for the width and height
 * adjust the viewing frustrum also
 *****************************************************************************/
void ResizeScene(int w, int h)
{
	float ratio, wd2;

	if (h == 0) h=1;

	wndWidth = w;
	wndHeight = h;
	glViewport(0,0,w,h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	ratio = w/(float)h;
	wd2 = near_clip * tan(0.5f * M_PI * aperture / 180.0f);
	glFrustum(-ratio*wd2,  // left
              ratio*wd2,   // right
              -wd2,        // bottom
              wd2,         // top
              near_clip, far_clip);

	glMatrixMode(GL_MODELVIEW);
}


/******************************************************************************
 ******************************************************************************
 * the algorithms below were outlined in OpenGL for CAD Applications
 * the implementation for the coaster curves were made using fourt control pts
 * the algorithms use the Bezier principle of parametrizing a curve
 * also the vector calculations were solved using Vector Dynamics for Engineers
 * i used that book for Engineering Dynamics last semester
 * since that course dealt with 3D objects and 3D rotations I was able to
 * better understand the principles of the motion of the roller coaster
 *****************************************************************************
 *****************************************************************************/


/******************************************************************************
 * calculate the length of a Bézier path by building it with 1000 segments
 * this algorithm was outlined in OpenGL for CAD Applications
 * this algorithm uses the Bezier principle of parametrizing a curve
 *****************************************************************************/
float BezierSegmentLegth(point *p1, point *p2, point *p3, point *p4)
{
	point a,b;
	float t,f;
	float d = 1.0f/1000.0f;
	a = *p1;
	f = 0.0f;
	for (t=d; t<=1.0f; t+=d)
	{
		// parametrize a curve with four control points
		// each step reprsents the coefficient of the Bezier equation
		Multiply(&b, (1-t)*(1-t)*(1-t), p1);
		AdditiveMultiply(&b, 3*t*(1-t)*(1-t), p2);
		AdditiveMultiply(&b, 3*t*t*(1-t), p3);
		AdditiveMultiply(&b, t*t*t, p4);
		Subtract(&a,&b,&a);
		f += Magnitude(&a);
		a = b;
	}
	return f;
}


/******************************************************************************
 * calculate the inverse radius of curvature of the Bézier curve I projected 
 * on the xy plan with x-coord T
 * this algorithm was outlined in OpenGL for CAD Applications
 * this algorithm outlines the projection of a curve onto a normal surface
 * the steps are as follows:
 *		project the vector acceleration on a unit vector orthoganal with the 
 *		trajectory using a scalar product
 *		divide the result by the square length of the flight path vector
 *		1 / R = ScalarProduct(&n,&g) / Magnitude(&v); where N is orthgnl to V
 *****************************************************************************/
float ProjectionOfCurve(int i, float t)
{
	point v,g;
	v.z = g.z = 0.0f;

	// n is projected orthoganly onto V
	Multiply(&v, -3*t*t+6*t-3, &pControl[3*i]);
	AdditiveMultiply(&v, 9*t*t-12*t+3, &pControl[3*i+1]);
	AdditiveMultiply(&v, -9*t*t+6*t, &pControl[3*i+2]);
	AdditiveMultiply(&v, 3*t*t, &pControl[3*i+3]);

	// g is the scaling factor
	Multiply(&g, -6*t+6, &pControl[3*i]);
	AdditiveMultiply(&g, 18*t-12, &pControl[3*i+1]);
	AdditiveMultiply(&g, -18*t+6, &pControl[3*i+2]);
	AdditiveMultiply(&g, 6*t, &pControl[3*i+3]);
	
	v.z = g.z = 0.0f;

	// formula : 1/R = det(v,g)/|v|^3 */
	// case: when trajectory is close to the vertical */
	if (Length(&v) < 0.5f) return 0.0f;
	
	return (v.x * g.y - v.y * g.x) / Length(&v) / Magnitude(&v);
}


/******************************************************************************
 * initialize the curve which represents the rails of the coaster
 * this involves initializing the segments that make each Bezier patch
 * and connecting the patches together
 *****************************************************************************/
void InitCurve(float segmentLength)
{
	int i,j,k;
	point a,b,c;
	float t,d;
	float *pnbsegment;

	// initialize pControl starting from pPointControl
	pControl = (point*) malloc((3*nbPointControl+1)*sizeof(point));
	for (i=0; i<nbPointControl; i++) 
	{
		pControl[3*i] = pPointControl[2*i];
		pControl[3*i+1] = pPointControl[2*i+1];
	}
	
	// calculate the intermediate control points
	// here the curve is created along the path
	// the control points are loaded
	// each control point is moved along
	// the necessary projection is calculated using a scalar
	pControl[3*nbPointControl] = pControl[0];
	for (i=0; i<3*nbPointControl; i+=3)
	{
		b = pControl[i+1];
		if (i>0) 
			j = i-3; 
		else 
			j = 3*nbPointControl-3;
		
		if (Magnitude(&b)<=0.1f) 
			Subtract(&b,&pControl[i+3],&pControl[j]);
		
		// point the curve in the correct direction
		// the vector is pointing inward
		NormalizePoint(&b);
		Subtract(&c,&pControl[i+3],&pControl[i]);
		if (ScalarProduct(&b,&c)<0.0f)
		{
			b.x = -b.x; 
			b.y = -b.y; 
			b.z = -b.z; 
		}

		// the vector is pointing inward again
		Multiply(&a,Magnitude(&c)*M_PI/5.5f/1.41421f,&b);
		Add(&pControl[i+1],&pControl[i],&a);
		b.x = -b.x; 
		b.y = -b.y; 
		b.z = -b.z;

		// finally adjust the # of ctrl points
		// and join the curve
		Subtract(&c,&pControl[i],&pControl[j]);
		Multiply(&a,Magnitude(&c)*M_PI/5.5f/1.41421f,&b);
		if (i>0) 
			k = i-1; 
		else 
			k = 3*nbPointControl-1;
		Add(&pControl[k],&pControl[i],&a);
	}

	// calculate the number of segments for each Bézier patch
	pnbsegment = (float*)malloc(nbPointControl*sizeof(float));
	nbLine = 0;
	for (i=0; i<nbPointControl; i++)
	{
		j = i * 3;
		t = BezierSegmentLegth(&pControl[j],&pControl[j+1],&pControl[j+2],&pControl[j+3]);
		pnbsegment[i] = t / segmentLength;
		if (!pnbsegment[i]) pnbsegment[i] = 1;
		nbLine += pnbsegment[i] + 1;
	}

	// calculate the unit of the Bezier patch
	// discretize the curve
	pLine = (point*)malloc(nbLine*sizeof(point));
	k = 0;
	for (i=0; i<nbPointControl; i++)
	{
		a = pControl[3*i];
		pLine[k++] = a;
		d = 1.0f / (pnbsegment[i]+1.0f);
		for (j=1,t=d; j<pnbsegment[i]; j++,t+=d)
		{
			// bezier equation for four points
			// parametrize the curve along these points
			Multiply(&b, (1-t)*(1-t)*(1-t), &pControl[3*i]);
			AdditiveMultiply(&b, 3*t*(1-t)*(1-t), &pControl[3*i+1]);
			AdditiveMultiply(&b, 3*t*t*(1-t), &pControl[3*i+2]);
			AdditiveMultiply(&b, t*t*t, &pControl[3*i+3]);
			pLine[k] = b;
			
			if (b.z > maxZ) 
				maxZ = b.z;
			if (b.x * b.x + b.y * b.y > maxDist) 
				maxDist = b.x * b.x + b.y * b.y;

			k++;
		}
	}
	// normalize the max distance of the patch from the origin
	maxDist = sqrt(maxDist);

	// free the allocated memory
	free(pnbsegment);
}


/******************************************************************************
 * intersection of a cylinder with axis v with the plane passing p of normal N
 * this will calculate a ring around a line/axis in the corresponding plane
 * basically a mesh is drawn around the line/axis which is filled with the
 * corresponding points of the cylinder
 *****************************************************************************/
void CalculateRing(point *p, point *z, point *v, int nbDimensions, float ray, float angle, point *n, point *pMesh)
{
	point r;
	int i;
	// six dimensions for smooth cylinder
	for (i=0; i<nbDimensions; i++)
	{
		// rotate r around v using diameter equation
		RotatePoint(&r, v, (angle+(float)i)*2.0f*M_PI/(float)nbDimensions, z);
		NormalizePoint(&r);
		Multiply(&r,ray,&r);
		Add(&r,p,&r);
		// calculate intersections of the ring drawn to the line/axis
		Intersection(&r,&r,v,p,n);
		// fill the mesh with the rotated/calculated ring
		pMesh[i] = r;
	}
}


/******************************************************************************
 * create a pipe around the line/axis
 * above the ring around the line/axis was made
 * use the function above to fill out the entire ring with points
 *****************************************************************************/
void CalculatePipe(point *pLine, int nbLine, int nbDimensions, float ray, float angle, point *pMesh)
{
	point a,b,c;
	point v,p;
	int i;
	a = pLine[nbLine-1];
	for (i=0; i<nbLine; i++)
	{
		b = pLine[i];
		c = pLine[(i+1)%nbLine];
		
		// normalize axis v
		Subtract(&v,&b,&a);
		NormalizePoint(&v);
		
		// normailize the rotational axis
		Subtract(&c,&c,&b);
		NormalizePoint(&c);
		
		Add(&p,&v,&c);
		c = pPos[i];
		// calcuate a ring with axis v with the plane passing b of normal p
		CalculateRing(&b, &c, &v, nbDimensions, ray, angle+pCurve[i], &p, &pMesh[i*nbDimensions]);
		// a is the same plane as b
		a = b;
	}
}


/******************************************************************************
 * calculate the running average of the number of elements
 * populate table with these averages
 * this average of point locations is performed to ensure a smooth curve
 * along any track especially when there are twists and tilts in the track
 * this was a late addition to the code to ensure the smooth track
 *****************************************************************************/
void Average(float *table, int nbElem, int nbPrev, int nbNext)
{
	float *buf;
	int i,j,k;
	float tot,f,m;

	// allocate some space
	buf = (float*)malloc(nbElem*sizeof(float));
	for (i=0; i<nbElem; i++) buf[i] = table[i];
	
	for (i=0; i<nbElem; i++)
	{
		// create a smooth curve that goes into the tilt
		k = i;
		tot = 0.0f, m = 0.0f;
		for (j=0; j<nbPrev; j++)
		{
			f = (float)j / nbPrev;
			f = 1.0f - f*f;
			m += f;
			tot += f * buf[k];
			k = (k-1+nbElem) % nbElem;
		}

		// create a smooth curve that exits the tilt
		k = (i+1) % nbElem;
		for (j=1; j<nbNext; j++)
		{
			f = (float)j / nbNext;
			f = 1.0f - f*f;
			m += f;
			tot += f * buf[k];
			k = (k+1) % nbElem;
		}
		// populate the average
		table[i] = tot / m;
	}
	free(buf);
}


/******************************************************************************
 * here the curve is converted to line segments so that the coaster can be
 * ridden in a way that makes the geometry mapping simpler that approximating
 * coordinates on a curve
 * the line segments generated from the curve help to place the columns,
 * bonds, and help in the banking of the track
 *****************************************************************************/
void InitLines()
{
	int i,j,k;
	float t;
	point a,b,c;
	point p,v;
	
	// initialization of pCurve and pNormal
	pCurve = (float*)malloc(nbLine*sizeof(float));
	pNormal = (float*)malloc(nbLine*sizeof(float));
	k=0;
	for (i=0; i<nbLine; i++)
	{
		Subtract(&p,&pLine[i],&pLine[(i-1+nbLine)%nbLine]);
		Subtract(&v,&pLine[(i+1)%nbLine],&pLine[i]);
		pNormal[i] = Magnitude(&p);
		NormalizePoint(&p);
		NormalizePoint(&v);
		p.z = v.z = 0;
		VectorProduct(&a,&v,&p);
		// tilt the track along the z-axis
		pCurve[i] = twistFactor * asin(a.z);
	}
	// calculating the next point by doing a running average
	// gives a smoother curve with continuous points in the
	// case of a turn or tilt as it happens on the track
	Average(pCurve,nbLine,7,5);
	
	// initialization of pPos
	pPos = (point*) malloc(nbLine*sizeof(point));
	c.x = 0.0f; 
	c.y = 0.0f; 
	c.z = 1.0f;
	
	// at point j, the curve's upward vector must be +z
	j = 0;
	for (i=0; i<nbLine; i++)
	{
		a = pLine[(j+i-1+nbLine)%nbLine];
		b = pLine[(j+i)%nbLine];
		Subtract(&p,&b,&a);
		NormalizePoint(&p);
		t = -ScalarProduct(&p,&c);
		AdditiveMultiply(&c,t,&p);
		NormalizePoint(&c);
		VectorProduct(&a,&p,&c);
		if (a.z != 0.0f)
		{
			a.z = 0.0f;
			VectorProduct(&c,&a,&p);
			NormalizePoint(&c);
		}
		pPos[(i+j)%nbLine] = c;
	}

	// initialization of pTrajectory and pTangent
	pTrajectory = (point*) malloc(nbLine*sizeof(point));
	pTangent = (point*) malloc(nbLine*sizeof(point));
	a = pLine[nbLine-1];
	for (i=0; i<nbLine; i++)
	{
		b = pLine[i];
		c = pLine[(i+1)%nbLine];
		Subtract(&p,&b,&a);
		NormalizePoint(&p);
		Subtract(&c,&c,&b);
		NormalizePoint(&c);
		Add(&v,&p,&c);
		pTangent[i] = v;
		NormalizePoint(&pTangent[i]);
		Add(&c,&pPos[i],&b);
		Intersection(&a,&c,&p,&b,&v);
		Subtract(&a,&a,&b);
		NormalizePoint(&a);
		pTrajectory[i] = a;
		a = b;
	}

	// initialization of pCylinder
	// calculate the pipe for this cylinder
	pCylinder = (point*) malloc(nbLine*4*sizeof(point));
	CalculatePipe(pLine,nbLine,4,lightCylinder,0.0f,pCylinder);
	
	// initialization of pLine1 and pLine2
	pLine1 = (point*) malloc(nbLine*sizeof(point));
	pLine2 = (point*) malloc(nbLine*sizeof(point));
	a = pLine[nbLine-1];
	for (i=0; i<nbLine; i++)
	{
		b = pLine[i];
		c = pLine[(i+1)%nbLine];
		Subtract(&p,&b,&a);
		NormalizePoint(&p);
		Subtract(&v,&c,&b);
		NormalizePoint(&v);
		t = pCurve[i];
		Add(&p,&v,&p);
		NormalizePoint(&p);
		v = pTrajectory[i];
		RotatePoint(&c,&p,t+1.0f*M_PI/3.0f,&v);
		Multiply(&c,lengthBond,&c);
		Add(&pLine1[i],&c,&b);
		RotatePoint(&c,&p,t-1.0f*M_PI/3.0f,&v);
		Multiply(&c,lengthBond,&c);
		Add(&pLine2[i],&c,&b);
		a = b;
	}

	// initialization of pRail1
	pRail1 = (point*) malloc(nbLine*nbDimensions*sizeof(point));
	CalculatePipe(pLine1,nbLine,nbDimensions,lightCylinder,0.0f,pRail1);
	
	// initialization of pRail2
	pRail2 = (point*) malloc(nbLine*nbDimensions*sizeof(point));
	CalculatePipe(pLine2,nbLine,nbDimensions,lightCylinder,0.0f,pRail2);
}


/******************************************************************************
 * initialize the bonds that connect the inner and outer rails
 * f is the distance between two consecutive  bonds
 *****************************************************************************/
void InitBonds(float f)
{ 
	point a,a1,a2,p,p1,p2,b,v,n;
	float length,l;
	int i;
	int ia,ib;
	float distBond;
	float t,t1,t2;
	length = 0.0f;
	for (i=0; i<nbLine; i++) length += pNormal[i];

	// initialization
	nbBonds = (int)floor(length / f);
	pBonds = (point*)malloc(nbBonds*12*sizeof(point));
	pnBonds = (point*)malloc(nbBonds*3*sizeof(point));
	distBond = length / (float)nbBonds;

	l = 0.0f;
	ia = nbLine-1;
	ib = 0;
	Subtract(&p,&pLine[ib],&pLine[ia]);
	Subtract(&p1,&pLine1[ib],&pLine1[ia]);
	Subtract(&p2,&pLine2[ib],&pLine2[ia]);
	
	for (i=0; i<nbBonds; i++)
	{
		while(l > pNormal[ib])
		{
			l -= pNormal[ib];
			ia = (ia + 1) % nbLine;
			ib = (ib + 1) % nbLine;
			Subtract(&p,&pLine[ib],&pLine[ia]);
			Subtract(&p1,&pLine1[ib],&pLine1[ia]);
			Subtract(&p2,&pLine2[ib],&pLine2[ia]);
		}

		// move along the bond geometry

		// first point
		t = l / pNormal[ib];
		a = pLine[ia];
		AdditiveMultiply(&a,t,&p);
		pnBonds[i*3] = a;
		
		// second point
		a1 = pLine1[ia];
		Subtract(&v,&a,&a1);
		t1 = ScalarProduct(&v,&p) / ScalarProduct(&p1,&p);
		AdditiveMultiply(&a1,t1,&p1);
		pnBonds[i*3+1] = a1;
		
		// third point
		a2 = pLine2[ia];
		Subtract(&v,&a,&a2);
		t2 = ScalarProduct(&v,&p) / ScalarProduct(&p2,&p);
		AdditiveMultiply(&a2,t2,&p2);
		pnBonds[i*3+2] = a2;
		
		// normalize the normal vector
		Subtract(&n,&a1,&a2);
		NormalizePoint(&n);
		
		// normalize the axis vector
		Subtract(&v,&a1,&a);
		NormalizePoint(&v);
		
		// normalize the bond vector
		VectorProduct(&b,&v,&p);
		NormalizePoint(&b);
		
		// calculate the centrail rail
		CalculateRing(&a,&b,&v,4,lightBond,0.0f,&n,&pBonds[i*12]);       // bond with the central rail
		CalculateRing(&a1,&b,&v,4,lightBond,0.0f,&n,&pBonds[i*12+4]);
		
		// calculate the outer ring to close the bond
		Subtract(&v,&a,&a2);
		NormalizePoint(&v);
		VectorProduct(&b,&v,&p);
		NormalizePoint(&b);
		CalculateRing(&a2,&b,&v,4,lightBond,0.0f,&n,&pBonds[i*12+8]);	// bond with the outer rail
		
		// next bond
		l += distBond;
	}
}


/******************************************************************************
 * initialize the normal vector tables for:
 *									two outer rails, the cylinder, the bonds
 * without initialization these tables cannot be used to store the necessary
 * data about the geometry of the rails, cylinders, and bonds
 *****************************************************************************/
void InitNormals()
{
	int i;
	
	// normals of the first outer rail
	pRail1n = (point*)malloc(nbLine*nbDimensions*sizeof(point));
	for (i=0; i<nbLine*nbDimensions; i++)
	{
		Subtract(&pRail1n[i],&pRail1[i],&pLine1[i/nbDimensions]);
		NormalizePoint(&pRail1n[i]);
	}
	
	// nomarls of the second outer rail
	pRail2n = (point*)malloc(nbLine*nbDimensions*sizeof(point));
	for (i=0; i<nbLine*nbDimensions; i++)
	{
		Subtract(&pRail2n[i],&pRail2[i],&pLine2[i/nbDimensions]);
		NormalizePoint(&pRail2n[i]);
	}
	
	// normals of the cylinders
	pCylindern = (point*)malloc(nbLine*4*sizeof(point));
	for (i=0; i<nbLine*4; i++)
	{
		Subtract(&pCylindern[i],&pCylinder[i],&pLine[i/4]);
		NormalizePoint(&pCylindern[i]);
	}
	
	// normals of the bonds
	pBondsn = (point*)malloc(nbBonds*12*sizeof(point));
	for (i=0; i<nbBonds*3; i++)
	{
		Subtract(&pBondsn[i*4],&pBonds[i*4],&pnBonds[i]);
		NormalizePoint(&pBondsn[i*4]);
		
		Subtract(&pBondsn[i*4+1],&pBonds[i*4+1],&pnBonds[i]);
		NormalizePoint(&pBondsn[i*4+1]);
		
		Subtract(&pBondsn[i*4+2],&pBonds[i*4+2],&pnBonds[i]);
		NormalizePoint(&pBondsn[i*4+2]);
		
		Subtract(&pBondsn[i*4+3],&pBonds[i*4+3],&pnBonds[i]);
		NormalizePoint(&pBondsn[i*4+3]);
	}
}


/******************************************************************************
 * initialize the colors of the components of the coaster
 * the coaster has yellow inner track and blue outer tracks
 * the bonds between the inner and outer tracks are also yellow but lighter
 *****************************************************************************/
void InitColors(point *lightDir)
{
	int i;
	float f;
	static GLfloat LightAmbient[] = {0.3f,0.3f,0.3f,1.0f};
	static GLfloat LightDiffuse[] = {1.1f,0.8f,1.1f,1.0f};

	NormalizePoint(lightDir);
	for (i=0; i<nbLine*nbDimensions; i++)
	{
		// blue right rail
		f = ScalarProduct(lightDir,&pRail1n[i]);
		if (f<=0.0f) f = 0.0f;
		pRail1n[i].x = LightDiffuse[0] * f + -1.0f * LightAmbient[0];
		pRail1n[i].y = LightDiffuse[1] * f + -1.0f * LightAmbient[1];
		pRail1n[i].z = LightDiffuse[2] * f + 2.0f * LightAmbient[2];
	
		// blue left rail
		f = ScalarProduct(lightDir,&pRail2n[i]);
		if (f<=0.0f) f = 0.0f;
		pRail2n[i].x = LightDiffuse[0] * f + -1.0f * LightAmbient[0];
		pRail2n[i].y = LightDiffuse[1] * f + -1.0f * LightAmbient[1];
		pRail2n[i].z = LightDiffuse[2] * f + 2.0f * LightAmbient[2];
	}

	// dark yellow inner tube
	for (i=0; i<nbLine*4; i++)
	{
		f = ScalarProduct(lightDir,&pCylindern[i]);
		if (f<=0.0f) f = 0.0f;
		pCylindern[i].x = LightDiffuse[0] * f + 2.0f * LightAmbient[0];
		pCylindern[i].y = LightDiffuse[1] * f + 2.0f * LightAmbient[1];
		pCylindern[i].z = LightDiffuse[2] * f + 0.0f * LightAmbient[2];
	}
	
	// yellow to white bonds
	for (i=0; i<nbBonds*12; i++)
	{
		f = ScalarProduct(lightDir,&pBondsn[i]);
		if (f<=0.0f) f = 0.0f;
		pBondsn[i].x = LightDiffuse[0] * f + 2.0f * LightAmbient[0];
		pBondsn[i].y = LightDiffuse[1] * f + 2.0f * LightAmbient[1];
		pBondsn[i].z = LightDiffuse[2] * f + 0.5f * LightAmbient[2];
	}
}


/******************************************************************************
 * draw the major rail
 * this rail will be closed onto itself
 *****************************************************************************/
void DrawMajorRail(point *pPipe, point *pNorm, int nbLine, int nbDimensions, float *mat, point *pLine)
{
	int i0,i,j,k;
	point v;
	int flag1,flag2;
	int index;
	VectorMultiply(&v,mat,&pLine[0]);
	flag1 = (v.x < -1.1f || v.x > 1.1f || v.y < -1.2f || v.y > 1.1f || v.z > 0.0f);
	i0 = 0;
	for (i=0; i<nbLine; i++)
	{
		VectorMultiply(&v,mat,&pLine[(i+1)%nbLine]);
		flag2 = (v.x < -1.1f || v.x > 1.1f || v.y < -1.2f || v.y > 1.1f || v.z > 0.0f);
		if (flag1 && flag2)
		{
			i0 = (i+1)%nbLine; 
			continue;
		}
		// draw the rail with strips
		glBegin(GL_QUAD_STRIP);
		for (j=0; j<=nbDimensions; j++)
		{
			k = j % nbDimensions;
			
			index = ((i+1)%nbLine)*nbDimensions+k;
			glColor3fv((float*)&pNorm[index]);
			glVertex3fv((float*)&pPipe[index]);
			
			index = i0*nbDimensions+k;
			glColor3fv((float*)&pNorm[index]);
			glVertex3fv((float*)&pPipe[index]);
		}
		glEnd();
		i0 = (i+1) % nbLine;
		flag1 = flag2;
		if ((v.z < -2.0f) && i<(nbLine-2) && (i%2)==0) i++;
	}
}


/******************************************************************************
 * draw the outer rail
 * like the major rail this rail will also be closed onto itself
 * this will draw a rail that is offset and translated wrt to the major rail
 *****************************************************************************/
void DrawOuterRail(point *pPipe, point *pNorm, int nbLine, int nbDimensions, float *mat, point *pLine)
{
	int i0,i,j,k;
	point v;
	int flag1,flag2;
	int shl,shl0;
	int index;
	
	VectorMultiply(&v,mat,&pLine[0]);
	flag1 = (v.x<-1.1f || v.x>1.1f || v.y<-1.2f || v.y>1.1f || v.z>0.0f);
	
	i0 = 0;
	shl0 = 0;
	for (i=0; i<nbLine; i++)
	{
		VectorMultiply(&v,mat,&pLine[(i+1)%nbLine]);
		flag2 = (v.x<-1.1f || v.x>1.1f || v.y<-1.2f || v.y>1.1f || v.z>0.0f);
		if ((v.z >= -2.0f)||(i==nbLine-1)) shl = 0; else shl = 1;
		if (flag1 && flag2)
		{
			i0 = (i+1)%nbLine;
			shl0 = shl;
			flag1 = flag2;
			continue;
		}
		
		// draws the rail around the normal
		if (!shl && shl0)
		{
			for (j=0; j<nbDimensions; j+=2)
			{
				glBegin(GL_TRIANGLE_STRIP);
				k = j;
				index = ((i+1)%nbLine)*nbDimensions+k;
				glColor3fv((float*)&pNorm[index]);
				glVertex3fv((float*)&pPipe[index]);
				index = i0*nbDimensions+k;
				glColor3fv((float*)&pNorm[index]);
				glVertex3fv((float*)&pPipe[index]);
				k = j+1;
				index = ((i+1)%nbLine)*nbDimensions+k;
				glColor3fv((float*)&pNorm[index]);
				glVertex3fv((float*)&pPipe[index]);
				k = (j+2) % nbDimensions;
				index = i0*nbDimensions+k;
				glColor3fv((float*)&pNorm[index]);
				glVertex3fv((float*)&pPipe[index]);
				index = ((i+1)%nbLine)*nbDimensions+k;
				glColor3fv((float*)&pNorm[index]);
				glVertex3fv((float*)&pPipe[index]);
				glEnd();
			}
		}
		// closes the rail
		else
		{
			glBegin(GL_QUAD_STRIP);
			for (j=0; j<=nbDimensions; j++)
			{
				k = j % nbDimensions;
				index = ((i+1)%nbLine)*nbDimensions+k;
				glColor3fv((float*)&pNorm[index]);
				glVertex3fv((float*)&pPipe[index]);
				index = i0*nbDimensions+k;
				glColor3fv((float*)&pNorm[index]);
				glVertex3fv((float*)&pPipe[index]);
			}
			glEnd();
		}
		
		i0 = (i+1)%nbLine;
		shl0 = shl;
		flag1 = flag2;
		if ((v.z < -2.0f) && i<(nbLine-2) && (i%2)==0) i++;
	}
}


/******************************************************************************
 * draw the bonds on the major rail to the outer rails
 * first the left bonds are drawn
 * tehn the right bonds are drawn
 *****************************************************************************/
void DrawBondsOnRail(point *pBonds, point *pBondsn, int nbBonds, float *mat, point *pnBonds)
{
	int i,j,k;
	int flag1,flag2,flag3,distance;
	point v;
	for (i=0; i<nbBonds; i++)
	{
		// set flags for distances between inner and outer rails
		VectorMultiply(&v,mat,&pnBonds[i*3]);
		flag1 = (v.x<-1.0f || v.x>1.0f || v.y<-1.0f || v.y>1.0f || v.z>0.0f);
		distance = (v.z < -2.0f);
		
		VectorMultiply(&v,mat,&pnBonds[i*3+1]);
		flag2 = (v.x<-1.0f || v.x>1.0f || v.y<-1.0f || v.y>1.0f || v.z>0.0f);
		
		VectorMultiply(&v,mat,&pnBonds[i*3+2]);
		flag3 = (v.x<-1.0f || v.x>1.0f || v.y<-1.0f || v.y>1.0f || v.z>0.0f);
		
		// draw the right bonds to the outer rails
		// bonds are smoothed out in the middle
		if (!(flag1&&flag2))
		{
			glBegin(GL_QUAD_STRIP);
			for (j=0; j<5; j++)
			{
				if (distance && j==3) j++;
				k = j % 4;
				glColor3fv((float*)&pBondsn[i*12+4+k]);
				glVertex3fv((float*)&pBonds[i*12+4+k]);
				glColor3fv((float*)&pBondsn[i*12+k]);
				glVertex3fv((float*)&pBonds[i*12+k]);
			}
			glEnd();
		}
		// draw the left bonds to the outer rails
		// bonds are smoothed out in the middle
		if (!(flag1&&flag3))
		{
			glBegin(GL_QUAD_STRIP);
			for (j=0; j<5; j++)
			{
				if (distance && j==3) j++;
				k = j % 4;
				glColor3fv((float*)&pBondsn[i*12+k]);
				glVertex3fv((float*)&pBonds[i*12+k]);
				glColor3fv((float*)&pBondsn[i*12+8+k]);
				glVertex3fv((float*)&pBonds[i*12+8+k]);
			}
			glEnd();
		}
	}
}


/******************************************************************************
 * draw sky and ground
 * first the sky is drawn
 * the sky is blended to show the changing of color in the sky
 * then the ground is drawn using a triangle fan instead of a single polygon
 * this step increased the frames as managing the entire ground is less intense
 *****************************************************************************/
void DrawSkyAndGround(int nbDimensions, float size, float texsize, float low, float high)
{ 
	int i;
	float x,y;
	float c,s;

	// draw the sky with strips
	glBegin(GL_QUAD_STRIP);
	glColor3f(0.3f,0.3f,0.7f);
	glVertex3f(size,0.0f,low);
	glColor3f(0.5f,0.5f,0.7f);
	glVertex3f(size,0.0f,high);
	// we are looking at a dome view
	// map the coordinates according to the angle
	for (i=1; i<nbDimensions; i++)
	{
		x = size * cos((float)i*2.0f*M_PI/(float)nbDimensions);
		y = size * sin((float)i*2.0f*M_PI/(float)nbDimensions);
		glColor3f(0.5f,0.5f,0.7f);
		glVertex3f(x,y,low);
		glColor3f(0.3f,0.3f,0.7f);
		glVertex3f(x,y,high);
	}
	glColor3f(0.3f,0.3f,0.7f);
	glVertex3f(size,0.0f,low);
	glColor3f(0.5f,0.5f,0.7f);
	glVertex3f(size,0.0f,high);
	glEnd();

	// draw the ground with many triangles instead of a single polygon
	glEnable(GL_TEXTURE_2D);
	glBmpBind(&grass);
	glColor3f(1.0f,1.0f,1.0f);
	glBegin(GL_TRIANGLE_FAN);
	glTexCoord2f(0.0f,0.0f);
	glVertex3f(0.0f,0.0f,low);
	// we are looking at a dome view
	// map the coordinates according to the angle
	for (i=0; i<nbDimensions; i++)
	{
		c = cos((float)i*2.0f*M_PI/(float)nbDimensions);
		s = sin((float)i*2.0f*M_PI/(float)nbDimensions);
		glTexCoord2f(texsize*c,texsize*s);
		glVertex3f(size*c,size*s,low);
	}
	glTexCoord2f(texsize,0.0f);
	glVertex3f(size,0.0f,low);
	glEnd();
	glDisable(GL_TEXTURE_2D);
}


/******************************************************************************
 * draw entire scene
 * this is the scene renderer
 * this is the flow of the application also
 *****************************************************************************/
void DrawEntireScene()
{
	// create a projection and model holders so that the scene can be drawn
	float m[20];
	float modelmat[20], projmat[20];

	glGetFloatv(GL_MODELVIEW_MATRIX,modelmat);
	glGetFloatv(GL_PROJECTION_MATRIX,projmat);
	MatrixMultiply(m,projmat,modelmat);
	m[2] = modelmat[2];
	m[6] = modelmat[6];
	m[10] = modelmat[10];
	m[14] = modelmat[14];
	m[18] = modelmat[18];

	// process of rendering the entire scene from start to finish
	DrawSkyAndGround(6,50.0f,60.0f,0.0f,10.0f);
	DrawTrees();
	DrawPlatform();
	DrawColumns();
	DrawMajorRail(pCylinder,pCylindern,nbLine,4,m,pLine);
	DrawBondsOnRail(pBonds,pBondsn,nbBonds,m,pnBonds);
	DrawOuterRail(pRail1,pRail1n,nbLine,nbDimensions,m,pLine1);
	DrawOuterRail(pRail2,pRail2n,nbLine,nbDimensions,m,pLine2);
}


/******************************************************************************
 * fade rendered scene based on opacity
 * the fade kicks in by using a timer
 *****************************************************************************/
void FadeScene(float opacity)
{
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColor4f(0.01f,0.01f,0.01f,opacity);
	glRectf(-1.0f, -1.0f, 1.0f, 1.0f);
	glDisable(GL_BLEND);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}


/******************************************************************************
 * draw the intial panoramic span of the rendered coaster
 * this involves rotating the camera about a circle in the sky
 * this is done using angles and a maximum distance away from the coaster
 * the angles to the coaster are calculated by mapping the distance to the 
 * angle and then rotating about that angle
 *****************************************************************************/
GLvoid DrawPanorama()
{
	// fps stuff
	static unsigned long begin_time;
	static unsigned long time, diftime, oldtimeframe;
	static int nbframe = 0;
	float fps;
	static int flag_firstcall = 1;
	static float angle = 0.0f;
	char wndtitle[BUFSIZ/2];
	
	float dist;

	// very first call to load time
	if (flag_firstcall)
	{
		oldtimeframe = timeGetTime();
		flag_firstcall = 0;
		begin_time = oldtimeframe;
	}
	time = timeGetTime();

	// calculate the necessary angles
	angle = (float)(timeGetTime() - begin_time) / 100; // angle in degree
	angle = (angle + 180.0f) / 180.0f * M_PI;          // angle in radian
	
	// compute maximum camera distance from the origin
	if (maxDist*tan(aperture*M_PI/360.0f)*2.0f < maxZ) 
		dist = maxZ / tan(aperture*M_PI/360.0f) / 2.0f;
	else 
		dist = 0.5f + maxDist;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	// the camera is placed along a circular track on the sky and
	// the panorama is completed
	// we always point up
	gluLookAt(dist*cos(angle), -dist*sin(angle),maxZ,				    // eye position
	        	0.0f, 0.0f, maxZ+1.0f-dist * tan(aperture*M_PI/360.0f),	// target
	        	0.0f, 0.0f, 1.0f);									    // up vector
	
	DrawEntireScene();
	
	// fade scene in at the beginning
	if (time - begin_time < 2000)
	{
		diftime = time - begin_time;
		FadeScene(1.0f-(float)diftime / 2000.0f);
	}
	// fade scene out at the end
	if (time - begin_time > 34000)
	{
		diftime = time - begin_time - 34000;
		FadeScene((float)diftime / 2000.0f);
	}
	// display fps
	if (time-oldtimeframe > 200.0f)
	{
		fps=(float)nbframe*1000.0f/(float)(time-oldtimeframe);
		oldtimeframe=time;
		nbframe=0;
		sprintf(wndtitle,"WIRED @ %1.2f FPS",fps);
		glutSetWindowTitle(wndtitle);
	}
	nbframe++;
}


/******************************************************************************
 * draw the ride from the 1st person viewpoint
 * the camera is placed where the first person riding it would see
 * there is no car unfortunately because of time constraints
 * i focused more on making the ride as smooth as possible
 * i also thought that the new coasters are floorless/carless where the rider
 * is seated and only sees the ride
 *
 * here the scene isn't moved through the eye position but the other way
 * one, i couldn't figure out how to move the entire scene
 * also, that is too intensive -- it will give low frame rates
 * second, moving the camera through the rendered scene requires mapping the
 * coordinate as the ride progresses
 *****************************************************************************/
GLvoid DrawRide()
{
	int i=0;
	static int index = 0;
	
	// mapping of the position of the camera as the ride goes on
	static float xcoord = 0.0f;
	
	// speed related
	static float oldz;
	static float speed = 0.0f;
	
	// set to zero when coaster doesn't rise and begins falling
	static int flag = 1;
	static int flag_firstcall = 1;
	static int flag_stopping = 0;

	// frames per seconds calculations
	float t, fps;
	char wndtitle[64];
	static unsigned long time, oldtime;
	static unsigned long begin_time;
	static unsigned long oldtimeframe;
	static unsigned long stoptime;
	static int nbframe = 0;

	// position coordinates
	// eye position
	// eye direction position (projection vectors)
	point p, position, eye;
	point ddv,ddv1,ddv2;     // front, up, side position vectors
	
	// time info on first call for frames/sec calc
	if (flag_firstcall)
	{
		oldtimeframe = timeGetTime();
		index = startSegment;
		time = begin_time = oldtimeframe;
	}

	// decrement the xcoord as moving along the coaster
	while(xcoord > pNormal[(index+1)%nbLine])
	{
		xcoord -= pNormal[(index+1)%nbLine];
		index++;
		index %= nbLine;
	}
	t = xcoord / pNormal[(index+1)%nbLine];
	position = pLine[index];
	Subtract(&p,&pLine[(index+1)%nbLine],&pLine[index]);
	AdditiveMultiply(&position,t,&p);
	
	// UP vector
	ddv1 = pTrajectory[index];
	Subtract(&p,&pTrajectory[(index+1)%nbLine],&pTrajectory[index]);
	AdditiveMultiply(&ddv1,t,&p);
	NormalizePoint(&ddv1);                // ddv1 vector points upward
	
	// FRONT vector
	ddv = pTangent[index];
	Subtract(&p,&pTangent[(index+1)%nbLine],&pTangent[index]);
	AdditiveMultiply(&ddv,t,&p);
	NormalizePoint(&ddv);                 // ddv vector points forward
	p = ddv1;
	
	// SIDE - RIGHT bector
	RotatePoint(&ddv1, &ddv, pCurve[index]+t*(pCurve[(index+1)%nbLine]-pCurve[index]), &p);
	VectorProduct(&ddv2,&ddv,&ddv1);     // ddv2 vector points rightward
	
	// compute camera position 
	Multiply(&eye,0.22f,&ddv1);	// elevate slightly off the track
	Add(&eye,&eye,&position);

	// rotate camera when on the coaster
	if (roll_angle!=0.0f) 
	{ 
		// rotate the front and side vectors around the up vector
	    p = ddv;		// Front vector
	    RotatePoint(&ddv, &ddv1, roll_angle, &p);
	    p = ddv2;		// Side vector
	    RotatePoint(&ddv2, &ddv1, roll_angle, &p);
	} 
	else if (pitch_angle!=0.0f) 
	{ 
		// rotate the up and front vectors around the side vector.
	    p = ddv;		// Front vector
	    RotatePoint(&ddv, &ddv2, pitch_angle, &p);
	    p = ddv1;		// Up vector
	    RotatePoint(&ddv1, &ddv2, pitch_angle, &p);
	}

	// position the eye at the front of the coaster
	glLoadIdentity();
	gluLookAt(eye.x, eye.y, eye.z,						// eye position
		    	eye.x+ddv.x, eye.y+ddv.y, eye.z+ddv.z,	// target
		    	ddv1.x, ddv1.y, ddv1.z);				// up vector
	// render the scene
	DrawEntireScene();
	
	// fade scene in at the beginning of the ride
	if (time - begin_time < 1500) FadeScene(1.0f-(float)(time - begin_time) / 1500.0f);

	time = timeGetTime();
	if (flag)
	{
		// going up the first hill
		if (oldz-position.z > 0.0f) flag = 0;
		else speed = 1.0f;
	}
	else
	{
		if (index >= brakeSegment || flag_stopping)
		{
			flag_stopping = 1;
			if (speed > 0.0f)
			{
				// time to brake do a hard stop
				// this simulates the physics of when a coaster comes to a stop
				speed = speed*speed - 50.0f * (float)(time - oldtime) / 1000.0f;
				if (speed <= 0.0f)
				{
					// stop it completely
					speed = 0.0f;
					stoptime = time;
				}
				// decrease the speed exponentially -- smooth
				else speed = sqrt(speed);
			}
			else
			{
				// it is almost time to stop
				if (time - stoptime > 5000)
				{
					flag_stopping = 0;
					flag = 1;
				}
			}
		}
		else
		{
			// we are at the top of a hill
			// the speed should pick up on the way down
			// the coaster halts at the very top and picks up speed
			speed = speed*speed + 4.2f*(oldz-position.z);
			if (speed <= 0.0f) speed = 0.0f;
			else speed = sqrt(speed);
		}
	}
	glFlush();
	if (flag_firstcall)
	{
		flag_firstcall = 0;
		oldtime = time-1;
	}
	// move the xcoord along the track wrt to speed
	xcoord += speed * (float)(time - oldtime) / 1000.0f;
	
	// calculate the fps and display
	// keep track of previous framerate calc and update
	if (time-oldtimeframe > 200.0f)
	{
		fps=(float)nbframe*1000.0f/(float)(time-oldtimeframe);
		oldtimeframe=time;
		nbframe=0;
		sprintf(wndtitle,"WIRED @ %1.2f FPS",fps);
		glutSetWindowTitle(wndtitle);
	}
	nbframe++;
	oldz = position.z;
	oldtime = time;
}


/******************************************************************************
 * initialize the roller coaster before it is rendered in the scene
 * initialize the scene itself
 * load the necessary textures
 * init curve
 * init segments
 * init bonds
 *****************************************************************************/
void InitializeRollerCoaster(int intro, int width, int height, float viewAngle, float eyeDist, float focalLength)
{
	static point lightDir={3.0f,6.0f,9.0f};
	doIntro = intro;
	aperture = viewAngle;
	eye_sep = eyeDist;
	focallength = focalLength;
	wndHeight=height;
	wndWidth=width;
	InitScene();
	ResizeScene(wndWidth,wndHeight);
	// parse the coaster definition file
	// this coaster will make you WIRED
	if (OpenInputFile("Data\\WIRED.coaster")) ParseFile(); CloseFile();
	InitCurve(avgSegmentLength);
	// saftey check
	if (startSegment < 0) startSegment = nbLine + startSegment;
	if (brakeSegment < 0) brakeSegment = nbLine + brakeSegment;
	LoadTextures();
	InitLines();
	InitBonds(0.45f);		// every three segments a bond is drawn
	InitNormals();
	InitColors(&lightDir);
	InitColumns(pLine,nbLine,pTangent);
}


/******************************************************************************
 * change roller coaster parameters such as viewing angle, 
 * panning, and focal length
 *****************************************************************************/
void ChangeRollerCoasterParameters(float viewAngle, float eyeDist, float focalLength, int left, int right)
{
	if (aperture != viewAngle) 
	{
		aperture = viewAngle;
		ResizeScene(wndWidth,wndHeight);
	}
	if (left) roll_angle += .025f;
	if (right) roll_angle -= .025f;
	aperture = viewAngle;
	eye_sep = eyeDist;
	focallength = focalLength;
}


/******************************************************************************
 * this section outlines functions for the drawing of columns and bonds
 * there is light information also
 *****************************************************************************/
point* pColumn;
point* pBond;
point* pColumnColor;
point* pBondColor;
int nbColumn;
int nbBond;

int nbColumnDimension = 6;			// smoother cylinders
int nbBondDimension = 4;			// smooth bonds notching the columns	

float rayColumn = 0.02f;
static float rayBond = 0.015f;

float distanceColumn = 0.4f;
float heightBond = 0.1f;

// light hitting the columns and bonds
point LightDirection = {3.0f,-1.5f,5.0f};
point LightAmbient = {0.4f,0.4f,0.4f};
point LightDiffuse = {1.0f,0.9f,1.0f};
point LightDefinition = {1.0f,0.5f,1.0f};


/******************************************************************************
 * calculate the light hitting an object (cylinder in most cases)
 * light is normalized so that it is constant at each length away from object
 *****************************************************************************/
void CreateCylinder(point* a, point* b, point* v, int nbCol, float ray, point *buf, point* bufn)
{
	int i,k;
	point p1,p2,p4;
	
	Subtract(&p1,a,b);
	NormalizePoint(&p1);
	
	p4 = *v;
	AdditiveMultiply(&p4,ScalarProduct(v,&p1)/Length(&p1),&p1);
	NormalizePoint(&p4);
	Multiply(&p4,ray,&p4);
	
	// create the column about the rotated point
	k = 0;
	for (i=0; i<nbCol; i++)
	{
		RotatePoint(&p2,&p1,(float)i*2.0f*M_PI/(float)nbCol,&p4);
		Add(&buf[k++],a,&p2);
		Add(&buf[k++],b,&p2);
	}
	
	// make columns shiny otherwise they look flat
	// don't reflect light
	k = 0;
	for (i=0; i<nbCol; i++)
	{
		Subtract(&p1,&buf[k],a);
		NormalizePoint(&p1);
		bufn[k++] = p1;
		Subtract(&p1,&buf[k],b);
		NormalizePoint(&p1);
		bufn[k++] = p1;
	}
}


/******************************************************************************
 * calculate the light hitting an object (cylinder in most cases)
 * light is normalized so that it is constant at each length away from object
 *****************************************************************************/
void CalculateLight(int nb, point* r, point* amb, point* dif, point* lig, point* normal)
{
	int i;
	float f;
	static point light = {0.4f,0.4f,0.4f};
	// the resulting light is normalized
	// by adjusting to the incident light
	for (i=0; i<nb; i++)
	{
		f = ScalarProduct(&LightDirection,&normal[i]);
		if (f<=0.0f) f = 0.0f;
		r[i].x = dif->x * f + light.x;
		r[i].y = dif->y * f + light.y;
		r[i].z = dif->z * f + light.z;
	}
}


/******************************************************************************
 * create a column specified by points p1 and p2
 * p is a point of the trajectory
 * v is the tangent with the trajectory in this point
 *****************************************************************************/
void CreateColumn(point* p1, point* p2, point* v, point* buf, point* bufc)
{
	point a,b;
	int i;
	a = *p1;
	a.z = 0.0f;
	b = *p2;
	b.z = 0.0f;
	if (v->z < 0.0f)
	{
		a.z = b.z = -v->z;
		v->z = 0.0f;
	}
	// create cylinder one for point one
	CreateCylinder(p1,&a,v,nbColumnDimension,rayColumn,buf,bufc);
	i = nbColumnDimension * 2;
	bufc[i].x = 0.0f; 
	bufc[i].y = 0.0f; 
	bufc[i].z = 1.0f;
	buf[i] = *p1;
	buf[i++].z += rayColumn;

	// create cylinder two for point two
	CreateCylinder(p2,&b,v,nbColumnDimension,rayColumn,&buf[i],&bufc[i]);
	i += nbColumnDimension * 2;
	bufc[i].x = 0.0f; 
	bufc[i].y = 0.0f; 
	bufc[i].z = 1.0f;
	buf[i] = *p2;
	buf[i++].z += rayColumn;

	// create a fading light effect
	CalculateLight(4*nbColumnDimension+2,bufc,&LightAmbient,&LightDiffuse,&LightDefinition,bufc);
}


/******************************************************************************
 * create a bond between two cylinders (rails or support columns)
 * store the bond information in the connection points p1 and p2
 *****************************************************************************/
void CreateBond(point* p, point* p1, point* p2, point* buf, point* bufc)
{
	point a;
	float z1, z2;
	z1 = p1->z;
	z2 = p2->z;
	p1->z = p->z - heightBond;
	p2->z = p->z - heightBond;
	a.x = 0.0f; a.y = 0.0f; a.z = 1.0f;
	// create the cylinder -- connection bond to point one
	CreateCylinder(p,p1,&a,nbBondDimension,rayBond,buf,bufc);
	// create the cylinder -- connection bond to point two
	CreateCylinder(p,p2,&a,nbBondDimension,rayBond,&buf[nbBondDimension*2],&bufc[nbBondDimension*2]);
	// create a light effect of fading light on the bond -- simulate a thinner bond
	CalculateLight(4*nbBondDimension,bufc,&LightAmbient,&LightDiffuse,&LightDefinition,bufc);
	p1->z = z1;
	p2->z = z2;
}


/******************************************************************************
 * draw the support columns
 * draw cylinders with triangles on top using GL_FAN
 *****************************************************************************/
void DrawColumn(int nb, point *buf, point *bufc)
{
	int i,j;
	nb *= 2;
	for (i=0; i<nb; i++)
	{
		// draw the columns up
		glBegin(GL_QUAD_STRIP);
		for (j=0; j<2*nbColumnDimension; j++)
		{
			glColor3fv((float*)&bufc[j+i*(2*nbColumnDimension+1)]);
			glVertex3fv((float*)&buf[j+i*(2*nbColumnDimension+1)]);
		}
		glColor3fv((float*)&bufc[i*(2*nbColumnDimension+1)]);
		glVertex3fv((float*)&buf[i*(2*nbColumnDimension+1)]);
		glColor3fv((float*)&bufc[i*(2*nbColumnDimension+1)+1]);
		glVertex3fv((float*)&buf[i*(2*nbColumnDimension+1)+1]);
		glEnd();
		
		// draw the triangular tips at the top of the columns
		// the indices use geometry of a triangle
		glBegin(GL_TRIANGLE_FAN);
		glColor3fv((float*)&bufc[(i+1)*(2*nbColumnDimension+1)-1]);
		glVertex3fv((float*)&buf[(i+1)*(2*nbColumnDimension+1)-1]);
		for (j=0; j<2*nbColumnDimension; j+=2)
		{
			glColor3fv((float*)&bufc[j+i*(2*nbColumnDimension+1)]);
			glVertex3fv((float*)&buf[j+i*(2*nbColumnDimension+1)]);
		}
		glColor3fv((float*)&bufc[i*(2*nbColumnDimension+1)]);
		glVertex3fv((float*)&buf[i*(2*nbColumnDimension+1)]);
		glEnd();
	}
}


/******************************************************************************
 * draw the bond at the top of the support columns
 * the bond joins the columns to the major rail
 *****************************************************************************/
void DrawColumnBond(int nb, point *buf, point *bufc)
{
	int i,j;
	nb *= 2;
	for (i=0; i<nb; i++)
	{
		glBegin(GL_QUAD_STRIP);
		// draw on the major rail
		for (j=0; j<2*nbBondDimension; j++)
		{
			glColor3fv((float*)&bufc[j+i*(2*nbBondDimension)]);
			glVertex3fv((float*)&buf[j+i*(2*nbBondDimension)]);
		}
		// draw on the support columns
		glColor3fv((float*)&bufc[i*(2*nbBondDimension)]);
		glVertex3fv((float*)&buf[i*(2*nbBondDimension)]);
		glColor3fv((float*)&bufc[i*(2*nbBondDimension)+1]);
		glVertex3fv((float*)&buf[i*(2*nbBondDimension)+1]);
		glEnd();
	}
}


/******************************************************************************
 * create a connection of the bonds from the support columns to the major rail
 * initialize the geometry
 * place the points at the top of the support columns
 *****************************************************************************/
int CreateConnection(point* r, point* p, point* p1, point* p2, int i, point* pLine, int nbLine)
{
	point v,n,z;
	float t;
	
	Subtract(&v,p1,p2);
	
	z.x = 0.0f; 
	z.y = 0.0f; 
	z.z = 1.0f;
	VectorProduct(&n,&v,&z);
	Subtract(&v,&pLine[(i+1)%nbLine],&pLine[i]);
	
	// check for intersection to the the major rail
	if (!Intersection(r,&pLine[i],&v,p1,&n)) return 0;
	
	// move lower
	Subtract(&z,p1,&pLine[i]);
	
	// check if the connection can be made
	t = ScalarProduct(&z,&n) / ScalarProduct(&v,&n);
	if (t<0.0f || t>=1.0f) return 0;
	
	Subtract(&v,r,p);
	v.z = 0.0f;
	
	// still too far
	if (Magnitude(&v)>0.3f) return 0;
	
	return 1;
}


/******************************************************************************
 * create the supporting columns for the track to rest on
 * without these columnds the track will look like it is suspended int in
 * free space and real coasters don't do that
 *
 * columns are drawn up to all rails so they have to be placed in exactly the
 * correct spot otherwise there will be collisions
 * i made a track which addresses this issue
 * as an extension, i might add the collision detection
 *****************************************************************************/
void CreateSupport(point* p, point* v, point* pLine, int nbLine)
{
	point a,b,r;
	point p1, p2;
	int i;
	int flag_last_i=0;
	float zmax;
	a = *v;
	a.z = 0.0f;
	NormalizePoint(&a);
	if (a.x*a.x>0.0001f)
	{
		a.x = -a.y / a.x;
		a.y = 1.0f;
	}
	else
	{
		a.y = -a.x / a.y;
		a.x = 1.0f;
	}
	NormalizePoint(&a);
	Multiply(&a,distanceColumn/2.0f,&a);
	b = *p;
	Add(&p1,&b,&a);
	Subtract(&p2,&b,&a);
	zmax = -1.0f;
	for (i=0; i<nbLine; i++)
	{
		// connect the columns to the major coaster rail
		// also connect them to each other at the middle
		if (CreateConnection(&r,p,&p1,&p2,i,pLine,nbLine))
		{
			// create a bond between the connection rail
			if (!flag_last_i && ((p->z*p->z)<0.0001f || p->z > r.z) && r.z-heightBond>-v->z)
			{
				// create a bond with the correct number of dimensions and color
				CreateBond(&r,&p1,&p2,&pBond[nbBond*4*nbBondDimension],&pBondColor[nbBond*4*nbBondDimension]);
				nbBond++;
				flag_last_i = 1;
				// check for the connection to the major rail
				if (r.z > zmax) zmax = r.z;
			}
			else flag_last_i = 0;
		}
		else flag_last_i = 0;
	}
	// create the support up to the last point on an incident rail
	if (zmax > -0.5f)
	{
		p1.z = p2.z = zmax;
		if (v->z > 0.001f);
		CreateColumn(&p1,&p2,v,&pColumn[nbColumn*(4*nbColumnDimension+2)],&pColumnColor[nbColumn*(4*nbColumnDimension+2)]);
		nbColumn++;
	}
}


point *pColumnCoord;		// coordinates where pillars have to be placed
int *pColumnCoordinate;		// table of x-coordinate where pillars have to be placed
int nbColumnCoord = 0;

/******************************************************************************
 * initialize the columns
 *****************************************************************************/
void InitColumns(point *pLine, int nbLine, point* pTangent)
{
	point p,v;
	point *table;
	int i, coulmn_table;
	int *table1;
	table = pColumnCoord;
	coulmn_table = nbColumnCoord;
	table1 = pColumnCoordinate;
	nbColumn = coulmn_table;
	nbBond = nbColumn*10;  // the maximum is 10 bonds for 1 column -- assumption
	
	// normalize the light in one direction
	NormalizePoint(&LightDirection);
	
	// initialize holders
	pColumn = (point*)malloc(nbColumn*(4*nbColumnDimension+2)*sizeof(point));
	pColumnColor = (point*)malloc(nbColumn*(4*nbColumnDimension+2)*sizeof(point));
	pBond = (point*)malloc(nbBond*(4*nbBondDimension)*sizeof(point));
	pBondColor = (point*)malloc(nbBond*(4*nbBondDimension)*sizeof(point));
	nbColumn = nbBond = 0;

	for (i=0; i<coulmn_table; i++)
	{
		p = pLine[table1[i]];
		v = pTangent[table1[i]];
		p.z = v.z = 0.0f;
		// create a support at each coordinate
		CreateSupport(&p,&v,pLine,nbLine);
	}
}


/******************************************************************************
 * draw all columns that will create the supports for the coaster rails
 *****************************************************************************/
void DrawColumns()
{
	DrawColumn(nbColumn,pColumn,pColumnColor);
	DrawColumnBond(nbBond,pBond,pBondColor);
}


// platform properties
float metalLength = 2.0f;
float metalAngle = 0.0f;
point metalPosition = {-1.0f,-1.0f,0.5f};

/******************************************************************************
 * draw the loading platform by loading the point coordinates
 * the platform is drawn in multiple steps
 * the sides are drawn
 * the top is drawn
 * the body is filled with the rest of the texture
 *****************************************************************************/
void DrawPlatform()
{
	// platform top 2D coordinates
	static coord2d tabt[12];
	// platform center coordinates
	static point tabc[] = {
		{-1.0f, 0.0f, 0.0f},
		{0.0f, 0.0f, 1.0f},
		{1.0f, 0.0f, 0.0f},
		{0.0f, 0.0f, 1.0f},
		{1.0f, 0.0f, 0.0f}
	};
	// platform vertical coordinates
	static point tab1v[] = {
		{-0.25f, 0.0f, -0.02f},
		{-0.25f, 0.0f, 0.04f},
		{-0.5f, 0.0f, 0.04f},
		{-0.5f, 0.0f, 0.0f},
		{0.25f, 0.0f, 0.0f},
		{0.25f, 0.0f, -0.02f}
	};
	
	static coord2d tab1t[6];
	static point tab1c = {0.0f, -1.0f, 0.0f};
	static point tab2v[6];
	static coord2d tab2t[6];
	static point tab2c = {0.0f, 1.0f, 0.0f};
	int i;
	static int first = 1;
	point p;
	point z = {0.0f, 0.0f, 1.0f};
	
	if (metalLength == 0.0f) return;
	if (first)
	{
		point LightAmbient = {0.4f,0.4f,0.4f};
		point LightDiffuse = {1.0f,1.0f,1.0f};
		point LightDefinition = {1.0f,1.0f,1.0f};
		float angle = metalAngle * M_PI / 180;

		// box is a 6 sided cube
		for (i=0; i<6; i++)
		{
			tab1t[i].x = tab1v[i].x * 20.0f;
			tab1t[i].y = (tab1v[i].z + (i!=3 && i!=4 ? metalPosition.z : 0.0f)) * 20.0f;
			// this creates the rising platform
			tabt[i*2].x = tabt[i*2+1].x = (tab1v[(9-i)%6].x + tab1v[(9-i)%6].z + (i!=0 && i!=5 ? metalPosition.z : 0.0f)) * 20.0f;
			tabt[i*2].y = 0.0f;
			tabt[i*2+1].y = metalLength * 20.0f;
		}
		// more initialization
		tab2t[0] = tab1t[0];
		for (i=1; i<6; i++) tab2t[i] = tab1t[6-i];
		tab2v[0] = tab1v[0];
		for (i=1; i<6; i++) tab2v[i] = tab1v[6-i];
		
		// define the coordinates of the platform
		// using the property of the platform
		for (i=0; i<6; i++)
		{
			tab1v[i].y -= metalLength / 2.0f;
			RotatePoint(&p,&z,angle,&tab1v[i]);
			Add(&tab1v[i],&metalPosition,&p);
			tab2v[i].y += metalLength / 2.0f;
			RotatePoint(&p,&z,angle,&tab2v[i]);
			Add(&tab2v[i],&metalPosition,&p);
		}
		tab1v[3].z = tab1v[4].z = 0.0f;
		tab2v[2].z = tab2v[3].z = 0.0f;
		
		// come around the corner of first side
		RotatePoint(&p,&z,angle,&tab1c);
		tab1c = p;
		
		// come around the corner of second side
		RotatePoint(&p,&z,angle,&tab2c);
		tab2c = p;
		
		// come around the center of the platform
		for (i=0; i<5; i++)
		{
			RotatePoint(&p,&z,angle,&tabc[i]);
			tabc[i] = p;
		}
		// calculate light hitting the platform
		CalculateLight(5,tabc,&LightAmbient,&LightDiffuse,&LightDefinition,tabc);
		CalculateLight(1,&tab1c,&LightAmbient,&LightDiffuse,&LightDefinition,&tab1c);
		CalculateLight(1,&tab2c,&LightAmbient,&LightDiffuse,&LightDefinition,&tab2c);
		first = 0;
	}

	// draw the top and the long sides of the platform
	glEnable(GL_TEXTURE_2D);
	glShadeModel(GL_FLAT);
	glBmpBind(&metal);
	glBegin(GL_QUAD_STRIP);
	for (i=0; i<6; i++)
	{
		glColor3fv((float*)&tabc[(i+4)%5]);
		glTexCoord2fv((float*)&tabt[2*i]);
		glVertex3fv((float*)&tab2v[(3+i)%6]);
		glTexCoord2fv((float*)&tabt[2*i+1]);
		glVertex3fv((float*)&tab1v[(9-i)%6]);
	}
	glEnd();
	
	// draw the upper end vertical of the platform
	glBegin(GL_TRIANGLE_FAN);
	glColor3fv((float*)&tab1c);
	for (i=0; i < sizeof(tab1v)/sizeof(point); i++)
	{
		glTexCoord2fv((float*)&tab1t[i]);
		glVertex3fv((float*)&tab1v[i]);
	}
	glEnd();
	
	// draw the lower end vertical of the platform
	glBegin(GL_TRIANGLE_FAN);
	glColor3fv((float*)&tab2c);
	for (i=0; i < sizeof(tab2v)/sizeof(point); i++)
	{
		glTexCoord2fv((float*)&tab2t[i]);
		glVertex3fv((float*)&tab2v[i]);
	}
	glEnd();
	
	glShadeModel(GL_SMOOTH);
	glDisable(GL_TEXTURE_2D);
}


/******************************************************************************
 * draw a tree by intersecting the loaded with itself at the midpoint at 90deg
 * the coordinates are mapped so that the tree is seen from all angles during
 * the initial panoramic scene
 *****************************************************************************/
point *pTree;		// table of coordinates of trees
int nbTree;			// number of trees

void DrawTree()
{
	// load the tree texture
	glEnable(GL_TEXTURE_2D);
	glBmpBind(&tree);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_NOTEQUAL,0.0f);
	glColor3f(1.0f,1.0f,1.0f);
	
	// map the tree texture
	glBegin(GL_QUADS);
		// draw once
		glTexCoord2f(1.0f, 1.0f); glVertex3f(-0.5f, 0.0f, 1.3f);
		glTexCoord2f(0.0f, 1.0f); glVertex3f(0.5f, 0.0f, 1.3f);
		glTexCoord2f(0.0f, 0.0f); glVertex3f(0.5f, 0.0f, 0.0f);
		glTexCoord2f(1.0f, 0.0f); glVertex3f(-0.5f, 0.0f, 0.0f);

		glTexCoord2f(1.0f, 0.0f); glVertex3f(-0.5f, 0.0f, 0.0f);
		glTexCoord2f(0.0f, 0.0f); glVertex3f(0.5f, 0.0f, 0.0f);
		glTexCoord2f(0.0f, 1.0f); glVertex3f(0.5f, 0.0f, 1.3f);
		glTexCoord2f(1.0f, 1.0f); glVertex3f(-0.5f, 0.0f, 1.3f);

		// intersect the tree with itself at the midpoint
		glTexCoord2f(1.0f, 1.0f); glVertex3f(0.0f, -0.5f, 1.3f);
		glTexCoord2f(0.0f, 1.0f); glVertex3f(0.0f, 0.5f, 1.3f);
		glTexCoord2f(0.0f, 0.0f); glVertex3f(0.0f, 0.5f, 0.0f);
		glTexCoord2f(1.0f, 0.0f); glVertex3f(0.0f, -0.5f, 0.0f);

		glTexCoord2f(1.0f, 0.0f); glVertex3f(0.0f, -0.5f, 0.0f);
		glTexCoord2f(0.0f, 0.0f); glVertex3f(0.0f, 0.5f, 0.0f);
		glTexCoord2f(0.0f, 1.0f); glVertex3f(0.0f, 0.5f, 1.3f);
		glTexCoord2f(1.0f, 1.0f); glVertex3f(0.0f, -0.5f, 1.3f);
	glEnd();
	
	glDisable(GL_ALPHA_TEST);
	glDisable(GL_TEXTURE_2D);
}


/******************************************************************************
 * draw trees specified by the coordinates in the *.coaster file
 * each tree is draw relative to the previous tree using a translation
 * the drawtree utility is used to place a tree at the corresponding location
 *****************************************************************************/
void DrawTrees()
{
	int i;
	point v;
	v.x = v.y = v.z = 0.0f;
	
	glPushMatrix();
	for (i=0; i<nbTree; i++)
	{
		glTranslatef(pTree[i].x - v.x, pTree[i].y - v.y, pTree[i].z - v.z);
		v = pTree[i];
		DrawTree();
	}
	glPopMatrix();
}


/******************************************************************************
 * draws the roller coaster
 * first the panorama then the ride itself
 *****************************************************************************/
void DrawRollerCoaster()
{
	static int flag_firstcall = 1;
	static unsigned long firsttime;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if (flag_firstcall)
	{
		flag_firstcall = 0;
		firsttime = timeGetTime();
	}

	// initially draw the panorama
	// roam around the roller coaster from the sky view
	if (timeGetTime()-firsttime < 36000 && doIntro)
	{
		glViewport(0,0,wndWidth,wndHeight);
		DrawPanorama();						// spans the entire ride
		glutSwapBuffers();
	}
	// once the panorama is complete draw the coaster ride
	else
	{
		glViewport(0,0,wndWidth,wndHeight);
		DrawRide();							// places camera on ride
		glutSwapBuffers();
	}
}