/******************************************************************************
 * Sukhchander Khanna [ khannas@rpi.edu ]
 * Computer Graphics Fall 04
 * Final Project Component One -- 'WIRED' --> My Roller Coaster Creation
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
 * File: RollerCoasterWIRED.c
 * Desc: everything starts here. glut is intialized. a window is created.
 *		 the intial panorama is shown then the ride begins.
 *		 the description of the coaster creation is in RollerCoaster.c
 *****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <GL/glut.h>
#include <windows.h>

#include "RollerCoaster.h"


int width = 1200, height = 800;

// scene information
float viewAngle = 60.0f;
static float eyeDist = 0.005f;
static float focalLength = 0.05f;
int intro = 1;


/******************************************************************************
 * displays windows dialog box with instructions for application
 * box will be shown with instructions
 * the scene will pause and resume
 *****************************************************************************/
void usage()
{
	char msg[BUFSIZ/2] = "";
	strcat(msg, "Coaster Cam \nOptions:\n'h'\tprint this help \n");
	strcat(msg, "'left'\tpan left when riding coaster \n");
	strcat(msg, "'right'\tpan right when riding coaster \n\n");
	MessageBox(NULL,msg,"Coaster Cam",MB_OK);
}


/******************************************************************************
 * standard keyboard operations
 *****************************************************************************/
void keyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
		case 'h':
			usage();
			break;
		case 'q':
		case 27 :
			exit(0);
	}
}


/******************************************************************************
 * special keyboard operations to interact with the scene
 *****************************************************************************/
void specialKeyboard(int key, int x, int y)
{
	switch (key) 
	{
		// change the viewing angle (zoom out)
		case GLUT_KEY_UP:
			if(viewAngle > 20.0f) viewAngle -= 0.5f;
			focalLength += 5.0f;
			ChangeRollerCoasterParameters(viewAngle,eyeDist,focalLength,0,0);
			break;
		// change the viewing angle (zoom in)
		case GLUT_KEY_DOWN:
			if(viewAngle < 90.0f) viewAngle += 0.5f;
			if(focalLength > 1.0f) focalLength -= 5.0f;
			ChangeRollerCoasterParameters(viewAngle,eyeDist,focalLength,0,0);
			break;
		// pan left
		case GLUT_KEY_LEFT:
			ChangeRollerCoasterParameters(viewAngle,eyeDist,focalLength,1,0);
			break;
		// pan left
		case GLUT_KEY_RIGHT:
			ChangeRollerCoasterParameters(viewAngle,eyeDist,focalLength,0,1);
			break;
	}
}


/******************************************************************************
 * display the roller coaster
 *****************************************************************************/
void display(void)
{
	DrawRollerCoaster();
}


/******************************************************************************
 * standard reshape
 *****************************************************************************/
void reshape(int w, int h)
{	
	ResizeScene(w,h);
	display();
}


/******************************************************************************
 * everything starts here
 *****************************************************************************/
int main(int argc, char *argv[])
{
	// initilization
	glutInit(&argc, argv);
	glutInitWindowSize(width, height);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(20,20);
	glutCreateWindow("WIRED");
	// function callbacks
	glutDisplayFunc(display);
	glutIdleFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(specialKeyboard);
	InitializeRollerCoaster(intro=1,width,height,viewAngle,eyeDist,focalLength);
	glutMainLoop();

	return 0;
}