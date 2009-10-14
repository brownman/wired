/*****************************************************************************
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
 *      SK = my SECTION
 *      WIRED = the new roller coaster ride that should be coming to an 
 *				amusement park near you in the future!
 *
 * File: RollerCoasterLoader.c
 * Desc: this is the parser for the coaster definition file
 *		 initially the data was static in a function
 *		 i created the parser that can allow for other coasters to be defined
 *		 the parser loads the file
 *				- validates all sections
 *				- parses each section
 *				- loads the corresponding data
 *				- please see COASTER.FORMAT.README.txt for details
 *
 * NOTE: i had written many parsers before for other courses and internships
 *		 used the same ideas in parsing the file
 *		 used the ANSI stream libraries and the vc++ stream libraries
 *****************************************************************************/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "Point.h"
#include "RollerCoaster.h"
#include "RollerCoasterLoader.h"

// file information
#define LF 0x0A  // windows line feed
#define CR 0x0D  // windows carriage return

FILE *file;
static int line_number = 1;


/******************************************************************************
 * print an error that occurred during the parsing
 * the line number will be specified as well the appropriate error
 * the use of varargs allows for specifying an error when it is encountered
 *****************************************************************************/
void PrintError(char *str, ...)
{
	char buf1[2*BUFSIZ];
	char buf2[BUFSIZ];
	sprintf(buf1,"Error line %d: ",line_number);
	sprintf(buf2,"%s\n",str);
	strcat(buf1,buf2);
	fprintf(stderr,buf1);
}


/******************************************************************************
 * this function makes many translations
 * upper case -> lower case : read lower case characters in the file
 * carriage return -> space : compatible with dos/windows file format
 *****************************************************************************/
int ReadCharacter()
{
	int c;
	c = getc(file);
	if (c >= 'A' && c<='Z') return c - 'A' + 'a';
	if (c == CR) return ' ';
	if (c == LF) line_number++;
	return c;
}


/******************************************************************************
 * unget character from stream
 *****************************************************************************/
int UngetCharacter(int c)
{
	if (c == LF) line_number--;
	if (c == EOF) return EOF;
	// pushes a character back onto the stream
	else return ungetc(c,file);
}


/******************************************************************************
 * skip whitespace
 *****************************************************************************/
void SkipWhitespace()
{
	int c;

	while((c = ReadCharacter()) != EOF)
	{
		if (c != ' ' && c != '\t')
		{
			UngetCharacter(c);
			return;
		}
	}
}


/******************************************************************************
 * skip comment
 * comment will start with '#'
 *****************************************************************************/
void SkipComment()
{
	int c;

	SkipWhitespace();
	if ((c = ReadCharacter()) == EOF) return;
	while(c == LF || c=='#')
	{
		if (c == '#')
		{
			c = ReadCharacter();
			while(c != EOF && c != LF) c = ReadCharacter();
		}
		SkipWhitespace();
		c = ReadCharacter();
	}
	UngetCharacter(c);
}


/******************************************************************************
 * read float
 * checks for positive and negative float
 * checks for a valid range of numbers
 *****************************************************************************/
int ReadFloat(float *f)
{
	int c;
	int minus = 0;
	float d = 0.0f;

	*f = 0.0f;
	
	if ((c = ReadCharacter()) == EOF) return 0;
	if (c == '-') minus = 1;
	else if (c != '+') UngetCharacter(c);

	c = ReadCharacter();
	if ((c < '0' || c > '9') && c != '.')
	{
		UngetCharacter(c);
		return 0;
	}
	UngetCharacter(c);
	while((c = ReadCharacter()) != EOF)
	{
		if (c == '.')
		{
			if (d == 0.0f)
			{
				d = 1.0f;
				continue;
			}
			else
			{
				UngetCharacter(c);
				return 1;
			}
		}
		if (c >= '0' && c <= '9')
		{
			if (d == 0.0f) *f = *f * 10.0f + (float)(c - '0');
			else *f += (float)(c - '0') / (d *= 10.0f);
		 	continue;
		}
		UngetCharacter(c);
		if (d != 1.0f)
		{
			if (minus) *f = -*f;
			return 1;
		}
		else return 0;
	}
	if (d != 1.0f)
	{
			if (minus) *f = -*f;
			return 1;
	}
	else return 0;
}


/******************************************************************************
 * rean an integer
 * validate the range of the integer
 *****************************************************************************/
int ReadInteger(int *i)
{
	int c;

	*i = 0;
	if ((c = ReadCharacter()) == EOF) return 0;
	if (c < '0' || c > '9')
	{
		UngetCharacter(c);
		return 0;
	}
	UngetCharacter(c);
	while((c = ReadCharacter()) != EOF)
	{
		if (c >= '0' && c <= '9')
		{
			*i = *i * 10 + c - '0';
		 	continue;
		}
		UngetCharacter(c);
		return 1;
	}
	return 1;
}


/******************************************************************************
 * read a section SECTION
 * construct the section SECTION character by character
 *****************************************************************************/
int ReadSectionName(char *str, int length)
{
	int c;
	int i = 0;

	str[i] = 0;
	if ((c = ReadCharacter()) == EOF) return 0;
	UngetCharacter(c);
	if (c < 'a' || c > 'z') return 0;
	while((c = ReadCharacter()) != EOF)
	{
		if (c >= 'a' && c <= 'z')
		{
			if (i < length - 1) str[i++] = c;
			else str[i] = 0;
			continue;
		}
		str[i] = 0;
		UngetCharacter(c);
		return 1;
	}
	str[i] = 0;
	return 1;
}


/******************************************************************************
 * skip newline
 *****************************************************************************/
int SkipNewline()
{
	int old = line_number;

	SkipComment();
	if (line_number - old > 0) return 1;
	else
	{
		PrintError("Newline expected");
		return 0;
	}
}


/******************************************************************************
 * read signed integer
 * used in the coaster definition file
 *****************************************************************************/
int ReadSignedInt(int *i)
{
	int c, minus = 0;

	if ((c = ReadCharacter()) == EOF)
	{
		PrintError("Integer expected");
		return 0;
	}
	switch(c)
	{
		case '-':
			minus = 1;
			break;
		case '+':
			break;
		default :
			UngetCharacter(c);
	}
	if (!ReadInteger(i))
	{
		PrintError("Integer expected");
		return 0;
	}
	if (minus) *i = -*i;
	return 1;
}


/******************************************************************************
 * read unsigned integer
 * this ensures that a positive integer is read
 *****************************************************************************/
int ReadUnsignedInt(int *i)
{
	if (!ReadInteger(i))
	{
		PrintError("Positive integer expected");
		return 0;
	}
	return 1;
}


/******************************************************************************
 * read point information and load it
 * point information will be read as x,y,z float coordinates
 *****************************************************************************/
int ReadPoint(point *p)
{
	if (!ReadFloat(&p->x))
	{
		PrintError("Number expected");
		return 0;
	}
	SkipWhitespace();
	if (!ReadFloat(&p->y))
	{
		PrintError("Number expected");
		return 0;
	}
	SkipWhitespace();
	if (!ReadFloat(&p->z))
	{
		PrintError("Number expected");
		return 0;
	}
	return 1;
}


/******************************************************************************
 * read separator
 * separator can be ';' or ','
 *****************************************************************************/
int ReadSeparator()
{
	int c;

	SkipWhitespace();
	c = ReadCharacter();
	if (c == ';' || c==',')
	{
		SkipWhitespace();
		return 1;
	}
	else
	{
		UngetCharacter(c);
		PrintError("Separator expected");
		return 0;
	}
}


/******************************************************************************
 * read the input file
 * parse each section as it occurs
 * load the data into the coaster variables
 * see COASTER.FORMAT.README.txt for details
 *****************************************************************************/
int ReadInputFile()
{
	char SECTION[BUFSIZ/2];
	int i, sectionOK;

	SkipComment();
	while(!feof(file))
	{
		if (!ReadSectionName(SECTION,BUFSIZ/2))
		{
			PrintError("Section SECTION expected");
			return 0;
		}
		// read each section and parse it
		sectionOK = 0;
		if (!strcmp(SECTION,"trackcontrolpoints"))
		{
			sectionOK = 1;
			if (!SkipNewline()) return 0;
			if (!ReadUnsignedInt(&nbPointControl)) return 0;
			
			pPointControl = (point*)malloc(nbPointControl*2*sizeof(point));
			for(i=0; i<nbPointControl; i++)
			{
				if (!SkipNewline()) return 0;
				if (!ReadPoint(&pPointControl[i*2])) return 0;
				if (!ReadSeparator()) return 0;
				if (!ReadPoint(&pPointControl[i*2+1])) return 0;
			}
		}
		if (!strcmp(SECTION,"supportscoordinate"))
		{
			sectionOK = 1;
			if (!SkipNewline()) return 0;
			if (!ReadUnsignedInt(&nbColumnCoord)) return 0;
			pColumnCoordinate = (int*)malloc(nbColumnCoord*sizeof(int));
			for(i=0; i<nbColumnCoord; i++)
			{
				if (!SkipNewline()) return 0;
				if (!ReadUnsignedInt(&pColumnCoordinate[i])) return 0;
			}
		}
		if (!strcmp(SECTION,"platform"))
		{
			sectionOK = 1;
			if (!SkipNewline()) return 0;
			if (!ReadPoint(&metalPosition)) return 0;
			if (!SkipNewline()) return 0;
			if (!ReadFloat(&metalLength))
			{
				PrintError("Number expected");
				return 0;
			}
			if (!SkipNewline()) return 0;
			if (!ReadFloat(&metalAngle))
			{
				PrintError("Number expected");
				return 0;
			}
		}
		if (!strcmp(SECTION,"trees"))
		{
			sectionOK = 1;
			if (!SkipNewline()) return 0;
			if (!ReadUnsignedInt(&nbTree)) return 0;
			pTree = (point*)malloc(nbTree*sizeof(point));
			for(i=0; i<nbTree; i++)
			{
				if (!SkipNewline()) return 0;
				if (!ReadPoint(&pTree[i])) return 0;
				pTree[i].z = -0.001f;
			}
		}
		if (!strcmp(SECTION,"startsegment"))
		{
			sectionOK = 1;
			if (!SkipNewline()) return 0;
			if (!ReadSignedInt(&startSegment)) return 0;
		}
		if (!strcmp(SECTION,"brakesegment"))
		{
			sectionOK = 1;
			if (!SkipNewline()) return 0;
			if (!ReadSignedInt(&brakeSegment)) return 0;
		}
		if (!strcmp(SECTION,"segmentlength"))
		{
			sectionOK = 1;
			if (!SkipNewline()) return 0;
			if (!ReadFloat(&avgSegmentLength))
			{
				PrintError("Number expected");
				return 0;
			}
		}
		if (!strcmp(SECTION,"bankfactor"))
		{
			sectionOK = 1;
			if (!SkipNewline()) return 0;
			if (!ReadFloat(&twistFactor))
			{
				PrintError("Number expected");
				return 0;
			}
		}
		if (!sectionOK) PrintError("Unknown section \"%s\"",SECTION);
		// don't print an error if there is no new line at the end of the file, 
		// don't directly call skipnewline
		i = line_number;
		SkipComment();
		if ((line_number - i == 0) && !feof(file)) PrintError("Newline expected");
	}
	return 1;
}


/******************************************************************************
 * initialize the variables imported from external source
 *****************************************************************************/
void Initialize()
{
	nbPointControl = 0;
	nbColumnCoord = 0;
	nbTree = 0;
	metalLength = 0.0f;
	startSegment = 7;
	brakeSegment = 0;
	twistFactor = 5.0f;
	avgSegmentLength = 0.15f;
}


/******************************************************************************
 * parse the file
 * call readfile which will parse the entire file
 *****************************************************************************/
int ParseFile()
{
	Initialize();
	if (ReadInputFile()) if (nbPointControl != 0) return 1;
	return 0;
}


/******************************************************************************
 * open file for reading
 * ensure file opened
 *****************************************************************************/
int OpenInputFile(char *str)
{
	file = fopen(str,"r");
	if (file == NULL) Initialize(); return 1;
	return 0;
}


/******************************************************************************
 * close the file after reading
 *****************************************************************************/
int CloseFile()
{
	return fclose(file);
}