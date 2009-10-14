/******************************************************************************
 * Sukhchander Khanna [ khannas@rpi.edu ]
 * Computer Graphics Fall 04
 * Final Project Component One -- 'WIRED' -- My Roller Coaster Creation
 *
 * File: ImageLoader.h
 * Desc: OpenGL C Port of Delhpi BMP and TGA Image Loading utility (SGI)
 *		 This came undocumented. Guaranteed to work.
 *****************************************************************************/

#ifndef IMAGE_LOADER_H
#define IMAGE_LOADER_H

#ifdef _WIN32
#include <windows.h>
#endif

#include <GL/gl.h>

#ifndef _WIN32
/*
 * the following definitions replaces windows headers definition
 * for BMP manipulation
 */
#define BI_RGB 0
#define BI_BITFIELDS 3

typedef unsigned char BYTE;
typedef unsigned short WORD;
typedef unsigned long DWORD;
typedef long LONG;

typedef struct
{
	BYTE rgbBlue;
	BYTE rgbGreen;
	BYTE rgbRed;
	BYTE rgbReserved;
} RGBQUAD;

typedef struct
{
	BYTE rgbtBlue;
	BYTE rgbtGreen;
	BYTE rgbtRed;
} RGBTRIPLE;

typedef struct
{
	WORD bfType __attribute__((packed));
	DWORD bfSize __attribute__((packed));
	WORD bfReserved1 __attribute__((packed));
	WORD bfReserved2 __attribute__((packed));
	DWORD bfOffBits __attribute__((packed));
} BITMAPFILEHEADER;

typedef struct
{
	DWORD biSize __attribute__((packed));
	LONG biWidth __attribute__((packed));
	LONG biHeight __attribute__((packed));
	WORD biPlanes __attribute__((packed));
	WORD biBitCount __attribute__((packed));
	DWORD biCompression __attribute__((packed));
	DWORD biSizeImage __attribute__((packed));
	LONG biXPelsPerMeter __attribute__((packed));
	LONG biYPelsPerMeter __attribute__((packed));
	DWORD biClrUsed __attribute__((packed));
	DWORD biClrImportant __attribute__((packed));
} BITMAPINFOHEADER;

typedef struct
{
	DWORD bcSize;
	WORD bcWidth;
	WORD bcHeight;
	WORD bcPlanes;
	WORD bcBitCount;
} BITMAPCOREHEADER;

#endif /* #ifndef _WIN32 */

typedef struct
{
	BYTE red;
	BYTE green;
	BYTE blue;
} GLRGBTRIPLE;

typedef struct
{
	BYTE red;
	BYTE green;
	BYTE blue;
	BYTE alpha;
} GLRGBQUAD;

typedef struct
{
	BYTE tfType;
	BYTE tfColorMapType;
	BYTE tfImageType;
	BYTE tfColorMapSpec[5];
	BYTE tfOrigX[2];
	BYTE tfOrigY[2];
	BYTE tfWidth[2];
	BYTE tfHeight[2];
	BYTE tfBpp;
	BYTE tfImageDes;
} TGAHEADER;


typedef struct
{
	GLuint texID;
	int width;
	int height;
	int colorDepth;

	GLRGBQUAD *pData;
	RGBQUAD colors[256];

	DWORD redMask;
	DWORD greenMask;
	DWORD blueMask;

	GLint texWrapS;
	GLint texWrapT;
	GLint minFilter;
	GLint magFilter;

	int imageSize;
} glBmpImage;

void glBmpInit(glBmpImage* img);

int glBmpInvert(glBmpImage* img);
int glBmpFlipVert(glBmpImage* img);
int glBmpFlipHorz(glBmpImage* img);
int glBmpRotate180(glBmpImage* img);

void glBmpGenTexture(glBmpImage* img);
void glBmpGenTextureMipMap(glBmpImage* img);
void glBmpSetTextureWrap(glBmpImage* img, GLint s, GLint y);
void glBmpSetFilter(glBmpImage* img, GLint min, GLint mag);
void glBmpBind(glBmpImage* img);

int glBmpGetWidth(glBmpImage* img);
int glBmpGetHeight(glBmpImage* img);
int glBmpGetColorDepth(glBmpImage* img);

int glBmpLoadImage(glBmpImage* img, char *szFileName);
int glBmpLoadBMP(glBmpImage* img, char *filename);
int glBmpLoadTGA(glBmpImage* img, char *filename);

int glBmpMemLoadBMP(glBmpImage* img, char *filePtr);
int glBmpMemLoadTGA(glBmpImage* img, char *filePtr);

int glBmpSaveBMP(glBmpImage* img, char *filename);
int glBmpSaveTGA(glBmpImage* img, char *filename);

int glBmpSaveScreen(glBmpImage* img);
int glBmpSaveGLBuffer(glBmpImage* img, char *szFileName);

#endif