/*
 *  SufaceGeometry.c
 *  Carbon OpenGL
 *
 *  Created by Geoff Stahl.
	Copyright:	Copyright © 2002-2003 Apple Computer, Inc., All Rights Reserved

	Disclaimer:	IMPORTANT:  This Apple software is supplied to you by Apple Computer, Inc.
			("Apple") in consideration of your agreement to the following terms, and your
			use, installation, modification or redistribution of this Apple software
			constitutes acceptance of these terms.  If you do not agree with these terms,
			please do not use, install, modify or redistribute this Apple software.

			In consideration of your agreement to abide by the following terms, and subject
			to these terms, Apple grants you a personal, non-exclusive license, under Apple’s
			copyrights in this original Apple software (the "Apple Software"), to use,
			reproduce, modify and redistribute the Apple Software, with or without
			modifications, in source and/or binary forms; provided that if you redistribute
			the Apple Software in its entirety and without modifications, you must retain
			this notice and the following text and disclaimers in all such redistributions of
			the Apple Software.  Neither the name, trademarks, service marks or logos of
			Apple Computer, Inc. may be used to endorse or promote products derived from the
			Apple Software without specific prior written permission from Apple.  Except as
			expressly stated in this notice, no other rights or licenses, express or implied,
			are granted by Apple herein, including but not limited to any patent rights that
			may be infringed by your derivative works or by other works in which the Apple
			Software may be incorporated.

			The Apple Software is provided by Apple on an "AS IS" basis.  APPLE MAKES NO
			WARRANTIES, EXPRESS OR IMPLIED, INCLUDING WITHOUT LIMITATION THE IMPLIED
			WARRANTIES OF NON-INFRINGEMENT, MERCHANTABILITY AND FITNESS FOR A PARTICULAR
			PURPOSE, REGARDING THE APPLE SOFTWARE OR ITS USE AND OPERATION ALONE OR IN
			COMBINATION WITH YOUR PRODUCTS.

			IN NO EVENT SHALL APPLE BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL OR
			CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
			GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
			ARISING IN ANY WAY OUT OF THE USE, REPRODUCTION, MODIFICATION AND/OR DISTRIBUTION
			OF THE APPLE SOFTWARE, HOWEVER CAUSED AND WHETHER UNDER THEORY OF CONTRACT, TORT
			(INCLUDING NEGLIGENCE), STRICT LIABILITY OR OTHERWISE, EVEN IF APPLE HAS BEEN
			ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

 /* Uses techniques described by Paul Bourke 1999 - 2002 */ 
 /* Tranguloid Trefoil and other example surfaces by Roger Bagula see <http://astronomy.swin.edu.au/~pbourke/surfaces/> */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "SurfaceGeometry.h"

#define TWOPI           6.283185307179586476925287
#define PI              3.141592653589793238462643

typedef struct {
   GLfloat x,y,z;
} recVec;

typedef struct {
   GLfloat r,g,b;
} recColor;

typedef struct {
   GLfloat s,t;
} recTexCoord;

/* Code based on work by Paul Bourke */

// globals for apps to use
// info
char gSurfName[256] = "";
char gSurfCredit[256] = "";
char gSurfX[256] = "";
char gSurfY[256] = "";
char gSurfZ[256] = "";
char gSurfRange[256] = "";


// simple cube data
GLint cube_num_vertices = 8;

GLfloat cube_vertices [8][3] = {
{1.0, 1.0, 1.0}, {1.0, -1.0, 1.0}, {-1.0, -1.0, 1.0}, {-1.0, 1.0, 1.0},
{1.0, 1.0, -1.0}, {1.0, -1.0, -1.0}, {-1.0, -1.0, -1.0}, {-1.0, 1.0, -1.0} };

GLfloat cube_vertex_colors [8][3] = {
{1.0, 1.0, 1.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 1.0, 1.0},
{1.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0} };

GLint cube_num_faces = 6;

short cube_faces [6][4] = {
{3, 2, 1, 0}, {2, 3, 7, 6}, {0, 1, 5, 4}, {3, 0, 4, 7}, {1, 2, 6, 5}, {4, 5, 6, 7} };

GLfloat cube_texCoords [2][4] = {
{0.0, 0.0, 1.0, 1.0}, {0.0, 1.0, 1.0, 0.0} };

recColor getColor(GLfloat v, GLfloat vmin, GLfloat vmax, int type)
{
   double dv,vmid;
	recColor c = {1.0,1.0,1.0};
	recColor c1,c2,c3;
	double ratio;

   if (v < vmin)
      v = vmin;
   if (v > vmax)
      v = vmax;
   dv = vmax - vmin;

	switch (type) {
	case 0:
		c.r = 1.0f;
		c.b = 1.0f;
		c.g = 1.0f;
		break;
	case 1:
   	if (v < (vmin + 0.25 * dv)) {
      	c.r = 0;
      	c.g = 4 * (v - vmin) / dv;
			c.b = 1;
   	} else if (v < (vmin + 0.5 * dv)) {
      	c.r = 0;
			c.g = 1;
      	c.b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
   	} else if (v < (vmin + 0.75 * dv)) {
      	c.r = 4 * (v - vmin - 0.5 * dv) / dv;
			c.g = 1;
      	c.b = 0;
   	} else {
			c.r = 1;
      	c.g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
      	c.b = 0;
   	}
		break;
	case 2:
		c.r = (v - vmin) / dv;
		c.g = 0;
		c.b = (vmax - v) / dv;
		break;
	case 3:
		c.r = (v - vmin) / dv;
		c.b = c.r;
		c.g = c.r;
		break;
	case 4:
      if (v < (vmin + dv / 6.0)) {
         c.r = 1; 
         c.g = 6 * (v - vmin) / dv;
         c.b = 0;
      } else if (v < (vmin + 2.0 * dv / 6.0)) {
         c.r = 1 + 6 * (vmin + dv / 6.0 - v) / dv;
         c.g = 1;
         c.b = 0;
      } else if (v < (vmin + 3.0 * dv / 6.0)) {
         c.r = 0;
         c.g = 1;
         c.b = 6 * (v - vmin - 2.0 * dv / 6.0) / dv;
      } else if (v < (vmin + 4.0 * dv / 6.0)) {
         c.r = 0;
         c.g = 1 + 6 * (vmin + 3.0 * dv / 6.0 - v) / dv;
         c.b = 1;
      } else if (v < (vmin + 5.0 * dv / 6.0)) {
         c.r = 6 * (v - vmin - 4.0 * dv / 6.0) / dv;
         c.g = 0;
         c.b = 1;
      } else {
         c.r = 1;
         c.g = 0;
         c.b = 1 + 6 * (vmin + 5.0 * dv / 6.0 - v) / dv;
      }
		break;
   case 5:
      c.r = (v - vmin) / (vmax - vmin);
      c.g = 1;
      c.b = 0;
		break;
   case 6:
      c.r = (v - vmin) / (vmax - vmin);
      c.g = (vmax - v) / (vmax - vmin);
      c.b = c.r;
		break;
   case 7:
      if (v < (vmin + 0.25 * dv)) {
         c.r = 0;
         c.g = 4 * (v - vmin) / dv;
         c.b = 1 - c.g;
      } else if (v < (vmin + 0.5 * dv)) {
			c.r = 4 * (v - vmin - 0.25 * dv) / dv;
         c.g = 1 - c.r;
         c.b = 0;
      } else if (v < (vmin + 0.75 * dv)) {
         c.g = 4 * (v - vmin - 0.5 * dv) / dv;
			c.r = 1 - c.g;
         c.b = 0;
      } else {
         c.r = 0;
         c.b = 4 * (v - vmin - 0.75 * dv) / dv;
			c.g = 1 - c.b;
      }
      break;
   case 8:
      if (v < (vmin + 0.5 * dv)) {
         c.r = 2 * (v - vmin) / dv;
         c.g = c.r;
         c.b = c.r;
      } else {
         c.r = 1 - 2 * (v - vmin - 0.5 * dv) / dv;
			c.g = c.r;
         c.b = c.r;
      }
      break;
   case 9:
      if (v < (vmin + dv / 3)) {
         c.b = 3 * (v - vmin) / dv;
			c.g = 0;
         c.r = 1 - c.b;
      } else if (v < (vmin + 2 * dv / 3)) {
         c.r = 0;
         c.g = 3 * (v - vmin - dv / 3) / dv;
         c.b = 1;
      } else {
         c.r = 3 * (v - vmin - 2 * dv / 3) / dv;
         c.g = 1 - c.r;
			c.b = 1;
      }
      break;
   case 10:
      if (v < (vmin + 0.2 * dv)) {
         c.r = 0;
         c.g = 5 * (v - vmin) / dv;
         c.b = 1;
      } else if (v < (vmin + 0.4 * dv)) {
         c.r = 0;
         c.g = 1;
         c.b = 1 + 5 * (vmin + 0.2 * dv - v) / dv;
      } else if (v < (vmin + 0.6 * dv)) {
         c.r = 5 * (v - vmin - 0.4 * dv) / dv;
         c.g = 1;
         c.b = 0;
      } else if (v < (vmin + 0.8 * dv)) {
         c.r = 1;
         c.g = 1 - 5 * (v - vmin - 0.6 * dv) / dv;
         c.b = 0;
      } else {
         c.r = 1;
         c.g = 5 * (v - vmin - 0.8 * dv) / dv;
         c.b = 5 * (v - vmin - 0.8 * dv) / dv;
      }
      break;
   case 11:
		c1.r = 200 / 255.0; c1.g =  60 / 255.0; c1.b =   0 / 255.0;
		c2.r = 250 / 255.0; c2.g = 160 / 255.0; c2.b = 110 / 255.0;
      c.r = (c2.r - c1.r) * (v - vmin) / dv + c1.r;
      c.g = (c2.g - c1.g) * (v - vmin) / dv + c1.g;
      c.b = (c2.b - c1.b) * (v - vmin) / dv + c1.b;
      break;
	case 12:
		c1.r =  55 / 255.0; c1.g =  55 / 255.0; c1.b =  45 / 255.0;
		c2.r = 200 / 255.0; c2.g =  60 / 255.0; c2.b =   0 / 255.0; 
//		c2.r = 235 / 255.0; c2.g =  90 / 255.0; c2.b =  30 / 255.0;
		c3.r = 250 / 255.0; c3.g = 160 / 255.0; c3.b = 110 / 255.0;
		ratio = 0.4;
		vmid = vmin + ratio * dv;
		if (v < vmid) {
      	c.r = (c2.r - c1.r) * (v - vmin) / (ratio*dv) + c1.r;
      	c.g = (c2.g - c1.g) * (v - vmin) / (ratio*dv) + c1.g;
      	c.b = (c2.b - c1.b) * (v - vmin) / (ratio*dv) + c1.b;
		} else {
         c.r = (c3.r - c2.r) * (v - vmid) / ((1-ratio)*dv) + c2.r;
         c.g = (c3.g - c2.g) * (v - vmid) / ((1-ratio)*dv) + c2.g;
         c.b = (c3.b - c2.b) * (v - vmid) / ((1-ratio)*dv) + c2.b;
		}
		break;
	case 13:
      c1.r =   0 / 255.0; c1.g = 255 / 255.0; c1.b =   0 / 255.0;
      c2.r = 255 / 255.0; c2.g = 150 / 255.0; c2.b =   0 / 255.0;
      c3.r = 255 / 255.0; c3.g = 250 / 255.0; c3.b = 240 / 255.0;
      ratio = 0.3;
      vmid = vmin + ratio * dv;
      if (v < vmid) {
         c.r = (c2.r - c1.r) * (v - vmin) / (ratio*dv) + c1.r;
         c.g = (c2.g - c1.g) * (v - vmin) / (ratio*dv) + c1.g;
         c.b = (c2.b - c1.b) * (v - vmin) / (ratio*dv) + c1.b;
      } else {
         c.r = (c3.r - c2.r) * (v - vmid) / ((1-ratio)*dv) + c2.r;
         c.g = (c3.g - c2.g) * (v - vmid) / ((1-ratio)*dv) + c2.g;
         c.b = (c3.b - c2.b) * (v - vmid) / ((1-ratio)*dv) + c2.b;
      }
		break;
   case 14:
      c.r = 1;
      c.g = 1 - (v - vmin) / dv;
      c.b = 0;
      break;
   case 15:
      if (v < (vmin + 0.25 * dv)) {
         c.r = 0;
         c.g = 4 * (v - vmin) / dv;
         c.b = 1;
      } else if (v < (vmin + 0.5 * dv)) {
         c.r = 0;
         c.g = 1;
         c.b = 1 - 4 * (v - vmin - 0.25 * dv) / dv;
      } else if (v < (vmin + 0.75 * dv)) {
			c.r = 4 * (v - vmin - 0.5 * dv) / dv;
         c.g = 1;
         c.b = 0;
      } else {
         c.r = 1;
			c.g = 1;
         c.b = 4 * (v - vmin - 0.75 * dv) / dv;
      }
      break;
   case 16:
      if (v < (vmin + 0.5 * dv)) {
         c.r = 0.0;
         c.g = 2 * (v - vmin) / dv;
         c.b = 1 - 2 * (v - vmin) / dv;
      } else {
         c.r = 2 * (v - vmin - 0.5 * dv) / dv;
         c.g = 1 - 2 * (v - vmin - 0.5 * dv) / dv;
         c.b = 0.0;
      }
      break;
   case 17:
      if (v < (vmin + 0.5 * dv)) {
         c.r = 1.0;
         c.g = 1 - 2 * (v - vmin) / dv;
         c.b = 2 * (v - vmin) / dv;
      } else {
         c.r = 1 - 2 * (v - vmin - 0.5 * dv) / dv;
         c.g = 2 * (v - vmin - 0.5 * dv) / dv;
         c.b = 1.0;
      }
      break;
   case 18:
      c.r = 0;
      c.g = (v - vmin) / (vmax - vmin);
      c.b = 1;
      break;
   case 19:
      c.r = (v - vmin) / (vmax - vmin);
      c.g = c.r;
      c.b = 1;
      break;
   case 20:
      c1.r =   0 / 255.0; c1.g = 160 / 255.0; c1.b =   0 / 255.0;
      c2.r = 180 / 255.0; c2.g = 220 / 255.0; c2.b =   0 / 255.0;
      c3.r = 250 / 255.0; c3.g = 220 / 255.0; c3.b = 170 / 255.0;
      ratio = 0.3;
      vmid = vmin + ratio * dv;
      if (v < vmid) {
         c.r = (c2.r - c1.r) * (v - vmin) / (ratio*dv) + c1.r;
         c.g = (c2.g - c1.g) * (v - vmin) / (ratio*dv) + c1.g;
         c.b = (c2.b - c1.b) * (v - vmin) / (ratio*dv) + c1.b;
      } else {
         c.r = (c3.r - c2.r) * (v - vmid) / ((1-ratio)*dv) + c2.r;
         c.g = (c3.g - c2.g) * (v - vmid) / ((1-ratio)*dv) + c2.g;
         c.b = (c3.b - c2.b) * (v - vmid) / ((1-ratio)*dv) + c2.b;
      }
      break;
	}
	return(c);
}

void normalise(recVec *p)
{
   double length;

   length = sqrt(p->x * p->x + p->y * p->y + p->z * p->z);
   if (length != 0) {
      p->x /= length;
      p->y /= length;
      p->z /= length;
   } else {
      p->x = 0;
      p->y = 0;
      p->z = 0;
   }
}

recVec CalcNormal(recVec p,recVec p1,recVec p2)
{
   recVec n,pa,pb;

   pa.x = p1.x - p.x;
   pa.y = p1.y - p.y;
   pa.z = p1.z - p.z;
   pb.x = p2.x - p.x;
   pb.y = p2.y - p.y;
   pb.z = p2.z - p.z;
   normalise(&pa);
   normalise(&pb);
  
   n.x = pa.y * pb.z - pa.z * pb.y;
   n.y = pa.z * pb.x - pa.x * pb.z;
   n.z = pa.x * pb.y - pa.y * pb.x;
   normalise(&n);

   return(n);
}

// expects u & v (-PI to PI)
recVec Eval(double u, double v, int type)
{
	recVec p;
	double temp;
   
	switch (type) {
		case kTranguloidTrefoil:
			p.x = sin(3*u) * 2 / (2 + cos(v));
			p.y = (sin(u) + 2 * sin(2*u)) * 2 / (2 + cos(v + TWOPI / 3));
			p.z = (cos(u) - 2 * cos(2*u)) * (2 + cos(v)) * (2 + cos(v + TWOPI/3))/4;
		break;
		case kTriaxialTritorus:
			p.x = 2.0 * sin (u) * (1 + cos (v));
			p.y = 2.0 * sin (u + 2 * PI / 3) * (1 + cos (v + 2 * PI / 3));
			p.z = 2.0 * sin (u + 4 * PI / 3) * (1 + cos (v + 4 * PI / 3));
		break;
		case kStilettoSurface:
			// reverse u and v for better distribution or points
			temp = u;
			u = v + PI; v = (temp + PI) / 2.0; // convert to: 0 <= u <= 2 pi, 0 <= v <= 2 pi 
			p.x = 4.0 *  (2.0 + cos(u)) * pow(cos(v), 3.0) * sin(v);
			p.y = 4.0 *  (2.0 + cos(u+TWOPI/3.0)) * pow (cos(v+TWOPI/3.0), 2.0) * pow (sin(v+TWOPI/3.0), 2.0);
			p.z = 4.0 * -(2.0 + cos(u-TWOPI/3.0)) * pow (cos(v+TWOPI/3.0), 2.0) * pow (sin(v+TWOPI/3.0), 2.0);
 		break;
		case kSlippersSurface:
			temp = u;
			u = v + PI * 2; v = temp + PI; // convert to: 0 <= u <= 4 pi, 0 <= v <= 2 pi 
			p.x = 4.0 *  (2 + cos (u)) * pow (cos (v), 3) * sin(v);
			p.y = 4.0 *  (2 + cos (u + TWOPI / 3)) * pow (cos (TWOPI / 3 + v), 2) * pow (sin (TWOPI / 3 + v), 2);
			p.z = 4.0 * -(2 + cos (u - TWOPI / 3)) * pow (cos (TWOPI / 3 - v), 2) * pow (sin (TWOPI / 3 - v), 3);
		break;
		case kMaedersOwl:
			u = (u + PI) * 2; v = (v + PI) / TWOPI; // convert to: 0 <= u <= 4 pi, 0 <= v <= 1 
			p.x = 3.0 *  v * cos(u) - 0.5 * v * v * cos(2 * u);
			p.y = 3.0 * -v * sin(u) - 0.5 * v * v * sin(2 * u);
			p.z = 3.0 *  4 * pow(v,1.5) * cos(1.5 * u) / 3;
		break;
		default:
			p.x = 0.0;
			p.y = 0.0;
			p.z = 0.0;
		break;
 	}
   return(p);
}

void GetStrings (unsigned int surface, char ** strName, char ** strAuthor, char ** strX, char ** strY, char ** strZ, char ** strRange)
{
	static char strings[6][6][256] = {{"Color Cube", 
									   " ", 
									   " ", 
									   " ", 
									   " ", 
									   " "},
									  {"Tranguloid Trefoil", 
									   "by Roger Bagula", 
									   "x = 2 sin(3 u) / (2 + cos(v))", 
									   "y = 2 (sin(u) + 2 sin(2 u)) / (2 + cos(v + 2 pi / 3))", 
									   "z = (cos(u) - 2 cos(2 u)) (2 + cos(v)) (2 + cos(v + 2 pi / 3)) / 4", 
									   "-pi <= u <= pi, -pi <= v <= pi"},
									  {"Triaxial Tritorus", 
									   "by Roger Bagula", 
									   "x = sin(u) (1 + cos(v))", 
									   "y = sin(u + 2pi / 3) (1 + cos(v + 2pi / 3))", 
									   "z = z = sin(u + 4pi / 3) (1 + cos(v + 4pi / 3))", 
									   "0 <= u <= 2 pi, 0 <= v <= 2 pi"},
									  {"Stiletto Surface", 
									   "by Roger Bagula", 
									   "x =  (2 + cos(u)) cos(v)^3 sin(v)", 
									   "y =  (2 + cos(u + 2pi /3)) cos(v + 2pi / 3)^2 sin(v + 2pi / 3)^2", 
									   "z = -(2 + cos(u - 2pi / 3)) cos(v + 2pi / 3)^2 sin(v + 2pi / 3)^2", 
									   "0 <= u <= 2 pi, 0 <= v <= 2 pi"},
									  {"Slippers Surface", 
									   "by Roger Bagula", 
									   "x =  (2 + cos(u)) cos(v)^3 sin(v)", 
									   "y =  (2 + cos(u + 2pi / 3)) cos(2pi / 3 + v)^2 sin(2pi / 3 + v)^2", 
									   "z = -(2 + cos(u - 2pi / 3)) cos(2pi / 3 - v)^2 sin(2pi / 3 - v)^3", 
									   "0 <= u <= 2 pi, 0 <= v <= 2 pi"},
									  {"Maeder's Owl", 
									   "by R. Maeder", 
									   "x = v cos(u) - 0.5 v^2 cos(2 u)", 
									   "y = - v sin(u) - 0.5 v^2 sin(2 u)", 
									   "z = 4 v^1.5 cos(3 u / 2) / 3", 
									   "0 <= u <= 4 pi, 0 <= v <= 1"}};
	*strName = strings[surface][0];
	*strAuthor = strings[surface][1];
	*strX = strings[surface][2];
	*strY = strings[surface][3];
	*strZ = strings[surface][4];
	*strRange = strings[surface][5];
}

// special case cube for test purposes
void BuildCube (GLuint * polyList, GLuint * lineList, GLuint * pointList, int colorScheme)
{
	float fSize = 2.0f;
	long f, i;
	*polyList = glGenLists (1);
	glNewList(*polyList, GL_COMPILE);
		
		glBegin (GL_QUADS);
		for (f = 0; f < cube_num_faces; f++)
			for (i = 0; i < 4; i++) {
				if (colorScheme)
					glColor3f (cube_vertex_colors[cube_faces[f][i]][0], cube_vertex_colors[cube_faces[f][i]][1], cube_vertex_colors[cube_faces[f][i]][2]);
				else
					glColor3f (1.0f, 1.0f, 1.0f);
				glTexCoord2f (cube_texCoords [0][i], cube_texCoords [1][i]);
				glVertex3f(cube_vertices[cube_faces[f][i]][0] * fSize, cube_vertices[cube_faces[f][i]][1] * fSize, cube_vertices[cube_faces[f][i]][2] * fSize);
			}
		glEnd ();
		glColor3f (0.0, 0.0, 0.0);
		for (f = 0; f < cube_num_faces; f++) {
			glBegin (GL_LINE_LOOP);
				for (i = 0; i < 4; i++)
					glVertex3f(cube_vertices[cube_faces[f][i]][0] * fSize, cube_vertices[cube_faces[f][i]][1] * fSize, cube_vertices[cube_faces[f][i]][2] * fSize);
			glEnd ();
		}
	glEndList ();
	
	*lineList = glGenLists (1);
	glNewList(*lineList, GL_COMPILE);
		glColor3f (1.0, 1.0, 1.0);
		for (f = 0; f < cube_num_faces; f++) {
			glBegin (GL_LINE_LOOP);
				for (i = 0; i < 4; i++)
					glVertex3f(cube_vertices[cube_faces[f][i]][0] * fSize, cube_vertices[cube_faces[f][i]][1] * fSize, cube_vertices[cube_faces[f][i]][2] * fSize);
			glEnd ();
		}
	glEndList ();
	
	*pointList = glGenLists (1);
	glNewList(*pointList, GL_COMPILE);
		glColor3f (1.0, 1.0, 1.0);
		for (f = 0; f < cube_num_vertices; f++) {
			glBegin (GL_POINTS);
					glVertex3f(cube_vertices[f][0] * fSize, cube_vertices[f][1] * fSize, cube_vertices[f][2] * fSize);
			glEnd ();
		}
	glEndList ();
}

void BuildGeometry (unsigned int surface, unsigned int colorScheme, unsigned int subdivisions, unsigned int xyRatio,
					GLuint * polyList, GLuint * lineList, GLuint * pointList)
{
	long i,j, index;
	long maxI = subdivisions * xyRatio, maxJ = subdivisions;
	double u, v, delta=0.001;
	recVec p1,p2;
	recVec *vertexPos = NULL,*vertexNormal = NULL;
	recColor *vertexColor = NULL;
	recTexCoord *vertexTexCoord = NULL;

	// set valid surface and color scheme
	surface %= kSurfaces;
	colorScheme %= kColorSchemes;

	// delete existing list
	if (*polyList)
		glDeleteLists (*polyList, 1);
	if (*lineList)
		glDeleteLists (*lineList, 1);
	if (*pointList)
		glDeleteLists (*pointList, 1);
	*polyList = *lineList = *pointList = 0;
	
	if (surface == kCube) // build the standard color cube (disregard color, subdivisions, and xyRatio)
		BuildCube (polyList, lineList, pointList, colorScheme);
	else {
		// build buffers
		vertexPos = (recVec*) malloc ((maxI) * (maxJ) * sizeof (recVec));
		if (vertexNormal)
			free (vertexNormal);
		vertexNormal = (recVec*) malloc ((maxI) * (maxJ) * sizeof (recVec));
		if (vertexColor)
			free (vertexColor);
		vertexColor = (recColor*) malloc ((maxI) * (maxJ) * sizeof (recColor));
		if (vertexTexCoord)
			free (vertexTexCoord);
		vertexTexCoord = (recTexCoord*) malloc ((maxI) * (maxJ) * sizeof (recTexCoord));
		if (!vertexPos || !vertexNormal || !vertexColor || !vertexTexCoord)
			return;
			
		// build surface
		for (i = 0; i < maxI; i++) {
			for (j = 0; j < maxJ; j++) {
				index = i * maxJ + j;
				u  = -PI + (i % maxI) * TWOPI / maxI;
				v  = -PI + (j % maxJ) * TWOPI / maxJ;
				vertexPos[index] = Eval(u,v, surface);
				p1 = Eval(u + delta, v, surface);
				p2 = Eval(u, v + delta, surface);
				vertexNormal[index] = CalcNormal(vertexPos[index],p1,p2);
				vertexColor[index] = getColor(u, -PI, PI, colorScheme);
				vertexTexCoord[index].s = (float) i * 5.0f / (float) maxI;
				vertexTexCoord[index].t = (float) j * 1.0f/ (float) maxJ;
			}
		}
		
		*polyList = glGenLists (1);
		glNewList(*polyList, GL_COMPILE);
			for (i=0; i< maxI; i++) {
				glBegin(GL_TRIANGLE_STRIP);
				for (j = 0; j <= maxJ; j++) {
					index = (i % maxI) * maxJ + (j % maxJ);
					glColor3fv (&vertexColor[index].r);
					glNormal3fv (&vertexNormal[index].x);
					glTexCoord2fv (&vertexTexCoord[index].s);
					glVertex3fv (&vertexPos[index].x);
		
					index = ((i + 1) % maxI) * maxJ + (j % maxJ);
					glColor3fv (&vertexColor[index].r);
					glNormal3fv (&vertexNormal[index].x);
					glTexCoord2fv (&vertexTexCoord[index].s);
					glVertex3fv (&vertexPos[index].x);
//					index = ((i - 1) % maxI) * maxJ + (j % maxJ);
				}
				glEnd ();
			}
		glEndList ();
	
		*lineList = glGenLists (1);
		glNewList(*lineList, GL_COMPILE);
			for (i=0; i< maxI; i++) {
				glBegin(GL_LINE_STRIP);
				for (j = 0; j < maxJ; j++) {
					index = i * maxJ + j;
					glColor3fv (&vertexColor[index].r);
					glVertex3fv (&vertexPos[index].x);
				}
				index = i * maxJ + 0;
				glColor3fv (&vertexColor[index].r);
				glVertex3fv (&vertexPos[index].x);
				glEnd ();
			}
			for (j=0; j< maxJ; j++) {
				glBegin(GL_LINE_STRIP);
				for (i = 0; i < maxI; i++) {
					index = i * maxJ + j;
					glColor3fv (&vertexColor[index].r);
					glVertex3fv (&vertexPos[index].x);
				}
				index = 0 + j;
				glColor3fv (&vertexColor[index].r);
				glVertex3fv (&vertexPos[index].x);
				glEnd ();
			}
		glEndList ();
	
		*pointList = glGenLists (1);
		glNewList(*pointList, GL_COMPILE);
			glBegin(GL_POINTS);
			for (i=0; i< maxI; i++) {
				for (j = 0; j < maxJ; j++) {
					index = i * maxJ + j;
					glColor3fv (&vertexColor[index].r);
					glVertex3fv (&vertexPos[index].x);
				}
			}
			glEnd ();
		glEndList ();
		free (vertexPos);
		free (vertexNormal);
		free (vertexColor);
	}
}
