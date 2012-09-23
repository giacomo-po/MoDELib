/* Copyright (c) Mark J. Kilgard, 1997. */

/* This program is freely distributable without licensing fees
 and is provided without guarantee or warrantee expressed or
 implied. This program is -not- in the public domain. */

/* Example showing how to use OpenGL's feedback mode to capture
 transformed vertices and output them as Encapsulated PostScript.
 Handles limited hidden surface removal by sorting and does
 smooth shading (albeit limited due to PostScript). */

/* Compile: cc -o rendereps rendereps.c -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm */

//#include <assert.h>
//#include <math.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <GLUT/glut.h>
//#include <GL/glut.h>

/* OpenGL's GL_3D_COLOR feedback vertex format. */
typedef struct _Feedback3Dcolor {
    GLfloat x;
    GLfloat y;
    GLfloat z;
    GLfloat red;
    GLfloat green;
    GLfloat blue;
    GLfloat alphaa;
} Feedback3Dcolor;

//int blackBackground = 0;  /* Initially use a white background. */
//int lighting = 0;       /* Initially disable lighting. */
//int polygonMode = 1;    /* Initially show wireframe. */
//int object = 1;         /* Initially show the torus. */

//GLfloat angle = 0.0;    /* Angle of rotation for object. */
//int moving, begin;      /* For interactive object rotation. */
int size = 1;           /* Size of lines and points. */

/* How many feedback buffer GLfloats each of the three objects need. */
int objectComplexity[3] =
{6000, 14000, 380000};  /* Teapot requires ~1.5 megabytes for
                         its feedback results! */

/* render gets called both by "display" (in OpenGL render mode)
 and by "outputEPS" (in OpenGL feedback mode). */
//void
//render(void)
//{
//    glPushMatrix();
//    glRotatef(angle, 0.0, 1.0, 0.0);
//    switch (object) {
//        case 0:
//            glutSolidSphere(1.0, 10, 10);
//            break;
//        case 1:
//            glutSolidTorus(0.5, 1.0, 15, 15);
//            break;
//        case 2:
//            glutSolidTeapot(1.0);
//            break;
//    }
//    glPopMatrix();
//}

//void
//display(void)
//{
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//    render();
//    glutSwapBuffers();
//}

//void
//updateBackground(void)
//{
//    if (blackBackground) {
//        /* Clear to black. */
//        glClearColor(0.0, 0.0, 0.0, 1.0);
//    } else {
//        /* Clear to white. */
//        glClearColor(1.0, 1.0, 1.0, 1.0);
//    }
//}

//void
//updateLighting(void)
//{
//    if (lighting) {
//        glEnable(GL_LIGHTING);
//    } else {
//        glDisable(GL_LIGHTING);
//    }
//}
//
//void
//updatePolygonMode(void)
//{
//    switch (polygonMode) {
//        case 0:
//            glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
//            break;
//        case 1:
//            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
//            break;
//        case 2:
//            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
//            break;
//    }
//}



GLfloat pointSize;

static char *gouraudtriangleEPS[] =
{
    "/bd{bind def}bind def /triangle { aload pop   setrgbcolor  aload pop 5 3",
    "roll 4 2 roll 3 2 roll exch moveto lineto lineto closepath fill } bd",
    "/computediff1 { 2 copy sub abs threshold ge {pop pop pop true} { exch 2",
    "index sub abs threshold ge { pop pop true} { sub abs threshold ge } ifelse",
    "} ifelse } bd /computediff3 { 3 copy 0 get 3 1 roll 0 get 3 1 roll 0 get",
    "computediff1 {true} { 3 copy 1 get 3 1 roll 1 get 3 1 roll 1 get",
    "computediff1 {true} { 3 copy 2 get 3 1 roll  2 get 3 1 roll 2 get",
    "computediff1 } ifelse } ifelse } bd /middlecolor { aload pop 4 -1 roll",
    "aload pop 4 -1 roll add 2 div 5 1 roll 3 -1 roll add 2 div 3 1 roll add 2",
    "div 3 1 roll exch 3 array astore } bd /gouraudtriangle { computediff3 { 4",
    "-1 roll aload 7 1 roll 6 -1 roll pop 3 -1 roll pop add 2 div 3 1 roll add",
    "2 div exch 3 -1 roll aload 7 1 roll exch pop 4 -1 roll pop add 2 div 3 1",
    "roll add 2 div exch 3 -1 roll aload 7 1 roll pop 3 -1 roll pop add 2 div 3",
    "1 roll add 2 div exch 7 3 roll 10 -3 roll dup 3 index middlecolor 4 1 roll",
    "2 copy middlecolor 4 1 roll 3 copy pop middlecolor 4 1 roll 13 -1 roll",
    "aload pop 17 index 6 index 15 index 19 index 6 index 17 index 6 array",
    "astore 10 index 10 index 14 index gouraudtriangle 17 index 5 index 17",
    "index 19 index 5 index 19 index 6 array astore 10 index 9 index 13 index",
    "gouraudtriangle 13 index 16 index 5 index 15 index 18 index 5 index 6",
    "array astore 12 index 12 index 9 index gouraudtriangle 17 index 16 index",
    "15 index 19 index 18 index 17 index 6 array astore 10 index 12 index 14",
    "index gouraudtriangle 18 {pop} repeat } { aload pop 5 3 roll aload pop 7 3",
    "roll aload pop 9 3 roll 4 index 6 index 4 index add add 3 div 10 1 roll 7",
    "index 5 index 3 index add add 3 div 10 1 roll 6 index 4 index 2 index add",
    "add 3 div 10 1 roll 9 {pop} repeat 3 array astore triangle } ifelse } bd",
    NULL
};

