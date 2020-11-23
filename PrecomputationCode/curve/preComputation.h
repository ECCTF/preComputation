#pragma once

#include "miracl.h"
 
long long cpucycles(void);

struct affinePoint
{//Points in Lambda coordinates. y is lambda
	big x;
	big y;
};


struct numberC
{//c[i]=g+h\tau\mu
	int g;
	int h;
};

struct affinePointLD
{//Points in Lambda coordinates. y is lambda
	big x;
	big y;
};

struct projectPointLD
{//Points in Lambda-projective coordinates. y is lambda
	big X;
	big Y;
	big Z;
};
 

struct projectPoint
{//Points in Lambda-projective coordinates. y is lambda
	big X;
	big Y;
	big Z;
};

struct affinePointu4
{//Points in mu4 coordinates.  
	big X0;
	big X1;
	//big X2;
	big X3;
};


struct projectPointu4
{//Points in mu4-projective coordinates.  if X2==1, it is a affine point
	big X0;
	big X1;
	big X2;
	big X3;
};
int tau(affinePoint P,  projectPoint &P1);
int tau(projectPoint P,  projectPoint &P1);

int tau(affinePointLD P, projectPointLD &P1);
int tau(projectPointLD P, projectPointLD &P1);

int tau(affinePointu4 P, projectPointu4 &P1);
int tau(projectPointu4 P, projectPointu4 &P1);

int mixedAddition(affinePoint P1, projectPoint P2,  projectPoint &P3);
int mixedAddition(affinePointLD P1, projectPointLD P2, projectPointLD &P3, int a=0);
int mixedAddition(affinePointu4 P1, projectPointu4 P2, projectPointu4 &P3, int a=0);

int mixedAddition(projectPoint P1, projectPoint P2,  projectPoint &P3);
int mixedAddition(projectPointLD P1, projectPointLD P2, projectPointLD &P3,int a=0);
int mixedAddition(projectPointu4 P1, projectPointu4 P2, projectPointu4 &P3, int a = 0);

int Addition(projectPoint P1, projectPoint P2,  projectPoint &P3);
int Addition(projectPointLD P1, projectPointLD P2, projectPointLD &P3, int a = 0);
int Addition(projectPointu4 P1, projectPointu4 P2, projectPointu4 &P3, int a = 0);


int PrecomputationNothing(affinePoint P0,  projectPoint *P, int w=4);
int PrecomputationNothing(affinePointu4 P0, projectPointu4 *P, int w = 4, int a = 0);
int Solinas(affinePoint P0,  projectPoint *P, int w=4);
int Solinas(affinePointu4 P0, projectPointu4 *P, long long &CPUcycles, int w = 4, int a = 0);
int Solinas(affinePointLD P0, projectPointLD *P, long long &CPUcycles, int w = 4, int a = 0);

int HMV(affinePoint P0,  projectPoint *P, int w=4);
int HMV(affinePointu4 P0, projectPointu4 *P, long long &CPUcycles, int w=4, int a = 0);
int HMV(affinePointLD P0, projectPointLD *P, long long &CPUcycles, int w = 4, int a = 0);

int xuPreComputation(affinePoint P0,  projectPoint *P, int w=4, int a=0);
int xuPreComputation(affinePointu4 P0, projectPointu4 *P, long long &CPUcycles, int w=4, int a = 0);
int xuPreComputation(affinePointLD P0, projectPointLD *P, long long &CPUcycles, int w = 4, int a = 0);

int OurPreComputation(affinePoint P0,  projectPoint *P, int w=4);
int OurPreComputation(affinePointu4 P0, projectPointu4 *P, long long &CPUcycles, int w=4, int a = 0);
int OurPreComputation(affinePointLD P0, projectPointLD *P, long long &CPUcycles, int w = 4, int a = 0);

int initialVauleC(numberC c[], int cases = 0, int w = 4, int a = 0);


int windowTauNAF(Big n1, Big n2, numberC c[], int *a, int &length, int w = 4, int parameterA = 0);
int windowNAF(Big n, int *a, int &length, int w=4);
int windowNAFregular(Big n, int *a, int &length, int w = 4);
int windowTauNAFregular(Big n1, Big n2, numberC c[], int *a, int &length, int w = 4, int parameterA = 0);
int MontgomeryTrick(projectPoint *P, int w = 4);
int MontgomeryTrick(projectPointLD *P, int w = 4);
int scalarMultiplication(affinePoint P0,affinePoint negative,  projectPoint *P, projectPoint negativePi, int *a, int aLength, projectPoint Q, int w=4, int parameterA = 0);
int scalarMultiplication(affinePointu4 P0, affinePointu4 negative, projectPointu4 *P, projectPointu4 negativePi, int *a, int aLength, projectPointu4 Q, int w = 4, int parameterA=0);
int scalarMultiplication(affinePointLD P0, affinePointLD negative, projectPointLD *P, projectPointLD negativePi, int *a, int aLength, projectPointLD Q, int w = 4, int parameterA = 0);
int scalarMultiplicationMontgomery(affinePoint P0,affinePoint negative,  projectPoint *P, projectPoint negativePi, int *a, int aLength, projectPoint Q, int w=4, int parameterA = 0);
int scalarMultiplicationMontgomery(affinePointLD P0, affinePointLD negative, projectPointLD *P, projectPointLD negativePi, int *a, int aLength, projectPointLD Q, int w = 4, int parameterA = 0);

Big rand(int n);