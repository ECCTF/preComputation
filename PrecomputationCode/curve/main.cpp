

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>

#include "poly.h"
#include<ctime>
#include "preComputation.h"
//#include "ec2.h"
//#include "ecn.h"
//#include "complex.h"
//#include "flpoly.h"
using namespace std;

miracl *mip;



int main()
{//The code of pre-computation and scalar multiplication is in "preComputation.h" and "preComputation.cpp"

	mip = mirsys(1024, 2);
	mip->IOBASE = 16;

	cout <<"int size:" << sizeof(mr_utype) << endl;
	long long CPUcycles;
	//y^2+xy=x^3+x^2+1
	Big x163("2fe13c0537bbc11acaa07d793de4e6d5e5c94eee8");
	Big y163("289070fb05d38ff58321f2e800536d538ccdaa3d9");

	Big xK1283("503213f78ca44883f1a3b8162f188e553cd265f23c1567a16876913b0c2ac2458492835");
	Big yK1283("5F0010016D6F6C2AEC9007618E60B28898C5ACD1845F077E4A278B01DE6F8B3EFC48EA4");

	//y^2+xy=x^3+1
	Big x233("17232ba853a7e731af129f22ff4149563a419c26bf50a4c9d6eefad6126");
	Big y233("1db537dece819b7f70f555a67c427a8cd9bf18aeb9b56e0c11056fae6a3");


	Big x283("503213f78ca44883f1a3b8162f188e553cd265f23c1567a16876913b0c2ac2458492836");
	Big y283("1ccda380f1c9e318d90f95d07e5426fe87e45c0e8184698e45962364e34116177dd2259");



	Big orderK1283("400000000000000000000000000000000002CA3A25F1511B3440100D775C3F3C3D3873F");



	Big x409("60f05f658f49c1ad3ab1890f7184210efd0987e307c84c27accfb8f9f67cc2c460189eb5aaaa62ee222eb1b35540cfe9023746");
	Big y409("1e369050b7c4e42acba1dacbf04299c3460782f918ea427e6325165e9ea10e3da5f6c42e9c55215aa9ca27a5863ec48d8e0286b");

	Big x571("26eb7a859923fbc82189631f8103fe4ac9ca2970012d5d46024804801841ca44370958493b205e647da304db4ceb08cbbd1ba39494776fb988b47174dca88c7e2945283a01c8972");
	Big y571("349dc807f4fbf374f4aeade3bca95314dd58cec9f307a54ffc61efc006d8a2c9d4979c0ac44aea74fbebbb9f772aedcb620b01a7ba7af1b320430c8591984f601cd4c143ef1c7a3");

	mip->IOBASE = 10;


	Big divisor("4374281369");
	Big Divisor1("185622216436668784049551");
	Big Divisor2("172585439722204322826079");



	mip->IOBASE = 16;


	bool check = 0;

	int secureParameter = 283;
	//if(prepare_basis(163,7,6,3,check))
		//if(prepare_basis(233,74,0,0,check))
	if (prepare_basis(283, 12, 7, 5, check))
		//if(prepare_basis(409,87,0,0,check))
		//if(prepare_basis(571,10,5,2,check)) 
	{
		//cout << "Basis is Ok!" << endl;
	}
	//Big x1283 = 0;



	Big basis = 0;
	basis = mip->modulus;
	//cout <<"basis:"<< basis <<endl;

	int parameterA = 0;//curve parameter,K1 283
	Big one = 1;

	numberC c[65];
	c[1].g = 1;
	c[1].h = 0;

	projectPointu4 P[65];//P[i]=iP for i\in I_w
	big X[65];
	big Y[65];
	big Z[65];
	big T[65];
	Big initialX, initialY;
	//srand(time(NULL));
	for (int i = 0; i < 65; i++)
	{
		X[i] = mirvar(0);
		P[i].X0 = X[i];
		Y[i] = mirvar(0);
		P[i].X1 = Y[i];
		Z[i] = mirvar(1);
		P[i].X2 = Z[i];
		T[i] = mirvar(1);
		P[i].X3 = T[i];
	}

	projectPoint PL[33];//P[i]=iP for i\in I_w
	big XL[33];
	big YL[33];
	big ZL[33];
	Big initialXL, initialYL;
	//srand(time(NULL));
	for (int i = 0; i<33; i++)
	{
		XL[i] = mirvar(0);
		PL[i].X = XL[i];
		YL[i] = mirvar(0);
		PL[i].Y = YL[i];
		ZL[i] = mirvar(1);
		PL[i].Z = ZL[i];
	}

	projectPointLD PLD[33];//P[i]=iP for i\in I_w
	big XLD[33];
	big YLD[33];
	big ZLD[33];
	Big initialXLD, initialYLD;
	//srand(time(NULL));
	for (int i = 0; i<33; i++)
	{
		XLD[i] = mirvar(0);
		PLD[i].X = XLD[i];
		YLD[i] = mirvar(0);
		PLD[i].Y = YLD[i];
		ZLD[i] = mirvar(1);
		PLD[i].Z = ZLD[i];
	}

	affinePoint PLambda;
	PLambda.x = mirvar(0);
	PLambda.y = mirvar(0);

	if (parameterA == 0)
	{
		copy(x571.getbig(), PLambda.x);//initial P0
		copy(y571.getbig(), PLambda.y);
	}
	else 
	{
		copy(xK1283.getbig(), PLambda.x);//initial P0
		copy(yK1283.getbig(), PLambda.y);
	}
	//change coordinates to \mu4
	affinePointu4 P0;
	P0.X0 = mirvar(0);
	P0.X1 = mirvar(0);
	//P0.X2 = mirvar(1);
	P0.X3 = mirvar(0);
	modsquare2(PLambda.x, P0.X0);
	add2(PLambda.y, P0.X0, P0.X1);
	add2(PLambda.x, P0.X1, P0.X3);
	
	copy(P0.X0, P[1].X0);
	copy(P0.X1, P[1].X1);
	copy(P0.X3, P[1].X3);
	//change coordinates to lambda coordinates.
	affinePoint PL0;
	PL0.x = mirvar(0);
	PL0.y = mirvar(0);
	inverse2(PLambda.x, PL0.x);
	modmult2(PL0.x, PLambda.y, PL0.y);
	copy(PLambda.x, PL0.x);

	affinePointLD PLD0;
	PLD0.x = mirvar(0);
	PLD0.y = mirvar(0);	 
	copy(PLambda.x, PLD0.x);
	copy(PLambda.y, PLD0.y);

	int aLength = 0;
	long long	start, end;


	int w = 6;
	int PreNothing = 0;
	const int CYCLES = 100;



	projectPointu4 Q;
	Q.X0 = mirvar(0);
	Q.X1 = mirvar(0);
	Q.X2 = mirvar(1);
	Q.X3 = mirvar(0);
	projectPoint QL;
	QL.X = mirvar(0);
	QL.Y = mirvar(0);
	QL.Z = mirvar(1);

	projectPointLD QLD;
	QLD.X = mirvar(0);
	QLD.Y = mirvar(0);
	QLD.Z = mirvar(1);

	Big  n = rand(283);
	int b[600];




	affinePointu4 negative;//negative=-P0
	negative.X0 = mirvar(0);
	negative.X1 = mirvar(0);
	negative.X3 = mirvar(0);
	copy(P0.X0, negative.X0);
	copy(P0.X1, negative.X3);
	copy(P0.X3, negative.X1);

	affinePoint negativeL;//negativeL=-PL0
	negativeL.x = mirvar(0);
	negativeL.y= mirvar(0);	
	copy(PL0.x, negativeL.x);
	add2(PL0.y, one.getbig(), negativeL.y);

	affinePointLD negativeLD;//negativeLD=-PLD0
	negativeL.x = mirvar(0);
	negativeL.y = mirvar(0);
	copy(PL0.x, negativeL.x);
	add2(PL0.y, PL0.x, negativeL.y);
	
	projectPointu4 negativePi;
	negativePi.X0 = mirvar(0);
	negativePi.X1 = mirvar(0);
	negativePi.X2 = mirvar(1);
	negativePi.X3 = mirvar(0);

	projectPoint negativePiL;
	negativePiL.X = mirvar(0);
	negativePiL.Y = mirvar(0);
	negativePiL.Z = mirvar(1);


	projectPointLD negativePiLD;
	negativePiLD.X = mirvar(0);
	negativePiLD.Y = mirvar(0);
	negativePiLD.Z = mirvar(1);

	Big n1 = 0;
	Big n2 = 0;

	double SolinasPreCost = 0;
	double HMVPreCost = 0;
	double xuPreCost = 0;
	double ourPreCost = 0;

	double SolinasCost = 0;
	double HMVCost = 0;
	double xuCost = 0;
	double xuOriginalLambdaCost = 0;
	double xuOriginalLDCost = 0;
	double ourCost = 0;
	int divide = 3300* CYCLES;
	xuPreComputation(PL0, PL, w, parameterA);

	for (int i = 0; i < CYCLES; i++)
	{
		if (w < 6)
		{//Solinas' precomputation only works when w<6
			Solinas(P0, P, CPUcycles, w, parameterA);
			SolinasPreCost += CPUcycles;
		}


		HMV(P0, P, CPUcycles, w, parameterA);
		HMVPreCost += CPUcycles;
		
		xuPreComputation(P0, P, CPUcycles, w, parameterA);
		xuPreCost += CPUcycles;
		
		OurPreComputation(P0, P, CPUcycles, w, parameterA);
		ourPreCost += CPUcycles;

	}
	if (w < 6)
	{//Solinas' precomputation only works when w<6		 
		cout << "The average time of Solinas' precomputation using mu4-coordinates  in microseconds on K-283 with w=4:" 
			<< SolinasPreCost / divide << endl;
	}
	//cout << "The average time of Hankerson, Menezes, and Vanstone's precomputation using mu4-coordinates  in microseconds   on K-283 with w=6:" 
	//	<< HMVPreCost / divide << endl;
	cout << "The average time of Trost and xu's precomputation using mu4-coordinates  in microseconds   on K-283 with w=6: " 
		<< xuPreCost / divide << endl;
	cout << "The average time of our precomputation using mu4-coordinates  in microseconds   on K-283 with w=6:" 
		 << ourPreCost/divide << endl;
 



	
	
	for (int i = 0; i < CYCLES; i++)
	{
		n1 = rand(secureParameter / 2 - 1), n2 = rand(secureParameter / 2 - 1);

		if (n1 % 2 == 0)
		{
			n1 = n1 + 1;
		}
		if (n2 % 2 == 0)
		{
			n2 = n2 + 1;
		}
		//n1 = 5; n2 = -1;
		mip->IOBASE = 10;
		//cout << "Input number:" << n1 << "+" << n2 << "\\tau" << endl;
		w = 6;
		initialVauleC(c, 2, w, parameterA);
		windowTauNAFregular(n1, n2, c, b, aLength, w, parameterA);
		start = cpucycles();
		HMV(P0, P, CPUcycles, w, parameterA);
		scalarMultiplication(P0, negative, P, negativePi, b, aLength, Q, w, parameterA);
		end = cpucycles();
		HMVCost += end - start;

		initialVauleC(c, 1, w, parameterA);
		windowTauNAFregular(n1, n2, c, b, aLength, w, parameterA);
		start = cpucycles();		 
		xuPreComputation(P0, P, CPUcycles, w, parameterA);
		scalarMultiplication(P0, negative, P, negativePi, b, aLength, Q, w, parameterA);
 		end = cpucycles();
		xuCost += end - start;
		
 
		start = cpucycles();
		xuPreComputation(PL0, PL, w, parameterA);
		scalarMultiplicationMontgomery(PL0, negativeL, PL, negativePiL, b, aLength, QL, w, parameterA);
		end = cpucycles();
		xuOriginalLambdaCost += end - start;

		start = cpucycles();
		xuPreComputation(PLD0, PLD, CPUcycles, w, parameterA);
		scalarMultiplicationMontgomery(PLD0, negativeLD, PLD, negativePiLD, b, aLength, QLD, w, parameterA);
		end = cpucycles();
		xuOriginalLDCost += end - start;

		w = 7;
		initialVauleC(c, 0, w, parameterA);
		windowTauNAFregular(n1, n2, c, b, aLength, w, parameterA);
		start = cpucycles();
		OurPreComputation(P0, P, CPUcycles, w, parameterA);
		scalarMultiplication(P0, negative, P, negativePi, b, aLength, Q, w, parameterA);
		end = cpucycles();
		ourCost += end - start;
	}
	//cout << "The time of scalar multiplication in microseconds using HMV's precomputation on K-283 in mu4-coordinates with w=6: " 
	//	<< HMVCost/divide << endl;	
	cout << "The time of scalar multiplication in microseconds using Trost and xu's precomputation in  lambda-coordinates on K-283 with w=6: "
		<< xuOriginalLambdaCost / divide << endl;
	cout << "The time of scalar multiplication in microseconds using Trost and xu's precomputation in  LD-coordinates on K-283 with w=6: " 
		<< xuOriginalLDCost / divide << endl;
	cout << "The time of scalar multiplication in microseconds using Trost and xu's precomputation in mu4-coordinates on K-283 with w=6: " 
		<< xuCost / divide << endl;
	cout << "The time of scalar multiplication in microseconds using our precomputation in mu4-coordinates on K-283 with w=7: " 
		<< ourCost / divide << endl;
	return 0;
}


