#include "poly.h"
#include "preComputation.h"
#include <cmath>

using namespace std;
//Solinas
//HMV
//xuPreComputation
///OurPreComputation

long long cpucycles(void)
{
	return __rdtsc();
}


int tauAffine(big x, big lambda, big &x1, big &lambda1)
{//x1=x^2, lambda1= lambda^2
	modsquare2(x, x1);
	modsquare2(lambda, lambda1);
	return 0;
}

int tauAffine(affinePoint P, affinePoint &P1)
{//P1=\tau P
	tauAffine(P.x, P.y, P1.x, P1.y);
	return 0;
}


int tau(affinePoint P, projectPoint &P1)
{//P1=\tau P in lambda-projective coordinate and the Z-component of P is 1
	big one = mirvar(1);
	tauAffine(P.x, P.y, P1.X, P1.Y);
	copy(one, P1.Z);
	return 0;
}

int tauProjective(big x, big lambda, big z, big &x1, big &lambda1, big &z1)
{//x1=x^2, lambda1= lambda^2, z1=z^2
	modsquare2(x, x1);
	modsquare2(lambda, lambda1);
	modsquare2(z, z1);
	return 0;
}

int tauProjective(projectPoint P, projectPoint &P1)
{//P1=\tau P
	tauProjective(P.X, P.Y, P.Z, P1.X, P1.Y, P1.Z);
	return 0;
}

int tau(projectPoint P, projectPoint &P1)
{//P1=\tau P
	tauProjective(P.X, P.Y, P.Z, P1.X, P1.Y, P1.Z);
	return 0;
}



int tau(projectPointLD P, projectPointLD &P1)
{//P1=\tau P
	tauProjective(P.X, P.Y, P.Z, P1.X, P1.Y, P1.Z);
	return 0;
}
int tau(affinePointLD P, projectPointLD &P1)
{//P1=\tau P in lambda-projective coordinate and the Z-component of P is 1
	big one = mirvar(1);
	tauAffine(P.x, P.y, P1.X, P1.Y);
	copy(one, P1.Z);
	return 0;
}

int tau(affinePointu4 P, projectPointu4 &P1)
{//P1=\tau P 
	big one = mirvar(1);
	modsquare2(P.X0, P1.X0);
	modsquare2(P.X1, P1.X1);
	modsquare2(P.X3, P1.X3);
	copy(one, P1.X2);
	return 0;
}

int tau(projectPointu4 P, projectPointu4 &P1)
{//P1=\tau P
	modsquare2(P.X0, P1.X0);
	modsquare2(P.X1, P1.X1);
	modsquare2(P.X2, P1.X2);
	modsquare2(P.X3, P1.X3);
	return 0;
}


int mixedAddition(big x1, big lambda1, big x2, big lambda2, big z2, big &x3, big &lambda3, big &z3)
{//(x1,lambda1)+(x2,lambda2,z2)=(x3,lambda3,z3)
	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0), A = mirvar(0), B = mirvar(0), one = mirvar(1);
	modmult2(lambda1, z2, R1);//R1=\lamba_PZ_Q

	add2(R1, lambda2, A);//A=\lamba_PZ_Q+lambda_Q

	modmult2(x1, z2, R1);//R1=x_PZ_Q
	add2(R1, x2, R2);//R2=(x_PZ_Q+X_Q)

	modsquare2(R2, B);//B=(x_PZ_Q+X_Q)^2


	modmult2(x2, A, R1);//R1=Ax_Q
	modmult2(z2, A, R2);//R2=AZ_Q
	modmult2(R2, R1, R3);//R3=Ax_QAZ_Q

	modmult2(R3, x1, x3);//x3=x_pAx_QAZ_Q
	modmult2(R2, B, z3);//z3=BAZ_Q


	add2(R1, B, R2);//R2= Ax_Q+B
	modsquare2(R2, R1);//R1=(AX_Q+B)^2
	add2(lambda1, one, R3);//R3=\lambda_p+1
	modmult2(z3, R3, R2);//R2=z_3(\lambda_p+1)

	add2(R1, R2, R3);
	copy(R3, lambda3);
	//add2(R1,R2, lambda3);//lambda3=(AX_Q+B)^2+z_3(\lambda_p+1)
	return 0;
}

int mixedAddition(affinePoint P1, projectPoint P2, projectPoint &P3)
{//P1+P2=P3

	mixedAddition(P1.x, P1.y, P2.X, P2.Y, P2.Z, P3.X, P3.Y, P3.Z);
	return 0;
}

int mixedAddition(affinePointLD P1, projectPointLD P2, projectPointLD &P3, int a)
{//P1+P2=P3

	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0),
		A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0);


	modsquare2(P2.Z, R3);//R3=(Z_Q)^2
	modmult2(R3, P1.y, R2);//R2=(Z_Q)^2Y_P
	add2(R2, P2.Y, A);//A=(Z_Q)^2Y_P+Y_Q

	modmult2(P2.Z, P1.x, R1);//R1=(Z_Q)X_P
	add2(R1, P2.X, B);//B=(Z_Q)X_P+X_Q

	modmult2(P2.Z, B, C);//C=(Z_Q)B

	modsquare2(C, P3.Z);//Z_{P+Q}=(C)^2
	modmult2(P3.Z, P1.x, D);//D=Z_{P+Q}X_{P}

	add2(P1.x, P1.y, E);//E=X_P+Y_P


	modsquare2(B, R2);//R2=(B)^2
	add2(R2, A, R1);//R1=A+B^2
	if (a == 1)
	{
		add2(R1, C, R1);//R1=A+B^2+aC
	}
	modmult2(R1, C, R2);//R2=C(A+B^2+aC)
	modsquare2(A, R1);//R1=(A)^2

	add2(R1, R2, P3.X);//X_{P+Q}=(A)^2+C(A+B^2+aC)

	modmult2(A, C, R1);//R1=AC
	add2(R1, P3.Z, R2);//R2=Z_{P+Q}+AC
	add2(D, P3.X, R3);//R3=D+X_{P+Q}
	modmult2(R2, R3, R1);//R1=(D+X_{P+Q})(Z_{P+Q}+AC)

	modsquare2(P3.Z, R3);//R3=(Z_{P+Q})^2	
	modmult2(R3, E, R2);//R2=E(Z_{P+Q})^2	

	add2(R1, R2, P3.Y);//Y_{P+Q}=(D+X_{P+Q})(Z_{P+Q}+AC)+E(Z_{P+Q})^2	
	return 0;
}


int mixedAddition(projectPointLD P1, projectPointLD P2, projectPointLD &P3, int a)
{//P1+P2=P3

	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0),
		A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), E = mirvar(0);

	modsquare2(P2.Z, R3);//R3=(Z_Q)^2
	modmult2(R3, P1.Y, R2);//R2=(Z_Q)^2Y_P
	add2(R2, P2.Y, A);//A=(Z_Q)^2Y_P+Y_Q

	modmult2(P2.Z, P1.X, R1);//R1=(Z_Q)X_P
	add2(R1, P2.X, B);//B=(Z_Q)X_P+X_Q

	modmult2(P2.Z, B, C);//C=(Z_Q)B


	modsquare2(B, R2);//R2=(B)^2
	if (a == 1)
	{
		add2(R3, C, R1);//R1=C+(Z_Q)^2
		modmult2(R1, R2, D);//D=B^2(C+a(Z_Q)^2)
	}
	else
	{
		modmult2(C, R2, D);//D=B^2(C+a(Z_Q)^2)
	}

	modmult2(A, C, E);//E=AC

	modsquare2(C, P3.Z);//Z_{P+Q}=(C)^2

	modsquare2(A, R1);//R1=(A)^2
	add2(R1, D, R2);//R2=(A)^2+D
	add2(R2, E, P3.X);//X_{P+Q}=(A)^2+D+E

	modmult2(P1.X, P3.Z, R1);//R1=X_P(Z_{P+Q})
	add2(R1, P3.X, R2);//R2=X_{P+Q}+X_P(Z_{P+Q})
	modmult2(R2, E, R3);//R3=E(X_{P+Q}+X_P(Z_{P+Q}))

	modmult2(P1.Y, P3.Z, R1);//R1=Y_P(Z_{P+Q})
	add2(R1, P3.X, R2);//R2=X_{P+Q}+Y_P(Z_{P+Q})
	modmult2(R2, P3.Z, R1);//R1=Z_{P+Q}(X_{P+Q}+Y_P(Z_{P+Q}))

	add2(R1, R3, P3.Y);//Y_{P+Q}=E(X_{P+Q}+X_P(Z_{P+Q}))+Z_{P+Q}(X_{P+Q}+Y_P(Z_{P+Q}))
	return 0;
}

int mixedAddition(affinePointu4 P1, projectPointu4 P2, projectPointu4 &P3, int a)
{//P1+P2=P3

	big U00 = mirvar(0), U11 = mirvar(0), U22 = mirvar(0), U33 = mirvar(0), R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0), H = mirvar(0);
	modmult2(P1.X0, P2.X0, U00);
	modmult2(P1.X1, P2.X1, U11);
	copy(P2.X2, U22);
	modmult2(P1.X3, P2.X3, U33);

	add2(U00, U22, R1); //R1 = U00 + U22;
	modsquare2(R1, P3.X0);//P3.X0

	add2(U11, U33, R2); //R2 = U11 + U33;
	modsquare2(R2, P3.X2);//P3.X2

	modmult2(R1, R2, R3);//R3=(U00 + U22)(U11 + U33)




	if (a == 1)
	{
		add2(P1.X1, P1.X3, R1);
		add2(P2.X1, P2.X3, R2);
		modmult2(R1, R2, H);

		add2(U11, H, U11);
		add2(U33, H, U33);
	}

	modmult2(U00, U11, R1);
	modmult2(U22, U33, R2);
	add2(R1, R2, P3.X1);
	add2(R3, P3.X1, P3.X3);

	return 0;
}

int mixedAddition(projectPoint P1, projectPoint P2, projectPoint &P3)
{//P1+P2=P3, P1.x=1

	mixedAddition(P1.X, P1.Y, P2.X, P2.Y, P2.Z, P3.X, P3.Y, P3.Z);
	return 0;
}

int mixedAddition(projectPointu4 P1, projectPointu4 P2, projectPointu4 &P3, int a)
{//P1+P2=P3

	big U00 = mirvar(0), U11 = mirvar(0), U22 = mirvar(0), U33 = mirvar(0), R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0), H = mirvar(0);
	modmult2(P1.X0, P2.X0, U00);
	modmult2(P1.X1, P2.X1, U11);
	copy(P2.X2, U22);
	modmult2(P1.X3, P2.X3, U33);

	add2(U00, U22, R1); //R1 = U00 + U22;
	modsquare2(R1, P3.X0);//P3.X0

	add2(U11, U33, R2); //R2 = U11 + U33;
	modsquare2(R2, P3.X2);//P3.X2

	modmult2(R1, R2, R3);//R3=(U00 + U22)(U11 + U33)




	if (a == 1)
	{
		add2(P1.X1, P1.X3, R1);
		add2(P2.X1, P2.X3, R2);
		modmult2(R1, R2, H);

		add2(U11, H, U11);
		add2(U33, H, U33);
	}

	modmult2(U00, U11, R1);
	modmult2(U22, U33, R2);
	add2(R1, R2, P3.X1);
	add2(R3, P3.X1, P3.X3);

	return 0;
}


int Addition(projectPointu4 P1, projectPointu4 P2, projectPointu4 &P3, int a)
{//P1+P2=P3

	big U00 = mirvar(0), U11 = mirvar(0), U22 = mirvar(0), U33 = mirvar(0), R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0), H = mirvar(0);
	modmult2(P1.X0, P2.X0, U00);
	modmult2(P1.X1, P2.X1, U11);
	modmult2(P1.X2, P2.X2, U22);
	modmult2(P1.X3, P2.X3, U33);

	add2(U00, U22, R1); //R1 = U00 + U22;
	modsquare2(R1, P3.X0);//P3.X0

	add2(U11, U33, R2); //R2 = U11 + U33;
	modsquare2(R2, P3.X2);//P3.X2

	modmult2(R1, R2, R3);//R3=(U00 + U22)(U11 + U33)




	if (a == 1)
	{
		add2(P1.X1, P1.X3, R1);
		add2(P2.X1, P2.X3, R2);
		modmult2(R1, R2, H);

		add2(U11, H, U11);
		add2(U33, H, U33);
	}

	modmult2(U00, U11, R1);
	modmult2(U22, U33, R2);
	add2(R1, R2, P3.X1);
	add2(R3, P3.X1, P3.X3);

	return 0;
}


int mixedAddition(big x1, big lambda1, big x2, big lambda2, big &x3, big &lambda3, big &z3)
{//(x1,lambda1)+(x2,lambda2,z2)=(x3,lambda3,z3)
	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0), A = mirvar(0), B = mirvar(0), one = mirvar(1);

	add2(lambda1, lambda2, A);//A=\lamba_PZ_Q+lambda_Q


	add2(x1, x2, R2);//R2=(x_PZ_Q+X_Q)

	modsquare2(R2, B);//B=(x_PZ_Q+X_Q)^2


	modmult2(x2, A, R1);//R1=Ax_Q

	modmult2(A, R1, R3);//R3=Ax_QA 

	modmult2(R3, x1, x3);//x3=x_pAx_QAZ_Q
	modmult2(A, B, z3);//z3=BA 


	add2(R1, B, R2);//R2= Ax_Q+B
	modsquare2(R2, R1);//R1=(AX_Q+B)^2
	add2(lambda1, one, R3);//R3=\lambda_p+1
	modmult2(z3, R3, R2);//R2=z_3(\lambda_p+1)

	add2(R1, R2, R3);
	copy(R3, lambda3);
	//add2(R1,R2, lambda3);//lambda3=(AX_Q+B)^2+z_3(\lambda_p+1)
	return 0;
}


int mixedAddition(affinePoint P1, affinePoint P2, projectPoint &P3)
{//P1+P2=P3

	mixedAddition(P1.x, P1.y, P2.x, P2.y, P3.X, P3.Y, P3.Z);
	return 0;
}


int Addition(projectPointLD P1, projectPointLD P2, projectPointLD &P3, int a)
{
	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0),
		A1 = mirvar(0), A2 = mirvar(0), B1 = mirvar(0), B2 = mirvar(0),
		C = mirvar(0), D = mirvar(0), E1 = mirvar(0), E2 = mirvar(0),
		F = mirvar(0), G = mirvar(0);


	modmult2(P1.X, P2.Z, A1);//A1=(X_P)Z_Q
	modmult2(P2.X, P1.Z, A2);//A2=(X_Q)Z_P
	add2(A1, A2, C);//C=A_1+A_2


	modsquare2(A1, B1);//B1=(A1)^2
	modsquare2(A2, B2);//B2=(A2)^2
	add2(B1, B2, D);//D=B1+B2

	modsquare2(P2.Z, R1);//R1=(Z_Q)^2
	modmult2(P1.Y, R1, E1);//E1=(Z_Q)^2Y_P
	modsquare2(P1.Z, R1);//R1=(Z_P)^2
	modmult2(P2.Y, R1, E2);//E2=(Z_P)^2Y_Q
	add2(E1, E2, F);//F=E1+E2

	modmult2(C, F, G);//G=CF

	modmult2(P1.Z, P2.Z, R1);//R1=Z_PZ_Q
	modmult2(R1, D, P3.Z);//Z_{P+Q}=Z_PZ_Q D

	add2(B2, E2, R1);//R1= B2+E2
	modmult2(R1, A1, R3);//R3=A1(B2+E2)
	add2(B1, E1, R1);//R1= B1+E1
	modmult2(R1, A2, R2);//R2=A2(B1+E1)

	add2(R2, R3, P3.X);//X_{P+Q}=A1(B2+E2)+A2(B1+E1)

	modmult2(A1, G, R1);//R1=A1G
	modmult2(E1, D, R2);//R2=E1D
	add2(R1, R2, R3);//R3=A1G+E1D
	modmult2(R3, D, R1);//R1=(A1G+E1D)D
	add2(G, P3.Z, R3);//R3=G+Z_{P+Q}
	modmult2(R3, P3.X, R2);//R2=X_{P+Q}(G+Z_{P+Q})

	add2(R1, R2, P3.Y);//Y_{P+Q}=(A1G+E1D)D+X_{P+Q}(G+Z_{P+Q})	
	return 0;
}


int Addition(big x1, big lambda1, big z1, big x2, big lambda2, big z2, big &x3, big &lambda3, big &z3)
{//(x1,lambda1,z1)+(x2,lambda2,z2)=(x3,lambda3,z3)
	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0), A = mirvar(0), B = mirvar(0), middle = mirvar(1);
	modmult2(lambda1, z2, R1);//R1=\lamba_PZ_Q

	modmult2(lambda2, z1, R2);//R2=\lamba_QZ_P
	add2(R1, R2, A);//A=\lamba_PZ_Q+lambda_QZ_P

	modmult2(x1, z2, R1);//R1=x_PZ_Q
	modmult2(x2, z1, R2);//R2=x_QZ_P
	add2(R1, R2, R3);//R3=(x_PZ_Q+X_QZ_P)

	modsquare2(R3, B);//B=(x_PZ_Q+X_QZ_P)^2


	modmult2(R2, A, middle);//middle=Ax_Qz_P
	modmult2(R1, A, R3);//R2=Ax_PZ_Q
	modmult2(middle, R3, x3);//x3=Ax_PZ_Q Ax_Qz_P

	modmult2(A, B, R1);//R1=AB
	modmult2(R1, z2, R3);//R3=ABZ_Q
	modmult2(R3, z1, z3);//z3=BAZ_QZ_P


	add2(middle, B, R2);//R2= Ax_QZ_P+B
	modsquare2(R2, R1);//R1=(AX_QZ_P+B)^2
	add2(lambda1, z1, A);//R3=\lambda_p+z_P
	modmult2(A, R3, R2);//R2=ABZ_Q(\lambda_p+z_P)

	add2(R1, R2, lambda3);//lambda3=(AX_Q+B)^2+z_3(\lambda_p+1)
	return 0;
}

int Addition(projectPoint P1, projectPoint P2, projectPoint &P3)
{//P1+P2=P3

	/*Addition(P1.X, P1.Y, P1.Z, P2.X, P2.Y, P2.Z, P3.X, P3.Y, P3.Z);
	return 0;*/


	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0), A = mirvar(0), B = mirvar(0), middle = mirvar(1);
	modmult2(P1.Y, P2.Z, R1);//R1=\lamba_PZ_Q

	modmult2(P2.Y, P1.Z, R2);//R2=\lamba_QZ_P
	add2(R1, R2, A);//A=\lamba_PZ_Q+lambda_QZ_P

	modmult2(P1.X, P2.Z, R1);//R1=x_PZ_Q
	modmult2(P2.X, P1.Z, R2);//R2=x_QZ_P
	add2(R1, R2, R3);//R3=(x_PZ_Q+X_QZ_P)

	modsquare2(R3, B);//B=(x_PZ_Q+X_QZ_P)^2


	modmult2(R2, A, middle);//middle=Ax_Qz_P
	modmult2(R1, A, R3);//R2=Ax_PZ_Q
	modmult2(middle, R3, P3.X);//x3=Ax_PZ_Q Ax_Qz_P

	modmult2(A, B, R1);//R1=AB
	modmult2(R1, P2.Z, R3);//R3=ABZ_Q
	modmult2(R3, P1.Z, P3.Z);//z3=BAZ_QZ_P


	add2(middle, B, R2);//R2= Ax_QZ_P+B
	modsquare2(R2, R1);//R1=(AX_QZ_P+B)^2
	add2(P1.Y, P1.Z, A);//R3=\lambda_p+z_P
	modmult2(A, R3, R2);//R2=ABZ_Q(\lambda_p+z_P)

	add2(R1, R2, P3.Y);//lambda3=(AX_Q+B)^2+z_3(\lambda_p+1)
	return 0;
}



int PpmQ(big x1, big lambda1, big x2, big lambda2, big z2, big &x3, big &lambda3, big &z3, big &x4, big &lambda4, big &z4)
{ //(x1,lambda1)+(x2,lambda2,z2)=(x3,lambda3,z3),(x1,lambda1)-(x2,lambda2,z2)=(x4,lambda4,z4)
	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0), A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0), one = mirvar(1);
	modmult2(lambda1, z2, R1);//R1=\lamba_PZ_Q

	add2(R1, lambda2, A);//A=\lamba_PZ_Q+lambda_Q

	modmult2(x1, z2, R1);//R1=x_PZ_Q
	add2(R1, x2, R2);//R2=(x_PZ_Q+X_Q)

	modsquare2(R2, B);//B=(x_PZ_Q+X_Q)^2


	modmult2(x2, z2, C);//C=X_QZ_Q
	modmult2(x1, C, D);//D=C x_P

	modsquare2(A, R1);//R1=A^2
	modmult2(R1, D, x3);//x3=A^2D
	modmult2(R2, B, z3);//z3=BAZ_Q


	add2(R1, B, R2);//R2= Ax_Q+B
	modsquare2(R2, R1);//R1=(AX_Q+B)^2
	add2(lambda1, one, R3);//R3=\lambda_p+1
	modmult2(z3, R3, R2);//R2=z_3(\lambda_p+1)

	add2(R1, R2, lambda3);//lambda3=(AX_Q+B)^2+z_3(\lambda_p+1)


	modsquare2(z2, R1);//R1=Z_Q^2


	modmult2(R1, D, A);//A=DZ_Q^2
	add2(x3, A, x4);//x4=x_{P+Q}+DZ_Q^2

	modmult2(R1, B, R2);//R2=BZ_Q^2
	add2(z3, R2, z4);//z4=BZ_Q^2+z_{p+Q}


	modmult2(R2, R3, one);//one=BZ_Q^2(\lambda_p+1)
	modsquare2(C, A);//A=C^2

	add2(A, one, R3);//R3=C^2+BZ_Q^2(\lambda_p+1)
	add2(R3, lambda3, lambda4);//lambda4=\lambda_{P+Q}+C^2+BZ_Q^2(\lambda_p+1)


	return 0;
}



int PpmQ(affinePoint P1, projectPoint P2, projectPoint &P3, projectPoint &P4)
{//P1+P2=P3,P1-P2=P4

	PpmQ(P1.x, P1.y, P2.X, P2.Y, P2.Z, P3.X, P3.Y, P3.Z, P4.X, P4.Y, P4.Z);
	return 0;
}

int PpmQ(projectPointu4 P1, projectPointu4 P2, projectPointu4 &P3, projectPointu4 &P4, int a = 0)
{//P1+P2=P3,P1-P2=P4, P1 affine point

	big U00 = mirvar(0), U11 = mirvar(0), U22 = mirvar(0), U33 = mirvar(0), U20 = mirvar(0), U02 = mirvar(0), R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0), R4 = mirvar(0), R5 = mirvar(0), H = mirvar(0);
	modmult2(P1.X0, P2.X0, U00);
	modmult2(P1.X1, P2.X1, U11);
	copy(P2.X2, U22);
	modmult2(P1.X3, P2.X3, U33);
	copy(P2.X0, U20);
	modmult2(P1.X0, P2.X2, U02);

	add2(U00, U22, R1); //R1 = U00 + U22;
	modsquare2(R1, P3.X0);//P3.X0

	add2(U11, U33, R5); //R5 = U11 + U33;
	modsquare2(R5, P3.X2);//P3.X2

	modmult2(R1, R2, R3);//R3=(U00 + U22)(U11 + U33)


	add2(U02, U20, R4); //R4 = U02 + U20;

	if (a == 1)
	{
		add2(P1.X1, P1.X3, R1);
		add2(P2.X1, P2.X3, R2);
		modmult2(R1, R2, H);

		add2(U11, H, U11);
		add2(U33, H, U33);
	}

	modmult2(U00, U11, R1);
	modmult2(U22, U33, R2);
	add2(R1, R2, P3.X1);
	add2(R3, P3.X1, P3.X3);

	copy(P3.X2, P4.X0); //P4.X0
	modsquare2(R4, P4.X2);//P4.X2

	modmult2(U02, U33, R1);
	modmult2(U20, U11, R2);
	add2(R1, R2, P4.X1);//P4.X1

	modmult2(R4, R5, R3);//R3=(U02 + U20)(U11 + U33)
	add2(R3, P3.X1, P4.X3);////P4.X3

	return 0;

}



int PpmQ(projectPointLD P1, projectPointLD P2, projectPointLD &P3, projectPointLD &P4, int a = 0)
{//P1+P2=P3,P1-P2=P4, P1 affine point

	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0), R4 = mirvar(0),
		A = mirvar(0), B = mirvar(0), C = mirvar(0), D = mirvar(0),
		E = mirvar(0), F = mirvar(0), G = mirvar(0);


	modsquare2(P2.Z, R3);//R3=(Z_Q)^2
	modmult2(R3, P1.Y, R2);//R2=(Z_Q)^2Y_P
	add2(R2, P2.Y, A);//A=(Z_Q)^2Y_P+Y_Q

	modmult2(P2.Z, P1.X, R1);//R1=(Z_Q)X_P
	add2(R1, P2.X, B);//B=(Z_Q)X_P+X_Q

	modmult2(P2.Z, B, C);//C=(Z_Q)B

	modsquare2(C, P3.Z);//Z_{P+Q}=(C)^2
	modmult2(P3.Z, P1.X, D);//D=Z_{P+Q}X_{P}

	add2(P1.X, P1.Y, E);//E=X_P+Y_P


	modsquare2(B, R2);//R2=(B)^2
	add2(R2, A, R1);//R1=A+B^2
	if (a == 1)
	{
		add2(R1, C, R1);//R1=A+B^2+aC
	}
	modmult2(R1, C, R2);//R2=C(A+B^2+aC)
	modsquare2(A, R1);//R1=(A)^2

	add2(R1, R2, P3.X);//X_{P+Q}=(A)^2+C(A+B^2+aC)

	modmult2(A, C, R1);//R1=AC
	add2(R1, P3.Z, R4);//R4=Z_{P+Q}+AC
	add2(D, P3.X, R3);//R3=D+X_{P+Q}
	modmult2(R4, R3, R1);//R1=(D+X_{P+Q})(Z_{P+Q}+AC)

	modsquare2(P3.Z, R3);//R3=(Z_{P+Q})^2	
	modmult2(R3, E, R2);//R2=E(Z_{P+Q})^2	

	add2(R1, R2, P3.Y);//Y_{P+Q}=(D+X_{P+Q})(Z_{P+Q}+AC)+E(Z_{P+Q})^2	

	copy(P3.Z, P4.Z);//Z_{P-Q}=Z_{P+Q}
	modmult2(P2.X, P2.Z, R1);//R1=X_QZ_Q
	modmult2(C, R1, F);//F=X_QZ_QC

	modsquare2(R1, R3);//R3=(X_QZ_Q)^2	
	add2(R3, F, G);//G=(X_QZ_Q)^2+F
	add2(P3.X, G, P4.X);//X_{P-Q}=X_{P+Q}+G

	add2(R4, P3.Z, R3);//G=AC+F+Z_{P+Q}
	modmult2(R3, G, R1);//R1=G(AC+F+Z_{P+Q})

	add2(D, P3.X, R3);//R3=D+X_{P+Q}
	modmult2(R3, F, R2);//R2=F(D+X_{P+Q})

	add2(R1, R2, R3);//R3=G(AC+F+Z_{P+Q})+F(D+X_{P+Q})
	add2(P3.Y, R3, P4.Y);//Y_{P-Q}=Y_{P+Q}+G(AC+F+Z_{P+Q})+F(D+X_{P+Q})

	return 0;

}


int barTau(projectPointu4 P, projectPointu4 &TP, int a = 0)
{//\mu\bar\tau P= TP
	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0);
	add2(P.X0, P.X2, R1);//R1=X0+X2
	add2(P.X1, P.X3, R2);//R2=X1+X3

	modsquare2(R1, TP.X0);//TP.X0
	modsquare2(R2, TP.X2);//TP.X2

	modmult2(R1, R2, R3);//R3=(X0 + X2)(X1 + X3)

	add2(P.X0, P.X1, R1);//R1=X0+X1
	add2(P.X2, P.X3, R2);//R2=X2+X3

	modmult2(R1, R2, R1);//R1=(X0 + X2)(X1 + X3)

	if (a == 0)
	{
		add2(R1, TP.X2, R1);
		//add2(R3, TP.X2, R3);
	}
	add2(R1, TP.X0, TP.X1);
	add2(R3, TP.X1, TP.X3);
	return 0;
}


int barTau(projectPointLD P, projectPointLD &TP, int a = 0)
{//\mu\bar\tau P= TP
	big R1 = mirvar(0), R2 = mirvar(0), R3 = mirvar(0);
	modmult2(P.X, P.Z, TP.Z);//Z_{TP}=X_P Z_P

	add2(P.X, P.Z, R1);//R1=X_P+Z_P
	modsquare2(R1, TP.X);//X_{TP}=(X_P+Z_P)^2

	if (a == 0)
	{
		add2(P.Y, TP.X, R1);//R1=Y_P+X_TP
		add2(P.Y, TP.Z, R2);//R2=Y_P+Z_TP
		modmult2(R1, R2, R3);//R3=(Y_P+X_TP)(Y_P+Z_TP)

		modsquare2(TP.Z, R1);//R1=(Z_{TP})^2

		add2(R1, R3, TP.Y);//Y_TP=(Y_P+X_TP)(Y_P+Z_TP)+(Z_{TP})^2
	}
	else if (a == 1)
	{
		add2(P.Y, TP.X, R1);//R1=Y_P+X_TP
		add2(R1, TP.Z, R2);//R2=Y_P+X_TP+Z_TP
		modmult2(P.Y, R2, TP.Y);//Y_TP=Y_P(Y_P+X_TP+Z_TP) 
	}
	return 0;
}

int barTau(big XP, big LambdaP, big ZP, big &XTP, big &LambdaTP, big &ZTP, big &alpha)
{//\mu\bar\tau(XP, LambdaP, ZP)=(XTP, LambdaTP, ZTP), alpha=x_PZ_P
	modmult2(XP, ZP, alpha);//alpha=x_PZ_P
	big R1 = mirvar(0), R2 = mirvar(0), A1 = mirvar(0);
	modsquare2(XP, R1);//R1=x_P^2
	modsquare2(ZP, R2);//R2=z_P^2
	add2(R1, R2, A1);//A1=x_P^2+z_P^2

	modsquare2(A1, XTP);//XTP=A_1^2
	modmult2(A1, alpha, ZTP);//ZTP=A_1\alpha

	modmult2(R1, alpha, R2);//R2=\alpha x_P^2
	modmult2(XP, A1, R1);//R1=x_PA1
	modmult2(R1, LambdaP, A1);//A1=x_P Lambda_P A1
	add2(R2, A1, LambdaTP);// LambdaTP=\alpha x_P^2+x_P Lambda_P A1




	return 0;
}





int barTau(projectPoint P, projectPoint &TP, big &alpha)
{//\mu\bar\tau P= TP, alpha=x_PZ_P
	barTau(P.X, P.Y, P.Z, TP.X, TP.Y, TP.Z, alpha);
	return 0;
}

////////////////////
int AffinePMinusTauP(big xP, big lambdaP, big &XTP, big &LambdaTP, big &ZTP)
{// (xP,lambdaP)-\tau(xP,lambdaP)=(XTP, LambdaTP, ZTP)
	big R1 = mirvar(0), R2 = mirvar(1), A1 = mirvar(0), one = mirvar(1);
	modsquare2(xP, R1);//R1=x_P^2
	modmult2(R1, xP, R2);//R2=x_P^3


	modsquare2(R1, A1);//A1=x_P^4
	add2(A1, one, XTP);//XTP=x_P^4+1

	add2(xP, R2, ZTP);//ZTP



	modmult2(ZTP, lambdaP, R1);//ZTP=A_1\alpha

	add2(R1, R2, LambdaTP);//ZTP

	return 0;
}


int AffinePMinusTauP(big xP, big lambdaP, big &XTP, big &LambdaTP, big &ZTP, big &X3)
{//(xP,lambdaP)-\tau(xP,lambdaP)=(XTP, LambdaTP, ZTP),X3=xP^3
	big R1 = mirvar(0), R2 = mirvar(1), A1 = mirvar(0), one = mirvar(1);
	modsquare2(xP, R1);//R1=x_P^2
	modmult2(R1, xP, R2);//R2=x_P^3
	copy(R2, X3);


	modsquare2(R1, A1);//A1=x_P^4
	add2(A1, one, XTP);//XTP=x_P^4+1

	add2(xP, R2, ZTP);//ZTP



	modmult2(ZTP, lambdaP, R1);//ZTP=A_1\alpha

	add2(R1, R2, LambdaTP);//ZTP

	return 0;
}

int AffinePMinusTauP(affinePoint P, projectPoint &TP)
{//P-\tauP=TP
	AffinePMinusTauP(P.x, P.y, TP.X, TP.Y, TP.Z);
	return 0;
}

int AffinePMinusTauP(affinePoint P, projectPoint &TP, big &X3)
{//P-\tauP=TP, X3= P.X^3
	AffinePMinusTauP(P.x, P.y, TP.X, TP.Y, TP.Z, X3);
	return 0;
}

////////////////////
int ProjectivePMinusTauP(big xP, big lambdaP, big zP, big &XTP, big &LambdaTP, big &ZTP)
{//(xP,lambdaP,zP)-\tau(xP,lambdaP,zP)=(XTP, LambdaTP, ZTP),X3=xP^3
	big R1 = mirvar(0), R2 = mirvar(1), A1 = mirvar(0), one = mirvar(1), R3 = mirvar(0);
	add2(xP, zP, R3);//R3=xP+zP
	modsquare2(R3, R1);//R1=x_P^2+z_P^2
	modmult2(R1, xP, R2);//R2=x_P(x_P^2+z_P^2)


	modsquare2(R1, XTP);//XTP=x_P^4+z_P^4

	modmult2(zP, R2, ZTP);//ZTP





	modmult2(ZTP, lambdaP, R1);//R1=ZTP lambdaP

	modsquare2(xP, R2);//R2=x_P^2

	modmult2(xP, R2, A1);//

	modmult2(zP, A1, R2);// 

	add2(R1, R2, LambdaTP);//LambdaTP

	return 0;
}


int ProjectivePMinusTauP(projectPoint P, projectPoint &TP)
{//P-\tauP=TP
	ProjectivePMinusTauP(P.X, P.Y, P.Z, TP.X, TP.Y, TP.Z);
	return 0;
}





int AffinePAddTauP(big xP, big lambdaP, big &XTP, big &LambdaTP, big &ZTP)
{////(xP,lambdaP)+\tau(xP,lambdaP)=(XTP, LambdaTP, ZTP),X3=xP^3
	big R1 = mirvar(0), R2 = mirvar(1), A1 = mirvar(0), one = mirvar(1), B = mirvar(0);
	modsquare2(xP, R1);//R1=x_P^2



	modsquare2(R1, R2);//R2=x_P^4
	add2(R1, one, A1);//A1=x_P^2+1
	add2(A1, R2, B);//B=x_P^4+x_P^2+1

	modsquare2(B, XTP);//XTP=B^2

	modmult2(R1, R2, one);//one=x_P^6

	modmult2(xP, A1, R1);//R1=x_P(x_P^2+1)
	modmult2(R1, B, ZTP);//R2=x_P(x_P^2+1)B
	modmult2(ZTP, lambdaP, R1);//R1=lambdaP x_P(x_P^2+1)B 
	modmult2(xP, one, R2);//R2=x_P^7
	add2(R1, R2, LambdaTP);


	return 0;
}

int AffinePAddTauP(big X3, big xP, big lambdaP, big &XTP, big &LambdaTP, big &ZTP)
{//(xP,lambdaP)+\tau(xP,lambdaP)=(XTP, LambdaTP, ZTP),X3=xP^3
	big R1 = mirvar(0), R2 = mirvar(1), A1 = mirvar(0), one = mirvar(1), B = mirvar(0);
	modsquare2(xP, R1);//R1=x_P^2



	modsquare2(R1, R2);//R2=x_P^4
	add2(R1, one, A1);//A1=x_P^2+1
	add2(A1, R2, B);//B=x_P^4+x_P^2+1

	modsquare2(B, XTP);//XTP=B^2

	modmult2(R1, R2, one);//one=x_P^6

	//modmult2(xP, A1,R1);//R1=x_P(x_P^2+1)



	add2(xP, X3, R1);//R1=x_P(x_P^2+1)


	modmult2(R1, B, ZTP);//R2=x_P(x_P^2+1)B
	modmult2(ZTP, lambdaP, R1);//R1=lambdaP x_P(x_P^2+1)B 
	modmult2(xP, one, R2);//R2=x_P^7
	add2(R1, R2, LambdaTP);


	return 0;
}



int AffinePAddTauP(big X3, affinePoint P, projectPoint &TP)
{//P+\tauP=TP
	AffinePAddTauP(X3, P.x, P.y, TP.X, TP.Y, TP.Z);
	return 0;
}

int AffinePAddTauP(affinePoint P, projectPoint &TP)
{//P+\tauP=TP
	AffinePAddTauP(P.x, P.y, TP.X, TP.Y, TP.Z);
	return 0;
}



int RebarTau(big XP, big LambdaP, big ZP, big alpha, big &XTP, big &LambdaTP, big &ZTP)
{// \mu\bar\tau (XP,LambdaP, ZP)=(XTP, LambdaTP, ZTP)
	big R1 = mirvar(0), R2 = mirvar(0), A2 = mirvar(0);
	modsquare2(alpha, R1);//R1=alpha^2
	//modsquare2(ZP, R2);//R2=z_P^2
	add2(R1, XP, A2);//A2

	modsquare2(A2, XTP);//XTP=A_2^2

	modmult2(A2, alpha, ZTP);//ZTP=A_2\alpha


	modmult2(XP, ZP, R1);// 
	modmult2(A2, LambdaP, R2);// 
	add2(R1, R2, LambdaTP);//  

	return 0;
}


int RebarTau(projectPoint P, big alpha, projectPoint &TP)
{//\bar\tau P=TP
	RebarTau(P.X, P.Y, P.Z, alpha, TP.X, TP.Y, TP.Z);
	return 0;
}


int ReRebarTau(big XP, big LambdaP, big ZP, big alpha, big &XTP, big &LambdaTP, big &ZTP)
{// \mu\bar\tau (XP,LambdaP, ZP)=(XTP, LambdaTP, ZTP)

	big R1 = mirvar(0), R2 = mirvar(0), A2 = mirvar(0);
	modsquare2(alpha, R1);//R1=alpha^2
	//modsquare2(ZP, R2);//R2=z_P^2
	add2(R1, XP, A2);//A2

	modsquare2(A2, XTP);//XTP=A_2^2

	modmult2(A2, alpha, ZTP);//ZTP=A_2\alpha


	modmult2(XP, ZP, R1);// 
	modmult2(A2, LambdaP, R2);// 
	add2(R1, R2, LambdaTP);//  

	return 0;
}



int ReRebarTau(projectPoint P, big alpha, projectPoint &TP)
{//\bar\tau P=TP
	ReRebarTau(P.X, P.Y, P.Z, alpha, TP.X, TP.Y, TP.Z);
	return 0;
}





int affineBarTau(big XP, big lambdaP, big &X1, big &Lambda1, big &Z1, big &X2, big &Lambda2, big &Z2)
{// \mu\bar\tau (XP,LambdaP)=(X1, Lambda1, Z1)  \mu\bar\tau^2 (XP,LambdaP)=(X2, Lambda2, Z2)
	big R1 = mirvar(0), R2 = mirvar(0), A2 = mirvar(0), beta = mirvar(0), one = mirvar(1);

	modsquare2(XP, beta);//beta=x_P^2
	modsquare2(beta, R1);//R1=beta^2
	add2(R1, one, X1);//A2


	modmult2(XP, beta, R2);//R2=x_p\beta
	add2(R2, XP, Z1);//z1=x_P\beta+x_P

	add2(lambdaP, one, R1);
	modmult2(R1, Z1, R2);//R2=(\lambda_p+1)Z1
	add2(R2, XP, Lambda1);//Lambda1=(\lambda_p+1)Z1+x_P

	add2(X1, beta, A2);
	modsquare2(A2, X2);//X2=A_2^2

	modmult2(A2, Z1, Z2);//Z2=A_2Z1


	modmult2(Z2, lambdaP, R2);// 
	add2(XP, R2, Lambda2);//  
	return 0;
}

int affineBarTau(affinePoint P, projectPoint &P1, projectPoint &P2)
{// \mu\bar\tau  P=P1, \mu\bar\tau^2  P=P2,
	affineBarTau(P.x, P.y, P1.X, P1.Y, P1.Z, P2.X, P2.Y, P2.Z);
	return 0;
}

int PrecomputationNothing(affinePoint P0, projectPoint *P, int w)
{//Solinas' pre-computation for w=4,5

	//big *x=(big *)mr_alloc(7,sizeof(big));
	big x[7], y[7], z[7], x0, y0;
	x0 = mirvar(0);
	y0 = mirvar(0);
	copy(P0.x, x0);
	copy(P0.y, y0);
	for (int i = 0; i < 7; i++)
	{
		x[i] = mirvar(0);
		y[i] = mirvar(0);
		z[i] = mirvar(0);
	}
	//big *y=(big *)mr_alloc(7,sizeof(big));
	//big *z=(big *)mr_alloc(7,sizeof(big));
	big middleNegative = mirvar(0);//the negative lambda of lambda coordinate 
	big middlex = mirvar(0), middley = mirvar(0), middlez = mirvar(1), midx = mirvar(0), midy = mirvar(0), midz = mirvar(1);

	big  one = mirvar(1);

	if (w == 5)
	{

	}

	return 0;

}


int Solinas(affinePoint P0, projectPoint *P, int w)
{//Solinas' pre-computation for w=4,5. Only works when a=0

	//big *x=(big *)mr_alloc(7,sizeof(big));
	big x[7], y[7], z[7], x0, y0;
	x0 = mirvar(0);
	y0 = mirvar(0);
	copy(P0.x, x0);
	copy(P0.y, y0);
	for (int i = 0; i < 7; i++)
	{
		x[i] = mirvar(0);
		y[i] = mirvar(0);
		z[i] = mirvar(0);
	}
	//big *y=(big *)mr_alloc(7,sizeof(big));
	//big *z=(big *)mr_alloc(7,sizeof(big));
	big middleNegative = mirvar(0);//the negative lambda of lambda coordinate 
	big middlex = mirvar(0), middley = mirvar(0), middlez = mirvar(1), midx = mirvar(0), midy = mirvar(0), midz = mirvar(1);

	big  one = mirvar(1);
	add2(y0, one, middleNegative);

	tauAffine(x0, y0, middlex, middley);//tau
	tauAffine(middlex, middley, midx, midy);//tau^2
	mixedAddition(x0, middleNegative, midx, midy, P[3].X, P[3].Y, P[3].Z);//Q3



	mixedAddition(x0, y0, midx, midy, P[5].X, P[5].Y, P[5].Z);//Q5

	copy(P[5].X, x[1]);
	copy(P[5].Y, y[1]);
	copy(P[5].Z, z[1]);

	tauAffine(midx, midy, middlex, middley);//tau^3
	add2(middley, one, middley);
	mixedAddition(x0, middleNegative, middlex, middley, P[7].X, P[7].Y, P[7].Z);//Q7


	if (w == 5)
	{
		tauAffine(middlex, middley, midx, midy);//tau^4
		mixedAddition(x0, middleNegative, midx, midy, P[15].X, P[15].Y, P[15].Z);//Q15
		tauProjective(x[1], y[1], z[1], middlex, middley, middlez);//tau Q_5
		tauProjective(middlex, middley, middlez, midx, midy, midz);//tau^2 Q_5
		add2(midy, midz, one);
		mixedAddition(x0, middleNegative, midx, one, midz, P[11].X, P[11].Y, P[11].Z);//Q11

		tauProjective(midx, one, midz, middlex, middley, middlez);//-tau^3 Q_5
		mixedAddition(x0, y0, middlex, middley, middlez, P[9].X, P[9].Y, P[9].Z);//Q9			
		mixedAddition(x0, y0, midx, one, midz, P[13].X, P[13].Y, P[13].Z);//Q13
	}

	return 0;

}

int negative(projectPointu4 P1, projectPointu4 &P2)
{
	copy(P1.X0, P2.X0);
	copy(P1.X1, P2.X3);
	copy(P1.X2, P2.X2);
	copy(P1.X3, P2.X1);

	return 0;
}

int negative(projectPointLD P1, projectPointLD &P2)
{
	big R = mirvar(0);
	copy(P1.X, P2.X);
	copy(P1.Z, P2.Z);
	modmult2(P1.X, P1.Z, R);
	add2(R, P1.Y, P2.Y);

	return 0;
}

int PrecomputationNothing(affinePointu4 P0, projectPointu4 *P, int w, int a)
{//  
	projectPointu4 P1[5];
	big X0[5], X1[5], X2[5], X3[5];
	for (int i = 0; i < 5; i++)
	{
		X0[i] = mirvar(0);
		X1[i] = mirvar(0);
		X2[i] = mirvar(1);
		X3[i] = mirvar(0);

		P1[i].X0 = X0[i];
		P1[i].X1 = X1[i];
		P1[i].X2 = X2[i];
		P1[i].X3 = X3[i];
	}
	copy(P0.X0, P1[0].X0);
	copy(P0.X1, P1[0].X0);
	copy(P0.X3, P1[0].X0);//P1[0]=P0

	tau(P1[0], P1[1]);//P1[1]=\tau P	
	negative(P1[1], P1[3]);//P1[3]=-\tau P 	
	negative(P1[0], P1[4]);//P1[4]=-P




	return 0;

}

int Solinas(affinePointu4 P0, projectPointu4 *P, long long &CPUcycles, int w, int a)
{//Solinas' pre-computation for w=4,5

 //big *x=(big *)mr_alloc(7,sizeof(big));
	projectPointu4 P1[7];
	big X0[7], X1[7], X2[7], X3[7];
	for (int i = 0; i < 7; i++)
	{
		X0[i] = mirvar(0);
		X1[i] = mirvar(0);
		X2[i] = mirvar(1);
		X3[i] = mirvar(0);

		P1[i].X0 = X0[i];
		P1[i].X1 = X1[i];
		P1[i].X2 = X2[i];
		P1[i].X3 = X3[i];
	}
	copy(P0.X0, P1[0].X0);
	copy(P0.X1, P1[0].X0);
	copy(P0.X3, P1[0].X0);//P1[0]=P0

	CPUcycles = cpucycles();
	tau(P1[0], P1[1]);//P1[1]=\tau P
	tau(P1[1], P1[2]);//P1[2]=\tau^2 P
	tau(P1[2], P1[3]);//P1[3]=\tau^3 P

	negative(P1[0], P1[6]);//P1[6]=-P

	Addition(P1[6], P1[2], P[3], a);//Q3=-P+\tau^2 P
	Addition(P1[0], P1[2], P[5], a);//Q5=P+\tau^2 P

	if (a == 0)
	{
		Addition(P1[0], P1[3], P[7], a);//Q7=P+\tau^3 P
	}
	else if (a == 1)
	{
		negative(P1[3], P1[4]);//P1[4]=-\tau^3 P
		Addition(P1[0], P1[4], P[7], a);//Q7=P-\tau^3 P
	}
	if (w == 5)
	{
		tau(P1[3], P1[4]);//P1[4]=\tau^4 P
		Addition(P1[6], P1[4], P[15], a);//Q15=-P+\tau^4 P

		tau(P[5], P1[1]);//P1[1]=\tau Q5
		tau(P1[1], P1[2]);//P1[2]=\tau^2 Q5
		tau(P1[2], P1[3]);//P1[3]=\tau^3 Q5

		negative(P1[2], P1[5]);//P1[5]=-\tau^2 Q5

		if (a == 0)
		{
			Addition(P1[0], P1[3], P[9], a);//Q9=P+\tau^3 Q5
		}
		else if (a == 1)
		{
			negative(P1[3], P1[4]);//P1[4]=-\tau^3 Q5
			Addition(P1[0], P1[4], P[9], a);//Q9=P-\tau^3 Q5 
		}

		Addition(P1[6], P1[5], P[11], a);//Q11=-P-\tau^2 Q5
		Addition(P1[0], P1[5], P[13], a);//Q13=P-\tau^2 Q5


	}
	CPUcycles = cpucycles() - CPUcycles;
	return 0;

}


int Solinas(affinePointLD P0, projectPointLD *P, long long &CPUcycles, int w, int a)
{//Solinas' pre-computation for w=4,5

 //big *x=(big *)mr_alloc(7,sizeof(big));
	projectPointLD P1[7];
	big X0[7], X1[7], X2[7];

	Big one = 1;
	for (int i = 0; i < 7; i++)
	{
		X0[i] = mirvar(0);
		X1[i] = mirvar(0);
		X2[i] = mirvar(1);


		P1[i].X = X0[i];
		P1[i].Y = X1[i];
		P1[i].Z = X2[i];

	}
	copy(P0.x, P1[0].X);
	copy(P0.y, P1[0].Y);
	//copy(one.getbig(), P1[0].Z);//P1[0]=P0

	CPUcycles = cpucycles();
	tau(P1[0], P1[1]);//P1[1]=\tau P
	tau(P1[1], P1[2]);//P1[2]=\tau^2 P
	tau(P1[2], P1[3]);//P1[3]=\tau^3 P

	negative(P1[0], P1[6]);//P1[6]=-P

	Addition(P1[6], P1[2], P[3], a);//Q3=-P+\tau^2 P
	Addition(P1[0], P1[2], P[5], a);//Q5=P+\tau^2 P

	if (a == 0)
	{
		Addition(P1[0], P1[3], P[7], a);//Q7=P+\tau^3 P
	}
	else if (a == 1)
	{
		negative(P1[3], P1[4]);//P1[4]=-\tau^3 P
		Addition(P1[0], P1[4], P[7], a);//Q7=P-\tau^3 P
	}
	if (w == 5)
	{
		tau(P1[3], P1[4]);//P1[4]=\tau^4 P
		Addition(P1[6], P1[4], P[15], a);//Q15=-P+\tau^4 P

		tau(P[5], P1[1]);//P1[1]=\tau Q5
		tau(P1[1], P1[2]);//P1[2]=\tau^2 Q5
		tau(P1[2], P1[3]);//P1[3]=\tau^3 Q5

		negative(P1[2], P1[5]);//P1[5]=-\tau^2 Q5

		if (a == 0)
		{
			Addition(P1[0], P1[3], P[9], a);//Q9=P+\tau^3 Q5
		}
		else if (a == 1)
		{
			negative(P1[3], P1[4]);//P1[4]=-\tau^3 Q5
			Addition(P1[0], P1[4], P[9], a);//Q9=P-\tau^3 Q5 
		}

		Addition(P1[6], P1[5], P[11], a);//Q11=-P-\tau^2 Q5
		Addition(P1[0], P1[5], P[13], a);//Q13=P-\tau^2 Q5


	}
	CPUcycles = cpucycles() - CPUcycles;
	return 0;

}


int HMV(affinePointu4 P0, projectPointu4 *P, long long &CPUcycles, int w, int a)
{//HMV' pre-computation for w=4,5,6

 //big *x=(big *)mr_alloc(7,sizeof(big));
	projectPointu4 P1[7];
	big X0[7], X1[7], X2[7], X3[7];
	for (int i = 0; i < 7; i++)
	{
		X0[i] = mirvar(0);
		X1[i] = mirvar(0);
		X2[i] = mirvar(1);
		X3[i] = mirvar(0);

		P1[i].X0 = X0[i];
		P1[i].X1 = X1[i];
		P1[i].X2 = X2[i];
		P1[i].X3 = X3[i];
	}
	copy(P0.X0, P1[0].X0);
	copy(P0.X1, P1[0].X1);
	copy(P0.X3, P1[0].X3);//P1[0]=P0

	CPUcycles = cpucycles();
	tau(P1[0], P1[1]);//P1[1]=\tau P
	tau(P1[1], P1[2]);//P1[2]=\tau^2 P
	tau(P1[2], P1[3]);//P1[3]=\tau^3 P

	negative(P1[0], P1[6]);//P1[6]=-P

	if (w < 6)
	{
		Addition(P1[6], P1[2], P[3], a);//Q3=-P+\tau^2 P
		Addition(P1[0], P1[2], P[5], a);//Q5=P+\tau^2 P

		if (a == 0)
		{
			Addition(P1[0], P1[3], P[7], a);//Q7=P+\tau^3 P
		}
		else if (a == 1)
		{
			negative(P1[3], P1[4]);//P1[4]=-\tau^3 P
			Addition(P1[0], P1[4], P[7], a);//Q7=P-\tau^3 P
		}
		if (w == 5)
		{


			tau(P[5], P1[1]);//P1[1]=\tau Q5
			tau(P1[1], P1[2]);//P1[2]=\tau^2 Q5
			tau(P1[2], P1[3]);//P1[3]=\tau^3 Q5

			negative(P1[2], P1[5]);//P1[5]=-\tau^2 Q5

			if (a == 0)
			{
				Addition(P1[0], P1[3], P[9], a);//Q9=P+\tau^3 Q5
			}
			else if (a == 1)
			{
				negative(P1[3], P1[4]);//P1[4]=-\tau^3 Q5
				Addition(P1[0], P1[4], P[9], a);//Q9=P-\tau^3 Q5 
			}

			Addition(P1[6], P1[5], P[11], a);//Q11=-P-\tau^2 Q5
			Addition(P1[0], P1[5], P[13], a);//Q13=P-\tau^2 Q5


			negative(P[5], P1[4]);//P1[4]=-\tau^3 Q5
			Addition(P1[2], P1[4], P[15], a);//Q15=-Q5+\tau^2 Q5

		}
	}
	else if (w == 6)
	{
		negative(P1[2], P1[2]);
		Addition(P1[6], P1[2], P[27], a);//Q27=-P-\tau^2 P
		Addition(P1[0], P1[2], P[29], a);//Q29=P-\tau^2 P

		negative(P1[3], P1[4]);//P1[4]=-\tau^3 P
		if (a == 0)
		{
			Addition(P1[0], P1[4], P[25], a);//Q25=P-\tau^3 P
			Addition(P1[6], P1[4], P[23], a);//Q23=-P-\tau^3 P
		}
		else if (a == 1)
		{
			Addition(P1[0], P1[4], P[23], a);//Q23=P-\tau^3 P
			Addition(P1[6], P1[4], P[25], a);//Q25=-P-\tau^3 P
		}
		tau(P[25], P1[1]);//P1[1]=\tau Q25
		tau(P1[1], P1[2]);//P1[2]=\tau^2 Q25

		Addition(P1[0], P1[2], P[5], a);//Q5=P+\tau^2 Q25
		Addition(P1[6], P1[2], P[3], a);//Q3=-P+\tau^2 Q25

		Addition(P[27], P1[2], P[31], a);//Q31=Q27+\tau^2 Q25

		tau(P[27], P1[1]);//P1[1]=\tau Q27
		tau(P1[1], P1[2]);//P1[2]=\tau^2 Q27
		negative(P1[2], P1[5]); //P1[5]=-\tau^2 Q27

		tau(P1[2], P1[3]);//P1[3]=\tau^3 Q27
		negative(P1[3], P1[4]); //P1[3]=-\tau^3 Q27
		if (a == 0)
		{
			Addition(P1[0], P1[4], P[9], a);//Q9=P-\tau^3 Q27
			Addition(P1[6], P1[4], P[7], a);//Q7=-P-\tau^3 Q27
		}
		else if (a == 1)
		{
			Addition(P1[0], P1[3], P[9], a);//Q9=P+\tau^3 Q27
			Addition(P1[6], P1[3], P[7], a);//Q7=-P+\tau^3 Q27
		}

		Addition(P1[0], P1[2], P[13], a);//Q13=P+\tau^2 Q27
		Addition(P1[6], P1[2], P[11], a);//Q11=-P+\tau^2 Q27

		Addition(P[27], P1[5], P[15], a);//Q15=Q27-\tau^2 Q27
		Addition(P[29], P1[5], P[17], a);//Q17=Q29-\tau^2 Q27

		tau(P[3], P1[1]);//P1[1]=\tau Q3
		tau(P1[1], P1[2]);//P1[2]=\tau^2 Q3
		negative(P1[2], P1[5]); //P1[5]=-\tau^2 Q3

		Addition(P1[6], P1[5], P[19], a);//Q19=-P-\tau^2 Q3

		tau(P[29], P1[1]);//P1[1]=\tau Q29
		tau(P1[1], P1[2]);//P1[2]=\tau^2 Q29
		Addition(P1[0], P1[2], P[21], a);//Q21=P+\tau^2 Q29

	}


	CPUcycles = cpucycles() - CPUcycles;
	return 0;

}


int HMV(affinePointLD P0, projectPointLD *P, long long &CPUcycles, int w, int a)
{//HMV' pre-computation for w=4,5,6

 //big *x=(big *)mr_alloc(7,sizeof(big));
	projectPointLD P1[7];
	big X0[7], X1[7], X2[7];
	for (int i = 0; i < 7; i++)
	{
		X0[i] = mirvar(0);
		X1[i] = mirvar(0);
		X2[i] = mirvar(1);

		P1[i].X = X0[i];
		P1[i].Y = X1[i];
		P1[i].Z = X2[i];
	}
	copy(P0.x, P1[0].X);
	copy(P0.y, P1[0].Y);
	//copy(P0.X3, P1[0].X0);//P1[0]=P0

	CPUcycles = cpucycles();
	tau(P1[0], P1[1]);//P1[1]=\tau P
	tau(P1[1], P1[2]);//P1[2]=\tau^2 P
	tau(P1[2], P1[3]);//P1[3]=\tau^3 P

	negative(P1[0], P1[6]);//P1[6]=-P

	if (w < 6)
	{
		Addition(P1[6], P1[2], P[3], a);//Q3=-P+\tau^2 P
		Addition(P1[0], P1[2], P[5], a);//Q5=P+\tau^2 P

		if (a == 0)
		{
			Addition(P1[0], P1[3], P[7], a);//Q7=P+\tau^3 P
		}
		else if (a == 1)
		{
			negative(P1[3], P1[4]);//P1[4]=-\tau^3 P
			Addition(P1[0], P1[4], P[7], a);//Q7=P-\tau^3 P
		}
		if (w == 5)
		{


			tau(P[5], P1[1]);//P1[1]=\tau Q5
			tau(P1[1], P1[2]);//P1[2]=\tau^2 Q5
			tau(P1[2], P1[3]);//P1[3]=\tau^3 Q5

			negative(P1[2], P1[5]);//P1[5]=-\tau^2 Q5

			if (a == 0)
			{
				Addition(P1[0], P1[3], P[9], a);//Q9=P+\tau^3 Q5
			}
			else if (a == 1)
			{
				negative(P1[3], P1[4]);//P1[4]=-\tau^3 Q5
				Addition(P1[0], P1[4], P[9], a);//Q9=P-\tau^3 Q5 
			}

			Addition(P1[6], P1[5], P[11], a);//Q11=-P-\tau^2 Q5
			Addition(P1[0], P1[5], P[13], a);//Q13=P-\tau^2 Q5


			negative(P[5], P1[4]);//P1[4]=-\tau^3 Q5
			Addition(P1[2], P1[4], P[15], a);//Q15=-Q5+\tau^2 Q5

		}
	}
	else if (w == 6)
	{
		negative(P1[2], P1[2]);
		Addition(P1[6], P1[2], P[27], a);//Q27=-P-\tau^2 P
		Addition(P1[0], P1[2], P[29], a);//Q29=P-\tau^2 P

		negative(P1[3], P1[4]);//P1[4]=-\tau^3 P
		if (a == 0)
		{
			Addition(P1[0], P1[4], P[25], a);//Q25=P-\tau^3 P
			Addition(P1[6], P1[4], P[23], a);//Q23=-P-\tau^3 P
		}
		else if (a == 1)
		{
			Addition(P1[0], P1[4], P[23], a);//Q23=P-\tau^3 P
			Addition(P1[6], P1[4], P[25], a);//Q25=-P-\tau^3 P
		}
		tau(P[25], P1[1]);//P1[1]=\tau Q25
		tau(P1[1], P1[2]);//P1[2]=\tau^2 Q25

		Addition(P1[0], P1[2], P[5], a);//Q5=P+\tau^2 Q25
		Addition(P1[6], P1[2], P[3], a);//Q3=-P+\tau^2 Q25

		Addition(P[27], P1[2], P[31], a);//Q31=Q27+\tau^2 Q25

		tau(P[27], P1[1]);//P1[1]=\tau Q27
		tau(P1[1], P1[2]);//P1[2]=\tau^2 Q27
		negative(P1[2], P1[5]); //P1[5]=-\tau^2 Q27

		tau(P1[2], P1[3]);//P1[3]=\tau^3 Q27
		negative(P1[3], P1[4]); //P1[3]=-\tau^3 Q27
		if (a == 0)
		{
			Addition(P1[0], P1[4], P[9], a);//Q9=P-\tau^3 Q27
			Addition(P1[6], P1[4], P[7], a);//Q7=-P-\tau^3 Q27
		}
		else if (a == 1)
		{
			Addition(P1[0], P1[3], P[9], a);//Q9=P+\tau^3 Q27
			Addition(P1[6], P1[3], P[7], a);//Q7=-P+\tau^3 Q27
		}

		Addition(P1[0], P1[2], P[13], a);//Q13=P+\tau^2 Q27
		Addition(P1[6], P1[2], P[11], a);//Q11=-P+\tau^2 Q27

		Addition(P[27], P1[5], P[15], a);//Q15=Q27-\tau^2 Q27
		Addition(P[29], P1[5], P[17], a);//Q17=Q29-\tau^2 Q27

		tau(P[3], P1[1]);//P1[1]=\tau Q3
		tau(P1[1], P1[2]);//P1[2]=\tau^2 Q3
		negative(P1[2], P1[5]); //P1[5]=-\tau^2 Q3

		Addition(P1[6], P1[5], P[19], a);//Q19=-P-\tau^2 Q3

		tau(P[29], P1[1]);//P1[1]=\tau Q29
		tau(P1[1], P1[2]);//P1[2]=\tau^2 Q29
		Addition(P1[0], P1[2], P[21], a);//Q21=P+\tau^2 Q29

	}


	CPUcycles = cpucycles() - CPUcycles;
	return 0;

}

int xuPreComputation(affinePointu4 P0, projectPointu4 *P, long long &CPUcycles, int w, int a)
{//xu' pre-computation for w=4,5,6

 //big *x=(big *)mr_alloc(7,sizeof(big));
	projectPointu4 P1[7];
	big X0[7], X1[7], X2[7], X3[7];
	for (int i = 0; i < 7; i++)
	{
		X0[i] = mirvar(0);
		X1[i] = mirvar(0);
		X2[i] = mirvar(1);
		X3[i] = mirvar(0);

		P1[i].X0 = X0[i];
		P1[i].X1 = X1[i];
		P1[i].X2 = X2[i];
		P1[i].X3 = X3[i];
	}
	copy(P0.X0, P1[0].X0);
	copy(P0.X1, P1[0].X1);
	copy(P0.X3, P1[0].X3);//P1[0]=P0

	CPUcycles = cpucycles();
	tau(P1[0], P1[1]);//P1[1]=\tau P
	tau(P1[1], P1[2]);//P1[2]=\tau^2 P
	negative(P1[1], P1[5]);//P1[5]=-\tau P 
	negative(P1[2], P1[4]);//P1[4]=-\tau^2 P
	negative(P1[0], P1[6]);//P1[6]=-P

	//if (a == 1)
	//{
	//	negative(P1[1], P1[3]);//P1[3]=-\tau P 
	//	negative(P1[3], P1[5]);//P1[5]=\tau P 
	//	negative(P1[5], P1[1]);//P1[1]=-\tau P 
	//}

	if (w < 6)
	{
		Addition(P1[6], P1[2], P[3], a);//Q3=-P+\tau^2 P


		if (a == 0)
		{
			Addition(P1[0], P1[5], P[7], a);//Q7=P-\tau P
			Addition(P1[6], P1[5], P[5], a);//Q5=-P-\tau P
		}
		else if (a == 1)
		{
			Addition(P1[0], P1[1], P[7], a);//Q7=P+\tau P
			Addition(P1[6], P1[1], P[5], a);//Q5=-P+\tau P
		}
		if (w == 5)
		{
			if (a == 0)
			{
				Addition(P1[5], P[3], P[9], a);//Q9=Q3-\tau P
				Addition(P1[5], P[5], P[11], a);//Q11=Q5-\tau P

				Addition(P1[5], P[7], P[13], a);//Q13=Q7-\tau P

				negative(P[11], P1[3]);
				Addition(P1[1], P1[3], P[11], a);//Q15=-Q11+\tau P
			}
			else if (a == 1)
			{
				Addition(P1[1], P[3], P[9], a);//Q9=Q3+\tau P
				Addition(P1[1], P[5], P[11], a);//Q11=Q5+\tau P

				Addition(P1[1], P[7], P[13], a);//Q13=Q7+\tau P

				negative(P[11], P1[3]);
				Addition(P1[5], P1[3], P[11], a);//Q15=-Q11-\tau P
			}


		}
	}
	else if (w == 6)
	{
		Addition(P1[0], P1[4], P[29], a);//Q29=P-\tau^2 P

		if (a == 1)
		{
			Addition(P1[0], P1[5], P[27], a);//Q27=P-\tau P
			Addition(P1[6], P1[5], P[25], a);//Q25=-P-\tau P

			Addition(P1[1], P[29], P[3], a);//Q3=Q29+\tau P
			negative(P[29], P1[3]);
			Addition(P1[1], P1[3], P[9], a);//Q9=-Q29+\tau P

			Addition(P1[4], P[3], P[31], a);//Q31=Q3-\tau^2 P

			Addition(P1[1], P[31], P[5], a);//Q5=Q31+\tau P

			negative(P[31], P1[3]);
			Addition(P1[1], P1[3], P[7], a);//Q7=-Q31+\tau P

			negative(P[27], P1[3]);
			Addition(P1[5], P1[3], P[11], a);//Q11=-Q27-\tau P

			negative(P[25], P1[3]);
			Addition(P1[1], P1[3], P[13], a);//Q13=-Q25+\tau P

			negative(P[11], P1[3]);
			Addition(P1[5], P1[3], P[15], a);//Q15=-Q11-\tau P

			negative(P[9], P1[3]);
			Addition(P1[5], P1[3], P[17], a);//Q17=-Q9-\tau P

			negative(P[7], P1[3]);
			Addition(P1[5], P1[3], P[19], a);//Q19=-Q7-\tau P

			negative(P[17], P1[3]);
			Addition(P1[1], P1[3], P[21], a);//Q21=-Q17+\tau P

			negative(P[3], P1[3]);
			Addition(P1[5], P1[3], P[23], a);//Q23=-Q3-\tau P
		}
		else if (a == 0)
		{
			Addition(P1[0], P1[1], P[27], a);//Q27=P+\tau P
			Addition(P1[6], P1[1], P[25], a);//Q25=-P+\tau P

			Addition(P1[5], P[29], P[3], a);//Q3=Q29-\tau P
			negative(P[29], P1[3]);
			Addition(P1[5], P1[3], P[9], a);//Q9=-Q29-\tau P

			Addition(P1[4], P[3], P[31], a);//Q31=Q3-\tau^2 P
	
			Addition(P1[5], P[31], P[5], a);//Q5=Q31-\tau P

			negative(P[31], P1[3]);
			Addition(P1[5], P1[3], P[7], a);//Q7=-Q31-\tau P

			negative(P[27], P1[3]);
			Addition(P1[5], P1[3], P[11], a);//Q11=-Q27-\tau P

			negative(P[25], P1[3]);
			Addition(P1[5], P1[3], P[13], a);//Q13=-Q25-\tau P

			negative(P[11], P1[3]);
			Addition(P1[1], P1[3], P[15], a);//Q15=-Q11+\tau P

			negative(P[9], P1[3]);
			Addition(P1[1], P1[3], P[17], a);//Q17=-Q9+\tau P

			negative(P[7], P1[3]);
			Addition(P1[1], P1[3], P[19], a);//Q19=-Q7+\tau P

			negative(P[17], P1[3]);
			Addition(P1[5], P1[3], P[21], a);//Q21=-Q17-\tau P

			negative(P[3], P1[3]);
			Addition(P1[1], P1[3], P[23], a);//Q23=-Q3+\tau P
		}

	}

	CPUcycles = cpucycles() - CPUcycles;
	return 0;

}


int xuPreComputation(affinePointLD P0, projectPointLD *P, long long &CPUcycles, int w, int a)
{//xu' pre-computation for w=4,5,6

 //big *x=(big *)mr_alloc(7,sizeof(big));
	projectPointLD P1[7];
	big X0[7], X1[7], X2[7];
	for (int i = 0; i < 7; i++)
	{
		X0[i] = mirvar(0);
		X1[i] = mirvar(0);
		X2[i] = mirvar(1);


		P1[i].X = X0[i];
		P1[i].Y = X1[i];
		P1[i].Z = X2[i];

	}
	copy(P0.x, P1[0].X);
	copy(P0.y, P1[0].Y);
	//copy(P0.X3, P1[0].X0);//P1[0]=P0

	CPUcycles = cpucycles();
	tau(P1[0], P1[1]);//P1[1]=\tau P
	tau(P1[1], P1[2]);//P1[2]=\tau^2 P
	negative(P1[1], P1[5]);//P1[5]=-\tau P 
	negative(P1[2], P1[4]);//P1[4]=-\tau^2 P
	negative(P1[0], P1[6]);//P1[6]=-P

						   //if (a == 1)
						   //{
						   //	negative(P1[1], P1[3]);//P1[3]=-\tau P 
						   //	negative(P1[3], P1[5]);//P1[5]=\tau P 
						   //	negative(P1[5], P1[1]);//P1[1]=-\tau P 
						   //}

	if (w < 6)
	{
		Addition(P1[6], P1[2], P[3], a);//Q3=-P+\tau^2 P


		if (a == 0)
		{
			Addition(P1[0], P1[5], P[7], a);//Q7=P-\tau P
			Addition(P1[6], P1[5], P[5], a);//Q5=-P-\tau P
		}
		else if (a == 1)
		{
			Addition(P1[0], P1[1], P[7], a);//Q7=P+\tau P
			Addition(P1[6], P1[1], P[5], a);//Q5=-P+\tau P
		}
		if (w == 5)
		{
			if (a == 0)
			{
				Addition(P1[5], P[3], P[9], a);//Q9=Q3-\tau P
				Addition(P1[5], P[5], P[11], a);//Q11=Q5-\tau P

				Addition(P1[5], P[7], P[13], a);//Q13=Q7-\tau P

				negative(P[11], P1[3]);
				Addition(P1[1], P1[3], P[11], a);//Q15=-Q11+\tau P
			}
			else if (a == 1)
			{
				Addition(P1[1], P[3], P[9], a);//Q9=Q3+\tau P
				Addition(P1[1], P[5], P[11], a);//Q11=Q5+\tau P

				Addition(P1[1], P[7], P[13], a);//Q13=Q7+\tau P

				negative(P[11], P1[3]);
				Addition(P1[5], P1[3], P[11], a);//Q15=-Q11-\tau P
			}


		}
	}
	else if (w == 6)
	{
		Addition(P1[0], P1[4], P[29], a);//Q29=P-\tau^2 P

		if (a == 1)
		{
			Addition(P1[0], P1[5], P[27], a);//Q27=P-\tau P
			Addition(P1[6], P1[5], P[25], a);//Q25=-P-\tau P

			Addition(P1[1], P[29], P[3], a);//Q3=Q29+\tau P
			negative(P[29], P1[3]);
			Addition(P1[1], P1[3], P[9], a);//Q9=-Q29+\tau P
		}
		else if (a == 0)
		{
			Addition(P1[0], P1[1], P[27], a);//Q27=P+\tau P
			Addition(P1[6], P1[1], P[25], a);//Q25=-P+\tau P

			Addition(P1[5], P[29], P[3], a);//Q3=Q29-\tau P
			negative(P[29], P1[3]);
			Addition(P1[5], P1[3], P[9], a);//Q9=-Q29-\tau P
		}
		Addition(P1[4], P[3], P[31], a);//Q31=Q3-\tau^2 P

		if (a == 0)
		{
			Addition(P1[5], P[31], P[5], a);//Q5=Q31-\tau P

			negative(P[31], P1[3]);
			Addition(P1[5], P1[3], P[7], a);//Q7=-Q31-\tau P

			negative(P[27], P1[3]);
			Addition(P1[5], P1[3], P[11], a);//Q11=-Q27-\tau P

			negative(P[25], P1[3]);
			Addition(P1[5], P1[3], P[13], a);//Q13=-Q25-\tau P

			negative(P[11], P1[3]);
			Addition(P1[1], P1[3], P[15], a);//Q15=-Q11+\tau P

			negative(P[9], P1[3]);
			Addition(P1[1], P1[3], P[17], a);//Q17=-Q9+\tau P

			negative(P[7], P1[3]);
			Addition(P1[1], P1[3], P[19], a);//Q19=-Q7+\tau P

			negative(P[17], P1[3]);
			Addition(P1[5], P1[3], P[21], a);//Q21=-Q17-\tau P

			negative(P[3], P1[3]);
			Addition(P1[1], P1[3], P[23], a);//Q23=-Q3+\tau P
		}
		else if (a == 1)
		{
			Addition(P1[1], P[31], P[5], a);//Q5=Q31+\tau P

			negative(P[31], P1[3]);
			Addition(P1[1], P1[3], P[7], a);//Q7=-Q31+\tau P

			negative(P[27], P1[3]);
			Addition(P1[5], P1[3], P[11], a);//Q11=-Q27-\tau P

			negative(P[25], P1[3]);
			Addition(P1[1], P1[3], P[13], a);//Q13=-Q25+\tau P

			negative(P[11], P1[3]);
			Addition(P1[5], P1[3], P[15], a);//Q15=-Q11-\tau P

			negative(P[9], P1[3]);
			Addition(P1[5], P1[3], P[17], a);//Q17=-Q9-\tau P

			negative(P[7], P1[3]);
			Addition(P1[5], P1[3], P[19], a);//Q19=-Q7-\tau P

			negative(P[17], P1[3]);
			Addition(P1[1], P1[3], P[21], a);//Q21=-Q17+\tau P

			negative(P[3], P1[3]);
			Addition(P1[5], P1[3], P[23], a);//Q23=-Q3-\tau P
		}
	}

	CPUcycles = cpucycles() - CPUcycles;
	return 0;

}

int OurPreComputation(affinePointu4 P0, projectPointu4 *P, long long &CPUcycles, int w, int a)
{//Our pre-computation for w=4,5,6,7

 //big *x=(big *)mr_alloc(7,sizeof(big));
	projectPointu4 P1[5];
	big X0[5], X1[5], X2[5], X3[5];
	for (int i = 0; i < 5; i++)
	{
		X0[i] = mirvar(0);
		X1[i] = mirvar(0);
		X2[i] = mirvar(1);
		X3[i] = mirvar(0);

		P1[i].X0 = X0[i];
		P1[i].X1 = X1[i];
		P1[i].X2 = X2[i];
		P1[i].X3 = X3[i];
	}
	copy(P0.X0, P1[0].X0);
	copy(P0.X1, P1[0].X1);
	copy(P0.X3, P1[0].X3);//P1[0]=P0


	tau(P1[0], P1[1]);//P1[1]=\tau P	
	negative(P1[1], P1[3]);//P1[3]=-\tau P 	
	negative(P1[0], P1[4]);//P1[4]=-P
	CPUcycles = cpucycles();


	if (w < 6)
	{
		barTau(P1[4], P[5], a);//P5=-\mu\bar\tau P
		barTau(P[5], P[7], a);//P7=-(\mu\bar\tau)^2 P
		barTau(P[7], P1[2], a);//P1[2]=-(\mu\bar\tau)^3 P

		negative(P1[2], P[3]);//P3=(\mu\bar\tau)^3 P



		if (w == 5)
		{
			barTau(P1[2], P[15], a);//P15=-(\mu\bar\tau)^4 P

			Addition(P1[1], P[5], P[11], a);//Q11=Q5-\mu\tau P

			barTau(P[11], P[9], a);//P9= (\mu\bar\tau) Q11
			barTau(P[9], P1[2], a);//P1[2]=(\mu\bar\tau)^2 Q11

			negative(P1[2], P[13]);//P13=-(\mu\bar\tau)^2 Q11
		}
	}
	else if (w == 6)
	{
		barTau(P1[0], P[27], a);//P27=\mu\bar\tau P
		barTau(P[27], P[25], a);//P25=(\mu\bar\tau)^2 P
		barTau(P[25], P1[2], a);//P1[2]=(\mu\bar\tau)^3 P

		negative(P1[2], P[29]);//P29=-(\mu\bar\tau)^3 P




		barTau(P[29], P[15], a);//P15=-(\mu\bar\tau)^4 P
		barTau(P[15], P[21], a);//P21=(\mu\bar\tau)^5 P

		PpmQ(P1[1], P[29], P[3], P[9], a);//Q3=\mu\tau P+Q29,Q9=\mu\tau P-Q29

		barTau(P[9], P1[2], a);//P1[2]=(\mu\bar\tau) Q9

		negative(P1[2], P[13]);//P13=-(\mu\bar\tau) Q9
		barTau(P[13], P[31], a);//P31=-(\mu\bar\tau)^2 Q9

		barTau(P[3], P[17], a);//P17=(\mu\bar\tau) Q3
		barTau(P[17], P[11]);//P11=(\mu\bar\tau)^2 Q3

		negative(P[15], P1[2]);
		mixedAddition(P1[2], P1[1], P[23], a);//$Q_{23}=\mu\tau P-Q_{15}$ 

		barTau(P[23], P1[2], a);//P1[2]=(\mu\bar\tau) Q23
		negative(P1[2], P[19]);//P19=-(\mu\bar\tau) Q23

		negative(P[21], P1[2]);
		mixedAddition(P1[2], P1[3], P[5], a);//$Q_{5}=-\mu\tau P -Q_{21}$ 

		barTau(P[5], P[7], a);//$Q_{7}=\mu\bar\tau Q_{5}$  


	}
	else if (w == 7)
	{
		barTau(P1[4], P[37], a);//$Q_{37}=-\mu\bar\tau  P$ 
		barTau(P[37], P[39], a);//$Q_{39}=-(\mu\bar\tau)^2 P $
		barTau(P[39], P1[2], a);//$Q_{35}=(\mu\bar\tau)^3 P$

		negative(P1[2], P[35]);//$Q_{35}=(\mu\bar\tau)^3 P$




		barTau(P1[2], P[15], a);//$Q_{15}=-(\mu\bar\tau)^4 P$ 
		barTau(P[15], P1[2], a);//$Q_{43}=(\mu\bar\tau)^5 P$
		negative(P1[2], P[43]);//$Q_{43}=(\mu\bar\tau)^5 P$


		PpmQ(P1[1], P[15], P[53], P[23], a);//$Q_{53}=\mu\tau P +Q_{15}$,$Q_{23}=\mu\tau P-Q_{15}$ 

		negative(P[53], P1[2]);
		//$Q_{41}=-\mu\bar\tau Q_{53}$
		barTau(P1[2], P[41], a);
		//$Q_{19}=-(\mu\bar\tau)^2 Q_{53}$  
		barTau(P[41], P[19], a);
		//$Q_{63}=(\mu\bar\tau)^3 Q_{53}$
		barTau(P[53], P1[2], a);
		negative(P1[2], P[63]);
		//$Q_{27}=-(\mu\bar\tau)^4 Q_{53}$ 
		barTau(P1[2], P[27], a);

		//$Q_{ 45 } = \mu\bar\tau Q_{ 23 }$
		barTau(P[23], P[45], a);

		//$Q_{ 3 } = \mu\tau P - Q_{ 35 }$, $Q_{ 55 } = -\mu\tau P - Q_{ 35 }$
		PpmQ(P1[1], P[35], P[3], P1[2], a);
		negative(P1[2], P[55]);


		//$Q_{ 17 } = \mu\bar\tau Q_{ 3 }$  
		barTau(P[3], P[17], a);
		//$Q_{ 11 } = (\mu\bar\tau) ^ 2 Q_{ 3 }$
		barTau(P[17], P[11], a);


		//$Q_{ 13 } = \mu\bar\tau Q_{ 55 }$
		barTau(P[55], P[13], a);
		//$Q_{ 31 } = (\mu\bar\tau) ^ 2 Q_{ 55 }$
		barTau(P[13], P[31], a);
		//$Q_{ 5 } = (\mu\bar\tau) ^ 3 Q_{ 55 }$
		barTau(P[31], P[5], a);

		//$Q_{ 51 } = \mu\tau P + Q_{ 13 }$,$Q_{25}=\mu\tau P-Q_{13}$
		PpmQ(P1[1], P[13], P[51], P[25], a);

		//$Q_{ 33 } = \mu\bar\tau Q_{ 51 }$
		barTau(P[51], P[33], a);
		//$Q_{ 59 } = (\mu\bar\tau) ^ 2 Q_{ 51 }$
		barTau(P[33], P[59], a);
		//$Q_{7}=-(\mu\bar\tau)^3 Q_{51}$
		barTau(P[59], P1[2], a);
		negative(P1[2], P[7]);



		//$Q_{29}=-\mu\bar\tau Q_{25}$ 
		barTau(P[25], P1[2], a);
		negative(P1[2], P[29]);

		//$Q_{ 49 } = -\mu\tau P - Q_{ 41 }$
		negative(P[41], P1[2]);
		mixedAddition(P1[2], P1[3], P[49], a);
		//$Q_{21}=-\mu\bar\tau Q_{49}$
		barTau(P[49], P1[2], a);
		negative(P1[2], P[21]);

		//$Q_{ 9 } = (\mu\bar\tau) ^ 2 Q_{ 49 }$
		barTau(P1[2], P[27], a);


		//$Q_{57}=-\mu\tau P -Q_{ 33}$ 
		negative(P[33], P1[2]);
		mixedAddition(P1[2], P1[3], P[57], a);
		//$Q_{61}=-\mu\bar\tau Q_{57}$ 
		barTau(P[57], P1[2], a);
		negative(P1[2], P[61]);

		//$Q_{47}=-(\mu\bar\tau)^2 Q_{57}$ 
		barTau(P[61], P[47], a);

	}

	CPUcycles = cpucycles() - CPUcycles;

	return 0;

}

int OurPreComputation(affinePointLD P0, projectPointLD *P, long long &CPUcycles, int w, int a)
{//Our pre-computation for w=4,5,6,7

 //big *x=(big *)mr_alloc(7,sizeof(big));
	projectPointLD P1[5];
	big X0[5], X1[5], X2[5];
	for (int i = 0; i < 5; i++)
	{
		X0[i] = mirvar(0);
		X1[i] = mirvar(0);
		X2[i] = mirvar(1);


		P1[i].X = X0[i];
		P1[i].Y = X1[i];
		P1[i].Z = X2[i];
	}
	copy(P0.x, P1[0].X);
	copy(P0.y, P1[0].Y);
	//copy(P0.X3, P1[0].X3);//P1[0]=P0


	tau(P1[0], P1[1]);//P1[1]=\tau P	
	negative(P1[1], P1[3]);//P1[3]=-\tau P 	
	negative(P1[0], P1[4]);//P1[4]=-P
	CPUcycles = cpucycles();


	if (w < 6)
	{
		barTau(P1[4], P[5], a);//P5=-\mu\bar\tau P
		barTau(P[5], P[7], a);//P7=-(\mu\bar\tau)^2 P
		barTau(P[7], P1[2], a);//P1[2]=-(\mu\bar\tau)^3 P

		negative(P1[2], P[3]);//P3=(\mu\bar\tau)^3 P



		if (w == 5)
		{
			barTau(P1[2], P[15], a);//P15=-(\mu\bar\tau)^4 P

			Addition(P1[1], P[5], P[11], a);//Q11=Q5-\mu\tau P

			barTau(P[11], P[9], a);//P9= (\mu\bar\tau) Q11
			barTau(P[9], P1[2], a);//P1[2]=(\mu\bar\tau)^2 Q11

			negative(P1[2], P[13]);//P13=-(\mu\bar\tau)^2 Q11
		}
	}
	else if (w == 6)
	{
		barTau(P1[0], P[27], a);//P27=\mu\bar\tau P
		barTau(P[27], P[25], a);//P25=(\mu\bar\tau)^2 P
		barTau(P[25], P1[2], a);//P1[2]=(\mu\bar\tau)^3 P

		negative(P1[2], P[29]);//P29=-(\mu\bar\tau)^3 P




		barTau(P[29], P[15], a);//P15=-(\mu\bar\tau)^4 P
		barTau(P[15], P[21], a);//P21=(\mu\bar\tau)^5 P

		PpmQ(P1[1], P[29], P[3], P[9], a);//Q3=\mu\tau P+Q29,Q9=\mu\tau P-Q29

		barTau(P[9], P1[2], a);//P1[2]=(\mu\bar\tau) Q9

		negative(P1[2], P[13]);//P13=-(\mu\bar\tau) Q9
		barTau(P[13], P[31], a);//P31=-(\mu\bar\tau)^2 Q9

		barTau(P[3], P[17], a);//P17=(\mu\bar\tau) Q3
		barTau(P[17], P[11]);//P11=(\mu\bar\tau)^2 Q3

		negative(P[15], P1[2]);
		mixedAddition(P1[2], P1[1], P[23], a);//$Q_{23}=\mu\tau P-Q_{15}$ 

		barTau(P[23], P1[2], a);//P1[2]=(\mu\bar\tau) Q23
		negative(P1[2], P[19]);//P19=-(\mu\bar\tau) Q23

		negative(P[21], P1[2]);
		mixedAddition(P1[2], P1[3], P[5], a);//$Q_{5}=-\mu\tau P -Q_{21}$ 

		barTau(P[5], P[7], a);//$Q_{7}=\mu\bar\tau Q_{5}$  


	}
	else if (w == 7)
	{
		barTau(P1[4], P[37], a);//$Q_{37}=-\mu\bar\tau  P$ 
		barTau(P[37], P[39], a);//$Q_{39}=-(\mu\bar\tau)^2 P $
		barTau(P[39], P1[2], a);//$Q_{35}=(\mu\bar\tau)^3 P$

		negative(P1[2], P[35]);//$Q_{35}=(\mu\bar\tau)^3 P$




		barTau(P1[2], P[15], a);//$Q_{15}=-(\mu\bar\tau)^4 P$ 
		barTau(P[15], P1[2], a);//$Q_{43}=(\mu\bar\tau)^5 P$
		negative(P1[2], P[43]);//$Q_{43}=(\mu\bar\tau)^5 P$


		PpmQ(P1[1], P[15], P[53], P[23], a);//$Q_{53}=\mu\tau P +Q_{15}$,$Q_{23}=\mu\tau P-Q_{15}$ 

		negative(P[53], P1[2]);
		//$Q_{41}=-\mu\bar\tau Q_{53}$
		barTau(P1[2], P[41], a);
		//$Q_{19}=-(\mu\bar\tau)^2 Q_{53}$  
		barTau(P[41], P[19], a);
		//$Q_{63}=(\mu\bar\tau)^3 Q_{53}$
		barTau(P[53], P1[2], a);
		negative(P1[2], P[63]);
		//$Q_{27}=-(\mu\bar\tau)^4 Q_{53}$ 
		barTau(P1[2], P[27], a);

		//$Q_{ 45 } = \mu\bar\tau Q_{ 23 }$
		barTau(P[23], P[45], a);

		//$Q_{ 3 } = \mu\tau P - Q_{ 35 }$, $Q_{ 55 } = -\mu\tau P - Q_{ 35 }$
		PpmQ(P1[1], P[35], P[3], P1[2], a);
		negative(P1[2], P[55]);


		//$Q_{ 17 } = \mu\bar\tau Q_{ 3 }$  
		barTau(P[3], P[17], a);
		//$Q_{ 11 } = (\mu\bar\tau) ^ 2 Q_{ 3 }$
		barTau(P[17], P[11], a);


		//$Q_{ 13 } = \mu\bar\tau Q_{ 55 }$
		barTau(P[55], P[13], a);
		//$Q_{ 31 } = (\mu\bar\tau) ^ 2 Q_{ 55 }$
		barTau(P[13], P[31], a);
		//$Q_{ 5 } = (\mu\bar\tau) ^ 3 Q_{ 55 }$
		barTau(P[31], P[5], a);

		//$Q_{ 51 } = \mu\tau P + Q_{ 13 }$,$Q_{25}=\mu\tau P-Q_{13}$
		PpmQ(P1[1], P[13], P[51], P[25], a);

		//$Q_{ 33 } = \mu\bar\tau Q_{ 51 }$
		barTau(P[51], P[33], a);
		//$Q_{ 59 } = (\mu\bar\tau) ^ 2 Q_{ 51 }$
		barTau(P[33], P[59], a);
		//$Q_{7}=-(\mu\bar\tau)^3 Q_{51}$
		barTau(P[59], P1[2], a);
		negative(P1[2], P[7]);



		//$Q_{29}=-\mu\bar\tau Q_{25}$ 
		barTau(P[25], P1[2], a);
		negative(P1[2], P[29]);

		//$Q_{ 49 } = -\mu\tau P - Q_{ 41 }$
		negative(P[41], P1[2]);
		mixedAddition(P1[2], P1[3], P[49], a);
		//$Q_{21}=-\mu\bar\tau Q_{49}$
		barTau(P[49], P1[2], a);
		negative(P1[2], P[21]);

		//$Q_{ 9 } = (\mu\bar\tau) ^ 2 Q_{ 49 }$
		barTau(P1[2], P[27], a);


		//$Q_{57}=-\mu\tau P -Q_{ 33}$ 
		negative(P[33], P1[2]);
		mixedAddition(P1[2], P1[3], P[57], a);
		//$Q_{61}=-\mu\bar\tau Q_{57}$ 
		barTau(P[57], P1[2], a);
		negative(P1[2], P[61]);

		//$Q_{47}=-(\mu\bar\tau)^2 Q_{57}$ 
		barTau(P[61], P[47], a);

	}

	CPUcycles = cpucycles() - CPUcycles;

	return 0;

}

int scalarMultiplication(affinePointu4 P0, affinePointu4 negative, projectPointu4 *P, projectPointu4 negativePi, int *a, int aLength, projectPointu4 Q, int w, int parameterA)
{ //using pre-computation in lambda-projective coordinates.
	tau(P0, Q); //Q=2P0;


	projectPointu4 middleP;
	middleP.X0 = mirvar(0);
	middleP.X1 = mirvar(0);
	middleP.X2 = mirvar(0);
	middleP.X3 = mirvar(0);

	int indexMiddle = a[aLength - 1];
	if (indexMiddle > 0)
	{
		copy(P[a[aLength - 1]].X0, Q.X0);
		copy(P[a[aLength - 1]].X1, Q.X1);

		copy(P[a[aLength - 1]].X2, Q.X2);
		copy(P[a[aLength - 1]].X3, Q.X3);
	}
	else
	{
		indexMiddle = -indexMiddle;
		copy(P[indexMiddle].X0, Q.X0);
		copy(P[indexMiddle].X1, Q.X3);

		copy(P[indexMiddle].X2, Q.X2);
		copy(P[indexMiddle].X3, Q.X1);
	}

	int b = parameterA;


	for (int i = aLength - 2; i > -1; i--)
	{
		if (a[i] > 0)
		{
			if (a[i] == 1)
			{
				mixedAddition(P0, Q, middleP, b);
			}
			else
			{
				Addition(P[a[i]], Q, middleP, b);
			}
		}
		else if (a[i] < 0)
		{
			if (a[i] == -1)
			{
				mixedAddition(negative, Q, middleP);//Q-P0
			}
			else
			{
				copy(P[-a[i]].X0, negativePi.X0);
				copy(P[-a[i]].X2, negativePi.X2);
				copy(P[-a[i]].X1, negativePi.X3);
				copy(P[-a[i]].X3, negativePi.X1);
				Addition(negativePi, Q, middleP, b);
			}
		}
		else
		{
			copy(Q.X0, middleP.X0);
			copy(Q.X1, middleP.X1);

			copy(Q.X2, middleP.X2);
			copy(Q.X3, middleP.X3);
		}

		tau(middleP, Q);
	}
	return 0;
}


int scalarMultiplication(affinePointLD P0, affinePointLD negative, projectPointLD *P, projectPointLD negativePi, int *a, int aLength, projectPointLD Q, int w, int parameterA)
{ //using pre-computation in LD coordinates.
	tau(P0, Q); //Q=2P0;


	projectPointLD middleP;
	middleP.X = mirvar(0);
	middleP.Y = mirvar(0);
	middleP.Z = mirvar(1);
	big R1=mirvar(0);
	 

	int indexMiddle = a[aLength - 1];
	if (indexMiddle > 0)
	{
		copy(P[a[aLength - 1]].X, Q.X);
		copy(P[a[aLength - 1]].Y, Q.Y);

		copy(P[a[aLength - 1]].Z, Q.Z);
		//copy(P[a[aLength - 1]].X3, Q.X3);
	}
	else
	{
		indexMiddle = -indexMiddle;
		copy(P[indexMiddle].X, Q.X);
		copy(P[indexMiddle].Z, Q.Z);
		modmult2(Q.X, Q.Z, R1);
		add2(P[indexMiddle].Y, R1, Q.Y);		 
	}

	int b = parameterA;


	for (int i = aLength - 2; i > -1; i--)
	{
		if (a[i] > 0)
		{
			if (a[i] == 1)
			{
				mixedAddition(P0, Q, middleP, b);
			}
			else
			{
				Addition(P[a[i]], Q, middleP, b);
			}
		}
		else if (a[i] < 0)
		{
			if (a[i] == -1)
			{
				mixedAddition(negative, Q, middleP);//Q-P0
			}
			else
			{
				copy(P[-a[i]].X, negativePi.X);
				copy(P[-a[i]].Z, negativePi.Z);
				modmult2(negativePi.X, negativePi.Z, R1);
				add2(P[-a[i]].Y, R1, negativePi.Y);
				//copy(P[-a[i]].X3, negativePi.X1);
				Addition(negativePi, Q, middleP, b);
			}
		}
		else
		{
			copy(Q.X, middleP.X);
			copy(Q.Y, middleP.Y);

			copy(Q.Z, middleP.Z);
			//copy(Q.X3, middleP.X3);
		}

		tau(middleP, Q);
	}
	return 0;
}


int HMV6(affinePoint P0, projectPoint *P, int w = 6)
{//HMV pre-computation for w=6. Only works when a=0

	//big *x=(big *)mr_alloc(7,sizeof(big));



	big x[15], y[15], z[15], x0, y0;

	x0 = mirvar(0);
	y0 = mirvar(0);
	copy(P0.x, x0);
	copy(P0.y, y0);
	for (int i = 0; i < 15; i++)
	{
		x[i] = mirvar(0);
		y[i] = mirvar(0);
		z[i] = mirvar(1);
	}
	//big *y=(big *)mr_alloc(7,sizeof(big));
	//big *z=(big *)mr_alloc(7,sizeof(big));
	big middleNegative = mirvar(0);//the negative lambda of lambda coordinate 
	big middlex = mirvar(0), middley = mirvar(0), middlez = mirvar(1), midx = mirvar(0), midy = mirvar(0), midz = mirvar(1);

	big  one = mirvar(1);
	add2(y0, one, middleNegative);//-P

	tauAffine(x0, middleNegative, middlex, middley);//-tau
	tauAffine(middlex, middley, midx, midy);//-tau^2
	mixedAddition(x0, middleNegative, midx, midy, x[12], y[12], z[12]);//Q27


	mixedAddition(x0, y0, midx, midy, x[13], y[13], z[13]);//Q29

	tauAffine(midx, midy, middlex, middley);//-tau^3

	mixedAddition(x0, middleNegative, middlex, middley, x[10], y[10], z[10]);//Q23
	mixedAddition(x0, y0, middlex, middley, x[11], y[11], z[11]);//Q25



	tauProjective(x[11], y[11], z[11], middlex, middley, middlez);//tau Q_25
	tauProjective(middlex, middley, middlez, midx, midy, midz);//tau^2 Q_25

	mixedAddition(x0, middleNegative, midx, midy, midz, x[0], y[0], z[0]);//Q3
	mixedAddition(x0, y0, midx, midy, midz, x[1], y[1], z[1]);//Q5


	tauProjective(x[12], y[12], z[12], middlex, middley, middlez);//tau Q_27
	tauProjective(middlex, middley, middlez, midx, midy, midz);//tau^2 Q_27

	mixedAddition(x0, middleNegative, midx, midy, midz, x[4], y[4], z[4]);//Q11
	mixedAddition(x0, y0, midx, midy, midz, x[5], y[5], z[5]);//Q13


	tauProjective(midx, midy, midz, middlex, middley, middlez);//tau^3 Q_27
	add2(midy, midz, one);


	mixedAddition(x0, middleNegative, middlex, one, middlez, x[2], y[2], z[2]);//Q7
	mixedAddition(x0, y0, middlex, one, middlez, x[3], y[3], z[3]);//Q9


	add2(midy, midz, one);
	Addition(x[12], y[12], z[12], midx, one, midz, x[6], y[6], z[6]);//Q15

	Addition(x[13], y[13], z[13], midx, one, midz, x[6], y[6], z[6]);//Q17



	tauProjective(x[13], y[13], z[13], middlex, middley, middlez);//tau Q_29
	tauProjective(middlex, middley, middlez, midx, midy, midz);//tau^2 Q_29

	mixedAddition(x0, y0, midx, midy, midz, x[9], y[9], z[9]);//Q21

	tauProjective(x[11], y[11], z[11], middlex, middley, middlez);//tau Q_25
	tauProjective(middlex, middley, middlez, midx, midy, midz);//tau^2 Q_25



	Addition(x[12], y[12], z[12], midx, midy, midz, x[14], y[14], z[14]);//Q31



	add2(y[0], z[0], midy);
	tauProjective(x[0], midy, z[0], middlex, middley, middlez);//-tau Q_3
	tauProjective(middlex, middley, middlez, midx, midy, midz);//-tau^2 Q_3



	mixedAddition(x0, middleNegative, midx, midy, midz, x[8], y[8], z[8]);//Q19


	for (int i = 0; i < 15; i++)
	{
		copy(x[i], P[2 * i + 3].X);
		copy(y[i], P[2 * i + 3].Y);
		copy(z[i], P[2 * i + 3].Z);
	}


	return 0;
}

int HMV(affinePoint P0, projectPoint *P, int w)
{//HMV pre-computation for w=4,5,6. Only works when a=0


	if (w == 6)
	{
		HMV6(P0, P);
		return 0;
	}

	big x[7], y[7], z[7], x0, y0;
	x0 = mirvar(0);
	y0 = mirvar(0);
	copy(P0.x, x0);
	copy(P0.y, y0);
	for (int i = 0; i < 7; i++)
	{
		x[i] = mirvar(0);
		y[i] = mirvar(0);
		z[i] = mirvar(0);
	}
	//big *y=(big *)mr_alloc(7,sizeof(big));
	//big *z=(big *)mr_alloc(7,sizeof(big));
	big middleNegative = mirvar(0);//the negative lambda of lambda coordinate 
	big middlex = mirvar(0), middley = mirvar(0), middlez = mirvar(1), midx = mirvar(0), midy = mirvar(0), midz = mirvar(1);

	big  one = mirvar(1);
	add2(y0, one, middleNegative);

	tauAffine(x0, y0, middlex, middley);//tau
	tauAffine(middlex, middley, midx, midy);//tau^2
	mixedAddition(x0, middleNegative, midx, midy, P[3].X, P[3].Y, P[3].Z);//Q3


	mixedAddition(x0, y0, midx, midy, P[5].X, P[5].Y, P[5].Z);//Q5

	copy(P[5].X, x[1]);
	copy(P[5].Y, y[1]);
	copy(P[5].Z, z[1]);


	tauAffine(midx, midy, middlex, middley);//tau^3
	add2(middley, one, middley);
	mixedAddition(x0, middleNegative, middlex, middley, P[7].X, P[7].Y, P[7].Z);//Q7


	if (w == 5)
	{

		tauProjective(x[1], y[1], z[1], middlex, middley, middlez);//tau Q_5
		tauProjective(middlex, middley, middlez, midx, midy, midz);//tau^2 Q_5

		add2(y[1], z[1], middley);

		Addition(x[1], y[1], z[1], midx, middley, midz, P[15].X, P[15].Y, P[15].Z);//Q15

		add2(midy, midz, one);
		mixedAddition(x0, middleNegative, midx, one, midz, P[11].X, P[11].Y, P[11].Z);//Q11

		tauProjective(midx, one, midz, middlex, middley, middlez);//-tau^3 Q_5
		mixedAddition(x0, y0, middlex, middley, middlez, P[9].X, P[9].Y, P[9].Z);//Q9			
		mixedAddition(x0, y0, midx, one, midz, P[13].X, P[13].Y, P[13].Z);//Q13
	}

	return 0;

}




int xu(big x0, big y0, int w = 4)
{//xu'pre-computation for w=4,5. Only works when a=0
	//big *x=(big *)mr_alloc(7,sizeof(big));
	big x[7], y[7], z[7];
	for (int i = 0; i < 7; i++)
	{
		x[i] = mirvar(0);
		y[i] = mirvar(0);
		z[i] = mirvar(0);
	}

	big middleNegative = mirvar(0);//the negative lambda of lambda coordinate 
	big middlex = mirvar(0), middley = mirvar(0), middlez = mirvar(1), midx = mirvar(0), midy = mirvar(0), midz = mirvar(1), middle = mirvar(0);

	big  one = mirvar(1);
	add2(y0, one, middleNegative);

	tauAffine(x0, y0, middlex, middley);//tau P
	tauAffine(middlex, middley, midx, midy);//tau^2 P


	AffinePMinusTauP(x0, y0, x[1], middle, z[1]);//Q5
	add2(middle, z[1], y[1]);
	AffinePAddTauP(x0, y0, x[2], y[2], z[2]);//Q7
	ProjectivePMinusTauP(x[2], y[2], z[2], x[0], y[0], z[0]);//Q3

	//mixedAddition(x0, middleNegative, midx,midy,midz, x[0], y[0], z[0]);//Q3
	//mixedAddition(x0, y0, middlex,middley,middlez, x[1], y[1], z[1]);//Q5
	//mixedAddition(x0, middleNegative, middlex,middley,middlez, x[2], y[2], z[2]);//Q7




	if (w == 5)
	{
		mixedAddition(middlex, middley, x[0], y[0], z[0], x[3], y[3], z[3]);//Q9
		mixedAddition(middlex, middley, x[1], y[1], z[1], x[4], y[4], z[4]);//Q11
		mixedAddition(middlex, middley, x[2], y[2], z[2], x[5], y[5], z[5]);//Q13
		mixedAddition(middlex, middley, x[4], y[4], z[4], x[6], middley, z[6]);//Q15
		add2(middley, z[6], y[6]);
	}

	return 0;

}




int xu(affinePoint P0, projectPoint *P, int w = 4, int a = 0)
{//Trost and xu's pre-computation. work for w=4,5.  

	//big *x=(big *)mr_alloc(7,sizeof(big));
	big x[7], y[7], z[7], x0, y0;
	big X3 = mirvar(0);
	x0 = mirvar(0);
	y0 = mirvar(0);
	copy(P0.x, x0);
	copy(P0.y, y0);
	for (int i = 0; i < 7; i++)
	{
		x[i] = mirvar(0);
		y[i] = mirvar(0);
		z[i] = mirvar(0);
	}

	big middleNegative = mirvar(0);//the negative lambda of lambda coordinate 
	big middlex = mirvar(0), middley = mirvar(0), middlez = mirvar(1), midx = mirvar(0), midy = mirvar(0), midz = mirvar(1), middle = mirvar(0);

	big  one = mirvar(1);
	//add2(y0, one, middleNegative);
	if (a == 0)
	{
		add2(y0, one, middleNegative);
	}
	tauAffine(x0, y0, middlex, middley);//tau P

	tauAffine(middlex, middley, midx, midy);//tau^2 P

	AffinePMinusTauP(x0, y0, x[1], middle, z[1], X3);//Q5

	//AffinePMinusTauP(x0, y0,   x[1], middle, z[1]);//Q5
	add2(middle, z[1], y[1]);


	AffinePAddTauP(X3, x0, y0, x[2], y[2], z[2]);//Q7
	//AffinePAddTauP(x0, y0,  x[2], y[2], z[2]);//Q7
	ProjectivePMinusTauP(x[2], y[2], z[2], x[0], y[0], z[0]);//Q3





	if (w == 5)
	{
		mixedAddition(middlex, middley, x[0], y[0], z[0], x[3], y[3], z[3]);//Q9
		mixedAddition(middlex, middley, x[1], y[1], z[1], x[4], y[4], z[4]);//Q11
		mixedAddition(middlex, middley, x[2], y[2], z[2], x[5], y[5], z[5]);//Q13
		mixedAddition(middlex, middley, x[4], y[4], z[4], x[6], middley, z[6]);//Q15
		add2(middley, z[6], y[6]);
	}

	for (int i = 0; i < 7; i++)
	{
		copy(x[i], P[2 * i + 3].X);
		copy(y[i], P[2 * i + 3].Y);
		copy(z[i], P[2 * i + 3].Z);
	}
	return 0;

}








int xu6(affinePoint P0, projectPoint *P, int a = 0)
{//Trost and xu's pre-computation. work for w=6.
	//big *x=(big *)mr_alloc(7,sizeof(big));
	big x0, y0;
	big X3 = mirvar(0);
	x0 = mirvar(0);
	y0 = mirvar(0);
	copy(P0.x, x0);
	copy(P0.y, y0);
	//big *x=(big *)mr_alloc(7,sizeof(big));



	big middleNegative = mirvar(0);//the negative lambda of lambda coordinate 
	big middlex = mirvar(0), middley = mirvar(0), middlez = mirvar(1), midx = mirvar(0), midy = mirvar(0), midz = mirvar(1), middle = mirvar(0);
	big middleyNegative = mirvar(0);
	big  one = mirvar(1);
	add2(y0, one, middleNegative);

	tauAffine(x0, y0, middlex, middley);//tau P
	add2(middley, one, middleyNegative);//-tau P

	if (a == 1)
	{//change tau P and -tau P
		tauAffine(x0, y0, middlex, middleyNegative);//tau P
		add2(middleyNegative, one, middley);//-tau P
	}
	tauAffine(middlex, middleyNegative, midx, midy);//-tau^2 P






	//AffinePMinusTauP(x0, y0,  x[12], y[12], z[12]);//Q27	//

	AffinePMinusTauP(x0, y0, P[27].X, P[27].Y, P[27].Z, X3);//Q27	//

	//AffinePAddTauP(x0, y0,  x[11], middle, z[11]);//Q25
	AffinePAddTauP(X3, x0, y0, P[25].X, middle, P[25].Z);//Q25
	add2(middle, P[25].Z, P[25].Y);	//
	ProjectivePMinusTauP(P[27].X, P[27].Y, P[27].Z, P[29].X, P[29].Y, P[29].Z);//Q29





	//mixedAddition(x0, y0, midx,midy,midz, x[13], y[13], z[13]);//Q29
	//mixedAddition(x0, y0, middlex,middleyNegative,middlez, x[12], y[12], z[12]);//Q27
	//mixedAddition(x0, middleNegative, middlex,middleyNegative,middlez, x[11], y[11], z[11]);//Q25


	mixedAddition(middlex, middley, P[29].X, P[29].Y, P[29].Z, P[3].X, P[3].Y, P[3].Z);//Q3
	add2(P[29].X, P[29].Y, one);
	mixedAddition(middlex, middley, P[29].X, one, P[29].Z, P[9].X, P[9].Y, P[9].Z);//Q9
	mixedAddition(midx, midy, P[3].X, P[3].Y, P[3].Z, P[31].X, P[31].Y, P[31].Z);//Q31


	mixedAddition(middlex, middley, P[31].X, P[31].Y, P[31].Z, P[5].X, P[5].Y, P[5].Z);//Q5
	add2(P[31].Y, P[31].Z, one);

	mixedAddition(middlex, middley, P[31].X, one, P[31].Z, P[7].X, P[7].Y, P[7].Z);//Q7

	add2(P[27].Y, P[27].Z, one);
	mixedAddition(middlex, middley, P[27].X, one, P[27].Z, P[11].X, P[11].Y, P[11].Z);//Q11
	add2(P[25].Y, P[25].Z, one);
	mixedAddition(middlex, middley, P[25].X, one, P[25].Z, P[13].X, P[13].Y, P[13].Z);//Q13
	add2(P[11].Y, P[11].Z, one);
	mixedAddition(middlex, middleyNegative, P[11].X, one, P[11].Z, P[15].X, P[15].Y, P[15].Z);//Q15

	add2(P[9].Y, P[9].Z, one);
	mixedAddition(middlex, middleyNegative, P[9].X, one, P[9].Z, P[17].X, P[17].Y, P[17].Z);//Q17
	add2(P[7].Y, P[7].Z, one);
	mixedAddition(middlex, middleyNegative, P[7].X, one, P[7].Z, P[19].X, P[19].Y, P[19].Z);//Q19

	add2(P[17].Y, P[17].Z, one);
	mixedAddition(middlex, middley, P[17].X, one, P[17].Z, P[21].X, P[21].Y, P[21].Z);//Q21

	add2(P[3].Y, P[3].Z, one);
	mixedAddition(middlex, middleyNegative, P[3].X, one, P[3].Z, P[23].X, P[23].Y, P[23].Z);//Q23



	return 0;

}




int xuPreComputation(affinePoint P0, projectPoint *P, int w, int a)
{//Trost and xu's pre-computation. work for w=4,5,6. Only works when a=0
	if (w == 6)
	{
		xu6(P0, P, a);
	}
	else if (w < 6)
	{
		if (w > 3)
		{
			xu(P0, P, w);
		}
	}
	return 0;
}


int optimal(affinePoint P0, projectPoint *P, int w = 4)
{//our pre-computation for w=4,5
	//big *x=(big *)mr_alloc(7,sizeof(big));
	big x[7], y[7], z[7], x0, y0;
	x0 = mirvar(0);
	y0 = mirvar(0);
	copy(P0.x, x0);
	copy(P0.y, y0);
	for (int i = 0; i < 7; i++)
	{
		x[i] = mirvar(0);
		y[i] = mirvar(0);
		z[i] = mirvar(0);
	}

	big middleNegative = mirvar(0);//the negative lambda of lambda coordinate 
	big middlex = mirvar(0), middley = mirvar(0), middlez = mirvar(1), midx = mirvar(0), midy = mirvar(0), midz = mirvar(1), middle = mirvar(0);

	big  one = mirvar(1);



	//barTau

	affineBarTau(x0, y0, x[1], middley, z[1], x[2], midy, z[2]);//Q5,Q7
	add2(middley, z[1], y[1]);
	add2(midy, z[2], y[2]);
	ReRebarTau(x[2], midy, z[2], z[1], x[0], y[0], z[0]);//Q3


	//RebarTau
	//ReRebarTau



	if (w == 5)
	{


		tauAffine(x0, y0, middlex, middley);//tau P

		ReRebarTau(x[0], y[0], z[0], z[2], x[6], middle, z[6]);//Q15		
		add2(middle, z[6], y[6]);
		mixedAddition(middlex, middley, x[1], y[1], z[1], x[4], y[4], z[4]);//Q11

		barTau(x[4], y[4], z[4], x[3], y[3], z[3], middle);//Q9
		RebarTau(x[3], y[3], z[3], middle, x[5], middley, z[5]);//Q13 

		add2(middley, z[5], y[5]);
	}


	for (int i = 0; i < 7; i++)
	{
		copy(x[i], P[2 * i + 3].X);
		copy(y[i], P[2 * i + 3].Y);
		copy(z[i], P[2 * i + 3].Z);
	}
	return 0;

}




int optimal6(affinePoint P0, projectPoint *P, int w = 4)
{//our pre-computation for w=6
	//big *x=(big *)mr_alloc(7,sizeof(big));
	big x[15], y[15], z[15], x0, y0;
	x0 = mirvar(0);
	y0 = mirvar(0);
	copy(P0.x, x0);
	copy(P0.y, y0);
	for (int i = 0; i < 15; i++)
	{
		x[i] = mirvar(0);
		y[i] = mirvar(0);
		z[i] = mirvar(0);
	}

	big middleNegative = mirvar(0);//the negative lambda of lambda coordinate 
	big middlex = mirvar(0), middley = mirvar(0), middlez = mirvar(1), midx = mirvar(0), midy = mirvar(0), midz = mirvar(1), middle = mirvar(0);

	big  one = mirvar(1);



	//barTau

	affineBarTau(x0, y0, x[12], y[12], z[12], x[11], y[11], z[11]);//Q27,Q25

	ReRebarTau(x[11], y[11], z[11], z[12], x[13], middle, z[13]);//Q29
	add2(middle, z[13], y[13]);

	ReRebarTau(x[13], y[13], z[13], z[11], x[6], y[6], z[6]);//Q15
	ReRebarTau(x[6], y[6], z[6], z[13], x[9], y[9], z[9]);//Q21




	//RebarTau
	//ReRebarTau
	tauAffine(x0, y0, middlex, middley);//tau P

	PpmQ(middlex, middley, x[13], y[13], z[13], x[0], y[0], z[0], x[3], y[3], z[3]);//Q3,Q9
	barTau(x[3], y[3], z[3], x[5], midy, z[5], middle);//Q13
	add2(midy, z[5], y[5]);
	RebarTau(x[5], y[5], z[5], middle, x[14], y[14], z[14]);//Q31 


	barTau(x[0], y[0], z[0], x[7], y[7], z[7], middle);//Q17

	RebarTau(x[7], y[7], z[7], middle, x[4], y[4], z[4]);//Q11 


	add2(y[6], z[6], middle);
	mixedAddition(middlex, middley, x[6], middle, z[6], x[10], y[10], z[10]);

	barTau(x[10], y[10], z[10], x[8], midy, z[8], middle);//Q19
	add2(midy, z[8], y[8]);



	PpmQ(middlex, middley, x[14], y[14], z[14], x[1], y[1], z[1], x[2], y[2], z[2]);//Q5,Q7




	for (int i = 0; i < 15; i++)
	{
		copy(x[i], P[2 * i + 3].X);
		copy(y[i], P[2 * i + 3].Y);
		copy(z[i], P[2 * i + 3].Z);
	}



	return 0;

}

int OurPreComputation(affinePoint P0, projectPoint *P, int w)
{
	if (w == 6)
	{
		optimal6(P0, P);
	}
	else if (w < 6)
	{
		if (w > 3)
		{
			optimal(P0, P, w);
		}
	}
	return 0;
}

//int generateRandomWindowTauNAF(int *a, int &length, int w)
//{//generate a random window tauNAF
//	int modNumber=(1<<w);
//	int judge=modNumber/2;
//	int value=rand()%judge;
//
//	return 0;
//
//}

int windowNAF(Big n, int *a, int &length, int w)
{//generate a window NAF. In fact, window tauNAF has a density of 1/(w+1) is the same as window NAF and {\pm 1,\pm 3, ...,\pm 2^{w-1}-1} are randomly chosen in the non-zero coefficients.
	int modNumber = (1 << w);
	int judge = modNumber / 2;

	int i = 0;
	int x = 0;

	while (n > 0)
	{
		if (n % 2 == 0)
		{
			a[i] = 0;
			i++;
			n = n / 2;
			continue;
		}
		a[i] = n%modNumber;
		//n=n/modNumber;
		if (a[i] > judge)
		{
			a[i] = a[i] - modNumber;

		}
		n = n - a[i];
		n = n / 2;
		i++;
	}
	length = i;
	return 0;
}

int initialVauleC(numberC c[], int cases, int w, int a)
{//cases 3: Solinas
 //case 2:	HMV
 //case 1: Xu
 //case 0:our

	if (cases == 0)
	{
		if (w < 6)
		{
			c[3].g = -3;
			c[3].h = -1;
			c[5].g = -1;
			c[5].h = -1;
			c[7].g = 1;
			c[7].h = -1;
			c[9].g = 3;
			c[9].h = -1;
			c[11].g = -1;
			c[11].h = -2;
			c[13].g = -5;
			c[13].h = -3;
			c[15].g = 1;
			c[15].h = 3;
		}
		else if (w == 6)
		{
			c[3].g = 3;
			c[3].h = 0;
			c[5].g = 5;
			c[5].h = 0;
			c[7].g = 5;
			c[7].h = -5;
			c[9].g = -3;
			c[9].h = -2;
			c[11].g = -3;
			c[11].h = 3;
			c[13].g = 1;
			c[13].h = -2;
			c[15].g = 1;
			c[15].h = 3;
			c[17].g = 3;
			c[17].h = 3;
			c[19].g = -7;
			c[19].h = 1;
			c[21].g = -5;
			c[21].h = 1;
			c[23].g = -1;
			c[23].h = -4;
			c[25].g = -1;
			c[25].h = 1;
			c[27].g = 1;
			c[27].h = 1;
			c[29].g = 3;
			c[29].h = 1;
			c[31].g = -7;
			c[31].h = -1;
		}
		else if (w == 7)
		{
			c[3].g = 3;
			c[3].h = 0;
			c[5].g = 5;
			c[5].h = -7;
			c[7].g = -7;
			c[7].h = 3;
			c[9].g = -1;
			c[9].h = -7;
			c[11].g = -3;
			c[11].h = 3;
			c[13].g = -1;
			c[13].h = 3;
			c[15].g = 1;
			c[15].h = 3;
			c[17].g = 3;
			c[17].h = 3;
			c[19].g = 5;
			c[19].h = 3;
			c[21].g = 7;
			c[21].h = 3;
			c[23].g = -1;
			c[23].h = -4;
			c[25].g = 1;
			c[25].h = -4;
			c[27].g = -11;
			c[27].h = -1;
			c[29].g = -9;
			c[29].h = -1;
			c[31].g = -7;
			c[31].h = -1;
			c[33].g = -5;
			c[33].h = -1;
			c[35].g = -3;
			c[35].h = -1;
			c[37].g = -1;
			c[37].h = -1;
			c[39].g = 1;
			c[39].h = -1;
			c[41].g = 3;
			c[41].h = -1;
			c[43].g = 5;
			c[43].h = -1;
			c[45].g = 7;
			c[45].h = -1;
			c[47].g = 9;
			c[47].h = -1;
			c[49].g = -3;
			c[49].h = 2;
			c[51].g = -1;
			c[51].h = 2;
			c[53].g = 1;
			c[53].h = 2;
			c[55].g = 3;
			c[55].h = 2;
			c[57].g = 5;
			c[57].h = 2;
			c[59].g = -3;
			c[59].h = -5;
			c[61].g = -1;
			c[61].h = -5;
			c[63].g = 1;
			c[63].h = -5;

		}
		if (a == 1)
		{
			for (int i = 1; i < 65; i++)
			{
				c[i].h = -c[i].h;
			}
		}
	}
	else if (cases == 1)
	{
		if (w < 6)
		{
			c[3].g = -3;
			c[3].h = -1;
			c[5].g = -1;
			c[5].h = -1;
			c[7].g = 1;
			c[7].h = -1;
			c[9].g = -3;
			c[9].h = -2;
			c[11].g = -1;
			c[11].h = -2;
			c[13].g = 1;
			c[13].h = -2;
			c[15].g = 1;
			c[15].h = 3;
		}
		else if (w == 6)
		{
			c[3].g = 3;
			c[3].h = 0;
			c[5].g = 5;
			c[5].h = 0;
			c[7].g = -5;
			c[7].h = -2;
			c[9].g = -3;
			c[9].h = -2;
			c[11].g = -1;
			c[11].h = -2;
			c[13].g = 1;
			c[13].h = -2;
			c[15].g = 1;
			c[15].h = 3;
			c[17].g = 3;
			c[17].h = 3;
			c[19].g = 5;
			c[19].h = 3;
			c[21].g = -3;
			c[21].h = -4;
			c[23].g = -3;
			c[23].h = 1;
			c[25].g = -1;
			c[25].h = 1;
			c[27].g = 1;
			c[27].h = 1;
			c[29].g = 3;
			c[29].h = 1;
			c[31].g = 5;
			c[31].h = 1;
		}
		if (a == 1)
		{
			for (int i = 1; i < 65; i++)
			{
				c[i].h = -c[i].h;
			}
		}
	}
	else if (cases == 2)
	{
		if (w < 6)
		{
			c[3].g = -3;
			c[3].h = -1;
			c[5].g = -1;
			c[5].h = -1;
			c[7].g = 1;
			c[7].h = -1;
			c[9].g = -3;
			c[9].h = -2;
			c[11].g = -1;
			c[11].h = -2;
			c[13].g = 1;
			c[13].h = -2;
			c[15].g = 1;
			c[15].h = 3;
		}
		else if (w == 6)
		{
			c[3].g = 3;
			c[3].h = 0;
			c[5].g = 5;
			c[5].h = 0;
			c[7].g = -5;
			c[7].h = -2;
			c[9].g = -3;
			c[9].h = -2;
			c[11].g = -1;
			c[11].h = -2;
			c[13].g = 1;
			c[13].h = -2;
			c[15].g = 1;
			c[15].h = 3;
			c[17].g = 3;
			c[17].h = 3;
			c[19].g = 5;
			c[19].h = 3;
			c[21].g = -3;
			c[21].h = -4;
			c[23].g = -3;
			c[23].h = 1;
			c[25].g = -1;
			c[25].h = 1;
			c[27].g = 1;
			c[27].h = 1;
			c[29].g = 3;
			c[29].h = 1;
			c[31].g = 5;
			c[31].h = 1;
		}
		if (a == 1)
		{
			for (int i = 1; i < 65; i++)
			{
				c[i].h = -c[i].h;
			}
		}
	}
	else if (cases == 3)
	{
		c[3].g = -3;
		c[3].h = -1;
		c[5].g = -1;
		c[5].h = -1;
		c[7].g = 1;
		c[7].h = -1;
		c[9].g = -3;
		c[9].h = -2;
		c[11].g = -1;
		c[11].h = -2;
		c[13].g = 1;
		c[13].h = -2;
		c[15].g = 1;
		c[15].h = 3;
		if (a == 1)
		{
			for (int i = 1; i < 65; i++)
			{
				c[i].h = -c[i].h;
			}
		}
	}
	return 0;
}



int windowTauNAF(Big n1, Big n2, numberC c[], int *a, int &length, int w, int parameterA)
{//generate a window NAF. In fact, window tauNAF has a density of 1/(w+1) is the same as window NAF and {\pm 1,\pm 3, ...,\pm 2^{w-1}-1} are randomly chosen in the non-zero coefficients.
	//n1+n2\tau,n1 and n2 are both odd
	int modNumber = (1 << w);
	int judge = modNumber / 2;
	Big middle = 0;
	int u = 1;
	if (parameterA == 0)
	{
		u = -1;
	}
	int s = 0;
	if (w < 2 || w>7)
	{
		cout << "Error" << endl;
		return 0;
	}
	if (w == 2)
	{
		s = 2;
	}
	else if (w < 6)
	{
		s = 6;
	}
	else
	{
		s = 38;
	}

	s = u*s;

	int i = 0;
	int x = 0;

	while ((n1 + n2*s) != 0)
	{
		if (n1 % 2 == 0)
		{
			a[i] = 0;
			i++;
			middle = n1 / 2;
			n1 = n2 + (middle)*u;
			n2 = -1 * middle;
			continue;
		}
		a[i] = (n1 + n2*s) % modNumber;

		if (a[i] < 0)
		{
			a[i] = a[i] + modNumber;

		}
		//n=n/modNumber;
		if (a[i] > judge)
		{
			a[i] = a[i] - modNumber;

		}
		int flag = 1;
		int indexMiddle = a[i];
		if (indexMiddle < 0)
		{
			flag = -1;
			indexMiddle = -indexMiddle;
		}
		n1 = n1 - flag*c[indexMiddle].g;
		n2 = n2 - flag*c[indexMiddle].h;
		middle = n1 / 2;
		n1 = n2 + (middle)*u;
		n2 = -1 * middle;
		i++;
	}
	length = i;
	return 0;
}

int windowNAFregular(Big n, int *a, int &length, int w)
{//return a regular  window NAF. n is odd, if n is even then n\right n+1
	int modNumber = (1 << (w));

	if (n % 2 == 0)
	{
		n = n + 1;
	}


	int i = 0;
	int x = 0;
	Big middle = 0;

	while (n > 0)
	{

		a[i] = n%modNumber;
		if (a[i] < 0)
		{
			a[i] = a[i] + modNumber;
		}
		middle = n - a[i];
		if (middle == 0)
		{
			break;
		}
		if ((middle % 2) == 0)
		{
			a[i] = a[i] - modNumber;
		}

		n = n - a[i];
		n = n / 2;
		i++;
	}
	length = i;
	return 0;
}

int windowTauNAFregular(Big n1, Big n2, numberC c[], int *a, int &length, int w, int parameterA)
{//generate a window NAF. In fact, window tauNAF has a density of 1/(w+1) is the same as window NAF and {\pm 1,\pm 3, ...,\pm 2^{w-1}-1} are randomly chosen in the non-zero coefficients.
 //n1+n2\tau,n1 and n2 are both odd
	int modNumber = (1 << (w));
	int judge = modNumber / 2;
	Big middle = 0;
	Big middleCompare = 0;
	Big middleCompare1 = 0;
	int u = 1;
	if (parameterA == 0)
	{
		u = -1;
	}
	int s = 0;
	if (w < 2 || w>7)
	{
		cout << "Error" << endl;
		return 0;
	}
	if (w == 2)
	{
		s = 2;
	}
	else if (w < 6)
	{
		s = 6;
	}
	else
	{
		s = 38;
	}

	s = u*s;

	int i = 0;
	int x = 0;

	while ((n1 + n2*s) != 0)
	{
		if (n1 % 2 == 0)
		{
			a[i] = 0;
			i++;
			middle = n1 / 2;
			n1 = n2 + (middle)*u;
			n2 = -1 * middle;
			continue;
		}
		middleCompare = n1 + n2*s;
		middleCompare1 = n1*n1 + 2 * n2*n2 + n1*n2*u;

		if ((middleCompare1<modNumber) && (middleCompare1>-1 * modNumber))
		{
			// cout << n1 << "+" << n2 << "\\tau" << endl;
			a[i] = middleCompare%modNumber;
			if (a[i] < 0)
			{
				a[i] = a[i] + modNumber;
			}
			if (a[i] > judge)
			{
				a[i] = a[i] - modNumber;
			}

		}
		else
		{
			a[i] = (middleCompare) % modNumber;
			//n=n/modNumber;
			if (a[i] < 0)
			{
				a[i] = a[i] + modNumber;

			}
			a[i] = a[i] - judge;
		}

		int flag = 1;
		int indexMiddle = a[i];
		if (indexMiddle < 0)
		{
			flag = -1;
			indexMiddle = -indexMiddle;
		}
		n1 = n1 - flag*c[indexMiddle].g;
		n2 = n2 - flag*c[indexMiddle].h;
		middle = n1 / 2;
		n1 = n2 + (middle)*u;
		n2 = -1 * middle;
		i++;
	}
	length = i;
	return 0;
}


int scalarMultiplication(affinePoint P0, affinePoint negative, projectPoint *P, projectPoint negativePi, int *a, int aLength, projectPoint Q, int w, int parameterA)
{ //using pre-computation in lambda-projective coordinates.
	tau(P0, Q); //Q=2P0;


	projectPoint middleP;
	middleP.X = mirvar(0);
	middleP.Y = mirvar(0);
	middleP.Z = mirvar(0);
	int indexMiddle = a[aLength - 1];
	if (indexMiddle >= 0)
	{
		copy(P[indexMiddle].X, Q.X);
		copy(P[indexMiddle].Y, Q.Y);

		copy(P[indexMiddle].Z, Q.Z);
	}
	else
	{
		indexMiddle = -indexMiddle;
		copy(P[indexMiddle].X, Q.X);
		copy(P[indexMiddle].Z, Q.Z);
		add2(P[indexMiddle].Y, Q.Z, Q.Y);



	}



	for (int i = aLength - 2; i > -1; i--)
	{
		if (a[i] > 0)
		{
			if (a[i] == 1)
			{
				mixedAddition(P0, Q, middleP);
			}
			else
			{
				Addition(P[a[i]], Q, middleP);
			}
		}
		else if (a[i] < 0)
		{
			if (a[i] == -1)
			{
				mixedAddition(negative, Q, middleP);//Q-P0
			}
			else
			{
				copy(P[-a[i]].X, negativePi.X);
				copy(P[-a[i]].Z, negativePi.Z);
				add2(P[-a[i]].Y, negativePi.Z, negativePi.Y);
				Addition(negativePi, Q, middleP);
			}
		}
		else
		{
			copy(Q.X, middleP.X);
			copy(Q.Y, middleP.Y);

			copy(Q.Z, middleP.Z);
		}

		//tau(middleP, Q);
		modsquare2(middleP.X, Q.X);
		modsquare2(middleP.Y, Q.Y);
		modsquare2(middleP.Z, Q.Z);
	}
	return 0;
}


int MontgomeryTrick(projectPoint *P, int w)
{//montgomery Trick

	big a[33];
	big b[33];
	big c[33];
	int n = (1 << (w - 2)) - 1;
	big one = mirvar(1), d = mirvar(1);
	for (int i = 0; i < 33; i++)
	{
		a[i] = mirvar(0);
		b[i] = mirvar(0);
		c[i] = mirvar(1);

	}
	for (int i = 1; i < n + 1; i++)
	{
		copy(P[2 * i + 1].Z, a[i]);
	}
	copy(a[1], c[1]);

	for (int i = 2; i < n + 1; i++)
	{
		modmult2(a[i], c[i - 1], c[i]);
	}

	inverse2(c[n], d);

	for (int i = n; i > 1; i--)
	{
		modmult2(d, c[i - 1], b[i]);
		modmult2(a[i], d, d);
	}
	copy(d, b[1]);


	for (int i = 1; i < n + 1; i++)
	{
		modmult2(b[i], P[2 * i + 1].X, P[2 * i + 1].X);
		modmult2(b[i], P[2 * i + 1].Y, P[2 * i + 1].Y);
		copy(one, P[2 * i + 1].Z);
	}

	//Big out=0;
	//for(int i=1;i<n+1;i++)
	//{	
	//	cout <<"Q_"<<2*i+1<<endl;
	//	out =P[2*i+1].X;
	//	cout <<"X:"<<out <<endl;
	//	out =P[2*i+1].Y;
	//	cout <<"Y:"<<out <<endl;
	//	out =P[2*i+1].Z;
	//	cout <<"Z:"<<out <<endl;
	//	 
	//}

	return 0;
}

int MontgomeryTrick(projectPointLD *P, int w)
{//montgomery Trick

	big a[33];
	big b[33];
	big c[33];
	int n = (1 << (w - 2)) - 1;
	big one = mirvar(1), d = mirvar(1), R1= mirvar(0);
	for (int i = 0; i < 33; i++)
	{
		a[i] = mirvar(0);
		b[i] = mirvar(0);
		c[i] = mirvar(1);

	}
	for (int i = 1; i < n + 1; i++)
	{
		copy(P[2 * i + 1].Z, a[i]);
	}
	copy(a[1], c[1]);

	for (int i = 2; i < n + 1; i++)
	{
		modmult2(a[i], c[i - 1], c[i]);
	}

	inverse2(c[n], d);

	for (int i = n; i > 1; i--)
	{
		modmult2(d, c[i - 1], b[i]);
		modmult2(a[i], d, d);
	}
	copy(d, b[1]);


	for (int i = 1; i < n + 1; i++)
	{
		modmult2(b[i], P[2 * i + 1].X, P[2 * i + 1].X);
		modsquare2(b[i], R1);
		modmult2(R1, P[2 * i + 1].Y, P[2 * i + 1].Y);
		copy(one, P[2 * i + 1].Z);
	}

	//Big out=0;
	//for(int i=1;i<n+1;i++)
	//{	
	//	cout <<"Q_"<<2*i+1<<endl;
	//	out =P[2*i+1].X;
	//	cout <<"X:"<<out <<endl;
	//	out =P[2*i+1].Y;
	//	cout <<"Y:"<<out <<endl;
	//	out =P[2*i+1].Z;
	//	cout <<"Z:"<<out <<endl;
	//	 
	//}

	return 0;
}

int scalarMultiplicationMontgomery(affinePoint P0, affinePoint negative, projectPoint *P, projectPoint negativePi, int *a, int aLength, projectPoint Q, int w, int parameterA)
{ //scalarMultiplicationMontgomery use pre-computations in lambda coordinates
	tau(P0, Q); //Q=2P0;


	MontgomeryTrick(P, w);
	projectPoint middleP;
	middleP.X = mirvar(0);
	middleP.Y = mirvar(0);
	middleP.Z = mirvar(0);

	int indexMiddle = a[aLength - 1];
	if (indexMiddle >= 0)
	{
		copy(P[indexMiddle].X, Q.X);
		copy(P[indexMiddle].Y, Q.Y);

		copy(P[indexMiddle].Z, Q.Z);
	}
	else
	{
		indexMiddle = -indexMiddle;
		copy(P[indexMiddle].X, Q.X);
		copy(P[indexMiddle].Z, Q.Z);
		add2(P[indexMiddle].Y, Q.Z, Q.Y);



	}



	for (int i = aLength - 2; i > -1; i--)
	{
		if (a[i] > 0)
		{
			//mixedAddition(P[a[i]], Q, middleP);
			if (a[i] == 1)
			{
				mixedAddition(P0, Q, middleP);
			}
			else
			{
				mixedAddition(P[a[i]], Q, middleP);
			}
		}
		else if (a[i] < 0)
		{
			if (a[i] == -1)
			{
				mixedAddition(negative, Q, middleP);//Q-P0
			}
			else
			{
				copy(P[-a[i]].X, negativePi.X);
				copy(P[-a[i]].Z, negativePi.Z);
				add2(P[-a[i]].Y, negativePi.Z, negativePi.Y);
				mixedAddition(negativePi, Q, middleP);
			}
		}
		else
		{
			copy(Q.X, middleP.X);
			copy(Q.Y, middleP.Y);

			copy(Q.Z, middleP.Z);
		}

		tau(middleP, Q);
	}
	return 0;
}

int scalarMultiplicationMontgomery(affinePointLD P0, affinePointLD negative, projectPointLD *P, projectPointLD negativePi, int *a, int aLength, projectPointLD Q, int w, int parameterA)
{ //scalarMultiplicationMontgomery use pre-computations in lambda coordinates
	tau(P0, Q); //Q=2P0;


	MontgomeryTrick(P, w);
	projectPointLD middleP;
	middleP.X = mirvar(0);
	middleP.Y = mirvar(0);
	middleP.Z = mirvar(0);

	big R1 = mirvar(0);

	int indexMiddle = a[aLength - 1];
	if (indexMiddle >= 0)
	{
		copy(P[indexMiddle].X, Q.X);
		copy(P[indexMiddle].Y, Q.Y);
		copy(P[indexMiddle].Z, Q.Z);
	}
	else
	{
		indexMiddle = -indexMiddle;
		copy(P[indexMiddle].X, Q.X);
		copy(P[indexMiddle].Z, Q.Z);
		modmult2(Q.X, Q.Z, R1);
		add2(P[indexMiddle].Y, R1, Q.Y);
	}



	for (int i = aLength - 2; i > -1; i--)
	{
		if (a[i] > 0)
		{
			//mixedAddition(P[a[i]], Q, middleP);
			if (a[i] == 1)
			{
				mixedAddition(P0, Q, middleP);
			}
			else
			{
				mixedAddition(P[a[i]], Q, middleP);
			}
		}
		else if (a[i] < 0)
		{
			if (a[i] == -1)
			{
				mixedAddition(negative, Q, middleP);//Q-P0
			}
			else
			{
				copy(P[-a[i]].X, negativePi.X);
				copy(P[-a[i]].Z, negativePi.Z);
				modmult2(negativePi.X, negativePi.Z, R1);
				add2(P[-a[i]].Y, R1, negativePi.Y);
				mixedAddition(negativePi, Q, middleP);
			}
		}
		else
		{
			copy(Q.X, middleP.X);
			copy(Q.Y, middleP.Y);
			copy(Q.Z, middleP.Z);
		}

		tau(middleP, Q);
	}
	return 0;
}

Big rand(int n)
{//n-bit random integer

	Big value = 1;
	for (int i = 1; i < n; i++)
	{
		value *= 2;
		value = value + (rand() % 2);
	}
	return value;
}