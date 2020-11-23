#pragma comment(lib,"miracl.lib")
//#define _MIPT_
//
//extern   "C"   
//{   
//#include "miracl.h"
//#include "mirdef.h"
//	//#pragma comment( lib, "ms32.lib") 
//}   
//#if _DEBUG
//#pragma comment(linker,"/NODEFAULTLIB:LIBC")
//#endif


#include <big.h>
#include <iostream>
using namespace std;


#include <miracl.h>





void swapr(int &a, int &b)
{
	int temp;
	temp = a;
	a = b;
	b = temp;
}






int main()
{


	int x = 3, y = 5;
	swapr(x, y);
	cout << x << y << endl;

	big ECC112, b, c;
	miracl *mip = mirsys(5000, 16);
	char str[1000] = "1FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF";
	mip->IOBASE=16;
	ECC112=mirvar(8);
	cinstr(ECC112, str);
	b=mirvar(3);
	c=mirvar(0);
	divide(ECC112, b, c);
	cotnum(ECC112, stdout);
	//cotnum(c, stdout);

	cinstr(ECC112, str);
	b=mirvar(4);
	divide(ECC112, b, c);
	cotnum(ECC112, stdout);

	//char s[100];
	//cotstr(ECC112,s);
	//cotnum(b, stdout);
	//cout << s<< endl;
	return 0;
}

//#include "miracl.h"
//#include <time.h>
//int main()
//{
//	int i;
//	big x, e, m, y;
//	FILE *fp;
//	clock_t tBegin, tEnd;
//	miracl *mip = mirsys(1000, 16);
//	x = mirvar(0);
//	e = mirvar(0);
//	m = mirvar(0);
//	y = mirvar(0);
//	fp = fopen("data.txt", "r+"); 
//	mip->IOBASE = 16;
//	cinnum(x, fp);
//	cinnum(e, fp);
//	cinnum(m, fp);
//	fclose(fp);
//	tBegin = clock();
//	for (i = 0; i < 100; i ++)
//		powmod(x, e, m, y);
//	tEnd = clock();
//	cotnum(x, stdout);
//	cotnum(e, stdout);
//	cotnum(m, stdout);
//	cotnum(y, stdout);
//	printf("\n\n进行100次1024比特的模指数运算所消耗的时间为:%ld ms\n\n", tEnd - tBegin);
//	return 0;
//}