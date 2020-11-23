#include "big.h"

#define MR_EPOINT_GENERAL    0
#define MR_EPOINT_NORMALIZED 1
#define MR_EPOINT_INFINITY   2


typedef struct {
	int marker;
	big X;
	big Y;
} eccPoint;

typedef struct{
	int intNumberI;
	int intNumberM;
	int intNumberS;
}scalarMultiplyCost;



typedef struct {
	/** The size of the curve in octets */
	int size;

	/** name of curve */
	char *name; 

	/** The prime that defines the field the curve is in (encoded in hex) */
	char *prime;

	/** The fields B param (hex) */
	char *B;

	/** The order of the curve (hex) */
	char *order;

	/** The x co-ordinate of the base point on the curve (hex) */
	char *Gx;

	/** The y co-ordinate of the base point on the curve (hex) */
	char *Gy;
} pkEccSetType;



class ECC
{
	eccPoint *P;
public:
	ECC(big x, big y);
	ECC(Big x, Big y);
	BOOL set(const Big& x,const Big& y)    {return  P->marker=0; P->X=x.getbig();P->Y=y.getbig();}
    BOOL iszero() const;
    int get(Big& x,Big& y) const;
	ECC & operator=(const ECC &b){ return *this;}


};

