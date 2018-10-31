#ifndef STRUCTS_DATA
#define STRUCTS_DATA


//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'
#define WORD_SIZE 32
#define POLY_DEG 12
#define NB_COEFF 13
#define NB_ADD_MAX 0

#define RHO_LOG2 25
//~ We will take : rho = 1 << RHO_LOG2.


typedef long long llong;
typedef unsigned long long ullong;


static uint amns_rho;

//~ will contain a representative of 'rho' in the amns
//~ important : this initial value is that of 'rho*phi'; correct value will be put during initialisation step
static int rho_rep[NB_COEFF] = {-353569, -1551726, -75969, -442411, 585906, -733358, -205736, -1223599, -861081, -806387, -312442, 257680, 117074};

//~ representatives of (RHO)^i (for i=2,4,...) in the amns (for convertions : int to amns)
static int RHO_POWS[(NB_COEFF - 2)][NB_COEFF];

static mpz_t modul_p;
static mpz_t gama_pow[POLY_DEG];

#endif

