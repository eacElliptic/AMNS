#ifndef STRUCTS_DATA
#define STRUCTS_DATA


//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'
#define WORD_SIZE 64
#define POLY_DEG 4
#define NB_COEFF 5
#define NB_ADD_MAX 4

#define RHO_LOG2 55
//~ We will take : rho = 1 << RHO_LOG2.


typedef __int128 int128;
typedef unsigned __int128 uint128;


static uint64_t amns_rho;

//~ will contain a representative of 'rho' in the amns
//~ important : this initial value is that of 'rho*phi'; correct value will be put during initialisation step
static int64_t rho_rep[NB_COEFF] = {1369907509596552, 745466846945416, 1077152682307042, -120724997098163, 504989843537910};

//~ representatives of (RHO)^i (for i=2,4,...) in the amns (for convertions : int to amns)
static int64_t RHO_POWS[(NB_COEFF - 2)][NB_COEFF];

static mpz_t modul_p;
static mpz_t gama_pow[POLY_DEG];

#endif

