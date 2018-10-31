#ifndef STRUCTS_DATA
#define STRUCTS_DATA


//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'
#define WORD_SIZE 64
#define POLY_DEG 7
#define NB_COEFF 8
#define NB_ADD_MAX 2

#define RHO_LOG2 54
//~ We will take : rho = 1 << RHO_LOG2.


typedef __int128 int128;
typedef unsigned __int128 uint128;


static uint64_t amns_rho;

//~ will contain a representative of 'rho' in the amns
//~ important : this initial value is that of 'rho*phi'; correct value will be put during initialisation step
static int64_t rho_rep[NB_COEFF] = {-645649518596818, -691385904768596, -77186075549493, 124971213534472, -271154449765646, 4341772222111, 720775245609795, -20344811836775};

//~ representatives of (RHO)^i (for i=2,4,...) in the amns (for convertions : int to amns)
static int64_t RHO_POWS[(NB_COEFF - 2)][NB_COEFF];

static mpz_t modul_p;
static mpz_t gama_pow[POLY_DEG];

#endif

