#include "add_mult_poly.h"


void add_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] + pb[j];
}

void sub_poly(int64_t *rop, int64_t *pa, int64_t *pb){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = pa[j] - pb[j];
}

void neg_poly(int64_t *rop, int64_t *op){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = -op[j];
}

//~ assumes 'scalar' and/or coeffs of 'op' small enough to avoid an overflow.
void scalar_mult_poly(int64_t *rop, int64_t *op, int64_t scalar){
	int j;
	for (j=0; j<NB_COEFF; j++)
		rop[j] = scalar * op[j];
}

//~ Computes pa(X)*pb(X) mod(X^n - c)
void mult_mod_poly(int64_t *rop, int64_t *pa, int64_t *pb){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17068803779512897489UL) + ((((uint64_t)op[1] * 3212326586507674772UL) + ((uint64_t)op[2] * 16268355636730401272UL) + ((uint64_t)op[3] * 2454018166354547410UL) + ((uint64_t)op[4] * 13544246269527365028UL) + ((uint64_t)op[5] * 11823279385442683832UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 11823279385442683832UL) + ((uint64_t)op[1] * 17068803779512897489UL) + ((((uint64_t)op[2] * 3212326586507674772UL) + ((uint64_t)op[3] * 16268355636730401272UL) + ((uint64_t)op[4] * 2454018166354547410UL) + ((uint64_t)op[5] * 13544246269527365028UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 13544246269527365028UL) + ((uint64_t)op[1] * 11823279385442683832UL) + ((uint64_t)op[2] * 17068803779512897489UL) + ((((uint64_t)op[3] * 3212326586507674772UL) + ((uint64_t)op[4] * 16268355636730401272UL) + ((uint64_t)op[5] * 2454018166354547410UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 2454018166354547410UL) + ((uint64_t)op[1] * 13544246269527365028UL) + ((uint64_t)op[2] * 11823279385442683832UL) + ((uint64_t)op[3] * 17068803779512897489UL) + ((((uint64_t)op[4] * 3212326586507674772UL) + ((uint64_t)op[5] * 16268355636730401272UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 16268355636730401272UL) + ((uint64_t)op[1] * 2454018166354547410UL) + ((uint64_t)op[2] * 13544246269527365028UL) + ((uint64_t)op[3] * 11823279385442683832UL) + ((uint64_t)op[4] * 17068803779512897489UL) + ((uint64_t)op[5] * 5597437727678852528UL);
	tmp_q[5] = ((uint64_t)op[0] * 3212326586507674772UL) + ((uint64_t)op[1] * 16268355636730401272UL) + ((uint64_t)op[2] * 2454018166354547410UL) + ((uint64_t)op[3] * 13544246269527365028UL) + ((uint64_t)op[4] * 11823279385442683832UL) + ((uint64_t)op[5] * 17068803779512897489UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 619731145439L) - ((-((int128)tmp_q[1] * 3965818997820L) - ((int128)tmp_q[2] * 272672001336L) + ((int128)tmp_q[3] * 2502311398578L) + ((int128)tmp_q[4] * 3830558876260L) + ((int128)tmp_q[5] * 1786630624440L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 1786630624440L) + ((int128)tmp_q[1] * 619731145439L) - ((-((int128)tmp_q[2] * 3965818997820L) - ((int128)tmp_q[3] * 272672001336L) + ((int128)tmp_q[4] * 2502311398578L) + ((int128)tmp_q[5] * 3830558876260L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 3830558876260L) + ((int128)tmp_q[1] * 1786630624440L) + ((int128)tmp_q[2] * 619731145439L) - ((-((int128)tmp_q[3] * 3965818997820L) - ((int128)tmp_q[4] * 272672001336L) + ((int128)tmp_q[5] * 2502311398578L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 2502311398578L) + ((int128)tmp_q[1] * 3830558876260L) + ((int128)tmp_q[2] * 1786630624440L) + ((int128)tmp_q[3] * 619731145439L) - ((-((int128)tmp_q[4] * 3965818997820L) - ((int128)tmp_q[5] * 272672001336L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 272672001336L) + ((int128)tmp_q[1] * 2502311398578L) + ((int128)tmp_q[2] * 3830558876260L) + ((int128)tmp_q[3] * 1786630624440L) + ((int128)tmp_q[4] * 619731145439L) + ((int128)tmp_q[5] * 15863275991280L);
	tmp_zero[5] = -((int128)tmp_q[0] * 3965818997820L) - ((int128)tmp_q[1] * 272672001336L) + ((int128)tmp_q[2] * 2502311398578L) + ((int128)tmp_q[3] * 3830558876260L) + ((int128)tmp_q[4] * 1786630624440L) + ((int128)tmp_q[5] * 619731145439L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

