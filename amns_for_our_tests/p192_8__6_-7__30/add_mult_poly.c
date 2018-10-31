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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15270589243618975668UL) + ((((uint64_t)op[1] * 16506100377991992447UL) + ((uint64_t)op[2] * 15981170267314657401UL) + ((uint64_t)op[3] * 1867366425783259618UL) + ((uint64_t)op[4] * 9569966564065088004UL) + ((uint64_t)op[5] * 9736386704657897575UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 9736386704657897575UL) + ((uint64_t)op[1] * 15270589243618975668UL) + ((((uint64_t)op[2] * 16506100377991992447UL) + ((uint64_t)op[3] * 15981170267314657401UL) + ((uint64_t)op[4] * 1867366425783259618UL) + ((uint64_t)op[5] * 9569966564065088004UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 9569966564065088004UL) + ((uint64_t)op[1] * 9736386704657897575UL) + ((uint64_t)op[2] * 15270589243618975668UL) + ((((uint64_t)op[3] * 16506100377991992447UL) + ((uint64_t)op[4] * 15981170267314657401UL) + ((uint64_t)op[5] * 1867366425783259618UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 1867366425783259618UL) + ((uint64_t)op[1] * 9569966564065088004UL) + ((uint64_t)op[2] * 9736386704657897575UL) + ((uint64_t)op[3] * 15270589243618975668UL) + ((((uint64_t)op[4] * 16506100377991992447UL) + ((uint64_t)op[5] * 15981170267314657401UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 15981170267314657401UL) + ((uint64_t)op[1] * 1867366425783259618UL) + ((uint64_t)op[2] * 9569966564065088004UL) + ((uint64_t)op[3] * 9736386704657897575UL) + ((uint64_t)op[4] * 15270589243618975668UL) + ((uint64_t)op[5] * 13584505870022914183UL);
	tmp_q[5] = ((uint64_t)op[0] * 16506100377991992447UL) + ((uint64_t)op[1] * 15981170267314657401UL) + ((uint64_t)op[2] * 1867366425783259618UL) + ((uint64_t)op[3] * 9569966564065088004UL) + ((uint64_t)op[4] * 9736386704657897575UL) + ((uint64_t)op[5] * 15270589243618975668UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3732018893L) - ((((int128)tmp_q[1] * 2412066920L) - ((int128)tmp_q[2] * 525652L) + ((int128)tmp_q[3] * 2462381359L) - ((int128)tmp_q[4] * 399414282L) + ((int128)tmp_q[5] * 891623001L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 891623001L) + ((int128)tmp_q[1] * 3732018893L) - ((((int128)tmp_q[2] * 2412066920L) - ((int128)tmp_q[3] * 525652L) + ((int128)tmp_q[4] * 2462381359L) - ((int128)tmp_q[5] * 399414282L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 399414282L) + ((int128)tmp_q[1] * 891623001L) + ((int128)tmp_q[2] * 3732018893L) - ((((int128)tmp_q[3] * 2412066920L) - ((int128)tmp_q[4] * 525652L) + ((int128)tmp_q[5] * 2462381359L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 2462381359L) - ((int128)tmp_q[1] * 399414282L) + ((int128)tmp_q[2] * 891623001L) + ((int128)tmp_q[3] * 3732018893L) - ((((int128)tmp_q[4] * 2412066920L) - ((int128)tmp_q[5] * 525652L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 525652L) + ((int128)tmp_q[1] * 2462381359L) - ((int128)tmp_q[2] * 399414282L) + ((int128)tmp_q[3] * 891623001L) + ((int128)tmp_q[4] * 3732018893L) - ((int128)tmp_q[5] * 16884468440L);
	tmp_zero[5] = ((int128)tmp_q[0] * 2412066920L) - ((int128)tmp_q[1] * 525652L) + ((int128)tmp_q[2] * 2462381359L) - ((int128)tmp_q[3] * 399414282L) + ((int128)tmp_q[4] * 891623001L) + ((int128)tmp_q[5] * 3732018893L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

