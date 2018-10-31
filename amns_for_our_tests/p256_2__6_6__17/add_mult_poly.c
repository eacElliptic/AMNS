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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7043238379221080421UL) + ((((uint64_t)op[1] * 8908091167324064332UL) + ((uint64_t)op[2] * 12495655176575685350UL) + ((uint64_t)op[3] * 72217953141739366UL) + ((uint64_t)op[4] * 10641461971690746741UL) + ((uint64_t)op[5] * 14511505162036972290UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 14511505162036972290UL) + ((uint64_t)op[1] * 7043238379221080421UL) + ((((uint64_t)op[2] * 8908091167324064332UL) + ((uint64_t)op[3] * 12495655176575685350UL) + ((uint64_t)op[4] * 72217953141739366UL) + ((uint64_t)op[5] * 10641461971690746741UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 10641461971690746741UL) + ((uint64_t)op[1] * 14511505162036972290UL) + ((uint64_t)op[2] * 7043238379221080421UL) + ((((uint64_t)op[3] * 8908091167324064332UL) + ((uint64_t)op[4] * 12495655176575685350UL) + ((uint64_t)op[5] * 72217953141739366UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 72217953141739366UL) + ((uint64_t)op[1] * 10641461971690746741UL) + ((uint64_t)op[2] * 14511505162036972290UL) + ((uint64_t)op[3] * 7043238379221080421UL) + ((((uint64_t)op[4] * 8908091167324064332UL) + ((uint64_t)op[5] * 12495655176575685350UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 12495655176575685350UL) + ((uint64_t)op[1] * 72217953141739366UL) + ((uint64_t)op[2] * 10641461971690746741UL) + ((uint64_t)op[3] * 14511505162036972290UL) + ((uint64_t)op[4] * 7043238379221080421UL) + ((uint64_t)op[5] * 16555058856525282760UL);
	tmp_q[5] = ((uint64_t)op[0] * 8908091167324064332UL) + ((uint64_t)op[1] * 12495655176575685350UL) + ((uint64_t)op[2] * 72217953141739366UL) + ((uint64_t)op[3] * 10641461971690746741UL) + ((uint64_t)op[4] * 14511505162036972290UL) + ((uint64_t)op[5] * 7043238379221080421UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2662612013203L) + ((-((int128)tmp_q[1] * 1964670033894L) + ((int128)tmp_q[2] * 3219182048919L) + ((int128)tmp_q[3] * 156285382022L) - ((int128)tmp_q[4] * 1822863574013L) + ((int128)tmp_q[5] * 3715089140630L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 3715089140630L) - ((int128)tmp_q[1] * 2662612013203L) + ((-((int128)tmp_q[2] * 1964670033894L) + ((int128)tmp_q[3] * 3219182048919L) + ((int128)tmp_q[4] * 156285382022L) - ((int128)tmp_q[5] * 1822863574013L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 1822863574013L) + ((int128)tmp_q[1] * 3715089140630L) - ((int128)tmp_q[2] * 2662612013203L) + ((-((int128)tmp_q[3] * 1964670033894L) + ((int128)tmp_q[4] * 3219182048919L) + ((int128)tmp_q[5] * 156285382022L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 156285382022L) - ((int128)tmp_q[1] * 1822863574013L) + ((int128)tmp_q[2] * 3715089140630L) - ((int128)tmp_q[3] * 2662612013203L) + ((-((int128)tmp_q[4] * 1964670033894L) + ((int128)tmp_q[5] * 3219182048919L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 3219182048919L) + ((int128)tmp_q[1] * 156285382022L) - ((int128)tmp_q[2] * 1822863574013L) + ((int128)tmp_q[3] * 3715089140630L) - ((int128)tmp_q[4] * 2662612013203L) - ((int128)tmp_q[5] * 11788020203364L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1964670033894L) + ((int128)tmp_q[1] * 3219182048919L) + ((int128)tmp_q[2] * 156285382022L) - ((int128)tmp_q[3] * 1822863574013L) + ((int128)tmp_q[4] * 3715089140630L) - ((int128)tmp_q[5] * 2662612013203L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

