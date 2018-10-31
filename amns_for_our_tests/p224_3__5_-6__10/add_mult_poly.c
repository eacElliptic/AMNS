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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3752868482141018997UL) + ((((uint64_t)op[1] * 5025900294417424582UL) + ((uint64_t)op[2] * 12458436945495784604UL) + ((uint64_t)op[3] * 3959996915376071351UL) + ((uint64_t)op[4] * 6245704865895534809UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 6245704865895534809UL) + ((uint64_t)op[1] * 3752868482141018997UL) + ((((uint64_t)op[2] * 5025900294417424582UL) + ((uint64_t)op[3] * 12458436945495784604UL) + ((uint64_t)op[4] * 3959996915376071351UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 3959996915376071351UL) + ((uint64_t)op[1] * 6245704865895534809UL) + ((uint64_t)op[2] * 3752868482141018997UL) + ((((uint64_t)op[3] * 5025900294417424582UL) + ((uint64_t)op[4] * 12458436945495784604UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 12458436945495784604UL) + ((uint64_t)op[1] * 3959996915376071351UL) + ((uint64_t)op[2] * 6245704865895534809UL) + ((uint64_t)op[3] * 3752868482141018997UL) + ((uint64_t)op[4] * 6738086380914555740UL);
	tmp_q[4] = ((uint64_t)op[0] * 5025900294417424582UL) + ((uint64_t)op[1] * 12458436945495784604UL) + ((uint64_t)op[2] * 3959996915376071351UL) + ((uint64_t)op[3] * 6245704865895534809UL) + ((uint64_t)op[4] * 3752868482141018997UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 5359592734607L) - ((-((int128)tmp_q[1] * 5545400517393L) + ((int128)tmp_q[2] * 17173757741031L) + ((int128)tmp_q[3] * 7815712861016L) + ((int128)tmp_q[4] * 13131233511067L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 13131233511067L) + ((int128)tmp_q[1] * 5359592734607L) - ((-((int128)tmp_q[2] * 5545400517393L) + ((int128)tmp_q[3] * 17173757741031L) + ((int128)tmp_q[4] * 7815712861016L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 7815712861016L) + ((int128)tmp_q[1] * 13131233511067L) + ((int128)tmp_q[2] * 5359592734607L) - ((-((int128)tmp_q[3] * 5545400517393L) + ((int128)tmp_q[4] * 17173757741031L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 17173757741031L) + ((int128)tmp_q[1] * 7815712861016L) + ((int128)tmp_q[2] * 13131233511067L) + ((int128)tmp_q[3] * 5359592734607L) + ((int128)tmp_q[4] * 33272403104358L);
	tmp_zero[4] = -((int128)tmp_q[0] * 5545400517393L) + ((int128)tmp_q[1] * 17173757741031L) + ((int128)tmp_q[2] * 7815712861016L) + ((int128)tmp_q[3] * 13131233511067L) + ((int128)tmp_q[4] * 5359592734607L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

