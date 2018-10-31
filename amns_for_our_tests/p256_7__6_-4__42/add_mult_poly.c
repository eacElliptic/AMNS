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
	tmp_q[0] = ((uint64_t)op[0] * 2711972735090405835UL) + ((((uint64_t)op[1] * 6714649536836943374UL) + ((uint64_t)op[2] * 1168840691177390540UL) + ((uint64_t)op[3] * 15378557961607468806UL) + ((uint64_t)op[4] * 3950081867525556912UL) + ((uint64_t)op[5] * 15790319750786532807UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 15790319750786532807UL) + ((uint64_t)op[1] * 2711972735090405835UL) + ((((uint64_t)op[2] * 6714649536836943374UL) + ((uint64_t)op[3] * 1168840691177390540UL) + ((uint64_t)op[4] * 15378557961607468806UL) + ((uint64_t)op[5] * 3950081867525556912UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 3950081867525556912UL) + ((uint64_t)op[1] * 15790319750786532807UL) + ((uint64_t)op[2] * 2711972735090405835UL) + ((((uint64_t)op[3] * 6714649536836943374UL) + ((uint64_t)op[4] * 1168840691177390540UL) + ((uint64_t)op[5] * 15378557961607468806UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 15378557961607468806UL) + ((uint64_t)op[1] * 3950081867525556912UL) + ((uint64_t)op[2] * 15790319750786532807UL) + ((uint64_t)op[3] * 2711972735090405835UL) + ((((uint64_t)op[4] * 6714649536836943374UL) + ((uint64_t)op[5] * 1168840691177390540UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 1168840691177390540UL) + ((uint64_t)op[1] * 15378557961607468806UL) + ((uint64_t)op[2] * 3950081867525556912UL) + ((uint64_t)op[3] * 15790319750786532807UL) + ((uint64_t)op[4] * 2711972735090405835UL) + ((uint64_t)op[5] * 10034890000071329736UL);
	tmp_q[5] = ((uint64_t)op[0] * 6714649536836943374UL) + ((uint64_t)op[1] * 1168840691177390540UL) + ((uint64_t)op[2] * 15378557961607468806UL) + ((uint64_t)op[3] * 3950081867525556912UL) + ((uint64_t)op[4] * 15790319750786532807UL) + ((uint64_t)op[5] * 2711972735090405835UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3166367786633L) - ((-((int128)tmp_q[1] * 3625933064965L) + ((int128)tmp_q[2] * 1501458471465L) - ((int128)tmp_q[3] * 1166172842479L) + ((int128)tmp_q[4] * 1721849513345L) - ((int128)tmp_q[5] * 2562401826589L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 2562401826589L) + ((int128)tmp_q[1] * 3166367786633L) - ((-((int128)tmp_q[2] * 3625933064965L) + ((int128)tmp_q[3] * 1501458471465L) - ((int128)tmp_q[4] * 1166172842479L) + ((int128)tmp_q[5] * 1721849513345L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 1721849513345L) - ((int128)tmp_q[1] * 2562401826589L) + ((int128)tmp_q[2] * 3166367786633L) - ((-((int128)tmp_q[3] * 3625933064965L) + ((int128)tmp_q[4] * 1501458471465L) - ((int128)tmp_q[5] * 1166172842479L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 1166172842479L) + ((int128)tmp_q[1] * 1721849513345L) - ((int128)tmp_q[2] * 2562401826589L) + ((int128)tmp_q[3] * 3166367786633L) - ((-((int128)tmp_q[4] * 3625933064965L) + ((int128)tmp_q[5] * 1501458471465L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 1501458471465L) - ((int128)tmp_q[1] * 1166172842479L) + ((int128)tmp_q[2] * 1721849513345L) - ((int128)tmp_q[3] * 2562401826589L) + ((int128)tmp_q[4] * 3166367786633L) + ((int128)tmp_q[5] * 14503732259860L);
	tmp_zero[5] = -((int128)tmp_q[0] * 3625933064965L) + ((int128)tmp_q[1] * 1501458471465L) - ((int128)tmp_q[2] * 1166172842479L) + ((int128)tmp_q[3] * 1721849513345L) - ((int128)tmp_q[4] * 2562401826589L) + ((int128)tmp_q[5] * 3166367786633L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

