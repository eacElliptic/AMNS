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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[3] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8891500004071297081UL) + ((((uint64_t)op[1] * 737184773870705192UL) + ((uint64_t)op[2] * 12691627108205306234UL) + ((uint64_t)op[3] * 3553242335259920712UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 3553242335259920712UL) + ((uint64_t)op[1] * 8891500004071297081UL) + ((((uint64_t)op[2] * 737184773870705192UL) + ((uint64_t)op[3] * 12691627108205306234UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 12691627108205306234UL) + ((uint64_t)op[1] * 3553242335259920712UL) + ((uint64_t)op[2] * 8891500004071297081UL) + ((uint64_t)op[3] * 16972374525968141232UL);
	tmp_q[3] = ((uint64_t)op[0] * 737184773870705192UL) + ((uint64_t)op[1] * 12691627108205306234UL) + ((uint64_t)op[2] * 3553242335259920712UL) + ((uint64_t)op[3] * 8891500004071297081UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 222440658892097L) - ((-((int128)tmp_q[1] * 11615600457912L) + ((int128)tmp_q[2] * 30144750289418L) - ((int128)tmp_q[3] * 95005445007544L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 95005445007544L) - ((int128)tmp_q[1] * 222440658892097L) - ((-((int128)tmp_q[2] * 11615600457912L) + ((int128)tmp_q[3] * 30144750289418L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 30144750289418L) - ((int128)tmp_q[1] * 95005445007544L) - ((int128)tmp_q[2] * 222440658892097L) + ((int128)tmp_q[3] * 23231200915824L);
	tmp_zero[3] = -((int128)tmp_q[0] * 11615600457912L) + ((int128)tmp_q[1] * 30144750289418L) - ((int128)tmp_q[2] * 95005445007544L) - ((int128)tmp_q[3] * 222440658892097L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

