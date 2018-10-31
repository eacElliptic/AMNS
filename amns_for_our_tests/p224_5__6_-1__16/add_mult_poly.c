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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - ((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - ((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - ((int128)pa[5] * pa[5]);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14158592410482891515UL) + ((((uint64_t)op[1] * 13506930773857554739UL) + ((uint64_t)op[2] * 1102931065813550284UL) + ((uint64_t)op[3] * 17273684843574626692UL) + ((uint64_t)op[4] * 4447914881363256369UL) + ((uint64_t)op[5] * 4402848692948435802UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 4402848692948435802UL) + ((uint64_t)op[1] * 14158592410482891515UL) + ((((uint64_t)op[2] * 13506930773857554739UL) + ((uint64_t)op[3] * 1102931065813550284UL) + ((uint64_t)op[4] * 17273684843574626692UL) + ((uint64_t)op[5] * 4447914881363256369UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 4447914881363256369UL) + ((uint64_t)op[1] * 4402848692948435802UL) + ((uint64_t)op[2] * 14158592410482891515UL) + ((((uint64_t)op[3] * 13506930773857554739UL) + ((uint64_t)op[4] * 1102931065813550284UL) + ((uint64_t)op[5] * 17273684843574626692UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 17273684843574626692UL) + ((uint64_t)op[1] * 4447914881363256369UL) + ((uint64_t)op[2] * 4402848692948435802UL) + ((uint64_t)op[3] * 14158592410482891515UL) + ((((uint64_t)op[4] * 13506930773857554739UL) + ((uint64_t)op[5] * 1102931065813550284UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 1102931065813550284UL) + ((uint64_t)op[1] * 17273684843574626692UL) + ((uint64_t)op[2] * 4447914881363256369UL) + ((uint64_t)op[3] * 4402848692948435802UL) + ((uint64_t)op[4] * 14158592410482891515UL) + ((uint64_t)op[5] * 4939813299851996877UL);
	tmp_q[5] = ((uint64_t)op[0] * 13506930773857554739UL) + ((uint64_t)op[1] * 1102931065813550284UL) + ((uint64_t)op[2] * 17273684843574626692UL) + ((uint64_t)op[3] * 4447914881363256369UL) + ((uint64_t)op[4] * 4402848692948435802UL) + ((uint64_t)op[5] * 14158592410482891515UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3076395023784179L) - (-((int128)tmp_q[1] * 29042231223318029L) - ((int128)tmp_q[2] * 27503065626573396L) - ((int128)tmp_q[3] * 6938303930640774L) - ((int128)tmp_q[4] * 24426670602789219L) + ((int128)tmp_q[5] * 22103927292677260L));
	tmp_zero[1] = ((int128)tmp_q[0] * 22103927292677260L) + ((int128)tmp_q[1] * 3076395023784179L) - (-((int128)tmp_q[2] * 29042231223318029L) - ((int128)tmp_q[3] * 27503065626573396L) - ((int128)tmp_q[4] * 6938303930640774L) - ((int128)tmp_q[5] * 24426670602789219L));
	tmp_zero[2] = -((int128)tmp_q[0] * 24426670602789219L) + ((int128)tmp_q[1] * 22103927292677260L) + ((int128)tmp_q[2] * 3076395023784179L) - (-((int128)tmp_q[3] * 29042231223318029L) - ((int128)tmp_q[4] * 27503065626573396L) - ((int128)tmp_q[5] * 6938303930640774L));
	tmp_zero[3] = -((int128)tmp_q[0] * 6938303930640774L) - ((int128)tmp_q[1] * 24426670602789219L) + ((int128)tmp_q[2] * 22103927292677260L) + ((int128)tmp_q[3] * 3076395023784179L) - (-((int128)tmp_q[4] * 29042231223318029L) - ((int128)tmp_q[5] * 27503065626573396L));
	tmp_zero[4] = -((int128)tmp_q[0] * 27503065626573396L) - ((int128)tmp_q[1] * 6938303930640774L) - ((int128)tmp_q[2] * 24426670602789219L) + ((int128)tmp_q[3] * 22103927292677260L) + ((int128)tmp_q[4] * 3076395023784179L) + ((int128)tmp_q[5] * 29042231223318029L);
	tmp_zero[5] = -((int128)tmp_q[0] * 29042231223318029L) - ((int128)tmp_q[1] * 27503065626573396L) - ((int128)tmp_q[2] * 6938303930640774L) - ((int128)tmp_q[3] * 24426670602789219L) + ((int128)tmp_q[4] * 22103927292677260L) + ((int128)tmp_q[5] * 3076395023784179L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

