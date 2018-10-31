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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13145711564329269325UL) + ((((uint64_t)op[1] * 15322048736200204296UL) + ((uint64_t)op[2] * 15494556806292442872UL) + ((uint64_t)op[3] * 4488228915165159448UL) + ((uint64_t)op[4] * 12020688878985909955UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 12020688878985909955UL) + ((uint64_t)op[1] * 13145711564329269325UL) + ((((uint64_t)op[2] * 15322048736200204296UL) + ((uint64_t)op[3] * 15494556806292442872UL) + ((uint64_t)op[4] * 4488228915165159448UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 4488228915165159448UL) + ((uint64_t)op[1] * 12020688878985909955UL) + ((uint64_t)op[2] * 13145711564329269325UL) + ((((uint64_t)op[3] * 15322048736200204296UL) + ((uint64_t)op[4] * 15494556806292442872UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 15494556806292442872UL) + ((uint64_t)op[1] * 4488228915165159448UL) + ((uint64_t)op[2] * 12020688878985909955UL) + ((uint64_t)op[3] * 13145711564329269325UL) + ((uint64_t)op[4] * 6249390675018694640UL);
	tmp_q[4] = ((uint64_t)op[0] * 15322048736200204296UL) + ((uint64_t)op[1] * 15494556806292442872UL) + ((uint64_t)op[2] * 4488228915165159448UL) + ((uint64_t)op[3] * 12020688878985909955UL) + ((uint64_t)op[4] * 13145711564329269325UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1413121304385753L) - ((-((int128)tmp_q[1] * 846469478082087L) + ((int128)tmp_q[2] * 881362564704801L) - ((int128)tmp_q[3] * 306817082595727L) + ((int128)tmp_q[4] * 604417652046969L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 604417652046969L) + ((int128)tmp_q[1] * 1413121304385753L) - ((-((int128)tmp_q[2] * 846469478082087L) + ((int128)tmp_q[3] * 881362564704801L) - ((int128)tmp_q[4] * 306817082595727L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 306817082595727L) + ((int128)tmp_q[1] * 604417652046969L) + ((int128)tmp_q[2] * 1413121304385753L) - ((-((int128)tmp_q[3] * 846469478082087L) + ((int128)tmp_q[4] * 881362564704801L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 881362564704801L) - ((int128)tmp_q[1] * 306817082595727L) + ((int128)tmp_q[2] * 604417652046969L) + ((int128)tmp_q[3] * 1413121304385753L) + ((int128)tmp_q[4] * 1692938956164174L);
	tmp_zero[4] = -((int128)tmp_q[0] * 846469478082087L) + ((int128)tmp_q[1] * 881362564704801L) - ((int128)tmp_q[2] * 306817082595727L) + ((int128)tmp_q[3] * 604417652046969L) + ((int128)tmp_q[4] * 1413121304385753L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

