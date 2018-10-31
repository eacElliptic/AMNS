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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10890567120602195863UL) + ((((uint64_t)op[1] * 5956619470078782671UL) + ((uint64_t)op[2] * 16656720969382742268UL) + ((uint64_t)op[3] * 13481373541443158714UL) + ((uint64_t)op[4] * 10246587891793302331UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 10246587891793302331UL) + ((uint64_t)op[1] * 10890567120602195863UL) + ((((uint64_t)op[2] * 5956619470078782671UL) + ((uint64_t)op[3] * 16656720969382742268UL) + ((uint64_t)op[4] * 13481373541443158714UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 13481373541443158714UL) + ((uint64_t)op[1] * 10246587891793302331UL) + ((uint64_t)op[2] * 10890567120602195863UL) + ((((uint64_t)op[3] * 5956619470078782671UL) + ((uint64_t)op[4] * 16656720969382742268UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 16656720969382742268UL) + ((uint64_t)op[1] * 13481373541443158714UL) + ((uint64_t)op[2] * 10246587891793302331UL) + ((uint64_t)op[3] * 10890567120602195863UL) + ((uint64_t)op[4] * 7687276460498393480UL);
	tmp_q[4] = ((uint64_t)op[0] * 5956619470078782671UL) + ((uint64_t)op[1] * 16656720969382742268UL) + ((uint64_t)op[2] * 13481373541443158714UL) + ((uint64_t)op[3] * 10246587891793302331UL) + ((uint64_t)op[4] * 10890567120602195863UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 22784123343L) - ((((int128)tmp_q[1] * 224317550178L) + ((int128)tmp_q[2] * 140813049699L) + ((int128)tmp_q[3] * 102958276387L) - ((int128)tmp_q[4] * 146551621845L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 146551621845L) - ((int128)tmp_q[1] * 22784123343L) - ((((int128)tmp_q[2] * 224317550178L) + ((int128)tmp_q[3] * 140813049699L) + ((int128)tmp_q[4] * 102958276387L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 102958276387L) - ((int128)tmp_q[1] * 146551621845L) - ((int128)tmp_q[2] * 22784123343L) - ((((int128)tmp_q[3] * 224317550178L) + ((int128)tmp_q[4] * 140813049699L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 140813049699L) + ((int128)tmp_q[1] * 102958276387L) - ((int128)tmp_q[2] * 146551621845L) - ((int128)tmp_q[3] * 22784123343L) - ((int128)tmp_q[4] * 1794540401424L);
	tmp_zero[4] = ((int128)tmp_q[0] * 224317550178L) + ((int128)tmp_q[1] * 140813049699L) + ((int128)tmp_q[2] * 102958276387L) - ((int128)tmp_q[3] * 146551621845L) - ((int128)tmp_q[4] * 22784123343L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

