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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11979896403154442959UL) + ((((uint64_t)op[1] * 9174583614823881405UL) + ((uint64_t)op[2] * 2763802635581170756UL) + ((uint64_t)op[3] * 575052460958758352UL) + ((uint64_t)op[4] * 13397466484509005391UL) + ((uint64_t)op[5] * 2467776485102580924UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 2467776485102580924UL) + ((uint64_t)op[1] * 11979896403154442959UL) + ((((uint64_t)op[2] * 9174583614823881405UL) + ((uint64_t)op[3] * 2763802635581170756UL) + ((uint64_t)op[4] * 575052460958758352UL) + ((uint64_t)op[5] * 13397466484509005391UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 13397466484509005391UL) + ((uint64_t)op[1] * 2467776485102580924UL) + ((uint64_t)op[2] * 11979896403154442959UL) + ((((uint64_t)op[3] * 9174583614823881405UL) + ((uint64_t)op[4] * 2763802635581170756UL) + ((uint64_t)op[5] * 575052460958758352UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 575052460958758352UL) + ((uint64_t)op[1] * 13397466484509005391UL) + ((uint64_t)op[2] * 2467776485102580924UL) + ((uint64_t)op[3] * 11979896403154442959UL) + ((((uint64_t)op[4] * 9174583614823881405UL) + ((uint64_t)op[5] * 2763802635581170756UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 2763802635581170756UL) + ((uint64_t)op[1] * 575052460958758352UL) + ((uint64_t)op[2] * 13397466484509005391UL) + ((uint64_t)op[3] * 2467776485102580924UL) + ((uint64_t)op[4] * 11979896403154442959UL) + ((uint64_t)op[5] * 18349167229647762810UL);
	tmp_q[5] = ((uint64_t)op[0] * 9174583614823881405UL) + ((uint64_t)op[1] * 2763802635581170756UL) + ((uint64_t)op[2] * 575052460958758352UL) + ((uint64_t)op[3] * 13397466484509005391UL) + ((uint64_t)op[4] * 2467776485102580924UL) + ((uint64_t)op[5] * 11979896403154442959UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 999816985L) + ((((int128)tmp_q[1] * 329683973L) + ((int128)tmp_q[2] * 2577349105L) - ((int128)tmp_q[3] * 967639438L) + ((int128)tmp_q[4] * 117559205L) + ((int128)tmp_q[5] * 644329036L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 644329036L) - ((int128)tmp_q[1] * 999816985L) + ((((int128)tmp_q[2] * 329683973L) + ((int128)tmp_q[3] * 2577349105L) - ((int128)tmp_q[4] * 967639438L) + ((int128)tmp_q[5] * 117559205L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 117559205L) + ((int128)tmp_q[1] * 644329036L) - ((int128)tmp_q[2] * 999816985L) + ((((int128)tmp_q[3] * 329683973L) + ((int128)tmp_q[4] * 2577349105L) - ((int128)tmp_q[5] * 967639438L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 967639438L) + ((int128)tmp_q[1] * 117559205L) + ((int128)tmp_q[2] * 644329036L) - ((int128)tmp_q[3] * 999816985L) + ((((int128)tmp_q[4] * 329683973L) + ((int128)tmp_q[5] * 2577349105L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 2577349105L) - ((int128)tmp_q[1] * 967639438L) + ((int128)tmp_q[2] * 117559205L) + ((int128)tmp_q[3] * 644329036L) - ((int128)tmp_q[4] * 999816985L) + ((int128)tmp_q[5] * 659367946L);
	tmp_zero[5] = ((int128)tmp_q[0] * 329683973L) + ((int128)tmp_q[1] * 2577349105L) - ((int128)tmp_q[2] * 967639438L) + ((int128)tmp_q[3] * 117559205L) + ((int128)tmp_q[4] * 644329036L) - ((int128)tmp_q[5] * 999816985L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

