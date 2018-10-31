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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 707893020676625611UL) + ((((uint64_t)op[1] * 11310691771883019790UL) + ((uint64_t)op[2] * 9276457577357394214UL) + ((uint64_t)op[3] * 15528741753804519637UL) + ((uint64_t)op[4] * 2666941202106167462UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 2666941202106167462UL) + ((uint64_t)op[1] * 707893020676625611UL) + ((((uint64_t)op[2] * 11310691771883019790UL) + ((uint64_t)op[3] * 9276457577357394214UL) + ((uint64_t)op[4] * 15528741753804519637UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 15528741753804519637UL) + ((uint64_t)op[1] * 2666941202106167462UL) + ((uint64_t)op[2] * 707893020676625611UL) + ((((uint64_t)op[3] * 11310691771883019790UL) + ((uint64_t)op[4] * 9276457577357394214UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 9276457577357394214UL) + ((uint64_t)op[1] * 15528741753804519637UL) + ((uint64_t)op[2] * 2666941202106167462UL) + ((uint64_t)op[3] * 707893020676625611UL) + ((uint64_t)op[4] * 10097465133596575688UL);
	tmp_q[4] = ((uint64_t)op[0] * 11310691771883019790UL) + ((uint64_t)op[1] * 9276457577357394214UL) + ((uint64_t)op[2] * 15528741753804519637UL) + ((uint64_t)op[3] * 2666941202106167462UL) + ((uint64_t)op[4] * 707893020676625611UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 193143438709L) - ((-((int128)tmp_q[1] * 97326717945L) - ((int128)tmp_q[2] * 29859193090L) - ((int128)tmp_q[3] * 127430039015L) - ((int128)tmp_q[4] * 138692936606L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 138692936606L) + ((int128)tmp_q[1] * 193143438709L) - ((-((int128)tmp_q[2] * 97326717945L) - ((int128)tmp_q[3] * 29859193090L) - ((int128)tmp_q[4] * 127430039015L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 127430039015L) - ((int128)tmp_q[1] * 138692936606L) + ((int128)tmp_q[2] * 193143438709L) - ((-((int128)tmp_q[3] * 97326717945L) - ((int128)tmp_q[4] * 29859193090L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 29859193090L) - ((int128)tmp_q[1] * 127430039015L) - ((int128)tmp_q[2] * 138692936606L) + ((int128)tmp_q[3] * 193143438709L) + ((int128)tmp_q[4] * 389306871780L);
	tmp_zero[4] = -((int128)tmp_q[0] * 97326717945L) - ((int128)tmp_q[1] * 29859193090L) - ((int128)tmp_q[2] * 127430039015L) - ((int128)tmp_q[3] * 138692936606L) + ((int128)tmp_q[4] * 193143438709L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

