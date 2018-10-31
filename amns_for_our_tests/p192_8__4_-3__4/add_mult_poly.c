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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[3] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17540045803519390407UL) + ((((uint64_t)op[1] * 6174151344266166023UL) + ((uint64_t)op[2] * 14693380158962934595UL) + ((uint64_t)op[3] * 1173026699617184038UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 1173026699617184038UL) + ((uint64_t)op[1] * 17540045803519390407UL) + ((((uint64_t)op[2] * 6174151344266166023UL) + ((uint64_t)op[3] * 14693380158962934595UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 14693380158962934595UL) + ((uint64_t)op[1] * 1173026699617184038UL) + ((uint64_t)op[2] * 17540045803519390407UL) + ((uint64_t)op[3] * 18371034114620605163UL);
	tmp_q[3] = ((uint64_t)op[0] * 6174151344266166023UL) + ((uint64_t)op[1] * 14693380158962934595UL) + ((uint64_t)op[2] * 1173026699617184038UL) + ((uint64_t)op[3] * 17540045803519390407UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 30287775527047L) - ((-((int128)tmp_q[1] * 146070463976520L) - ((int128)tmp_q[2] * 174527251981329L) - ((int128)tmp_q[3] * 97238988816209L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 97238988816209L) + ((int128)tmp_q[1] * 30287775527047L) - ((-((int128)tmp_q[2] * 146070463976520L) - ((int128)tmp_q[3] * 174527251981329L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 174527251981329L) - ((int128)tmp_q[1] * 97238988816209L) + ((int128)tmp_q[2] * 30287775527047L) + ((int128)tmp_q[3] * 438211391929560L);
	tmp_zero[3] = -((int128)tmp_q[0] * 146070463976520L) - ((int128)tmp_q[1] * 174527251981329L) - ((int128)tmp_q[2] * 97238988816209L) + ((int128)tmp_q[3] * 30287775527047L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

