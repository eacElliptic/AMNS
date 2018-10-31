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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17253757146333624087UL) + ((((uint64_t)op[1] * 18347504255334709693UL) + ((uint64_t)op[2] * 6479984006563990855UL) + ((uint64_t)op[3] * 5036921686197915088UL) + ((uint64_t)op[4] * 9027911415125877338UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 9027911415125877338UL) + ((uint64_t)op[1] * 17253757146333624087UL) + ((((uint64_t)op[2] * 18347504255334709693UL) + ((uint64_t)op[3] * 6479984006563990855UL) + ((uint64_t)op[4] * 5036921686197915088UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 5036921686197915088UL) + ((uint64_t)op[1] * 9027911415125877338UL) + ((uint64_t)op[2] * 17253757146333624087UL) + ((((uint64_t)op[3] * 18347504255334709693UL) + ((uint64_t)op[4] * 6479984006563990855UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 6479984006563990855UL) + ((uint64_t)op[1] * 5036921686197915088UL) + ((uint64_t)op[2] * 9027911415125877338UL) + ((uint64_t)op[3] * 17253757146333624087UL) + ((uint64_t)op[4] * 496199091874209615UL);
	tmp_q[4] = ((uint64_t)op[0] * 18347504255334709693UL) + ((uint64_t)op[1] * 6479984006563990855UL) + ((uint64_t)op[2] * 5036921686197915088UL) + ((uint64_t)op[3] * 9027911415125877338UL) + ((uint64_t)op[4] * 17253757146333624087UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 22892099954L) - ((((int128)tmp_q[1] * 139302163089L) - ((int128)tmp_q[2] * 166209480267L) - ((int128)tmp_q[3] * 136023699186L) - ((int128)tmp_q[4] * 49653684393L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 49653684393L) - ((int128)tmp_q[1] * 22892099954L) - ((((int128)tmp_q[2] * 139302163089L) - ((int128)tmp_q[3] * 166209480267L) - ((int128)tmp_q[4] * 136023699186L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 136023699186L) - ((int128)tmp_q[1] * 49653684393L) - ((int128)tmp_q[2] * 22892099954L) - ((((int128)tmp_q[3] * 139302163089L) - ((int128)tmp_q[4] * 166209480267L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 166209480267L) - ((int128)tmp_q[1] * 136023699186L) - ((int128)tmp_q[2] * 49653684393L) - ((int128)tmp_q[3] * 22892099954L) - ((int128)tmp_q[4] * 696510815445L);
	tmp_zero[4] = ((int128)tmp_q[0] * 139302163089L) - ((int128)tmp_q[1] * 166209480267L) - ((int128)tmp_q[2] * 136023699186L) - ((int128)tmp_q[3] * 49653684393L) - ((int128)tmp_q[4] * 22892099954L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

