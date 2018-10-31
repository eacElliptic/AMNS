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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1332852251530242303UL) + ((((uint64_t)op[1] * 7400953639656990962UL) + ((uint64_t)op[2] * 9788944336353377862UL) + ((uint64_t)op[3] * 5830034090827453182UL) + ((uint64_t)op[4] * 15054180281085851043UL) + ((uint64_t)op[5] * 7698363859669426020UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 7698363859669426020UL) + ((uint64_t)op[1] * 1332852251530242303UL) + ((((uint64_t)op[2] * 7400953639656990962UL) + ((uint64_t)op[3] * 9788944336353377862UL) + ((uint64_t)op[4] * 5830034090827453182UL) + ((uint64_t)op[5] * 15054180281085851043UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 15054180281085851043UL) + ((uint64_t)op[1] * 7698363859669426020UL) + ((uint64_t)op[2] * 1332852251530242303UL) + ((((uint64_t)op[3] * 7400953639656990962UL) + ((uint64_t)op[4] * 9788944336353377862UL) + ((uint64_t)op[5] * 5830034090827453182UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 5830034090827453182UL) + ((uint64_t)op[1] * 15054180281085851043UL) + ((uint64_t)op[2] * 7698363859669426020UL) + ((uint64_t)op[3] * 1332852251530242303UL) + ((((uint64_t)op[4] * 7400953639656990962UL) + ((uint64_t)op[5] * 9788944336353377862UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 9788944336353377862UL) + ((uint64_t)op[1] * 5830034090827453182UL) + ((uint64_t)op[2] * 15054180281085851043UL) + ((uint64_t)op[3] * 7698363859669426020UL) + ((uint64_t)op[4] * 1332852251530242303UL) + ((uint64_t)op[5] * 3644836794395569692UL);
	tmp_q[5] = ((uint64_t)op[0] * 7400953639656990962UL) + ((uint64_t)op[1] * 9788944336353377862UL) + ((uint64_t)op[2] * 5830034090827453182UL) + ((uint64_t)op[3] * 15054180281085851043UL) + ((uint64_t)op[4] * 7698363859669426020UL) + ((uint64_t)op[5] * 1332852251530242303UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 458352449L) - ((((int128)tmp_q[1] * 989701694L) - ((int128)tmp_q[2] * 200230531L) + ((int128)tmp_q[3] * 2815744618L) - ((int128)tmp_q[4] * 1166050351L) - ((int128)tmp_q[5] * 439610496L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 439610496L) - ((int128)tmp_q[1] * 458352449L) - ((((int128)tmp_q[2] * 989701694L) - ((int128)tmp_q[3] * 200230531L) + ((int128)tmp_q[4] * 2815744618L) - ((int128)tmp_q[5] * 1166050351L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 1166050351L) - ((int128)tmp_q[1] * 439610496L) - ((int128)tmp_q[2] * 458352449L) - ((((int128)tmp_q[3] * 989701694L) - ((int128)tmp_q[4] * 200230531L) + ((int128)tmp_q[5] * 2815744618L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 2815744618L) - ((int128)tmp_q[1] * 1166050351L) - ((int128)tmp_q[2] * 439610496L) - ((int128)tmp_q[3] * 458352449L) - ((((int128)tmp_q[4] * 989701694L) - ((int128)tmp_q[5] * 200230531L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 200230531L) + ((int128)tmp_q[1] * 2815744618L) - ((int128)tmp_q[2] * 1166050351L) - ((int128)tmp_q[3] * 439610496L) - ((int128)tmp_q[4] * 458352449L) - ((int128)tmp_q[5] * 1979403388L);
	tmp_zero[5] = ((int128)tmp_q[0] * 989701694L) - ((int128)tmp_q[1] * 200230531L) + ((int128)tmp_q[2] * 2815744618L) - ((int128)tmp_q[3] * 1166050351L) - ((int128)tmp_q[4] * 439610496L) - ((int128)tmp_q[5] * 458352449L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

