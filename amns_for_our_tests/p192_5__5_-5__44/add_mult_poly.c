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
	tmp_q[0] = ((uint64_t)op[0] * 1944507284282689804UL) + ((((uint64_t)op[1] * 109462678467589714UL) + ((uint64_t)op[2] * 724849353384994048UL) + ((uint64_t)op[3] * 9509411241765944216UL) + ((uint64_t)op[4] * 16055563047416542191UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 16055563047416542191UL) + ((uint64_t)op[1] * 1944507284282689804UL) + ((((uint64_t)op[2] * 109462678467589714UL) + ((uint64_t)op[3] * 724849353384994048UL) + ((uint64_t)op[4] * 9509411241765944216UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 9509411241765944216UL) + ((uint64_t)op[1] * 16055563047416542191UL) + ((uint64_t)op[2] * 1944507284282689804UL) + ((((uint64_t)op[3] * 109462678467589714UL) + ((uint64_t)op[4] * 724849353384994048UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 724849353384994048UL) + ((uint64_t)op[1] * 9509411241765944216UL) + ((uint64_t)op[2] * 16055563047416542191UL) + ((uint64_t)op[3] * 1944507284282689804UL) + ((uint64_t)op[4] * 17899430681371603046UL);
	tmp_q[4] = ((uint64_t)op[0] * 109462678467589714UL) + ((uint64_t)op[1] * 724849353384994048UL) + ((uint64_t)op[2] * 9509411241765944216UL) + ((uint64_t)op[3] * 16055563047416542191UL) + ((uint64_t)op[4] * 1944507284282689804UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 182546839972L) - ((-((int128)tmp_q[1] * 76811744013L) - ((int128)tmp_q[2] * 60799347108L) + ((int128)tmp_q[3] * 133505619970L) + ((int128)tmp_q[4] * 57753965344L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 57753965344L) + ((int128)tmp_q[1] * 182546839972L) - ((-((int128)tmp_q[2] * 76811744013L) - ((int128)tmp_q[3] * 60799347108L) + ((int128)tmp_q[4] * 133505619970L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 133505619970L) + ((int128)tmp_q[1] * 57753965344L) + ((int128)tmp_q[2] * 182546839972L) - ((-((int128)tmp_q[3] * 76811744013L) - ((int128)tmp_q[4] * 60799347108L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 60799347108L) + ((int128)tmp_q[1] * 133505619970L) + ((int128)tmp_q[2] * 57753965344L) + ((int128)tmp_q[3] * 182546839972L) + ((int128)tmp_q[4] * 384058720065L);
	tmp_zero[4] = -((int128)tmp_q[0] * 76811744013L) - ((int128)tmp_q[1] * 60799347108L) + ((int128)tmp_q[2] * 133505619970L) + ((int128)tmp_q[3] * 57753965344L) + ((int128)tmp_q[4] * 182546839972L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

