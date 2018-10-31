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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14097315190055142593UL) + ((((uint64_t)op[1] * 1510578797264576295UL) + ((uint64_t)op[2] * 15811374148742718950UL) + ((uint64_t)op[3] * 13266087004496956037UL) + ((uint64_t)op[4] * 9173293526995768907UL) + ((uint64_t)op[5] * 16625677194410474715UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 16625677194410474715UL) + ((uint64_t)op[1] * 14097315190055142593UL) + ((((uint64_t)op[2] * 1510578797264576295UL) + ((uint64_t)op[3] * 15811374148742718950UL) + ((uint64_t)op[4] * 13266087004496956037UL) + ((uint64_t)op[5] * 9173293526995768907UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 9173293526995768907UL) + ((uint64_t)op[1] * 16625677194410474715UL) + ((uint64_t)op[2] * 14097315190055142593UL) + ((((uint64_t)op[3] * 1510578797264576295UL) + ((uint64_t)op[4] * 15811374148742718950UL) + ((uint64_t)op[5] * 13266087004496956037UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 13266087004496956037UL) + ((uint64_t)op[1] * 9173293526995768907UL) + ((uint64_t)op[2] * 16625677194410474715UL) + ((uint64_t)op[3] * 14097315190055142593UL) + ((((uint64_t)op[4] * 1510578797264576295UL) + ((uint64_t)op[5] * 15811374148742718950UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 15811374148742718950UL) + ((uint64_t)op[1] * 13266087004496956037UL) + ((uint64_t)op[2] * 9173293526995768907UL) + ((uint64_t)op[3] * 16625677194410474715UL) + ((uint64_t)op[4] * 14097315190055142593UL) + ((uint64_t)op[5] * 4531736391793728885UL);
	tmp_q[5] = ((uint64_t)op[0] * 1510578797264576295UL) + ((uint64_t)op[1] * 15811374148742718950UL) + ((uint64_t)op[2] * 13266087004496956037UL) + ((uint64_t)op[3] * 9173293526995768907UL) + ((uint64_t)op[4] * 16625677194410474715UL) + ((uint64_t)op[5] * 14097315190055142593UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 98448110239L) + ((-((int128)tmp_q[1] * 14165170153L) - ((int128)tmp_q[2] * 32224917265L) + ((int128)tmp_q[3] * 57104147465L) + ((int128)tmp_q[4] * 52516208446L) + ((int128)tmp_q[5] * 88910914667L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 88910914667L) + ((int128)tmp_q[1] * 98448110239L) + ((-((int128)tmp_q[2] * 14165170153L) - ((int128)tmp_q[3] * 32224917265L) + ((int128)tmp_q[4] * 57104147465L) + ((int128)tmp_q[5] * 52516208446L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 52516208446L) + ((int128)tmp_q[1] * 88910914667L) + ((int128)tmp_q[2] * 98448110239L) + ((-((int128)tmp_q[3] * 14165170153L) - ((int128)tmp_q[4] * 32224917265L) + ((int128)tmp_q[5] * 57104147465L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 57104147465L) + ((int128)tmp_q[1] * 52516208446L) + ((int128)tmp_q[2] * 88910914667L) + ((int128)tmp_q[3] * 98448110239L) + ((-((int128)tmp_q[4] * 14165170153L) - ((int128)tmp_q[5] * 32224917265L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 32224917265L) + ((int128)tmp_q[1] * 57104147465L) + ((int128)tmp_q[2] * 52516208446L) + ((int128)tmp_q[3] * 88910914667L) + ((int128)tmp_q[4] * 98448110239L) - ((int128)tmp_q[5] * 42495510459L);
	tmp_zero[5] = -((int128)tmp_q[0] * 14165170153L) - ((int128)tmp_q[1] * 32224917265L) + ((int128)tmp_q[2] * 57104147465L) + ((int128)tmp_q[3] * 52516208446L) + ((int128)tmp_q[4] * 88910914667L) + ((int128)tmp_q[5] * 98448110239L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

