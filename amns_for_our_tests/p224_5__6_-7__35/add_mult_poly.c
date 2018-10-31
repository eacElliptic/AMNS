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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12604627612719650387UL) + ((((uint64_t)op[1] * 13154888970187876582UL) + ((uint64_t)op[2] * 4778688337829736160UL) + ((uint64_t)op[3] * 814332132363204109UL) + ((uint64_t)op[4] * 13346388622619900489UL) + ((uint64_t)op[5] * 11747671193618827936UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 11747671193618827936UL) + ((uint64_t)op[1] * 12604627612719650387UL) + ((((uint64_t)op[2] * 13154888970187876582UL) + ((uint64_t)op[3] * 4778688337829736160UL) + ((uint64_t)op[4] * 814332132363204109UL) + ((uint64_t)op[5] * 13346388622619900489UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 13346388622619900489UL) + ((uint64_t)op[1] * 11747671193618827936UL) + ((uint64_t)op[2] * 12604627612719650387UL) + ((((uint64_t)op[3] * 13154888970187876582UL) + ((uint64_t)op[4] * 4778688337829736160UL) + ((uint64_t)op[5] * 814332132363204109UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 814332132363204109UL) + ((uint64_t)op[1] * 13346388622619900489UL) + ((uint64_t)op[2] * 11747671193618827936UL) + ((uint64_t)op[3] * 12604627612719650387UL) + ((((uint64_t)op[4] * 13154888970187876582UL) + ((uint64_t)op[5] * 4778688337829736160UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 4778688337829736160UL) + ((uint64_t)op[1] * 814332132363204109UL) + ((uint64_t)op[2] * 13346388622619900489UL) + ((uint64_t)op[3] * 11747671193618827936UL) + ((uint64_t)op[4] * 12604627612719650387UL) + ((uint64_t)op[5] * 149497577232622006UL);
	tmp_q[5] = ((uint64_t)op[0] * 13154888970187876582UL) + ((uint64_t)op[1] * 4778688337829736160UL) + ((uint64_t)op[2] * 814332132363204109UL) + ((uint64_t)op[3] * 13346388622619900489UL) + ((uint64_t)op[4] * 11747671193618827936UL) + ((uint64_t)op[5] * 12604627612719650387UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 69546247214L) - ((((int128)tmp_q[1] * 42207188773L) + ((int128)tmp_q[2] * 77667022831L) - ((int128)tmp_q[3] * 21681169658L) + ((int128)tmp_q[4] * 108628675259L) + ((int128)tmp_q[5] * 791115882L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 791115882L) - ((int128)tmp_q[1] * 69546247214L) - ((((int128)tmp_q[2] * 42207188773L) + ((int128)tmp_q[3] * 77667022831L) - ((int128)tmp_q[4] * 21681169658L) + ((int128)tmp_q[5] * 108628675259L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 108628675259L) + ((int128)tmp_q[1] * 791115882L) - ((int128)tmp_q[2] * 69546247214L) - ((((int128)tmp_q[3] * 42207188773L) + ((int128)tmp_q[4] * 77667022831L) - ((int128)tmp_q[5] * 21681169658L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 21681169658L) + ((int128)tmp_q[1] * 108628675259L) + ((int128)tmp_q[2] * 791115882L) - ((int128)tmp_q[3] * 69546247214L) - ((((int128)tmp_q[4] * 42207188773L) + ((int128)tmp_q[5] * 77667022831L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 77667022831L) - ((int128)tmp_q[1] * 21681169658L) + ((int128)tmp_q[2] * 108628675259L) + ((int128)tmp_q[3] * 791115882L) - ((int128)tmp_q[4] * 69546247214L) - ((int128)tmp_q[5] * 295450321411L);
	tmp_zero[5] = ((int128)tmp_q[0] * 42207188773L) + ((int128)tmp_q[1] * 77667022831L) - ((int128)tmp_q[2] * 21681169658L) + ((int128)tmp_q[3] * 108628675259L) + ((int128)tmp_q[4] * 791115882L) - ((int128)tmp_q[5] * 69546247214L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

