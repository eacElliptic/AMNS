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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5394687908824669695UL) + ((((uint64_t)op[1] * 13390047099021006630UL) + ((uint64_t)op[2] * 248928498935051015UL) + ((uint64_t)op[3] * 13095035533029186881UL) + ((uint64_t)op[4] * 6907916571583774306UL) + ((uint64_t)op[5] * 672129539317006598UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 672129539317006598UL) + ((uint64_t)op[1] * 5394687908824669695UL) + ((((uint64_t)op[2] * 13390047099021006630UL) + ((uint64_t)op[3] * 248928498935051015UL) + ((uint64_t)op[4] * 13095035533029186881UL) + ((uint64_t)op[5] * 6907916571583774306UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 6907916571583774306UL) + ((uint64_t)op[1] * 672129539317006598UL) + ((uint64_t)op[2] * 5394687908824669695UL) + ((((uint64_t)op[3] * 13390047099021006630UL) + ((uint64_t)op[4] * 248928498935051015UL) + ((uint64_t)op[5] * 13095035533029186881UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 13095035533029186881UL) + ((uint64_t)op[1] * 6907916571583774306UL) + ((uint64_t)op[2] * 672129539317006598UL) + ((uint64_t)op[3] * 5394687908824669695UL) + ((((uint64_t)op[4] * 13390047099021006630UL) + ((uint64_t)op[5] * 248928498935051015UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 248928498935051015UL) + ((uint64_t)op[1] * 13095035533029186881UL) + ((uint64_t)op[2] * 6907916571583774306UL) + ((uint64_t)op[3] * 672129539317006598UL) + ((uint64_t)op[4] * 5394687908824669695UL) + ((uint64_t)op[5] * 6836740799733173314UL);
	tmp_q[5] = ((uint64_t)op[0] * 13390047099021006630UL) + ((uint64_t)op[1] * 248928498935051015UL) + ((uint64_t)op[2] * 13095035533029186881UL) + ((uint64_t)op[3] * 6907916571583774306UL) + ((uint64_t)op[4] * 672129539317006598UL) + ((uint64_t)op[5] * 5394687908824669695UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 201929110742L) - ((-((int128)tmp_q[1] * 38894995520L) + ((int128)tmp_q[2] * 1946430409799L) - ((int128)tmp_q[3] * 632824306484L) - ((int128)tmp_q[4] * 1180869633919L) - ((int128)tmp_q[5] * 1075432817225L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 1075432817225L) - ((int128)tmp_q[1] * 201929110742L) - ((-((int128)tmp_q[2] * 38894995520L) + ((int128)tmp_q[3] * 1946430409799L) - ((int128)tmp_q[4] * 632824306484L) - ((int128)tmp_q[5] * 1180869633919L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 1180869633919L) - ((int128)tmp_q[1] * 1075432817225L) - ((int128)tmp_q[2] * 201929110742L) - ((-((int128)tmp_q[3] * 38894995520L) + ((int128)tmp_q[4] * 1946430409799L) - ((int128)tmp_q[5] * 632824306484L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 632824306484L) - ((int128)tmp_q[1] * 1180869633919L) - ((int128)tmp_q[2] * 1075432817225L) - ((int128)tmp_q[3] * 201929110742L) - ((-((int128)tmp_q[4] * 38894995520L) + ((int128)tmp_q[5] * 1946430409799L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 1946430409799L) - ((int128)tmp_q[1] * 632824306484L) - ((int128)tmp_q[2] * 1180869633919L) - ((int128)tmp_q[3] * 1075432817225L) - ((int128)tmp_q[4] * 201929110742L) + ((int128)tmp_q[5] * 194474977600L);
	tmp_zero[5] = -((int128)tmp_q[0] * 38894995520L) + ((int128)tmp_q[1] * 1946430409799L) - ((int128)tmp_q[2] * 632824306484L) - ((int128)tmp_q[3] * 1180869633919L) - ((int128)tmp_q[4] * 1075432817225L) - ((int128)tmp_q[5] * 201929110742L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

