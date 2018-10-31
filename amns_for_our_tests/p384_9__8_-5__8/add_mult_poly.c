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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 10);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4264656876667743041UL) + ((((uint64_t)op[1] * 3516525646781885433UL) + ((uint64_t)op[2] * 16517192504157867571UL) + ((uint64_t)op[3] * 3396433451930822368UL) + ((uint64_t)op[4] * 14862203410043898454UL) + ((uint64_t)op[5] * 16516201265497537573UL) + ((uint64_t)op[6] * 14302250310894396953UL) + ((uint64_t)op[7] * 2291464230428305948UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 2291464230428305948UL) + ((uint64_t)op[1] * 4264656876667743041UL) + ((((uint64_t)op[2] * 3516525646781885433UL) + ((uint64_t)op[3] * 16517192504157867571UL) + ((uint64_t)op[4] * 3396433451930822368UL) + ((uint64_t)op[5] * 14862203410043898454UL) + ((uint64_t)op[6] * 16516201265497537573UL) + ((uint64_t)op[7] * 14302250310894396953UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 14302250310894396953UL) + ((uint64_t)op[1] * 2291464230428305948UL) + ((uint64_t)op[2] * 4264656876667743041UL) + ((((uint64_t)op[3] * 3516525646781885433UL) + ((uint64_t)op[4] * 16517192504157867571UL) + ((uint64_t)op[5] * 3396433451930822368UL) + ((uint64_t)op[6] * 14862203410043898454UL) + ((uint64_t)op[7] * 16516201265497537573UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 16516201265497537573UL) + ((uint64_t)op[1] * 14302250310894396953UL) + ((uint64_t)op[2] * 2291464230428305948UL) + ((uint64_t)op[3] * 4264656876667743041UL) + ((((uint64_t)op[4] * 3516525646781885433UL) + ((uint64_t)op[5] * 16517192504157867571UL) + ((uint64_t)op[6] * 3396433451930822368UL) + ((uint64_t)op[7] * 14862203410043898454UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 14862203410043898454UL) + ((uint64_t)op[1] * 16516201265497537573UL) + ((uint64_t)op[2] * 14302250310894396953UL) + ((uint64_t)op[3] * 2291464230428305948UL) + ((uint64_t)op[4] * 4264656876667743041UL) + ((((uint64_t)op[5] * 3516525646781885433UL) + ((uint64_t)op[6] * 16517192504157867571UL) + ((uint64_t)op[7] * 3396433451930822368UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 3396433451930822368UL) + ((uint64_t)op[1] * 14862203410043898454UL) + ((uint64_t)op[2] * 16516201265497537573UL) + ((uint64_t)op[3] * 14302250310894396953UL) + ((uint64_t)op[4] * 2291464230428305948UL) + ((uint64_t)op[5] * 4264656876667743041UL) + ((((uint64_t)op[6] * 3516525646781885433UL) + ((uint64_t)op[7] * 16517192504157867571UL)) * 18446744073709551611);
	tmp_q[6] = ((uint64_t)op[0] * 16517192504157867571UL) + ((uint64_t)op[1] * 3396433451930822368UL) + ((uint64_t)op[2] * 14862203410043898454UL) + ((uint64_t)op[3] * 16516201265497537573UL) + ((uint64_t)op[4] * 14302250310894396953UL) + ((uint64_t)op[5] * 2291464230428305948UL) + ((uint64_t)op[6] * 4264656876667743041UL) + ((uint64_t)op[7] * 864115839800124451UL);
	tmp_q[7] = ((uint64_t)op[0] * 3516525646781885433UL) + ((uint64_t)op[1] * 16517192504157867571UL) + ((uint64_t)op[2] * 3396433451930822368UL) + ((uint64_t)op[3] * 14862203410043898454UL) + ((uint64_t)op[4] * 16516201265497537573UL) + ((uint64_t)op[5] * 14302250310894396953UL) + ((uint64_t)op[6] * 2291464230428305948UL) + ((uint64_t)op[7] * 4264656876667743041UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 20998768949603L) - ((-((int128)tmp_q[1] * 171180078951421L) - ((int128)tmp_q[2] * 127718445788275L) + ((int128)tmp_q[3] * 163975002531550L) - ((int128)tmp_q[4] * 89853312953312L) + ((int128)tmp_q[5] * 105083557347175L) - ((int128)tmp_q[6] * 132621320272631L) - ((int128)tmp_q[7] * 163278761810614L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 163278761810614L) + ((int128)tmp_q[1] * 20998768949603L) - ((-((int128)tmp_q[2] * 171180078951421L) - ((int128)tmp_q[3] * 127718445788275L) + ((int128)tmp_q[4] * 163975002531550L) - ((int128)tmp_q[5] * 89853312953312L) + ((int128)tmp_q[6] * 105083557347175L) - ((int128)tmp_q[7] * 132621320272631L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 132621320272631L) - ((int128)tmp_q[1] * 163278761810614L) + ((int128)tmp_q[2] * 20998768949603L) - ((-((int128)tmp_q[3] * 171180078951421L) - ((int128)tmp_q[4] * 127718445788275L) + ((int128)tmp_q[5] * 163975002531550L) - ((int128)tmp_q[6] * 89853312953312L) + ((int128)tmp_q[7] * 105083557347175L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 105083557347175L) - ((int128)tmp_q[1] * 132621320272631L) - ((int128)tmp_q[2] * 163278761810614L) + ((int128)tmp_q[3] * 20998768949603L) - ((-((int128)tmp_q[4] * 171180078951421L) - ((int128)tmp_q[5] * 127718445788275L) + ((int128)tmp_q[6] * 163975002531550L) - ((int128)tmp_q[7] * 89853312953312L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 89853312953312L) + ((int128)tmp_q[1] * 105083557347175L) - ((int128)tmp_q[2] * 132621320272631L) - ((int128)tmp_q[3] * 163278761810614L) + ((int128)tmp_q[4] * 20998768949603L) - ((-((int128)tmp_q[5] * 171180078951421L) - ((int128)tmp_q[6] * 127718445788275L) + ((int128)tmp_q[7] * 163975002531550L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 163975002531550L) - ((int128)tmp_q[1] * 89853312953312L) + ((int128)tmp_q[2] * 105083557347175L) - ((int128)tmp_q[3] * 132621320272631L) - ((int128)tmp_q[4] * 163278761810614L) + ((int128)tmp_q[5] * 20998768949603L) - ((-((int128)tmp_q[6] * 171180078951421L) - ((int128)tmp_q[7] * 127718445788275L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 127718445788275L) + ((int128)tmp_q[1] * 163975002531550L) - ((int128)tmp_q[2] * 89853312953312L) + ((int128)tmp_q[3] * 105083557347175L) - ((int128)tmp_q[4] * 132621320272631L) - ((int128)tmp_q[5] * 163278761810614L) + ((int128)tmp_q[6] * 20998768949603L) + ((int128)tmp_q[7] * 855900394757105L);
	tmp_zero[7] = -((int128)tmp_q[0] * 171180078951421L) - ((int128)tmp_q[1] * 127718445788275L) + ((int128)tmp_q[2] * 163975002531550L) - ((int128)tmp_q[3] * 89853312953312L) + ((int128)tmp_q[4] * 105083557347175L) - ((int128)tmp_q[5] * 132621320272631L) - ((int128)tmp_q[6] * 163278761810614L) + ((int128)tmp_q[7] * 20998768949603L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

