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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9586530944354143139UL) + ((((uint64_t)op[1] * 15285586328573300460UL) + ((uint64_t)op[2] * 12208075818940603434UL) + ((uint64_t)op[3] * 10026890669994123268UL) + ((uint64_t)op[4] * 9539525055937119270UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 9539525055937119270UL) + ((uint64_t)op[1] * 9586530944354143139UL) + ((((uint64_t)op[2] * 15285586328573300460UL) + ((uint64_t)op[3] * 12208075818940603434UL) + ((uint64_t)op[4] * 10026890669994123268UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 10026890669994123268UL) + ((uint64_t)op[1] * 9539525055937119270UL) + ((uint64_t)op[2] * 9586530944354143139UL) + ((((uint64_t)op[3] * 15285586328573300460UL) + ((uint64_t)op[4] * 12208075818940603434UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 12208075818940603434UL) + ((uint64_t)op[1] * 10026890669994123268UL) + ((uint64_t)op[2] * 9539525055937119270UL) + ((uint64_t)op[3] * 9586530944354143139UL) + ((uint64_t)op[4] * 14765383931465345140UL);
	tmp_q[4] = ((uint64_t)op[0] * 15285586328573300460UL) + ((uint64_t)op[1] * 12208075818940603434UL) + ((uint64_t)op[2] * 10026890669994123268UL) + ((uint64_t)op[3] * 9539525055937119270UL) + ((uint64_t)op[4] * 9586530944354143139UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1857514939137075L) + ((-((int128)tmp_q[1] * 925429728102356L) + ((int128)tmp_q[2] * 215354248904322L) + ((int128)tmp_q[3] * 734701406510864L) + ((int128)tmp_q[4] * 1385447535388754L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 1385447535388754L) - ((int128)tmp_q[1] * 1857514939137075L) + ((-((int128)tmp_q[2] * 925429728102356L) + ((int128)tmp_q[3] * 215354248904322L) + ((int128)tmp_q[4] * 734701406510864L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 734701406510864L) + ((int128)tmp_q[1] * 1385447535388754L) - ((int128)tmp_q[2] * 1857514939137075L) + ((-((int128)tmp_q[3] * 925429728102356L) + ((int128)tmp_q[4] * 215354248904322L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 215354248904322L) + ((int128)tmp_q[1] * 734701406510864L) + ((int128)tmp_q[2] * 1385447535388754L) - ((int128)tmp_q[3] * 1857514939137075L) - ((int128)tmp_q[4] * 6478008096716492L);
	tmp_zero[4] = -((int128)tmp_q[0] * 925429728102356L) + ((int128)tmp_q[1] * 215354248904322L) + ((int128)tmp_q[2] * 734701406510864L) + ((int128)tmp_q[3] * 1385447535388754L) - ((int128)tmp_q[4] * 1857514939137075L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

