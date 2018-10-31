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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18005679977992584999UL) + ((((uint64_t)op[1] * 6732816038746592516UL) + ((uint64_t)op[2] * 2595083198489157826UL) + ((uint64_t)op[3] * 2527677030169184442UL) + ((uint64_t)op[4] * 1915797302469859288UL) + ((uint64_t)op[5] * 4947528299278130297UL) + ((uint64_t)op[6] * 10356253790886250974UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 10356253790886250974UL) + ((uint64_t)op[1] * 18005679977992584999UL) + ((((uint64_t)op[2] * 6732816038746592516UL) + ((uint64_t)op[3] * 2595083198489157826UL) + ((uint64_t)op[4] * 2527677030169184442UL) + ((uint64_t)op[5] * 1915797302469859288UL) + ((uint64_t)op[6] * 4947528299278130297UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 4947528299278130297UL) + ((uint64_t)op[1] * 10356253790886250974UL) + ((uint64_t)op[2] * 18005679977992584999UL) + ((((uint64_t)op[3] * 6732816038746592516UL) + ((uint64_t)op[4] * 2595083198489157826UL) + ((uint64_t)op[5] * 2527677030169184442UL) + ((uint64_t)op[6] * 1915797302469859288UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 1915797302469859288UL) + ((uint64_t)op[1] * 4947528299278130297UL) + ((uint64_t)op[2] * 10356253790886250974UL) + ((uint64_t)op[3] * 18005679977992584999UL) + ((((uint64_t)op[4] * 6732816038746592516UL) + ((uint64_t)op[5] * 2595083198489157826UL) + ((uint64_t)op[6] * 2527677030169184442UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 2527677030169184442UL) + ((uint64_t)op[1] * 1915797302469859288UL) + ((uint64_t)op[2] * 4947528299278130297UL) + ((uint64_t)op[3] * 10356253790886250974UL) + ((uint64_t)op[4] * 18005679977992584999UL) + ((((uint64_t)op[5] * 6732816038746592516UL) + ((uint64_t)op[6] * 2595083198489157826UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 2595083198489157826UL) + ((uint64_t)op[1] * 2527677030169184442UL) + ((uint64_t)op[2] * 1915797302469859288UL) + ((uint64_t)op[3] * 4947528299278130297UL) + ((uint64_t)op[4] * 10356253790886250974UL) + ((uint64_t)op[5] * 18005679977992584999UL) + ((uint64_t)op[6] * 13465632077493185032UL);
	tmp_q[6] = ((uint64_t)op[0] * 6732816038746592516UL) + ((uint64_t)op[1] * 2595083198489157826UL) + ((uint64_t)op[2] * 2527677030169184442UL) + ((uint64_t)op[3] * 1915797302469859288UL) + ((uint64_t)op[4] * 4947528299278130297UL) + ((uint64_t)op[5] * 10356253790886250974UL) + ((uint64_t)op[6] * 18005679977992584999UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2118248133L) + ((((int128)tmp_q[1] * 54087457229L) - ((int128)tmp_q[2] * 15604995278L) + ((int128)tmp_q[3] * 40214619883L) + ((int128)tmp_q[4] * 39432410918L) + ((int128)tmp_q[5] * 29558820289L) - ((int128)tmp_q[6] * 60823143812L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 60823143812L) + ((int128)tmp_q[1] * 2118248133L) + ((((int128)tmp_q[2] * 54087457229L) - ((int128)tmp_q[3] * 15604995278L) + ((int128)tmp_q[4] * 40214619883L) + ((int128)tmp_q[5] * 39432410918L) + ((int128)tmp_q[6] * 29558820289L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 29558820289L) - ((int128)tmp_q[1] * 60823143812L) + ((int128)tmp_q[2] * 2118248133L) + ((((int128)tmp_q[3] * 54087457229L) - ((int128)tmp_q[4] * 15604995278L) + ((int128)tmp_q[5] * 40214619883L) + ((int128)tmp_q[6] * 39432410918L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 39432410918L) + ((int128)tmp_q[1] * 29558820289L) - ((int128)tmp_q[2] * 60823143812L) + ((int128)tmp_q[3] * 2118248133L) + ((((int128)tmp_q[4] * 54087457229L) - ((int128)tmp_q[5] * 15604995278L) + ((int128)tmp_q[6] * 40214619883L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 40214619883L) + ((int128)tmp_q[1] * 39432410918L) + ((int128)tmp_q[2] * 29558820289L) - ((int128)tmp_q[3] * 60823143812L) + ((int128)tmp_q[4] * 2118248133L) + ((((int128)tmp_q[5] * 54087457229L) - ((int128)tmp_q[6] * 15604995278L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 15604995278L) + ((int128)tmp_q[1] * 40214619883L) + ((int128)tmp_q[2] * 39432410918L) + ((int128)tmp_q[3] * 29558820289L) - ((int128)tmp_q[4] * 60823143812L) + ((int128)tmp_q[5] * 2118248133L) + ((int128)tmp_q[6] * 108174914458L);
	tmp_zero[6] = ((int128)tmp_q[0] * 54087457229L) - ((int128)tmp_q[1] * 15604995278L) + ((int128)tmp_q[2] * 40214619883L) + ((int128)tmp_q[3] * 39432410918L) + ((int128)tmp_q[4] * 29558820289L) - ((int128)tmp_q[5] * 60823143812L) + ((int128)tmp_q[6] * 2118248133L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

