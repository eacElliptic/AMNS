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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 12);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17470229705395733405UL) + ((((uint64_t)op[1] * 15532517842973838724UL) + ((uint64_t)op[2] * 14752980371075837620UL) + ((uint64_t)op[3] * 4399058119748821886UL) + ((uint64_t)op[4] * 7508204202181736297UL) + ((uint64_t)op[5] * 1428179629751378636UL) + ((uint64_t)op[6] * 2203435527871082494UL) + ((uint64_t)op[7] * 7657961395254077391UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 7657961395254077391UL) + ((uint64_t)op[1] * 17470229705395733405UL) + ((((uint64_t)op[2] * 15532517842973838724UL) + ((uint64_t)op[3] * 14752980371075837620UL) + ((uint64_t)op[4] * 4399058119748821886UL) + ((uint64_t)op[5] * 7508204202181736297UL) + ((uint64_t)op[6] * 1428179629751378636UL) + ((uint64_t)op[7] * 2203435527871082494UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 2203435527871082494UL) + ((uint64_t)op[1] * 7657961395254077391UL) + ((uint64_t)op[2] * 17470229705395733405UL) + ((((uint64_t)op[3] * 15532517842973838724UL) + ((uint64_t)op[4] * 14752980371075837620UL) + ((uint64_t)op[5] * 4399058119748821886UL) + ((uint64_t)op[6] * 7508204202181736297UL) + ((uint64_t)op[7] * 1428179629751378636UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 1428179629751378636UL) + ((uint64_t)op[1] * 2203435527871082494UL) + ((uint64_t)op[2] * 7657961395254077391UL) + ((uint64_t)op[3] * 17470229705395733405UL) + ((((uint64_t)op[4] * 15532517842973838724UL) + ((uint64_t)op[5] * 14752980371075837620UL) + ((uint64_t)op[6] * 4399058119748821886UL) + ((uint64_t)op[7] * 7508204202181736297UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 7508204202181736297UL) + ((uint64_t)op[1] * 1428179629751378636UL) + ((uint64_t)op[2] * 2203435527871082494UL) + ((uint64_t)op[3] * 7657961395254077391UL) + ((uint64_t)op[4] * 17470229705395733405UL) + ((((uint64_t)op[5] * 15532517842973838724UL) + ((uint64_t)op[6] * 14752980371075837620UL) + ((uint64_t)op[7] * 4399058119748821886UL)) * 6);
	tmp_q[5] = ((uint64_t)op[0] * 4399058119748821886UL) + ((uint64_t)op[1] * 7508204202181736297UL) + ((uint64_t)op[2] * 1428179629751378636UL) + ((uint64_t)op[3] * 2203435527871082494UL) + ((uint64_t)op[4] * 7657961395254077391UL) + ((uint64_t)op[5] * 17470229705395733405UL) + ((((uint64_t)op[6] * 15532517842973838724UL) + ((uint64_t)op[7] * 14752980371075837620UL)) * 6);
	tmp_q[6] = ((uint64_t)op[0] * 14752980371075837620UL) + ((uint64_t)op[1] * 4399058119748821886UL) + ((uint64_t)op[2] * 7508204202181736297UL) + ((uint64_t)op[3] * 1428179629751378636UL) + ((uint64_t)op[4] * 2203435527871082494UL) + ((uint64_t)op[5] * 7657961395254077391UL) + ((uint64_t)op[6] * 17470229705395733405UL) + ((uint64_t)op[7] * 961386689295274264UL);
	tmp_q[7] = ((uint64_t)op[0] * 15532517842973838724UL) + ((uint64_t)op[1] * 14752980371075837620UL) + ((uint64_t)op[2] * 4399058119748821886UL) + ((uint64_t)op[3] * 7508204202181736297UL) + ((uint64_t)op[4] * 1428179629751378636UL) + ((uint64_t)op[5] * 2203435527871082494UL) + ((uint64_t)op[6] * 7657961395254077391UL) + ((uint64_t)op[7] * 17470229705395733405UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 16585483307485L) + ((-((int128)tmp_q[1] * 62602581045897L) + ((int128)tmp_q[2] * 3168492922356L) - ((int128)tmp_q[3] * 18602018818473L) + ((int128)tmp_q[4] * 109207065768134L) + ((int128)tmp_q[5] * 164004945202577L) - ((int128)tmp_q[6] * 65636387480059L) + ((int128)tmp_q[7] * 114081497853579L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 114081497853579L) + ((int128)tmp_q[1] * 16585483307485L) + ((-((int128)tmp_q[2] * 62602581045897L) + ((int128)tmp_q[3] * 3168492922356L) - ((int128)tmp_q[4] * 18602018818473L) + ((int128)tmp_q[5] * 109207065768134L) + ((int128)tmp_q[6] * 164004945202577L) - ((int128)tmp_q[7] * 65636387480059L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 65636387480059L) + ((int128)tmp_q[1] * 114081497853579L) + ((int128)tmp_q[2] * 16585483307485L) + ((-((int128)tmp_q[3] * 62602581045897L) + ((int128)tmp_q[4] * 3168492922356L) - ((int128)tmp_q[5] * 18602018818473L) + ((int128)tmp_q[6] * 109207065768134L) + ((int128)tmp_q[7] * 164004945202577L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 164004945202577L) - ((int128)tmp_q[1] * 65636387480059L) + ((int128)tmp_q[2] * 114081497853579L) + ((int128)tmp_q[3] * 16585483307485L) + ((-((int128)tmp_q[4] * 62602581045897L) + ((int128)tmp_q[5] * 3168492922356L) - ((int128)tmp_q[6] * 18602018818473L) + ((int128)tmp_q[7] * 109207065768134L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 109207065768134L) + ((int128)tmp_q[1] * 164004945202577L) - ((int128)tmp_q[2] * 65636387480059L) + ((int128)tmp_q[3] * 114081497853579L) + ((int128)tmp_q[4] * 16585483307485L) + ((-((int128)tmp_q[5] * 62602581045897L) + ((int128)tmp_q[6] * 3168492922356L) - ((int128)tmp_q[7] * 18602018818473L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 18602018818473L) + ((int128)tmp_q[1] * 109207065768134L) + ((int128)tmp_q[2] * 164004945202577L) - ((int128)tmp_q[3] * 65636387480059L) + ((int128)tmp_q[4] * 114081497853579L) + ((int128)tmp_q[5] * 16585483307485L) + ((-((int128)tmp_q[6] * 62602581045897L) + ((int128)tmp_q[7] * 3168492922356L)) * 6);
	tmp_zero[6] = ((int128)tmp_q[0] * 3168492922356L) - ((int128)tmp_q[1] * 18602018818473L) + ((int128)tmp_q[2] * 109207065768134L) + ((int128)tmp_q[3] * 164004945202577L) - ((int128)tmp_q[4] * 65636387480059L) + ((int128)tmp_q[5] * 114081497853579L) + ((int128)tmp_q[6] * 16585483307485L) - ((int128)tmp_q[7] * 375615486275382L);
	tmp_zero[7] = -((int128)tmp_q[0] * 62602581045897L) + ((int128)tmp_q[1] * 3168492922356L) - ((int128)tmp_q[2] * 18602018818473L) + ((int128)tmp_q[3] * 109207065768134L) + ((int128)tmp_q[4] * 164004945202577L) - ((int128)tmp_q[5] * 65636387480059L) + ((int128)tmp_q[6] * 114081497853579L) + ((int128)tmp_q[7] * 16585483307485L);

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

