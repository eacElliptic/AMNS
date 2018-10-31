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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12049267256870079623UL) + ((((uint64_t)op[1] * 1864190256744609498UL) + ((uint64_t)op[2] * 3027664733266812212UL) + ((uint64_t)op[3] * 16041485510208712861UL) + ((uint64_t)op[4] * 11130272110782902774UL) + ((uint64_t)op[5] * 4703473291293941965UL) + ((uint64_t)op[6] * 3203456267554366155UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 3203456267554366155UL) + ((uint64_t)op[1] * 12049267256870079623UL) + ((((uint64_t)op[2] * 1864190256744609498UL) + ((uint64_t)op[3] * 3027664733266812212UL) + ((uint64_t)op[4] * 16041485510208712861UL) + ((uint64_t)op[5] * 11130272110782902774UL) + ((uint64_t)op[6] * 4703473291293941965UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 4703473291293941965UL) + ((uint64_t)op[1] * 3203456267554366155UL) + ((uint64_t)op[2] * 12049267256870079623UL) + ((((uint64_t)op[3] * 1864190256744609498UL) + ((uint64_t)op[4] * 3027664733266812212UL) + ((uint64_t)op[5] * 16041485510208712861UL) + ((uint64_t)op[6] * 11130272110782902774UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 11130272110782902774UL) + ((uint64_t)op[1] * 4703473291293941965UL) + ((uint64_t)op[2] * 3203456267554366155UL) + ((uint64_t)op[3] * 12049267256870079623UL) + ((((uint64_t)op[4] * 1864190256744609498UL) + ((uint64_t)op[5] * 3027664733266812212UL) + ((uint64_t)op[6] * 16041485510208712861UL)) * 4);
	tmp_q[4] = ((uint64_t)op[0] * 16041485510208712861UL) + ((uint64_t)op[1] * 11130272110782902774UL) + ((uint64_t)op[2] * 4703473291293941965UL) + ((uint64_t)op[3] * 3203456267554366155UL) + ((uint64_t)op[4] * 12049267256870079623UL) + ((((uint64_t)op[5] * 1864190256744609498UL) + ((uint64_t)op[6] * 3027664733266812212UL)) * 4);
	tmp_q[5] = ((uint64_t)op[0] * 3027664733266812212UL) + ((uint64_t)op[1] * 16041485510208712861UL) + ((uint64_t)op[2] * 11130272110782902774UL) + ((uint64_t)op[3] * 4703473291293941965UL) + ((uint64_t)op[4] * 3203456267554366155UL) + ((uint64_t)op[5] * 12049267256870079623UL) + ((uint64_t)op[6] * 7456761026978437992UL);
	tmp_q[6] = ((uint64_t)op[0] * 1864190256744609498UL) + ((uint64_t)op[1] * 3027664733266812212UL) + ((uint64_t)op[2] * 16041485510208712861UL) + ((uint64_t)op[3] * 11130272110782902774UL) + ((uint64_t)op[4] * 4703473291293941965UL) + ((uint64_t)op[5] * 3203456267554366155UL) + ((uint64_t)op[6] * 12049267256870079623UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3599983437L) + ((-((int128)tmp_q[1] * 8274232024L) - ((int128)tmp_q[2] * 32264893192L) - ((int128)tmp_q[3] * 57690655590L) + ((int128)tmp_q[4] * 42840821531L) - ((int128)tmp_q[5] * 33881238154L) + ((int128)tmp_q[6] * 59991359247L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 59991359247L) + ((int128)tmp_q[1] * 3599983437L) + ((-((int128)tmp_q[2] * 8274232024L) - ((int128)tmp_q[3] * 32264893192L) - ((int128)tmp_q[4] * 57690655590L) + ((int128)tmp_q[5] * 42840821531L) - ((int128)tmp_q[6] * 33881238154L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 33881238154L) + ((int128)tmp_q[1] * 59991359247L) + ((int128)tmp_q[2] * 3599983437L) + ((-((int128)tmp_q[3] * 8274232024L) - ((int128)tmp_q[4] * 32264893192L) - ((int128)tmp_q[5] * 57690655590L) + ((int128)tmp_q[6] * 42840821531L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 42840821531L) - ((int128)tmp_q[1] * 33881238154L) + ((int128)tmp_q[2] * 59991359247L) + ((int128)tmp_q[3] * 3599983437L) + ((-((int128)tmp_q[4] * 8274232024L) - ((int128)tmp_q[5] * 32264893192L) - ((int128)tmp_q[6] * 57690655590L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 57690655590L) + ((int128)tmp_q[1] * 42840821531L) - ((int128)tmp_q[2] * 33881238154L) + ((int128)tmp_q[3] * 59991359247L) + ((int128)tmp_q[4] * 3599983437L) + ((-((int128)tmp_q[5] * 8274232024L) - ((int128)tmp_q[6] * 32264893192L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 32264893192L) - ((int128)tmp_q[1] * 57690655590L) + ((int128)tmp_q[2] * 42840821531L) - ((int128)tmp_q[3] * 33881238154L) + ((int128)tmp_q[4] * 59991359247L) + ((int128)tmp_q[5] * 3599983437L) - ((int128)tmp_q[6] * 33096928096L);
	tmp_zero[6] = -((int128)tmp_q[0] * 8274232024L) - ((int128)tmp_q[1] * 32264893192L) - ((int128)tmp_q[2] * 57690655590L) + ((int128)tmp_q[3] * 42840821531L) - ((int128)tmp_q[4] * 33881238154L) + ((int128)tmp_q[5] * 59991359247L) + ((int128)tmp_q[6] * 3599983437L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

