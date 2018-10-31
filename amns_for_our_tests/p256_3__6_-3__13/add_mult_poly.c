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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 400482793683452116UL) + ((((uint64_t)op[1] * 363292294635378755UL) + ((uint64_t)op[2] * 376339896896434186UL) + ((uint64_t)op[3] * 10899005823500871748UL) + ((uint64_t)op[4] * 10958077207646776032UL) + ((uint64_t)op[5] * 14553881571584550236UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 14553881571584550236UL) + ((uint64_t)op[1] * 400482793683452116UL) + ((((uint64_t)op[2] * 363292294635378755UL) + ((uint64_t)op[3] * 376339896896434186UL) + ((uint64_t)op[4] * 10899005823500871748UL) + ((uint64_t)op[5] * 10958077207646776032UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 10958077207646776032UL) + ((uint64_t)op[1] * 14553881571584550236UL) + ((uint64_t)op[2] * 400482793683452116UL) + ((((uint64_t)op[3] * 363292294635378755UL) + ((uint64_t)op[4] * 376339896896434186UL) + ((uint64_t)op[5] * 10899005823500871748UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 10899005823500871748UL) + ((uint64_t)op[1] * 10958077207646776032UL) + ((uint64_t)op[2] * 14553881571584550236UL) + ((uint64_t)op[3] * 400482793683452116UL) + ((((uint64_t)op[4] * 363292294635378755UL) + ((uint64_t)op[5] * 376339896896434186UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 376339896896434186UL) + ((uint64_t)op[1] * 10899005823500871748UL) + ((uint64_t)op[2] * 10958077207646776032UL) + ((uint64_t)op[3] * 14553881571584550236UL) + ((uint64_t)op[4] * 400482793683452116UL) + ((uint64_t)op[5] * 17356867189803415351UL);
	tmp_q[5] = ((uint64_t)op[0] * 363292294635378755UL) + ((uint64_t)op[1] * 376339896896434186UL) + ((uint64_t)op[2] * 10899005823500871748UL) + ((uint64_t)op[3] * 10958077207646776032UL) + ((uint64_t)op[4] * 14553881571584550236UL) + ((uint64_t)op[5] * 400482793683452116UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2806306640302L) - ((-((int128)tmp_q[1] * 1998694738936L) - ((int128)tmp_q[2] * 909092519784L) - ((int128)tmp_q[3] * 2089849352964L) - ((int128)tmp_q[4] * 562736584476L) + ((int128)tmp_q[5] * 1108392432697L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 1108392432697L) - ((int128)tmp_q[1] * 2806306640302L) - ((-((int128)tmp_q[2] * 1998694738936L) - ((int128)tmp_q[3] * 909092519784L) - ((int128)tmp_q[4] * 2089849352964L) - ((int128)tmp_q[5] * 562736584476L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 562736584476L) + ((int128)tmp_q[1] * 1108392432697L) - ((int128)tmp_q[2] * 2806306640302L) - ((-((int128)tmp_q[3] * 1998694738936L) - ((int128)tmp_q[4] * 909092519784L) - ((int128)tmp_q[5] * 2089849352964L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 2089849352964L) - ((int128)tmp_q[1] * 562736584476L) + ((int128)tmp_q[2] * 1108392432697L) - ((int128)tmp_q[3] * 2806306640302L) - ((-((int128)tmp_q[4] * 1998694738936L) - ((int128)tmp_q[5] * 909092519784L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 909092519784L) - ((int128)tmp_q[1] * 2089849352964L) - ((int128)tmp_q[2] * 562736584476L) + ((int128)tmp_q[3] * 1108392432697L) - ((int128)tmp_q[4] * 2806306640302L) + ((int128)tmp_q[5] * 5996084216808L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1998694738936L) - ((int128)tmp_q[1] * 909092519784L) - ((int128)tmp_q[2] * 2089849352964L) - ((int128)tmp_q[3] * 562736584476L) + ((int128)tmp_q[4] * 1108392432697L) - ((int128)tmp_q[5] * 2806306640302L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

