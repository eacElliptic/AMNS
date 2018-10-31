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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11888472257379557390UL) + ((((uint64_t)op[1] * 874620960336474022UL) + ((uint64_t)op[2] * 9609871584561477466UL) + ((uint64_t)op[3] * 14315182500685834588UL) + ((uint64_t)op[4] * 3415890777584446117UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 3415890777584446117UL) + ((uint64_t)op[1] * 11888472257379557390UL) + ((((uint64_t)op[2] * 874620960336474022UL) + ((uint64_t)op[3] * 9609871584561477466UL) + ((uint64_t)op[4] * 14315182500685834588UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 14315182500685834588UL) + ((uint64_t)op[1] * 3415890777584446117UL) + ((uint64_t)op[2] * 11888472257379557390UL) + ((((uint64_t)op[3] * 874620960336474022UL) + ((uint64_t)op[4] * 9609871584561477466UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 9609871584561477466UL) + ((uint64_t)op[1] * 14315182500685834588UL) + ((uint64_t)op[2] * 3415890777584446117UL) + ((uint64_t)op[3] * 11888472257379557390UL) + ((uint64_t)op[4] * 2623862881009422066UL);
	tmp_q[4] = ((uint64_t)op[0] * 874620960336474022UL) + ((uint64_t)op[1] * 9609871584561477466UL) + ((uint64_t)op[2] * 14315182500685834588UL) + ((uint64_t)op[3] * 3415890777584446117UL) + ((uint64_t)op[4] * 11888472257379557390UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 103456836088L) + ((((int128)tmp_q[1] * 180400154089L) - ((int128)tmp_q[2] * 259046304242L) - ((int128)tmp_q[3] * 1172992150L) + ((int128)tmp_q[4] * 143323924370L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 143323924370L) + ((int128)tmp_q[1] * 103456836088L) + ((((int128)tmp_q[2] * 180400154089L) - ((int128)tmp_q[3] * 259046304242L) - ((int128)tmp_q[4] * 1172992150L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 1172992150L) + ((int128)tmp_q[1] * 143323924370L) + ((int128)tmp_q[2] * 103456836088L) + ((((int128)tmp_q[3] * 180400154089L) - ((int128)tmp_q[4] * 259046304242L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 259046304242L) - ((int128)tmp_q[1] * 1172992150L) + ((int128)tmp_q[2] * 143323924370L) + ((int128)tmp_q[3] * 103456836088L) + ((int128)tmp_q[4] * 541200462267L);
	tmp_zero[4] = ((int128)tmp_q[0] * 180400154089L) - ((int128)tmp_q[1] * 259046304242L) - ((int128)tmp_q[2] * 1172992150L) + ((int128)tmp_q[3] * 143323924370L) + ((int128)tmp_q[4] * 103456836088L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

