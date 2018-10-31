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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10807802383227392861UL) + ((((uint64_t)op[1] * 16653579646578055840UL) + ((uint64_t)op[2] * 11620478823915146104UL) + ((uint64_t)op[3] * 848713169091623505UL) + ((uint64_t)op[4] * 704100969385986626UL) + ((uint64_t)op[5] * 14673003289287387071UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 14673003289287387071UL) + ((uint64_t)op[1] * 10807802383227392861UL) + ((((uint64_t)op[2] * 16653579646578055840UL) + ((uint64_t)op[3] * 11620478823915146104UL) + ((uint64_t)op[4] * 848713169091623505UL) + ((uint64_t)op[5] * 704100969385986626UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 704100969385986626UL) + ((uint64_t)op[1] * 14673003289287387071UL) + ((uint64_t)op[2] * 10807802383227392861UL) + ((((uint64_t)op[3] * 16653579646578055840UL) + ((uint64_t)op[4] * 11620478823915146104UL) + ((uint64_t)op[5] * 848713169091623505UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 848713169091623505UL) + ((uint64_t)op[1] * 704100969385986626UL) + ((uint64_t)op[2] * 14673003289287387071UL) + ((uint64_t)op[3] * 10807802383227392861UL) + ((((uint64_t)op[4] * 16653579646578055840UL) + ((uint64_t)op[5] * 11620478823915146104UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 11620478823915146104UL) + ((uint64_t)op[1] * 848713169091623505UL) + ((uint64_t)op[2] * 704100969385986626UL) + ((uint64_t)op[3] * 14673003289287387071UL) + ((uint64_t)op[4] * 10807802383227392861UL) + ((uint64_t)op[5] * 9480921938052072736UL);
	tmp_q[5] = ((uint64_t)op[0] * 16653579646578055840UL) + ((uint64_t)op[1] * 11620478823915146104UL) + ((uint64_t)op[2] * 848713169091623505UL) + ((uint64_t)op[3] * 704100969385986626UL) + ((uint64_t)op[4] * 14673003289287387071UL) + ((uint64_t)op[5] * 10807802383227392861UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1290609199624L) + ((((int128)tmp_q[1] * 877183084321L) + ((int128)tmp_q[2] * 2991466486949L) - ((int128)tmp_q[3] * 4923945313038L) - ((int128)tmp_q[4] * 3844698520308L) + ((int128)tmp_q[5] * 212459643961L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 212459643961L) - ((int128)tmp_q[1] * 1290609199624L) + ((((int128)tmp_q[2] * 877183084321L) + ((int128)tmp_q[3] * 2991466486949L) - ((int128)tmp_q[4] * 4923945313038L) - ((int128)tmp_q[5] * 3844698520308L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 3844698520308L) + ((int128)tmp_q[1] * 212459643961L) - ((int128)tmp_q[2] * 1290609199624L) + ((((int128)tmp_q[3] * 877183084321L) + ((int128)tmp_q[4] * 2991466486949L) - ((int128)tmp_q[5] * 4923945313038L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 4923945313038L) - ((int128)tmp_q[1] * 3844698520308L) + ((int128)tmp_q[2] * 212459643961L) - ((int128)tmp_q[3] * 1290609199624L) + ((((int128)tmp_q[4] * 877183084321L) + ((int128)tmp_q[5] * 2991466486949L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 2991466486949L) - ((int128)tmp_q[1] * 4923945313038L) - ((int128)tmp_q[2] * 3844698520308L) + ((int128)tmp_q[3] * 212459643961L) - ((int128)tmp_q[4] * 1290609199624L) + ((int128)tmp_q[5] * 4385915421605L);
	tmp_zero[5] = ((int128)tmp_q[0] * 877183084321L) + ((int128)tmp_q[1] * 2991466486949L) - ((int128)tmp_q[2] * 4923945313038L) - ((int128)tmp_q[3] * 3844698520308L) + ((int128)tmp_q[4] * 212459643961L) - ((int128)tmp_q[5] * 1290609199624L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

