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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 984514166811453756UL) + ((((uint64_t)op[1] * 18170997079008614544UL) + ((uint64_t)op[2] * 4384447832457150726UL) + ((uint64_t)op[3] * 74571241062417764UL) + ((uint64_t)op[4] * 707475017012242111UL) + ((uint64_t)op[5] * 1149181324167426293UL) + ((uint64_t)op[6] * 14681101852772758639UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 14681101852772758639UL) + ((uint64_t)op[1] * 984514166811453756UL) + ((((uint64_t)op[2] * 18170997079008614544UL) + ((uint64_t)op[3] * 4384447832457150726UL) + ((uint64_t)op[4] * 74571241062417764UL) + ((uint64_t)op[5] * 707475017012242111UL) + ((uint64_t)op[6] * 1149181324167426293UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 1149181324167426293UL) + ((uint64_t)op[1] * 14681101852772758639UL) + ((uint64_t)op[2] * 984514166811453756UL) + ((((uint64_t)op[3] * 18170997079008614544UL) + ((uint64_t)op[4] * 4384447832457150726UL) + ((uint64_t)op[5] * 74571241062417764UL) + ((uint64_t)op[6] * 707475017012242111UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 707475017012242111UL) + ((uint64_t)op[1] * 1149181324167426293UL) + ((uint64_t)op[2] * 14681101852772758639UL) + ((uint64_t)op[3] * 984514166811453756UL) + ((((uint64_t)op[4] * 18170997079008614544UL) + ((uint64_t)op[5] * 4384447832457150726UL) + ((uint64_t)op[6] * 74571241062417764UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 74571241062417764UL) + ((uint64_t)op[1] * 707475017012242111UL) + ((uint64_t)op[2] * 1149181324167426293UL) + ((uint64_t)op[3] * 14681101852772758639UL) + ((uint64_t)op[4] * 984514166811453756UL) + ((((uint64_t)op[5] * 18170997079008614544UL) + ((uint64_t)op[6] * 4384447832457150726UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 4384447832457150726UL) + ((uint64_t)op[1] * 74571241062417764UL) + ((uint64_t)op[2] * 707475017012242111UL) + ((uint64_t)op[3] * 1149181324167426293UL) + ((uint64_t)op[4] * 14681101852772758639UL) + ((uint64_t)op[5] * 984514166811453756UL) + ((uint64_t)op[6] * 17068009100204866256UL);
	tmp_q[6] = ((uint64_t)op[0] * 18170997079008614544UL) + ((uint64_t)op[1] * 4384447832457150726UL) + ((uint64_t)op[2] * 74571241062417764UL) + ((uint64_t)op[3] * 707475017012242111UL) + ((uint64_t)op[4] * 1149181324167426293UL) + ((uint64_t)op[5] * 14681101852772758639UL) + ((uint64_t)op[6] * 984514166811453756UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1583169782L) + ((((int128)tmp_q[1] * 2840771027L) - ((int128)tmp_q[2] * 14084275535L) + ((int128)tmp_q[3] * 27585588845L) + ((int128)tmp_q[4] * 37108610894L) + ((int128)tmp_q[5] * 59691156181L) + ((int128)tmp_q[6] * 34817783145L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 34817783145L) - ((int128)tmp_q[1] * 1583169782L) + ((((int128)tmp_q[2] * 2840771027L) - ((int128)tmp_q[3] * 14084275535L) + ((int128)tmp_q[4] * 27585588845L) + ((int128)tmp_q[5] * 37108610894L) + ((int128)tmp_q[6] * 59691156181L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 59691156181L) + ((int128)tmp_q[1] * 34817783145L) - ((int128)tmp_q[2] * 1583169782L) + ((((int128)tmp_q[3] * 2840771027L) - ((int128)tmp_q[4] * 14084275535L) + ((int128)tmp_q[5] * 27585588845L) + ((int128)tmp_q[6] * 37108610894L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 37108610894L) + ((int128)tmp_q[1] * 59691156181L) + ((int128)tmp_q[2] * 34817783145L) - ((int128)tmp_q[3] * 1583169782L) + ((((int128)tmp_q[4] * 2840771027L) - ((int128)tmp_q[5] * 14084275535L) + ((int128)tmp_q[6] * 27585588845L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 27585588845L) + ((int128)tmp_q[1] * 37108610894L) + ((int128)tmp_q[2] * 59691156181L) + ((int128)tmp_q[3] * 34817783145L) - ((int128)tmp_q[4] * 1583169782L) + ((((int128)tmp_q[5] * 2840771027L) - ((int128)tmp_q[6] * 14084275535L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 14084275535L) + ((int128)tmp_q[1] * 27585588845L) + ((int128)tmp_q[2] * 37108610894L) + ((int128)tmp_q[3] * 59691156181L) + ((int128)tmp_q[4] * 34817783145L) - ((int128)tmp_q[5] * 1583169782L) + ((int128)tmp_q[6] * 14203855135L);
	tmp_zero[6] = ((int128)tmp_q[0] * 2840771027L) - ((int128)tmp_q[1] * 14084275535L) + ((int128)tmp_q[2] * 27585588845L) + ((int128)tmp_q[3] * 37108610894L) + ((int128)tmp_q[4] * 59691156181L) + ((int128)tmp_q[5] * 34817783145L) - ((int128)tmp_q[6] * 1583169782L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

