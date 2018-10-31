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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13975973438994428237UL) + ((((uint64_t)op[1] * 17760231537470501828UL) + ((uint64_t)op[2] * 18164281117135308293UL) + ((uint64_t)op[3] * 1047890448969067876UL) + ((uint64_t)op[4] * 5884837336007042774UL) + ((uint64_t)op[5] * 10310243852062552014UL) + ((uint64_t)op[6] * 17755450115243187172UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 17755450115243187172UL) + ((uint64_t)op[1] * 13975973438994428237UL) + ((((uint64_t)op[2] * 17760231537470501828UL) + ((uint64_t)op[3] * 18164281117135308293UL) + ((uint64_t)op[4] * 1047890448969067876UL) + ((uint64_t)op[5] * 5884837336007042774UL) + ((uint64_t)op[6] * 10310243852062552014UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 10310243852062552014UL) + ((uint64_t)op[1] * 17755450115243187172UL) + ((uint64_t)op[2] * 13975973438994428237UL) + ((((uint64_t)op[3] * 17760231537470501828UL) + ((uint64_t)op[4] * 18164281117135308293UL) + ((uint64_t)op[5] * 1047890448969067876UL) + ((uint64_t)op[6] * 5884837336007042774UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 5884837336007042774UL) + ((uint64_t)op[1] * 10310243852062552014UL) + ((uint64_t)op[2] * 17755450115243187172UL) + ((uint64_t)op[3] * 13975973438994428237UL) + ((((uint64_t)op[4] * 17760231537470501828UL) + ((uint64_t)op[5] * 18164281117135308293UL) + ((uint64_t)op[6] * 1047890448969067876UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 1047890448969067876UL) + ((uint64_t)op[1] * 5884837336007042774UL) + ((uint64_t)op[2] * 10310243852062552014UL) + ((uint64_t)op[3] * 17755450115243187172UL) + ((uint64_t)op[4] * 13975973438994428237UL) + ((((uint64_t)op[5] * 17760231537470501828UL) + ((uint64_t)op[6] * 18164281117135308293UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 18164281117135308293UL) + ((uint64_t)op[1] * 1047890448969067876UL) + ((uint64_t)op[2] * 5884837336007042774UL) + ((uint64_t)op[3] * 10310243852062552014UL) + ((uint64_t)op[4] * 17755450115243187172UL) + ((uint64_t)op[5] * 13975973438994428237UL) + ((uint64_t)op[6] * 12954643783797153312UL);
	tmp_q[6] = ((uint64_t)op[0] * 17760231537470501828UL) + ((uint64_t)op[1] * 18164281117135308293UL) + ((uint64_t)op[2] * 1047890448969067876UL) + ((uint64_t)op[3] * 5884837336007042774UL) + ((uint64_t)op[4] * 10310243852062552014UL) + ((uint64_t)op[5] * 17755450115243187172UL) + ((uint64_t)op[6] * 13975973438994428237UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 31792746843L) + ((-((int128)tmp_q[1] * 13283295872L) - ((int128)tmp_q[2] * 67160073803L) + ((int128)tmp_q[3] * 58844436672L) + ((int128)tmp_q[4] * 12884990734L) + ((int128)tmp_q[5] * 27262962094L) - ((int128)tmp_q[6] * 2115311644L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 2115311644L) + ((int128)tmp_q[1] * 31792746843L) + ((-((int128)tmp_q[2] * 13283295872L) - ((int128)tmp_q[3] * 67160073803L) + ((int128)tmp_q[4] * 58844436672L) + ((int128)tmp_q[5] * 12884990734L) + ((int128)tmp_q[6] * 27262962094L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 27262962094L) - ((int128)tmp_q[1] * 2115311644L) + ((int128)tmp_q[2] * 31792746843L) + ((-((int128)tmp_q[3] * 13283295872L) - ((int128)tmp_q[4] * 67160073803L) + ((int128)tmp_q[5] * 58844436672L) + ((int128)tmp_q[6] * 12884990734L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 12884990734L) + ((int128)tmp_q[1] * 27262962094L) - ((int128)tmp_q[2] * 2115311644L) + ((int128)tmp_q[3] * 31792746843L) + ((-((int128)tmp_q[4] * 13283295872L) - ((int128)tmp_q[5] * 67160073803L) + ((int128)tmp_q[6] * 58844436672L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 58844436672L) + ((int128)tmp_q[1] * 12884990734L) + ((int128)tmp_q[2] * 27262962094L) - ((int128)tmp_q[3] * 2115311644L) + ((int128)tmp_q[4] * 31792746843L) + ((-((int128)tmp_q[5] * 13283295872L) - ((int128)tmp_q[6] * 67160073803L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 67160073803L) + ((int128)tmp_q[1] * 58844436672L) + ((int128)tmp_q[2] * 12884990734L) + ((int128)tmp_q[3] * 27262962094L) - ((int128)tmp_q[4] * 2115311644L) + ((int128)tmp_q[5] * 31792746843L) - ((int128)tmp_q[6] * 106266366976L);
	tmp_zero[6] = -((int128)tmp_q[0] * 13283295872L) - ((int128)tmp_q[1] * 67160073803L) + ((int128)tmp_q[2] * 58844436672L) + ((int128)tmp_q[3] * 12884990734L) + ((int128)tmp_q[4] * 27262962094L) - ((int128)tmp_q[5] * 2115311644L) + ((int128)tmp_q[6] * 31792746843L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

