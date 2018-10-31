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
	tmp_q[0] = ((uint64_t)op[0] * 18054503162246253377UL) + ((((uint64_t)op[1] * 6445158496448159872UL) + ((uint64_t)op[2] * 10891103429788260653UL) + ((uint64_t)op[3] * 9390622438804773131UL) + ((uint64_t)op[4] * 15343301489650747173UL) + ((uint64_t)op[5] * 6671969937432533122UL) + ((uint64_t)op[6] * 5238730448298173726UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 5238730448298173726UL) + ((uint64_t)op[1] * 18054503162246253377UL) + ((((uint64_t)op[2] * 6445158496448159872UL) + ((uint64_t)op[3] * 10891103429788260653UL) + ((uint64_t)op[4] * 9390622438804773131UL) + ((uint64_t)op[5] * 15343301489650747173UL) + ((uint64_t)op[6] * 6671969937432533122UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 6671969937432533122UL) + ((uint64_t)op[1] * 5238730448298173726UL) + ((uint64_t)op[2] * 18054503162246253377UL) + ((((uint64_t)op[3] * 6445158496448159872UL) + ((uint64_t)op[4] * 10891103429788260653UL) + ((uint64_t)op[5] * 9390622438804773131UL) + ((uint64_t)op[6] * 15343301489650747173UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 15343301489650747173UL) + ((uint64_t)op[1] * 6671969937432533122UL) + ((uint64_t)op[2] * 5238730448298173726UL) + ((uint64_t)op[3] * 18054503162246253377UL) + ((((uint64_t)op[4] * 6445158496448159872UL) + ((uint64_t)op[5] * 10891103429788260653UL) + ((uint64_t)op[6] * 9390622438804773131UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 9390622438804773131UL) + ((uint64_t)op[1] * 15343301489650747173UL) + ((uint64_t)op[2] * 6671969937432533122UL) + ((uint64_t)op[3] * 5238730448298173726UL) + ((uint64_t)op[4] * 18054503162246253377UL) + ((((uint64_t)op[5] * 6445158496448159872UL) + ((uint64_t)op[6] * 10891103429788260653UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 10891103429788260653UL) + ((uint64_t)op[1] * 9390622438804773131UL) + ((uint64_t)op[2] * 15343301489650747173UL) + ((uint64_t)op[3] * 6671969937432533122UL) + ((uint64_t)op[4] * 5238730448298173726UL) + ((uint64_t)op[5] * 18054503162246253377UL) + ((uint64_t)op[6] * 14667779824166175744UL);
	tmp_q[6] = ((uint64_t)op[0] * 6445158496448159872UL) + ((uint64_t)op[1] * 10891103429788260653UL) + ((uint64_t)op[2] * 9390622438804773131UL) + ((uint64_t)op[3] * 15343301489650747173UL) + ((uint64_t)op[4] * 6671969937432533122UL) + ((uint64_t)op[5] * 5238730448298173726UL) + ((uint64_t)op[6] * 18054503162246253377UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 6886841023L) + ((-((int128)tmp_q[1] * 7093349213L) + ((int128)tmp_q[2] * 9720874073L) + ((int128)tmp_q[3] * 4926907011L) - ((int128)tmp_q[4] * 50622270443L) - ((int128)tmp_q[5] * 10188283290L) - ((int128)tmp_q[6] * 24573965450L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 24573965450L) + ((int128)tmp_q[1] * 6886841023L) + ((-((int128)tmp_q[2] * 7093349213L) + ((int128)tmp_q[3] * 9720874073L) + ((int128)tmp_q[4] * 4926907011L) - ((int128)tmp_q[5] * 50622270443L) - ((int128)tmp_q[6] * 10188283290L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 10188283290L) - ((int128)tmp_q[1] * 24573965450L) + ((int128)tmp_q[2] * 6886841023L) + ((-((int128)tmp_q[3] * 7093349213L) + ((int128)tmp_q[4] * 9720874073L) + ((int128)tmp_q[5] * 4926907011L) - ((int128)tmp_q[6] * 50622270443L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 50622270443L) - ((int128)tmp_q[1] * 10188283290L) - ((int128)tmp_q[2] * 24573965450L) + ((int128)tmp_q[3] * 6886841023L) + ((-((int128)tmp_q[4] * 7093349213L) + ((int128)tmp_q[5] * 9720874073L) + ((int128)tmp_q[6] * 4926907011L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 4926907011L) - ((int128)tmp_q[1] * 50622270443L) - ((int128)tmp_q[2] * 10188283290L) - ((int128)tmp_q[3] * 24573965450L) + ((int128)tmp_q[4] * 6886841023L) + ((-((int128)tmp_q[5] * 7093349213L) + ((int128)tmp_q[6] * 9720874073L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 9720874073L) + ((int128)tmp_q[1] * 4926907011L) - ((int128)tmp_q[2] * 50622270443L) - ((int128)tmp_q[3] * 10188283290L) - ((int128)tmp_q[4] * 24573965450L) + ((int128)tmp_q[5] * 6886841023L) - ((int128)tmp_q[6] * 56746793704L);
	tmp_zero[6] = -((int128)tmp_q[0] * 7093349213L) + ((int128)tmp_q[1] * 9720874073L) + ((int128)tmp_q[2] * 4926907011L) - ((int128)tmp_q[3] * 50622270443L) - ((int128)tmp_q[4] * 10188283290L) - ((int128)tmp_q[5] * 24573965450L) + ((int128)tmp_q[6] * 6886841023L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

