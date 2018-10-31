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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1313071582598928273UL) + ((((uint64_t)op[1] * 18358500035359575787UL) + ((uint64_t)op[2] * 15518936604089108241UL) + ((uint64_t)op[3] * 4755934093075828093UL) + ((uint64_t)op[4] * 15922085050982225161UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 15922085050982225161UL) + ((uint64_t)op[1] * 1313071582598928273UL) + ((((uint64_t)op[2] * 18358500035359575787UL) + ((uint64_t)op[3] * 15518936604089108241UL) + ((uint64_t)op[4] * 4755934093075828093UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 4755934093075828093UL) + ((uint64_t)op[1] * 15922085050982225161UL) + ((uint64_t)op[2] * 1313071582598928273UL) + ((((uint64_t)op[3] * 18358500035359575787UL) + ((uint64_t)op[4] * 15518936604089108241UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 15518936604089108241UL) + ((uint64_t)op[1] * 4755934093075828093UL) + ((uint64_t)op[2] * 15922085050982225161UL) + ((uint64_t)op[3] * 1313071582598928273UL) + ((uint64_t)op[4] * 18093767920309648300UL);
	tmp_q[4] = ((uint64_t)op[0] * 18358500035359575787UL) + ((uint64_t)op[1] * 15518936604089108241UL) + ((uint64_t)op[2] * 4755934093075828093UL) + ((uint64_t)op[3] * 15922085050982225161UL) + ((uint64_t)op[4] * 1313071582598928273UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3619258415701L) + ((-((int128)tmp_q[1] * 1510778318922L) + ((int128)tmp_q[2] * 13574053266960L) - ((int128)tmp_q[3] * 21931371610916L) - ((int128)tmp_q[4] * 12882329642331L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 12882329642331L) - ((int128)tmp_q[1] * 3619258415701L) + ((-((int128)tmp_q[2] * 1510778318922L) + ((int128)tmp_q[3] * 13574053266960L) - ((int128)tmp_q[4] * 21931371610916L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 21931371610916L) - ((int128)tmp_q[1] * 12882329642331L) - ((int128)tmp_q[2] * 3619258415701L) + ((-((int128)tmp_q[3] * 1510778318922L) + ((int128)tmp_q[4] * 13574053266960L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 13574053266960L) - ((int128)tmp_q[1] * 21931371610916L) - ((int128)tmp_q[2] * 12882329642331L) - ((int128)tmp_q[3] * 3619258415701L) - ((int128)tmp_q[4] * 6043113275688L);
	tmp_zero[4] = -((int128)tmp_q[0] * 1510778318922L) + ((int128)tmp_q[1] * 13574053266960L) - ((int128)tmp_q[2] * 21931371610916L) - ((int128)tmp_q[3] * 12882329642331L) - ((int128)tmp_q[4] * 3619258415701L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

