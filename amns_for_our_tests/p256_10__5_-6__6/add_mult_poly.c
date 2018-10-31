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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12030442742899727789UL) + ((((uint64_t)op[1] * 1594089333880806781UL) + ((uint64_t)op[2] * 13874035965640586246UL) + ((uint64_t)op[3] * 6552569562813210010UL) + ((uint64_t)op[4] * 16188312538540940950UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 16188312538540940950UL) + ((uint64_t)op[1] * 12030442742899727789UL) + ((((uint64_t)op[2] * 1594089333880806781UL) + ((uint64_t)op[3] * 13874035965640586246UL) + ((uint64_t)op[4] * 6552569562813210010UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 6552569562813210010UL) + ((uint64_t)op[1] * 16188312538540940950UL) + ((uint64_t)op[2] * 12030442742899727789UL) + ((((uint64_t)op[3] * 1594089333880806781UL) + ((uint64_t)op[4] * 13874035965640586246UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 13874035965640586246UL) + ((uint64_t)op[1] * 6552569562813210010UL) + ((uint64_t)op[2] * 16188312538540940950UL) + ((uint64_t)op[3] * 12030442742899727789UL) + ((uint64_t)op[4] * 8882208070424710930UL);
	tmp_q[4] = ((uint64_t)op[0] * 1594089333880806781UL) + ((uint64_t)op[1] * 13874035965640586246UL) + ((uint64_t)op[2] * 6552569562813210010UL) + ((uint64_t)op[3] * 16188312538540940950UL) + ((uint64_t)op[4] * 12030442742899727789UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 558641914211093L) - ((((int128)tmp_q[1] * 1402914204515701L) + ((int128)tmp_q[2] * 1794965304409884L) + ((int128)tmp_q[3] * 547775021085810L) - ((int128)tmp_q[4] * 1015036615127714L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 1015036615127714L) - ((int128)tmp_q[1] * 558641914211093L) - ((((int128)tmp_q[2] * 1402914204515701L) + ((int128)tmp_q[3] * 1794965304409884L) + ((int128)tmp_q[4] * 547775021085810L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 547775021085810L) - ((int128)tmp_q[1] * 1015036615127714L) - ((int128)tmp_q[2] * 558641914211093L) - ((((int128)tmp_q[3] * 1402914204515701L) + ((int128)tmp_q[4] * 1794965304409884L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 1794965304409884L) + ((int128)tmp_q[1] * 547775021085810L) - ((int128)tmp_q[2] * 1015036615127714L) - ((int128)tmp_q[3] * 558641914211093L) - ((int128)tmp_q[4] * 8417485227094206L);
	tmp_zero[4] = ((int128)tmp_q[0] * 1402914204515701L) + ((int128)tmp_q[1] * 1794965304409884L) + ((int128)tmp_q[2] * 547775021085810L) - ((int128)tmp_q[3] * 1015036615127714L) - ((int128)tmp_q[4] * 558641914211093L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

