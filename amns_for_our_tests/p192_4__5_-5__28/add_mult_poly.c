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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8936213679215798929UL) + ((((uint64_t)op[1] * 4237868564112465652UL) + ((uint64_t)op[2] * 111157696738214677UL) + ((uint64_t)op[3] * 18274258011644698120UL) + ((uint64_t)op[4] * 11359027464177122745UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 11359027464177122745UL) + ((uint64_t)op[1] * 8936213679215798929UL) + ((((uint64_t)op[2] * 4237868564112465652UL) + ((uint64_t)op[3] * 111157696738214677UL) + ((uint64_t)op[4] * 18274258011644698120UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 18274258011644698120UL) + ((uint64_t)op[1] * 11359027464177122745UL) + ((uint64_t)op[2] * 8936213679215798929UL) + ((((uint64_t)op[3] * 4237868564112465652UL) + ((uint64_t)op[4] * 111157696738214677UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 111157696738214677UL) + ((uint64_t)op[1] * 18274258011644698120UL) + ((uint64_t)op[2] * 11359027464177122745UL) + ((uint64_t)op[3] * 8936213679215798929UL) + ((uint64_t)op[4] * 15704145326856774972UL);
	tmp_q[4] = ((uint64_t)op[0] * 4237868564112465652UL) + ((uint64_t)op[1] * 111157696738214677UL) + ((uint64_t)op[2] * 18274258011644698120UL) + ((uint64_t)op[3] * 11359027464177122745UL) + ((uint64_t)op[4] * 8936213679215798929UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 166458493942L) - ((((int128)tmp_q[1] * 26098668900L) + ((int128)tmp_q[2] * 259169195L) + ((int128)tmp_q[3] * 105199138579L) - ((int128)tmp_q[4] * 220409532347L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 220409532347L) - ((int128)tmp_q[1] * 166458493942L) - ((((int128)tmp_q[2] * 26098668900L) + ((int128)tmp_q[3] * 259169195L) + ((int128)tmp_q[4] * 105199138579L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 105199138579L) - ((int128)tmp_q[1] * 220409532347L) - ((int128)tmp_q[2] * 166458493942L) - ((((int128)tmp_q[3] * 26098668900L) + ((int128)tmp_q[4] * 259169195L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 259169195L) + ((int128)tmp_q[1] * 105199138579L) - ((int128)tmp_q[2] * 220409532347L) - ((int128)tmp_q[3] * 166458493942L) - ((int128)tmp_q[4] * 130493344500L);
	tmp_zero[4] = ((int128)tmp_q[0] * 26098668900L) + ((int128)tmp_q[1] * 259169195L) + ((int128)tmp_q[2] * 105199138579L) - ((int128)tmp_q[3] * 220409532347L) - ((int128)tmp_q[4] * 166458493942L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

