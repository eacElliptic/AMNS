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
	tmp_q[0] = ((uint64_t)op[0] * 6378122719252998468UL) + ((((uint64_t)op[1] * 13062113999861344062UL) + ((uint64_t)op[2] * 15363476072483424375UL) + ((uint64_t)op[3] * 1935930233447910015UL) + ((uint64_t)op[4] * 13675834284163184066UL) + ((uint64_t)op[5] * 5345057877985595377UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 5345057877985595377UL) + ((uint64_t)op[1] * 6378122719252998468UL) + ((((uint64_t)op[2] * 13062113999861344062UL) + ((uint64_t)op[3] * 15363476072483424375UL) + ((uint64_t)op[4] * 1935930233447910015UL) + ((uint64_t)op[5] * 13675834284163184066UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 13675834284163184066UL) + ((uint64_t)op[1] * 5345057877985595377UL) + ((uint64_t)op[2] * 6378122719252998468UL) + ((((uint64_t)op[3] * 13062113999861344062UL) + ((uint64_t)op[4] * 15363476072483424375UL) + ((uint64_t)op[5] * 1935930233447910015UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 1935930233447910015UL) + ((uint64_t)op[1] * 13675834284163184066UL) + ((uint64_t)op[2] * 5345057877985595377UL) + ((uint64_t)op[3] * 6378122719252998468UL) + ((((uint64_t)op[4] * 13062113999861344062UL) + ((uint64_t)op[5] * 15363476072483424375UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 15363476072483424375UL) + ((uint64_t)op[1] * 1935930233447910015UL) + ((uint64_t)op[2] * 13675834284163184066UL) + ((uint64_t)op[3] * 5345057877985595377UL) + ((uint64_t)op[4] * 6378122719252998468UL) + ((uint64_t)op[5] * 9970337778178065462UL);
	tmp_q[5] = ((uint64_t)op[0] * 13062113999861344062UL) + ((uint64_t)op[1] * 15363476072483424375UL) + ((uint64_t)op[2] * 1935930233447910015UL) + ((uint64_t)op[3] * 13675834284163184066UL) + ((uint64_t)op[4] * 5345057877985595377UL) + ((uint64_t)op[5] * 6378122719252998468UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1595028748L) + ((((int128)tmp_q[1] * 1998306132L) + ((int128)tmp_q[2] * 1018163203L) - ((int128)tmp_q[3] * 784614155L) + ((int128)tmp_q[4] * 1342124984L) - ((int128)tmp_q[5] * 579288443L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 579288443L) - ((int128)tmp_q[1] * 1595028748L) + ((((int128)tmp_q[2] * 1998306132L) + ((int128)tmp_q[3] * 1018163203L) - ((int128)tmp_q[4] * 784614155L) + ((int128)tmp_q[5] * 1342124984L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 1342124984L) - ((int128)tmp_q[1] * 579288443L) - ((int128)tmp_q[2] * 1595028748L) + ((((int128)tmp_q[3] * 1998306132L) + ((int128)tmp_q[4] * 1018163203L) - ((int128)tmp_q[5] * 784614155L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 784614155L) + ((int128)tmp_q[1] * 1342124984L) - ((int128)tmp_q[2] * 579288443L) - ((int128)tmp_q[3] * 1595028748L) + ((((int128)tmp_q[4] * 1998306132L) + ((int128)tmp_q[5] * 1018163203L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 1018163203L) - ((int128)tmp_q[1] * 784614155L) + ((int128)tmp_q[2] * 1342124984L) - ((int128)tmp_q[3] * 579288443L) - ((int128)tmp_q[4] * 1595028748L) + ((int128)tmp_q[5] * 9991530660L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1998306132L) + ((int128)tmp_q[1] * 1018163203L) - ((int128)tmp_q[2] * 784614155L) + ((int128)tmp_q[3] * 1342124984L) - ((int128)tmp_q[4] * 579288443L) - ((int128)tmp_q[5] * 1595028748L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

