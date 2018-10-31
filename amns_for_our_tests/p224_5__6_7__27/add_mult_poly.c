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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18256110393538185054UL) + ((((uint64_t)op[1] * 904444724886311174UL) + ((uint64_t)op[2] * 232097748770364067UL) + ((uint64_t)op[3] * 17855972179381968168UL) + ((uint64_t)op[4] * 14182573698899986417UL) + ((uint64_t)op[5] * 14228080845406982669UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 14228080845406982669UL) + ((uint64_t)op[1] * 18256110393538185054UL) + ((((uint64_t)op[2] * 904444724886311174UL) + ((uint64_t)op[3] * 232097748770364067UL) + ((uint64_t)op[4] * 17855972179381968168UL) + ((uint64_t)op[5] * 14182573698899986417UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 14182573698899986417UL) + ((uint64_t)op[1] * 14228080845406982669UL) + ((uint64_t)op[2] * 18256110393538185054UL) + ((((uint64_t)op[3] * 904444724886311174UL) + ((uint64_t)op[4] * 232097748770364067UL) + ((uint64_t)op[5] * 17855972179381968168UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 17855972179381968168UL) + ((uint64_t)op[1] * 14182573698899986417UL) + ((uint64_t)op[2] * 14228080845406982669UL) + ((uint64_t)op[3] * 18256110393538185054UL) + ((((uint64_t)op[4] * 904444724886311174UL) + ((uint64_t)op[5] * 232097748770364067UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 232097748770364067UL) + ((uint64_t)op[1] * 17855972179381968168UL) + ((uint64_t)op[2] * 14182573698899986417UL) + ((uint64_t)op[3] * 14228080845406982669UL) + ((uint64_t)op[4] * 18256110393538185054UL) + ((uint64_t)op[5] * 6331113074204178218UL);
	tmp_q[5] = ((uint64_t)op[0] * 904444724886311174UL) + ((uint64_t)op[1] * 232097748770364067UL) + ((uint64_t)op[2] * 17855972179381968168UL) + ((uint64_t)op[3] * 14182573698899986417UL) + ((uint64_t)op[4] * 14228080845406982669UL) + ((uint64_t)op[5] * 18256110393538185054UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 68926294181L) + ((((int128)tmp_q[1] * 53147867398L) - ((int128)tmp_q[2] * 84488905997L) + ((int128)tmp_q[3] * 11004449919L) - ((int128)tmp_q[4] * 9510698650L) - ((int128)tmp_q[5] * 41778421220L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 41778421220L) - ((int128)tmp_q[1] * 68926294181L) + ((((int128)tmp_q[2] * 53147867398L) - ((int128)tmp_q[3] * 84488905997L) + ((int128)tmp_q[4] * 11004449919L) - ((int128)tmp_q[5] * 9510698650L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 9510698650L) - ((int128)tmp_q[1] * 41778421220L) - ((int128)tmp_q[2] * 68926294181L) + ((((int128)tmp_q[3] * 53147867398L) - ((int128)tmp_q[4] * 84488905997L) + ((int128)tmp_q[5] * 11004449919L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 11004449919L) - ((int128)tmp_q[1] * 9510698650L) - ((int128)tmp_q[2] * 41778421220L) - ((int128)tmp_q[3] * 68926294181L) + ((((int128)tmp_q[4] * 53147867398L) - ((int128)tmp_q[5] * 84488905997L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 84488905997L) + ((int128)tmp_q[1] * 11004449919L) - ((int128)tmp_q[2] * 9510698650L) - ((int128)tmp_q[3] * 41778421220L) - ((int128)tmp_q[4] * 68926294181L) + ((int128)tmp_q[5] * 372035071786L);
	tmp_zero[5] = ((int128)tmp_q[0] * 53147867398L) - ((int128)tmp_q[1] * 84488905997L) + ((int128)tmp_q[2] * 11004449919L) - ((int128)tmp_q[3] * 9510698650L) - ((int128)tmp_q[4] * 41778421220L) - ((int128)tmp_q[5] * 68926294181L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

