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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8959746220229806233UL) + ((((uint64_t)op[1] * 14311312274313856304UL) + ((uint64_t)op[2] * 17251584217284999877UL) + ((uint64_t)op[3] * 18386379201540816375UL) + ((uint64_t)op[4] * 4171926639446308307UL) + ((uint64_t)op[5] * 9041631466732722768UL) + ((uint64_t)op[6] * 9027328966508184044UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 9027328966508184044UL) + ((uint64_t)op[1] * 8959746220229806233UL) + ((((uint64_t)op[2] * 14311312274313856304UL) + ((uint64_t)op[3] * 17251584217284999877UL) + ((uint64_t)op[4] * 18386379201540816375UL) + ((uint64_t)op[5] * 4171926639446308307UL) + ((uint64_t)op[6] * 9041631466732722768UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 9041631466732722768UL) + ((uint64_t)op[1] * 9027328966508184044UL) + ((uint64_t)op[2] * 8959746220229806233UL) + ((((uint64_t)op[3] * 14311312274313856304UL) + ((uint64_t)op[4] * 17251584217284999877UL) + ((uint64_t)op[5] * 18386379201540816375UL) + ((uint64_t)op[6] * 4171926639446308307UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 4171926639446308307UL) + ((uint64_t)op[1] * 9041631466732722768UL) + ((uint64_t)op[2] * 9027328966508184044UL) + ((uint64_t)op[3] * 8959746220229806233UL) + ((((uint64_t)op[4] * 14311312274313856304UL) + ((uint64_t)op[5] * 17251584217284999877UL) + ((uint64_t)op[6] * 18386379201540816375UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 18386379201540816375UL) + ((uint64_t)op[1] * 4171926639446308307UL) + ((uint64_t)op[2] * 9041631466732722768UL) + ((uint64_t)op[3] * 9027328966508184044UL) + ((uint64_t)op[4] * 8959746220229806233UL) + ((((uint64_t)op[5] * 14311312274313856304UL) + ((uint64_t)op[6] * 17251584217284999877UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 17251584217284999877UL) + ((uint64_t)op[1] * 18386379201540816375UL) + ((uint64_t)op[2] * 4171926639446308307UL) + ((uint64_t)op[3] * 9041631466732722768UL) + ((uint64_t)op[4] * 9027328966508184044UL) + ((uint64_t)op[5] * 8959746220229806233UL) + ((uint64_t)op[6] * 6365846722664620256UL);
	tmp_q[6] = ((uint64_t)op[0] * 14311312274313856304UL) + ((uint64_t)op[1] * 17251584217284999877UL) + ((uint64_t)op[2] * 18386379201540816375UL) + ((uint64_t)op[3] * 4171926639446308307UL) + ((uint64_t)op[4] * 9041631466732722768UL) + ((uint64_t)op[5] * 9027328966508184044UL) + ((uint64_t)op[6] * 8959746220229806233UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 6849970297L) - ((((int128)tmp_q[1] * 36079418355L) + ((int128)tmp_q[2] * 3419489097L) - ((int128)tmp_q[3] * 47095985581L) - ((int128)tmp_q[4] * 40279304185L) + ((int128)tmp_q[5] * 48271186850L) - ((int128)tmp_q[6] * 29734172434L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 29734172434L) - ((int128)tmp_q[1] * 6849970297L) - ((((int128)tmp_q[2] * 36079418355L) + ((int128)tmp_q[3] * 3419489097L) - ((int128)tmp_q[4] * 47095985581L) - ((int128)tmp_q[5] * 40279304185L) + ((int128)tmp_q[6] * 48271186850L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 48271186850L) - ((int128)tmp_q[1] * 29734172434L) - ((int128)tmp_q[2] * 6849970297L) - ((((int128)tmp_q[3] * 36079418355L) + ((int128)tmp_q[4] * 3419489097L) - ((int128)tmp_q[5] * 47095985581L) - ((int128)tmp_q[6] * 40279304185L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 40279304185L) + ((int128)tmp_q[1] * 48271186850L) - ((int128)tmp_q[2] * 29734172434L) - ((int128)tmp_q[3] * 6849970297L) - ((((int128)tmp_q[4] * 36079418355L) + ((int128)tmp_q[5] * 3419489097L) - ((int128)tmp_q[6] * 47095985581L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 47095985581L) - ((int128)tmp_q[1] * 40279304185L) + ((int128)tmp_q[2] * 48271186850L) - ((int128)tmp_q[3] * 29734172434L) - ((int128)tmp_q[4] * 6849970297L) - ((((int128)tmp_q[5] * 36079418355L) + ((int128)tmp_q[6] * 3419489097L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 3419489097L) - ((int128)tmp_q[1] * 47095985581L) - ((int128)tmp_q[2] * 40279304185L) + ((int128)tmp_q[3] * 48271186850L) - ((int128)tmp_q[4] * 29734172434L) - ((int128)tmp_q[5] * 6849970297L) - ((int128)tmp_q[6] * 216476510130L);
	tmp_zero[6] = ((int128)tmp_q[0] * 36079418355L) + ((int128)tmp_q[1] * 3419489097L) - ((int128)tmp_q[2] * 47095985581L) - ((int128)tmp_q[3] * 40279304185L) + ((int128)tmp_q[4] * 48271186850L) - ((int128)tmp_q[5] * 29734172434L) - ((int128)tmp_q[6] * 6849970297L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

