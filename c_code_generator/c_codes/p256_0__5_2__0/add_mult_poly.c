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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14356718804419984039UL) + ((((uint64_t)op[1] * 10537774877680181295UL) + ((uint64_t)op[2] * 5592215177109344109UL) + ((uint64_t)op[3] * 6508683138508426415UL) + ((uint64_t)op[4] * 1794657206921851374UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 1794657206921851374UL) + ((uint64_t)op[1] * 14356718804419984039UL) + ((((uint64_t)op[2] * 10537774877680181295UL) + ((uint64_t)op[3] * 5592215177109344109UL) + ((uint64_t)op[4] * 6508683138508426415UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 6508683138508426415UL) + ((uint64_t)op[1] * 1794657206921851374UL) + ((uint64_t)op[2] * 14356718804419984039UL) + ((((uint64_t)op[3] * 10537774877680181295UL) + ((uint64_t)op[4] * 5592215177109344109UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 5592215177109344109UL) + ((uint64_t)op[1] * 6508683138508426415UL) + ((uint64_t)op[2] * 1794657206921851374UL) + ((uint64_t)op[3] * 14356718804419984039UL) + ((uint64_t)op[4] * 2628805681650810974UL);
	tmp_q[4] = ((uint64_t)op[0] * 10537774877680181295UL) + ((uint64_t)op[1] * 5592215177109344109UL) + ((uint64_t)op[2] * 6508683138508426415UL) + ((uint64_t)op[3] * 1794657206921851374UL) + ((uint64_t)op[4] * 14356718804419984039UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1626219128380213L) + ((((int128)tmp_q[1] * 292601798855450L) - ((int128)tmp_q[2] * 1345411614681303L) - ((int128)tmp_q[3] * 1663301451987807L) - ((int128)tmp_q[4] * 982952150664854L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 982952150664854L) + ((int128)tmp_q[1] * 1626219128380213L) + ((((int128)tmp_q[2] * 292601798855450L) - ((int128)tmp_q[3] * 1345411614681303L) - ((int128)tmp_q[4] * 1663301451987807L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 1663301451987807L) - ((int128)tmp_q[1] * 982952150664854L) + ((int128)tmp_q[2] * 1626219128380213L) + ((((int128)tmp_q[3] * 292601798855450L) - ((int128)tmp_q[4] * 1345411614681303L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 1345411614681303L) - ((int128)tmp_q[1] * 1663301451987807L) - ((int128)tmp_q[2] * 982952150664854L) + ((int128)tmp_q[3] * 1626219128380213L) + ((int128)tmp_q[4] * 585203597710900L);
	tmp_zero[4] = ((int128)tmp_q[0] * 292601798855450L) - ((int128)tmp_q[1] * 1345411614681303L) - ((int128)tmp_q[2] * 1663301451987807L) - ((int128)tmp_q[3] * 982952150664854L) + ((int128)tmp_q[4] * 1626219128380213L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

