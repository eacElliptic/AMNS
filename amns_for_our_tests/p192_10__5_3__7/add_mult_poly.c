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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16189722089169904876UL) + ((((uint64_t)op[1] * 12619629207834115536UL) + ((uint64_t)op[2] * 7667584631139487150UL) + ((uint64_t)op[3] * 394077755286943512UL) + ((uint64_t)op[4] * 12238765475845336765UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 12238765475845336765UL) + ((uint64_t)op[1] * 16189722089169904876UL) + ((((uint64_t)op[2] * 12619629207834115536UL) + ((uint64_t)op[3] * 7667584631139487150UL) + ((uint64_t)op[4] * 394077755286943512UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 394077755286943512UL) + ((uint64_t)op[1] * 12238765475845336765UL) + ((uint64_t)op[2] * 16189722089169904876UL) + ((((uint64_t)op[3] * 12619629207834115536UL) + ((uint64_t)op[4] * 7667584631139487150UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 7667584631139487150UL) + ((uint64_t)op[1] * 394077755286943512UL) + ((uint64_t)op[2] * 12238765475845336765UL) + ((uint64_t)op[3] * 16189722089169904876UL) + ((uint64_t)op[4] * 965399476083243376UL);
	tmp_q[4] = ((uint64_t)op[0] * 12619629207834115536UL) + ((uint64_t)op[1] * 7667584631139487150UL) + ((uint64_t)op[2] * 394077755286943512UL) + ((uint64_t)op[3] * 12238765475845336765UL) + ((uint64_t)op[4] * 16189722089169904876UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 130070011376L) + ((((int128)tmp_q[1] * 87348731673L) - ((int128)tmp_q[2] * 87573080688L) + ((int128)tmp_q[3] * 112541027296L) + ((int128)tmp_q[4] * 116466877310L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 116466877310L) + ((int128)tmp_q[1] * 130070011376L) + ((((int128)tmp_q[2] * 87348731673L) - ((int128)tmp_q[3] * 87573080688L) + ((int128)tmp_q[4] * 112541027296L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 112541027296L) + ((int128)tmp_q[1] * 116466877310L) + ((int128)tmp_q[2] * 130070011376L) + ((((int128)tmp_q[3] * 87348731673L) - ((int128)tmp_q[4] * 87573080688L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 87573080688L) + ((int128)tmp_q[1] * 112541027296L) + ((int128)tmp_q[2] * 116466877310L) + ((int128)tmp_q[3] * 130070011376L) + ((int128)tmp_q[4] * 262046195019L);
	tmp_zero[4] = ((int128)tmp_q[0] * 87348731673L) - ((int128)tmp_q[1] * 87573080688L) + ((int128)tmp_q[2] * 112541027296L) + ((int128)tmp_q[3] * 116466877310L) + ((int128)tmp_q[4] * 130070011376L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

