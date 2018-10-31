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
	tmp_q[0] = ((uint64_t)op[0] * 81582114544656335UL) + ((((uint64_t)op[1] * 12509319367620255923UL) + ((uint64_t)op[2] * 15130455437319501117UL) + ((uint64_t)op[3] * 7596186924194773562UL) + ((uint64_t)op[4] * 15820959919153275706UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 15820959919153275706UL) + ((uint64_t)op[1] * 81582114544656335UL) + ((((uint64_t)op[2] * 12509319367620255923UL) + ((uint64_t)op[3] * 15130455437319501117UL) + ((uint64_t)op[4] * 7596186924194773562UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 7596186924194773562UL) + ((uint64_t)op[1] * 15820959919153275706UL) + ((uint64_t)op[2] * 81582114544656335UL) + ((((uint64_t)op[3] * 12509319367620255923UL) + ((uint64_t)op[4] * 15130455437319501117UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 15130455437319501117UL) + ((uint64_t)op[1] * 7596186924194773562UL) + ((uint64_t)op[2] * 15820959919153275706UL) + ((uint64_t)op[3] * 81582114544656335UL) + ((uint64_t)op[4] * 6571894661530960230UL);
	tmp_q[4] = ((uint64_t)op[0] * 12509319367620255923UL) + ((uint64_t)op[1] * 15130455437319501117UL) + ((uint64_t)op[2] * 7596186924194773562UL) + ((uint64_t)op[3] * 15820959919153275706UL) + ((uint64_t)op[4] * 81582114544656335UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 662029473251229L) + ((-((int128)tmp_q[1] * 96065912122983L) + ((int128)tmp_q[2] * 1431297676867947L) + ((int128)tmp_q[3] * 804327015687470L) - ((int128)tmp_q[4] * 571292729887560L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 571292729887560L) + ((int128)tmp_q[1] * 662029473251229L) + ((-((int128)tmp_q[2] * 96065912122983L) + ((int128)tmp_q[3] * 1431297676867947L) + ((int128)tmp_q[4] * 804327015687470L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 804327015687470L) - ((int128)tmp_q[1] * 571292729887560L) + ((int128)tmp_q[2] * 662029473251229L) + ((-((int128)tmp_q[3] * 96065912122983L) + ((int128)tmp_q[4] * 1431297676867947L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 1431297676867947L) + ((int128)tmp_q[1] * 804327015687470L) - ((int128)tmp_q[2] * 571292729887560L) + ((int128)tmp_q[3] * 662029473251229L) - ((int128)tmp_q[4] * 192131824245966L);
	tmp_zero[4] = -((int128)tmp_q[0] * 96065912122983L) + ((int128)tmp_q[1] * 1431297676867947L) + ((int128)tmp_q[2] * 804327015687470L) - ((int128)tmp_q[3] * 571292729887560L) + ((int128)tmp_q[4] * 662029473251229L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

