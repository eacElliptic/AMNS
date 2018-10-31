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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2985493728568643965UL) + ((((uint64_t)op[1] * 2650229916051806235UL) + ((uint64_t)op[2] * 4576961498302186248UL) + ((uint64_t)op[3] * 8847239199743534341UL) + ((uint64_t)op[4] * 5571098869991723768UL) + ((uint64_t)op[5] * 16272522481725653513UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 16272522481725653513UL) + ((uint64_t)op[1] * 2985493728568643965UL) + ((((uint64_t)op[2] * 2650229916051806235UL) + ((uint64_t)op[3] * 4576961498302186248UL) + ((uint64_t)op[4] * 8847239199743534341UL) + ((uint64_t)op[5] * 5571098869991723768UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 5571098869991723768UL) + ((uint64_t)op[1] * 16272522481725653513UL) + ((uint64_t)op[2] * 2985493728568643965UL) + ((((uint64_t)op[3] * 2650229916051806235UL) + ((uint64_t)op[4] * 4576961498302186248UL) + ((uint64_t)op[5] * 8847239199743534341UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 8847239199743534341UL) + ((uint64_t)op[1] * 5571098869991723768UL) + ((uint64_t)op[2] * 16272522481725653513UL) + ((uint64_t)op[3] * 2985493728568643965UL) + ((((uint64_t)op[4] * 2650229916051806235UL) + ((uint64_t)op[5] * 4576961498302186248UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 4576961498302186248UL) + ((uint64_t)op[1] * 8847239199743534341UL) + ((uint64_t)op[2] * 5571098869991723768UL) + ((uint64_t)op[3] * 16272522481725653513UL) + ((uint64_t)op[4] * 2985493728568643965UL) + ((uint64_t)op[5] * 7845824409502326676UL);
	tmp_q[5] = ((uint64_t)op[0] * 2650229916051806235UL) + ((uint64_t)op[1] * 4576961498302186248UL) + ((uint64_t)op[2] * 8847239199743534341UL) + ((uint64_t)op[3] * 5571098869991723768UL) + ((uint64_t)op[4] * 16272522481725653513UL) + ((uint64_t)op[5] * 2985493728568643965UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 13327935675L) - ((-((int128)tmp_q[1] * 48167222389L) + ((int128)tmp_q[2] * 28625061989L) - ((int128)tmp_q[3] * 31281689578L) - ((int128)tmp_q[4] * 84226787161L) + ((int128)tmp_q[5] * 40155664585L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 40155664585L) + ((int128)tmp_q[1] * 13327935675L) - ((-((int128)tmp_q[2] * 48167222389L) + ((int128)tmp_q[3] * 28625061989L) - ((int128)tmp_q[4] * 31281689578L) - ((int128)tmp_q[5] * 84226787161L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 84226787161L) + ((int128)tmp_q[1] * 40155664585L) + ((int128)tmp_q[2] * 13327935675L) - ((-((int128)tmp_q[3] * 48167222389L) + ((int128)tmp_q[4] * 28625061989L) - ((int128)tmp_q[5] * 31281689578L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 31281689578L) - ((int128)tmp_q[1] * 84226787161L) + ((int128)tmp_q[2] * 40155664585L) + ((int128)tmp_q[3] * 13327935675L) - ((-((int128)tmp_q[4] * 48167222389L) + ((int128)tmp_q[5] * 28625061989L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 28625061989L) - ((int128)tmp_q[1] * 31281689578L) - ((int128)tmp_q[2] * 84226787161L) + ((int128)tmp_q[3] * 40155664585L) + ((int128)tmp_q[4] * 13327935675L) + ((int128)tmp_q[5] * 192668889556L);
	tmp_zero[5] = -((int128)tmp_q[0] * 48167222389L) + ((int128)tmp_q[1] * 28625061989L) - ((int128)tmp_q[2] * 31281689578L) - ((int128)tmp_q[3] * 84226787161L) + ((int128)tmp_q[4] * 40155664585L) + ((int128)tmp_q[5] * 13327935675L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

