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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1115206617461371825UL) + ((((uint64_t)op[1] * 4254450486739099141UL) + ((uint64_t)op[2] * 15817923304991519584UL) + ((uint64_t)op[3] * 3698370697641263762UL) + ((uint64_t)op[4] * 13906107359886561899UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 13906107359886561899UL) + ((uint64_t)op[1] * 1115206617461371825UL) + ((((uint64_t)op[2] * 4254450486739099141UL) + ((uint64_t)op[3] * 15817923304991519584UL) + ((uint64_t)op[4] * 3698370697641263762UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 3698370697641263762UL) + ((uint64_t)op[1] * 13906107359886561899UL) + ((uint64_t)op[2] * 1115206617461371825UL) + ((((uint64_t)op[3] * 4254450486739099141UL) + ((uint64_t)op[4] * 15817923304991519584UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 15817923304991519584UL) + ((uint64_t)op[1] * 3698370697641263762UL) + ((uint64_t)op[2] * 13906107359886561899UL) + ((uint64_t)op[3] * 1115206617461371825UL) + ((uint64_t)op[4] * 2825508359985944089UL);
	tmp_q[4] = ((uint64_t)op[0] * 4254450486739099141UL) + ((uint64_t)op[1] * 15817923304991519584UL) + ((uint64_t)op[2] * 3698370697641263762UL) + ((uint64_t)op[3] * 13906107359886561899UL) + ((uint64_t)op[4] * 1115206617461371825UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 59985229713L) + ((-((int128)tmp_q[1] * 91730827718L) - ((int128)tmp_q[2] * 43146785215L) + ((int128)tmp_q[3] * 15808178235L) - ((int128)tmp_q[4] * 131500914242L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 131500914242L) + ((int128)tmp_q[1] * 59985229713L) + ((-((int128)tmp_q[2] * 91730827718L) - ((int128)tmp_q[3] * 43146785215L) + ((int128)tmp_q[4] * 15808178235L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 15808178235L) - ((int128)tmp_q[1] * 131500914242L) + ((int128)tmp_q[2] * 59985229713L) + ((-((int128)tmp_q[3] * 91730827718L) - ((int128)tmp_q[4] * 43146785215L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 43146785215L) + ((int128)tmp_q[1] * 15808178235L) - ((int128)tmp_q[2] * 131500914242L) + ((int128)tmp_q[3] * 59985229713L) - ((int128)tmp_q[4] * 458654138590L);
	tmp_zero[4] = -((int128)tmp_q[0] * 91730827718L) - ((int128)tmp_q[1] * 43146785215L) + ((int128)tmp_q[2] * 15808178235L) - ((int128)tmp_q[3] * 131500914242L) + ((int128)tmp_q[4] * 59985229713L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

