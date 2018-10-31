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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18313103972657734961UL) + ((((uint64_t)op[1] * 12224908551704187103UL) + ((uint64_t)op[2] * 9039832678153265410UL) + ((uint64_t)op[3] * 11009378302419516115UL) + ((uint64_t)op[4] * 17596241560430020058UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 17596241560430020058UL) + ((uint64_t)op[1] * 18313103972657734961UL) + ((((uint64_t)op[2] * 12224908551704187103UL) + ((uint64_t)op[3] * 9039832678153265410UL) + ((uint64_t)op[4] * 11009378302419516115UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 11009378302419516115UL) + ((uint64_t)op[1] * 17596241560430020058UL) + ((uint64_t)op[2] * 18313103972657734961UL) + ((((uint64_t)op[3] * 12224908551704187103UL) + ((uint64_t)op[4] * 9039832678153265410UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 9039832678153265410UL) + ((uint64_t)op[1] * 11009378302419516115UL) + ((uint64_t)op[2] * 17596241560430020058UL) + ((uint64_t)op[3] * 18313103972657734961UL) + ((uint64_t)op[4] * 5565548045085738744UL);
	tmp_q[4] = ((uint64_t)op[0] * 12224908551704187103UL) + ((uint64_t)op[1] * 9039832678153265410UL) + ((uint64_t)op[2] * 11009378302419516115UL) + ((uint64_t)op[3] * 17596241560430020058UL) + ((uint64_t)op[4] * 18313103972657734961UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 751005797345L) + ((-((int128)tmp_q[1] * 2915474161454L) + ((int128)tmp_q[2] * 12436052620422L) - ((int128)tmp_q[3] * 1096807264321L) - ((int128)tmp_q[4] * 328932706782L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 328932706782L) - ((int128)tmp_q[1] * 751005797345L) + ((-((int128)tmp_q[2] * 2915474161454L) + ((int128)tmp_q[3] * 12436052620422L) - ((int128)tmp_q[4] * 1096807264321L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 1096807264321L) - ((int128)tmp_q[1] * 328932706782L) - ((int128)tmp_q[2] * 751005797345L) + ((-((int128)tmp_q[3] * 2915474161454L) + ((int128)tmp_q[4] * 12436052620422L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 12436052620422L) - ((int128)tmp_q[1] * 1096807264321L) - ((int128)tmp_q[2] * 328932706782L) - ((int128)tmp_q[3] * 751005797345L) - ((int128)tmp_q[4] * 23323793291632L);
	tmp_zero[4] = -((int128)tmp_q[0] * 2915474161454L) + ((int128)tmp_q[1] * 12436052620422L) - ((int128)tmp_q[2] * 1096807264321L) - ((int128)tmp_q[3] * 328932706782L) - ((int128)tmp_q[4] * 751005797345L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

