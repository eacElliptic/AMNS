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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6469446716430804343UL) + ((((uint64_t)op[1] * 12201639717172045805UL) + ((uint64_t)op[2] * 14205360562393552975UL) + ((uint64_t)op[3] * 13047828430016713343UL) + ((uint64_t)op[4] * 10400004969933974769UL) + ((uint64_t)op[5] * 1600996880961475564UL) + ((uint64_t)op[6] * 10568795448901603256UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 10568795448901603256UL) + ((uint64_t)op[1] * 6469446716430804343UL) + ((((uint64_t)op[2] * 12201639717172045805UL) + ((uint64_t)op[3] * 14205360562393552975UL) + ((uint64_t)op[4] * 13047828430016713343UL) + ((uint64_t)op[5] * 10400004969933974769UL) + ((uint64_t)op[6] * 1600996880961475564UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 1600996880961475564UL) + ((uint64_t)op[1] * 10568795448901603256UL) + ((uint64_t)op[2] * 6469446716430804343UL) + ((((uint64_t)op[3] * 12201639717172045805UL) + ((uint64_t)op[4] * 14205360562393552975UL) + ((uint64_t)op[5] * 13047828430016713343UL) + ((uint64_t)op[6] * 10400004969933974769UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 10400004969933974769UL) + ((uint64_t)op[1] * 1600996880961475564UL) + ((uint64_t)op[2] * 10568795448901603256UL) + ((uint64_t)op[3] * 6469446716430804343UL) + ((((uint64_t)op[4] * 12201639717172045805UL) + ((uint64_t)op[5] * 14205360562393552975UL) + ((uint64_t)op[6] * 13047828430016713343UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 13047828430016713343UL) + ((uint64_t)op[1] * 10400004969933974769UL) + ((uint64_t)op[2] * 1600996880961475564UL) + ((uint64_t)op[3] * 10568795448901603256UL) + ((uint64_t)op[4] * 6469446716430804343UL) + ((((uint64_t)op[5] * 12201639717172045805UL) + ((uint64_t)op[6] * 14205360562393552975UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 14205360562393552975UL) + ((uint64_t)op[1] * 13047828430016713343UL) + ((uint64_t)op[2] * 10400004969933974769UL) + ((uint64_t)op[3] * 1600996880961475564UL) + ((uint64_t)op[4] * 10568795448901603256UL) + ((uint64_t)op[5] * 6469446716430804343UL) + ((uint64_t)op[6] * 12778777708977977439UL);
	tmp_q[6] = ((uint64_t)op[0] * 12201639717172045805UL) + ((uint64_t)op[1] * 14205360562393552975UL) + ((uint64_t)op[2] * 13047828430016713343UL) + ((uint64_t)op[3] * 10400004969933974769UL) + ((uint64_t)op[4] * 1600996880961475564UL) + ((uint64_t)op[5] * 10568795448901603256UL) + ((uint64_t)op[6] * 6469446716430804343UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 47997668399L) - ((-((int128)tmp_q[1] * 6965071870L) + ((int128)tmp_q[2] * 52234942448L) + ((int128)tmp_q[3] * 4749388067L) + ((int128)tmp_q[4] * 58827936760L) + ((int128)tmp_q[5] * 17142729153L) - ((int128)tmp_q[6] * 42713997498L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 42713997498L) - ((int128)tmp_q[1] * 47997668399L) - ((-((int128)tmp_q[2] * 6965071870L) + ((int128)tmp_q[3] * 52234942448L) + ((int128)tmp_q[4] * 4749388067L) + ((int128)tmp_q[5] * 58827936760L) + ((int128)tmp_q[6] * 17142729153L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 17142729153L) - ((int128)tmp_q[1] * 42713997498L) - ((int128)tmp_q[2] * 47997668399L) - ((-((int128)tmp_q[3] * 6965071870L) + ((int128)tmp_q[4] * 52234942448L) + ((int128)tmp_q[5] * 4749388067L) + ((int128)tmp_q[6] * 58827936760L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 58827936760L) + ((int128)tmp_q[1] * 17142729153L) - ((int128)tmp_q[2] * 42713997498L) - ((int128)tmp_q[3] * 47997668399L) - ((-((int128)tmp_q[4] * 6965071870L) + ((int128)tmp_q[5] * 52234942448L) + ((int128)tmp_q[6] * 4749388067L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 4749388067L) + ((int128)tmp_q[1] * 58827936760L) + ((int128)tmp_q[2] * 17142729153L) - ((int128)tmp_q[3] * 42713997498L) - ((int128)tmp_q[4] * 47997668399L) - ((-((int128)tmp_q[5] * 6965071870L) + ((int128)tmp_q[6] * 52234942448L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 52234942448L) + ((int128)tmp_q[1] * 4749388067L) + ((int128)tmp_q[2] * 58827936760L) + ((int128)tmp_q[3] * 17142729153L) - ((int128)tmp_q[4] * 42713997498L) - ((int128)tmp_q[5] * 47997668399L) + ((int128)tmp_q[6] * 34825359350L);
	tmp_zero[6] = -((int128)tmp_q[0] * 6965071870L) + ((int128)tmp_q[1] * 52234942448L) + ((int128)tmp_q[2] * 4749388067L) + ((int128)tmp_q[3] * 58827936760L) + ((int128)tmp_q[4] * 17142729153L) - ((int128)tmp_q[5] * 42713997498L) - ((int128)tmp_q[6] * 47997668399L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

