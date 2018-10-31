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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 542479798435667425UL) + ((((uint64_t)op[1] * 8656453118512067968UL) + ((uint64_t)op[2] * 11381030008431161371UL) + ((uint64_t)op[3] * 6658199394523247080UL) + ((uint64_t)op[4] * 2043390167181852647UL) + ((uint64_t)op[5] * 16525031290983770810UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 16525031290983770810UL) + ((uint64_t)op[1] * 542479798435667425UL) + ((((uint64_t)op[2] * 8656453118512067968UL) + ((uint64_t)op[3] * 11381030008431161371UL) + ((uint64_t)op[4] * 6658199394523247080UL) + ((uint64_t)op[5] * 2043390167181852647UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 2043390167181852647UL) + ((uint64_t)op[1] * 16525031290983770810UL) + ((uint64_t)op[2] * 542479798435667425UL) + ((((uint64_t)op[3] * 8656453118512067968UL) + ((uint64_t)op[4] * 11381030008431161371UL) + ((uint64_t)op[5] * 6658199394523247080UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 6658199394523247080UL) + ((uint64_t)op[1] * 2043390167181852647UL) + ((uint64_t)op[2] * 16525031290983770810UL) + ((uint64_t)op[3] * 542479798435667425UL) + ((((uint64_t)op[4] * 8656453118512067968UL) + ((uint64_t)op[5] * 11381030008431161371UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 11381030008431161371UL) + ((uint64_t)op[1] * 6658199394523247080UL) + ((uint64_t)op[2] * 2043390167181852647UL) + ((uint64_t)op[3] * 16525031290983770810UL) + ((uint64_t)op[4] * 542479798435667425UL) + ((uint64_t)op[5] * 3401513510056247040UL);
	tmp_q[5] = ((uint64_t)op[0] * 8656453118512067968UL) + ((uint64_t)op[1] * 11381030008431161371UL) + ((uint64_t)op[2] * 6658199394523247080UL) + ((uint64_t)op[3] * 2043390167181852647UL) + ((uint64_t)op[4] * 16525031290983770810UL) + ((uint64_t)op[5] * 542479798435667425UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4644547944003L) - ((((int128)tmp_q[1] * 1872997286338L) + ((int128)tmp_q[2] * 2010844898686L) - ((int128)tmp_q[3] * 4609314329884L) + ((int128)tmp_q[4] * 3672711338029L) - ((int128)tmp_q[5] * 3377709171030L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 3377709171030L) - ((int128)tmp_q[1] * 4644547944003L) - ((((int128)tmp_q[2] * 1872997286338L) + ((int128)tmp_q[3] * 2010844898686L) - ((int128)tmp_q[4] * 4609314329884L) + ((int128)tmp_q[5] * 3672711338029L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 3672711338029L) - ((int128)tmp_q[1] * 3377709171030L) - ((int128)tmp_q[2] * 4644547944003L) - ((((int128)tmp_q[3] * 1872997286338L) + ((int128)tmp_q[4] * 2010844898686L) - ((int128)tmp_q[5] * 4609314329884L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 4609314329884L) + ((int128)tmp_q[1] * 3672711338029L) - ((int128)tmp_q[2] * 3377709171030L) - ((int128)tmp_q[3] * 4644547944003L) - ((((int128)tmp_q[4] * 1872997286338L) + ((int128)tmp_q[5] * 2010844898686L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 2010844898686L) - ((int128)tmp_q[1] * 4609314329884L) + ((int128)tmp_q[2] * 3672711338029L) - ((int128)tmp_q[3] * 3377709171030L) - ((int128)tmp_q[4] * 4644547944003L) - ((int128)tmp_q[5] * 11237983718028L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1872997286338L) + ((int128)tmp_q[1] * 2010844898686L) - ((int128)tmp_q[2] * 4609314329884L) + ((int128)tmp_q[3] * 3672711338029L) - ((int128)tmp_q[4] * 3377709171030L) - ((int128)tmp_q[5] * 4644547944003L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

