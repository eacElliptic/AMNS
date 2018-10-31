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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15975358528394620521UL) + ((((uint64_t)op[1] * 6451265420875091333UL) + ((uint64_t)op[2] * 2421498852260413766UL) + ((uint64_t)op[3] * 13432822657732623350UL) + ((uint64_t)op[4] * 2177146259598475085UL) + ((uint64_t)op[5] * 9031564410494283146UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 9031564410494283146UL) + ((uint64_t)op[1] * 15975358528394620521UL) + ((((uint64_t)op[2] * 6451265420875091333UL) + ((uint64_t)op[3] * 2421498852260413766UL) + ((uint64_t)op[4] * 13432822657732623350UL) + ((uint64_t)op[5] * 2177146259598475085UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 2177146259598475085UL) + ((uint64_t)op[1] * 9031564410494283146UL) + ((uint64_t)op[2] * 15975358528394620521UL) + ((((uint64_t)op[3] * 6451265420875091333UL) + ((uint64_t)op[4] * 2421498852260413766UL) + ((uint64_t)op[5] * 13432822657732623350UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 13432822657732623350UL) + ((uint64_t)op[1] * 2177146259598475085UL) + ((uint64_t)op[2] * 9031564410494283146UL) + ((uint64_t)op[3] * 15975358528394620521UL) + ((((uint64_t)op[4] * 6451265420875091333UL) + ((uint64_t)op[5] * 2421498852260413766UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 2421498852260413766UL) + ((uint64_t)op[1] * 13432822657732623350UL) + ((uint64_t)op[2] * 2177146259598475085UL) + ((uint64_t)op[3] * 9031564410494283146UL) + ((uint64_t)op[4] * 15975358528394620521UL) + ((uint64_t)op[5] * 4637161043043646567UL);
	tmp_q[5] = ((uint64_t)op[0] * 6451265420875091333UL) + ((uint64_t)op[1] * 2421498852260413766UL) + ((uint64_t)op[2] * 13432822657732623350UL) + ((uint64_t)op[3] * 2177146259598475085UL) + ((uint64_t)op[4] * 9031564410494283146UL) + ((uint64_t)op[5] * 15975358528394620521UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 221882031L) - ((-((int128)tmp_q[1] * 1780603447L) + ((int128)tmp_q[2] * 11851842L) + ((int128)tmp_q[3] * 1553534248L) - ((int128)tmp_q[4] * 717975683L) - ((int128)tmp_q[5] * 2523875552L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 2523875552L) - ((int128)tmp_q[1] * 221882031L) - ((-((int128)tmp_q[2] * 1780603447L) + ((int128)tmp_q[3] * 11851842L) + ((int128)tmp_q[4] * 1553534248L) - ((int128)tmp_q[5] * 717975683L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 717975683L) - ((int128)tmp_q[1] * 2523875552L) - ((int128)tmp_q[2] * 221882031L) - ((-((int128)tmp_q[3] * 1780603447L) + ((int128)tmp_q[4] * 11851842L) + ((int128)tmp_q[5] * 1553534248L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 1553534248L) - ((int128)tmp_q[1] * 717975683L) - ((int128)tmp_q[2] * 2523875552L) - ((int128)tmp_q[3] * 221882031L) - ((-((int128)tmp_q[4] * 1780603447L) + ((int128)tmp_q[5] * 11851842L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 11851842L) + ((int128)tmp_q[1] * 1553534248L) - ((int128)tmp_q[2] * 717975683L) - ((int128)tmp_q[3] * 2523875552L) - ((int128)tmp_q[4] * 221882031L) + ((int128)tmp_q[5] * 8903017235L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1780603447L) + ((int128)tmp_q[1] * 11851842L) + ((int128)tmp_q[2] * 1553534248L) - ((int128)tmp_q[3] * 717975683L) - ((int128)tmp_q[4] * 2523875552L) - ((int128)tmp_q[5] * 221882031L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

