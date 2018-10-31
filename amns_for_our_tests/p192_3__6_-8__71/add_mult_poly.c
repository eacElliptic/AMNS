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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10474533078825869455UL) + ((((uint64_t)op[1] * 15145395330266705061UL) + ((uint64_t)op[2] * 8367450681941160190UL) + ((uint64_t)op[3] * 1734796831262786229UL) + ((uint64_t)op[4] * 256977571243125980UL) + ((uint64_t)op[5] * 18067408700023312420UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 18067408700023312420UL) + ((uint64_t)op[1] * 10474533078825869455UL) + ((((uint64_t)op[2] * 15145395330266705061UL) + ((uint64_t)op[3] * 8367450681941160190UL) + ((uint64_t)op[4] * 1734796831262786229UL) + ((uint64_t)op[5] * 256977571243125980UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 256977571243125980UL) + ((uint64_t)op[1] * 18067408700023312420UL) + ((uint64_t)op[2] * 10474533078825869455UL) + ((((uint64_t)op[3] * 15145395330266705061UL) + ((uint64_t)op[4] * 8367450681941160190UL) + ((uint64_t)op[5] * 1734796831262786229UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 1734796831262786229UL) + ((uint64_t)op[1] * 256977571243125980UL) + ((uint64_t)op[2] * 18067408700023312420UL) + ((uint64_t)op[3] * 10474533078825869455UL) + ((((uint64_t)op[4] * 15145395330266705061UL) + ((uint64_t)op[5] * 8367450681941160190UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 8367450681941160190UL) + ((uint64_t)op[1] * 1734796831262786229UL) + ((uint64_t)op[2] * 256977571243125980UL) + ((uint64_t)op[3] * 18067408700023312420UL) + ((uint64_t)op[4] * 10474533078825869455UL) + ((uint64_t)op[5] * 7964045873833220824UL);
	tmp_q[5] = ((uint64_t)op[0] * 15145395330266705061UL) + ((uint64_t)op[1] * 8367450681941160190UL) + ((uint64_t)op[2] * 1734796831262786229UL) + ((uint64_t)op[3] * 256977571243125980UL) + ((uint64_t)op[4] * 18067408700023312420UL) + ((uint64_t)op[5] * 10474533078825869455UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 27471823004361L) - ((((int128)tmp_q[1] * 82383079949893L) - ((int128)tmp_q[2] * 5820854906594L) - ((int128)tmp_q[3] * 39517259251475L) + ((int128)tmp_q[4] * 52033218555420L) - ((int128)tmp_q[5] * 40354399613212L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 40354399613212L) + ((int128)tmp_q[1] * 27471823004361L) - ((((int128)tmp_q[2] * 82383079949893L) - ((int128)tmp_q[3] * 5820854906594L) - ((int128)tmp_q[4] * 39517259251475L) + ((int128)tmp_q[5] * 52033218555420L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 52033218555420L) - ((int128)tmp_q[1] * 40354399613212L) + ((int128)tmp_q[2] * 27471823004361L) - ((((int128)tmp_q[3] * 82383079949893L) - ((int128)tmp_q[4] * 5820854906594L) - ((int128)tmp_q[5] * 39517259251475L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 39517259251475L) + ((int128)tmp_q[1] * 52033218555420L) - ((int128)tmp_q[2] * 40354399613212L) + ((int128)tmp_q[3] * 27471823004361L) - ((((int128)tmp_q[4] * 82383079949893L) - ((int128)tmp_q[5] * 5820854906594L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 5820854906594L) - ((int128)tmp_q[1] * 39517259251475L) + ((int128)tmp_q[2] * 52033218555420L) - ((int128)tmp_q[3] * 40354399613212L) + ((int128)tmp_q[4] * 27471823004361L) - ((int128)tmp_q[5] * 659064639599144L);
	tmp_zero[5] = ((int128)tmp_q[0] * 82383079949893L) - ((int128)tmp_q[1] * 5820854906594L) - ((int128)tmp_q[2] * 39517259251475L) + ((int128)tmp_q[3] * 52033218555420L) - ((int128)tmp_q[4] * 40354399613212L) + ((int128)tmp_q[5] * 27471823004361L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

