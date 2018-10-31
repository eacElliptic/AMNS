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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11497820160928870070UL) + ((((uint64_t)op[1] * 6056848936982990146UL) + ((uint64_t)op[2] * 15797286168697067516UL) + ((uint64_t)op[3] * 12793909848839756840UL) + ((uint64_t)op[4] * 11129678055050142041UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 11129678055050142041UL) + ((uint64_t)op[1] * 11497820160928870070UL) + ((((uint64_t)op[2] * 6056848936982990146UL) + ((uint64_t)op[3] * 15797286168697067516UL) + ((uint64_t)op[4] * 12793909848839756840UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 12793909848839756840UL) + ((uint64_t)op[1] * 11129678055050142041UL) + ((uint64_t)op[2] * 11497820160928870070UL) + ((((uint64_t)op[3] * 6056848936982990146UL) + ((uint64_t)op[4] * 15797286168697067516UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 15797286168697067516UL) + ((uint64_t)op[1] * 12793909848839756840UL) + ((uint64_t)op[2] * 11129678055050142041UL) + ((uint64_t)op[3] * 11497820160928870070UL) + ((uint64_t)op[4] * 276197262760581178UL);
	tmp_q[4] = ((uint64_t)op[0] * 6056848936982990146UL) + ((uint64_t)op[1] * 15797286168697067516UL) + ((uint64_t)op[2] * 12793909848839756840UL) + ((uint64_t)op[3] * 11129678055050142041UL) + ((uint64_t)op[4] * 11497820160928870070UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 11242198948180L) - ((((int128)tmp_q[1] * 11219809329387L) - ((int128)tmp_q[2] * 14602682223866L) + ((int128)tmp_q[3] * 16040950789118L) - ((int128)tmp_q[4] * 5216603543540L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 5216603543540L) - ((int128)tmp_q[1] * 11242198948180L) - ((((int128)tmp_q[2] * 11219809329387L) - ((int128)tmp_q[3] * 14602682223866L) + ((int128)tmp_q[4] * 16040950789118L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 16040950789118L) - ((int128)tmp_q[1] * 5216603543540L) - ((int128)tmp_q[2] * 11242198948180L) - ((((int128)tmp_q[3] * 11219809329387L) - ((int128)tmp_q[4] * 14602682223866L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 14602682223866L) + ((int128)tmp_q[1] * 16040950789118L) - ((int128)tmp_q[2] * 5216603543540L) - ((int128)tmp_q[3] * 11242198948180L) - ((int128)tmp_q[4] * 33659427988161L);
	tmp_zero[4] = ((int128)tmp_q[0] * 11219809329387L) - ((int128)tmp_q[1] * 14602682223866L) + ((int128)tmp_q[2] * 16040950789118L) - ((int128)tmp_q[3] * 5216603543540L) - ((int128)tmp_q[4] * 11242198948180L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

