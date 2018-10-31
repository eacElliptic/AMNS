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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7131364654903183533UL) + ((((uint64_t)op[1] * 17484520504772810064UL) + ((uint64_t)op[2] * 13693451780488243106UL) + ((uint64_t)op[3] * 18216924798689426524UL) + ((uint64_t)op[4] * 4051513795887237651UL) + ((uint64_t)op[5] * 9207144138968656777UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 9207144138968656777UL) + ((uint64_t)op[1] * 7131364654903183533UL) + ((((uint64_t)op[2] * 17484520504772810064UL) + ((uint64_t)op[3] * 13693451780488243106UL) + ((uint64_t)op[4] * 18216924798689426524UL) + ((uint64_t)op[5] * 4051513795887237651UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 4051513795887237651UL) + ((uint64_t)op[1] * 9207144138968656777UL) + ((uint64_t)op[2] * 7131364654903183533UL) + ((((uint64_t)op[3] * 17484520504772810064UL) + ((uint64_t)op[4] * 13693451780488243106UL) + ((uint64_t)op[5] * 18216924798689426524UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 18216924798689426524UL) + ((uint64_t)op[1] * 4051513795887237651UL) + ((uint64_t)op[2] * 9207144138968656777UL) + ((uint64_t)op[3] * 7131364654903183533UL) + ((((uint64_t)op[4] * 17484520504772810064UL) + ((uint64_t)op[5] * 13693451780488243106UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 13693451780488243106UL) + ((uint64_t)op[1] * 18216924798689426524UL) + ((uint64_t)op[2] * 4051513795887237651UL) + ((uint64_t)op[3] * 9207144138968656777UL) + ((uint64_t)op[4] * 7131364654903183533UL) + ((uint64_t)op[5] * 1924447137873483104UL);
	tmp_q[5] = ((uint64_t)op[0] * 17484520504772810064UL) + ((uint64_t)op[1] * 13693451780488243106UL) + ((uint64_t)op[2] * 18216924798689426524UL) + ((uint64_t)op[3] * 4051513795887237651UL) + ((uint64_t)op[4] * 9207144138968656777UL) + ((uint64_t)op[5] * 7131364654903183533UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1322781628165L) - ((((int128)tmp_q[1] * 4855639246000L) - ((int128)tmp_q[2] * 335845104773L) + ((int128)tmp_q[3] * 3620888444965L) + ((int128)tmp_q[4] * 2358577929202L) + ((int128)tmp_q[5] * 2328009936011L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 2328009936011L) + ((int128)tmp_q[1] * 1322781628165L) - ((((int128)tmp_q[2] * 4855639246000L) - ((int128)tmp_q[3] * 335845104773L) + ((int128)tmp_q[4] * 3620888444965L) + ((int128)tmp_q[5] * 2358577929202L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 2358577929202L) + ((int128)tmp_q[1] * 2328009936011L) + ((int128)tmp_q[2] * 1322781628165L) - ((((int128)tmp_q[3] * 4855639246000L) - ((int128)tmp_q[4] * 335845104773L) + ((int128)tmp_q[5] * 3620888444965L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 3620888444965L) + ((int128)tmp_q[1] * 2358577929202L) + ((int128)tmp_q[2] * 2328009936011L) + ((int128)tmp_q[3] * 1322781628165L) - ((((int128)tmp_q[4] * 4855639246000L) - ((int128)tmp_q[5] * 335845104773L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 335845104773L) + ((int128)tmp_q[1] * 3620888444965L) + ((int128)tmp_q[2] * 2358577929202L) + ((int128)tmp_q[3] * 2328009936011L) + ((int128)tmp_q[4] * 1322781628165L) - ((int128)tmp_q[5] * 9711278492000L);
	tmp_zero[5] = ((int128)tmp_q[0] * 4855639246000L) - ((int128)tmp_q[1] * 335845104773L) + ((int128)tmp_q[2] * 3620888444965L) + ((int128)tmp_q[3] * 2358577929202L) + ((int128)tmp_q[4] * 2328009936011L) + ((int128)tmp_q[5] * 1322781628165L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

