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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[3] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5176777346965077543UL) + ((((uint64_t)op[1] * 811001465455443571UL) + ((uint64_t)op[2] * 748709281601232991UL) + ((uint64_t)op[3] * 5041497157299708546UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 5041497157299708546UL) + ((uint64_t)op[1] * 5176777346965077543UL) + ((((uint64_t)op[2] * 811001465455443571UL) + ((uint64_t)op[3] * 748709281601232991UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 748709281601232991UL) + ((uint64_t)op[1] * 5041497157299708546UL) + ((uint64_t)op[2] * 5176777346965077543UL) + ((uint64_t)op[3] * 14391736746432333761UL);
	tmp_q[3] = ((uint64_t)op[0] * 811001465455443571UL) + ((uint64_t)op[1] * 748709281601232991UL) + ((uint64_t)op[2] * 5041497157299708546UL) + ((uint64_t)op[3] * 5176777346965077543UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 213611517574501L) - ((((int128)tmp_q[1] * 27798708696234L) - ((int128)tmp_q[2] * 18157284594551L) + ((int128)tmp_q[3] * 53358293763193L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 53358293763193L) - ((int128)tmp_q[1] * 213611517574501L) - ((((int128)tmp_q[2] * 27798708696234L) - ((int128)tmp_q[3] * 18157284594551L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 18157284594551L) + ((int128)tmp_q[1] * 53358293763193L) - ((int128)tmp_q[2] * 213611517574501L) - ((int128)tmp_q[3] * 138993543481170L);
	tmp_zero[3] = ((int128)tmp_q[0] * 27798708696234L) - ((int128)tmp_q[1] * 18157284594551L) + ((int128)tmp_q[2] * 53358293763193L) - ((int128)tmp_q[3] * 213611517574501L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

