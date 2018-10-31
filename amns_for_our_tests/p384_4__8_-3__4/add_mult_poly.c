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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11820736586833506044UL) + ((((uint64_t)op[1] * 18228292647586744779UL) + ((uint64_t)op[2] * 4997204453471229747UL) + ((uint64_t)op[3] * 14078546843282398177UL) + ((uint64_t)op[4] * 7671415241050607138UL) + ((uint64_t)op[5] * 12228041573823938292UL) + ((uint64_t)op[6] * 12287244903541561023UL) + ((uint64_t)op[7] * 16450833074681140551UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 16450833074681140551UL) + ((uint64_t)op[1] * 11820736586833506044UL) + ((((uint64_t)op[2] * 18228292647586744779UL) + ((uint64_t)op[3] * 4997204453471229747UL) + ((uint64_t)op[4] * 14078546843282398177UL) + ((uint64_t)op[5] * 7671415241050607138UL) + ((uint64_t)op[6] * 12228041573823938292UL) + ((uint64_t)op[7] * 12287244903541561023UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 12287244903541561023UL) + ((uint64_t)op[1] * 16450833074681140551UL) + ((uint64_t)op[2] * 11820736586833506044UL) + ((((uint64_t)op[3] * 18228292647586744779UL) + ((uint64_t)op[4] * 4997204453471229747UL) + ((uint64_t)op[5] * 14078546843282398177UL) + ((uint64_t)op[6] * 7671415241050607138UL) + ((uint64_t)op[7] * 12228041573823938292UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 12228041573823938292UL) + ((uint64_t)op[1] * 12287244903541561023UL) + ((uint64_t)op[2] * 16450833074681140551UL) + ((uint64_t)op[3] * 11820736586833506044UL) + ((((uint64_t)op[4] * 18228292647586744779UL) + ((uint64_t)op[5] * 4997204453471229747UL) + ((uint64_t)op[6] * 14078546843282398177UL) + ((uint64_t)op[7] * 7671415241050607138UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 7671415241050607138UL) + ((uint64_t)op[1] * 12228041573823938292UL) + ((uint64_t)op[2] * 12287244903541561023UL) + ((uint64_t)op[3] * 16450833074681140551UL) + ((uint64_t)op[4] * 11820736586833506044UL) + ((((uint64_t)op[5] * 18228292647586744779UL) + ((uint64_t)op[6] * 4997204453471229747UL) + ((uint64_t)op[7] * 14078546843282398177UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 14078546843282398177UL) + ((uint64_t)op[1] * 7671415241050607138UL) + ((uint64_t)op[2] * 12228041573823938292UL) + ((uint64_t)op[3] * 12287244903541561023UL) + ((uint64_t)op[4] * 16450833074681140551UL) + ((uint64_t)op[5] * 11820736586833506044UL) + ((((uint64_t)op[6] * 18228292647586744779UL) + ((uint64_t)op[7] * 4997204453471229747UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 4997204453471229747UL) + ((uint64_t)op[1] * 14078546843282398177UL) + ((uint64_t)op[2] * 7671415241050607138UL) + ((uint64_t)op[3] * 12228041573823938292UL) + ((uint64_t)op[4] * 12287244903541561023UL) + ((uint64_t)op[5] * 16450833074681140551UL) + ((uint64_t)op[6] * 11820736586833506044UL) + ((uint64_t)op[7] * 655354278368420511UL);
	tmp_q[7] = ((uint64_t)op[0] * 18228292647586744779UL) + ((uint64_t)op[1] * 4997204453471229747UL) + ((uint64_t)op[2] * 14078546843282398177UL) + ((uint64_t)op[3] * 7671415241050607138UL) + ((uint64_t)op[4] * 12228041573823938292UL) + ((uint64_t)op[5] * 12287244903541561023UL) + ((uint64_t)op[6] * 16450833074681140551UL) + ((uint64_t)op[7] * 11820736586833506044UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 23506084610933L) - ((-((int128)tmp_q[1] * 115814941855457L) - ((int128)tmp_q[2] * 191670083413542L) + ((int128)tmp_q[3] * 49949993018854L) + ((int128)tmp_q[4] * 52927152733635L) + ((int128)tmp_q[5] * 100924580009277L) - ((int128)tmp_q[6] * 28797846053500L) + ((int128)tmp_q[7] * 50049554693603L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 50049554693603L) - ((int128)tmp_q[1] * 23506084610933L) - ((-((int128)tmp_q[2] * 115814941855457L) - ((int128)tmp_q[3] * 191670083413542L) + ((int128)tmp_q[4] * 49949993018854L) + ((int128)tmp_q[5] * 52927152733635L) + ((int128)tmp_q[6] * 100924580009277L) - ((int128)tmp_q[7] * 28797846053500L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 28797846053500L) + ((int128)tmp_q[1] * 50049554693603L) - ((int128)tmp_q[2] * 23506084610933L) - ((-((int128)tmp_q[3] * 115814941855457L) - ((int128)tmp_q[4] * 191670083413542L) + ((int128)tmp_q[5] * 49949993018854L) + ((int128)tmp_q[6] * 52927152733635L) + ((int128)tmp_q[7] * 100924580009277L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 100924580009277L) - ((int128)tmp_q[1] * 28797846053500L) + ((int128)tmp_q[2] * 50049554693603L) - ((int128)tmp_q[3] * 23506084610933L) - ((-((int128)tmp_q[4] * 115814941855457L) - ((int128)tmp_q[5] * 191670083413542L) + ((int128)tmp_q[6] * 49949993018854L) + ((int128)tmp_q[7] * 52927152733635L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 52927152733635L) + ((int128)tmp_q[1] * 100924580009277L) - ((int128)tmp_q[2] * 28797846053500L) + ((int128)tmp_q[3] * 50049554693603L) - ((int128)tmp_q[4] * 23506084610933L) - ((-((int128)tmp_q[5] * 115814941855457L) - ((int128)tmp_q[6] * 191670083413542L) + ((int128)tmp_q[7] * 49949993018854L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 49949993018854L) + ((int128)tmp_q[1] * 52927152733635L) + ((int128)tmp_q[2] * 100924580009277L) - ((int128)tmp_q[3] * 28797846053500L) + ((int128)tmp_q[4] * 50049554693603L) - ((int128)tmp_q[5] * 23506084610933L) - ((-((int128)tmp_q[6] * 115814941855457L) - ((int128)tmp_q[7] * 191670083413542L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 191670083413542L) + ((int128)tmp_q[1] * 49949993018854L) + ((int128)tmp_q[2] * 52927152733635L) + ((int128)tmp_q[3] * 100924580009277L) - ((int128)tmp_q[4] * 28797846053500L) + ((int128)tmp_q[5] * 50049554693603L) - ((int128)tmp_q[6] * 23506084610933L) + ((int128)tmp_q[7] * 347444825566371L);
	tmp_zero[7] = -((int128)tmp_q[0] * 115814941855457L) - ((int128)tmp_q[1] * 191670083413542L) + ((int128)tmp_q[2] * 49949993018854L) + ((int128)tmp_q[3] * 52927152733635L) + ((int128)tmp_q[4] * 100924580009277L) - ((int128)tmp_q[5] * 28797846053500L) + ((int128)tmp_q[6] * 50049554693603L) - ((int128)tmp_q[7] * 23506084610933L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

