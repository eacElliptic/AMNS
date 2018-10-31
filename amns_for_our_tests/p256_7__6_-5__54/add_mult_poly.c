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
	tmp_q[0] = ((uint64_t)op[0] * 14511390028919607225UL) + ((((uint64_t)op[1] * 3857959091845338816UL) + ((uint64_t)op[2] * 15762030743693339941UL) + ((uint64_t)op[3] * 13660541380005041413UL) + ((uint64_t)op[4] * 14837632725230805438UL) + ((uint64_t)op[5] * 5623985744115279896UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 5623985744115279896UL) + ((uint64_t)op[1] * 14511390028919607225UL) + ((((uint64_t)op[2] * 3857959091845338816UL) + ((uint64_t)op[3] * 15762030743693339941UL) + ((uint64_t)op[4] * 13660541380005041413UL) + ((uint64_t)op[5] * 14837632725230805438UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 14837632725230805438UL) + ((uint64_t)op[1] * 5623985744115279896UL) + ((uint64_t)op[2] * 14511390028919607225UL) + ((((uint64_t)op[3] * 3857959091845338816UL) + ((uint64_t)op[4] * 15762030743693339941UL) + ((uint64_t)op[5] * 13660541380005041413UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 13660541380005041413UL) + ((uint64_t)op[1] * 14837632725230805438UL) + ((uint64_t)op[2] * 5623985744115279896UL) + ((uint64_t)op[3] * 14511390028919607225UL) + ((((uint64_t)op[4] * 3857959091845338816UL) + ((uint64_t)op[5] * 15762030743693339941UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 15762030743693339941UL) + ((uint64_t)op[1] * 13660541380005041413UL) + ((uint64_t)op[2] * 14837632725230805438UL) + ((uint64_t)op[3] * 5623985744115279896UL) + ((uint64_t)op[4] * 14511390028919607225UL) + ((uint64_t)op[5] * 17603692688192409152UL);
	tmp_q[5] = ((uint64_t)op[0] * 3857959091845338816UL) + ((uint64_t)op[1] * 15762030743693339941UL) + ((uint64_t)op[2] * 13660541380005041413UL) + ((uint64_t)op[3] * 14837632725230805438UL) + ((uint64_t)op[4] * 5623985744115279896UL) + ((uint64_t)op[5] * 14511390028919607225UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2406021202802L) - ((((int128)tmp_q[1] * 1051572001198L) + ((int128)tmp_q[2] * 572430178697L) + ((int128)tmp_q[3] * 1988284149382L) + ((int128)tmp_q[4] * 4123636896107L) + ((int128)tmp_q[5] * 3909736030771L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 3909736030771L) - ((int128)tmp_q[1] * 2406021202802L) - ((((int128)tmp_q[2] * 1051572001198L) + ((int128)tmp_q[3] * 572430178697L) + ((int128)tmp_q[4] * 1988284149382L) + ((int128)tmp_q[5] * 4123636896107L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 4123636896107L) + ((int128)tmp_q[1] * 3909736030771L) - ((int128)tmp_q[2] * 2406021202802L) - ((((int128)tmp_q[3] * 1051572001198L) + ((int128)tmp_q[4] * 572430178697L) + ((int128)tmp_q[5] * 1988284149382L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 1988284149382L) + ((int128)tmp_q[1] * 4123636896107L) + ((int128)tmp_q[2] * 3909736030771L) - ((int128)tmp_q[3] * 2406021202802L) - ((((int128)tmp_q[4] * 1051572001198L) + ((int128)tmp_q[5] * 572430178697L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 572430178697L) + ((int128)tmp_q[1] * 1988284149382L) + ((int128)tmp_q[2] * 4123636896107L) + ((int128)tmp_q[3] * 3909736030771L) - ((int128)tmp_q[4] * 2406021202802L) - ((int128)tmp_q[5] * 5257860005990L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1051572001198L) + ((int128)tmp_q[1] * 572430178697L) + ((int128)tmp_q[2] * 1988284149382L) + ((int128)tmp_q[3] * 4123636896107L) + ((int128)tmp_q[4] * 3909736030771L) - ((int128)tmp_q[5] * 2406021202802L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

