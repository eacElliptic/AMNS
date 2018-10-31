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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11668275959916591091UL) + ((((uint64_t)op[1] * 15521776687693108466UL) + ((uint64_t)op[2] * 680279334860060282UL) + ((uint64_t)op[3] * 1665862681096992536UL) + ((uint64_t)op[4] * 16370171215657301964UL) + ((uint64_t)op[5] * 10757411797891729406UL) + ((uint64_t)op[6] * 2060703917908691828UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 2060703917908691828UL) + ((uint64_t)op[1] * 11668275959916591091UL) + ((((uint64_t)op[2] * 15521776687693108466UL) + ((uint64_t)op[3] * 680279334860060282UL) + ((uint64_t)op[4] * 1665862681096992536UL) + ((uint64_t)op[5] * 16370171215657301964UL) + ((uint64_t)op[6] * 10757411797891729406UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 10757411797891729406UL) + ((uint64_t)op[1] * 2060703917908691828UL) + ((uint64_t)op[2] * 11668275959916591091UL) + ((((uint64_t)op[3] * 15521776687693108466UL) + ((uint64_t)op[4] * 680279334860060282UL) + ((uint64_t)op[5] * 1665862681096992536UL) + ((uint64_t)op[6] * 16370171215657301964UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 16370171215657301964UL) + ((uint64_t)op[1] * 10757411797891729406UL) + ((uint64_t)op[2] * 2060703917908691828UL) + ((uint64_t)op[3] * 11668275959916591091UL) + ((((uint64_t)op[4] * 15521776687693108466UL) + ((uint64_t)op[5] * 680279334860060282UL) + ((uint64_t)op[6] * 1665862681096992536UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 1665862681096992536UL) + ((uint64_t)op[1] * 16370171215657301964UL) + ((uint64_t)op[2] * 10757411797891729406UL) + ((uint64_t)op[3] * 2060703917908691828UL) + ((uint64_t)op[4] * 11668275959916591091UL) + ((((uint64_t)op[5] * 15521776687693108466UL) + ((uint64_t)op[6] * 680279334860060282UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 680279334860060282UL) + ((uint64_t)op[1] * 1665862681096992536UL) + ((uint64_t)op[2] * 16370171215657301964UL) + ((uint64_t)op[3] * 10757411797891729406UL) + ((uint64_t)op[4] * 2060703917908691828UL) + ((uint64_t)op[5] * 11668275959916591091UL) + ((uint64_t)op[6] * 4952995014421993584UL);
	tmp_q[6] = ((uint64_t)op[0] * 15521776687693108466UL) + ((uint64_t)op[1] * 680279334860060282UL) + ((uint64_t)op[2] * 1665862681096992536UL) + ((uint64_t)op[3] * 16370171215657301964UL) + ((uint64_t)op[4] * 10757411797891729406UL) + ((uint64_t)op[5] * 2060703917908691828UL) + ((uint64_t)op[6] * 11668275959916591091UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 30737191163L) - ((-((int128)tmp_q[1] * 22230314982L) - ((int128)tmp_q[2] * 36865841398L) - ((int128)tmp_q[3] * 57426008180L) + ((int128)tmp_q[4] * 65910221628L) + ((int128)tmp_q[5] * 2598467038L) + ((int128)tmp_q[6] * 46388043668L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 46388043668L) - ((int128)tmp_q[1] * 30737191163L) - ((-((int128)tmp_q[2] * 22230314982L) - ((int128)tmp_q[3] * 36865841398L) - ((int128)tmp_q[4] * 57426008180L) + ((int128)tmp_q[5] * 65910221628L) + ((int128)tmp_q[6] * 2598467038L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 2598467038L) + ((int128)tmp_q[1] * 46388043668L) - ((int128)tmp_q[2] * 30737191163L) - ((-((int128)tmp_q[3] * 22230314982L) - ((int128)tmp_q[4] * 36865841398L) - ((int128)tmp_q[5] * 57426008180L) + ((int128)tmp_q[6] * 65910221628L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 65910221628L) + ((int128)tmp_q[1] * 2598467038L) + ((int128)tmp_q[2] * 46388043668L) - ((int128)tmp_q[3] * 30737191163L) - ((-((int128)tmp_q[4] * 22230314982L) - ((int128)tmp_q[5] * 36865841398L) - ((int128)tmp_q[6] * 57426008180L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 57426008180L) + ((int128)tmp_q[1] * 65910221628L) + ((int128)tmp_q[2] * 2598467038L) + ((int128)tmp_q[3] * 46388043668L) - ((int128)tmp_q[4] * 30737191163L) - ((-((int128)tmp_q[5] * 22230314982L) - ((int128)tmp_q[6] * 36865841398L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 36865841398L) - ((int128)tmp_q[1] * 57426008180L) + ((int128)tmp_q[2] * 65910221628L) + ((int128)tmp_q[3] * 2598467038L) + ((int128)tmp_q[4] * 46388043668L) - ((int128)tmp_q[5] * 30737191163L) + ((int128)tmp_q[6] * 177842519856L);
	tmp_zero[6] = -((int128)tmp_q[0] * 22230314982L) - ((int128)tmp_q[1] * 36865841398L) - ((int128)tmp_q[2] * 57426008180L) + ((int128)tmp_q[3] * 65910221628L) + ((int128)tmp_q[4] * 2598467038L) + ((int128)tmp_q[5] * 46388043668L) - ((int128)tmp_q[6] * 30737191163L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

