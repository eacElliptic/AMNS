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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11337178660009888593UL) + ((((uint64_t)op[1] * 15349282852981703588UL) + ((uint64_t)op[2] * 15039242982457981184UL) + ((uint64_t)op[3] * 7307750830529502813UL) + ((uint64_t)op[4] * 17502397674372328731UL) + ((uint64_t)op[5] * 15165559855966190264UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 15165559855966190264UL) + ((uint64_t)op[1] * 11337178660009888593UL) + ((((uint64_t)op[2] * 15349282852981703588UL) + ((uint64_t)op[3] * 15039242982457981184UL) + ((uint64_t)op[4] * 7307750830529502813UL) + ((uint64_t)op[5] * 17502397674372328731UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 17502397674372328731UL) + ((uint64_t)op[1] * 15165559855966190264UL) + ((uint64_t)op[2] * 11337178660009888593UL) + ((((uint64_t)op[3] * 15349282852981703588UL) + ((uint64_t)op[4] * 15039242982457981184UL) + ((uint64_t)op[5] * 7307750830529502813UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 7307750830529502813UL) + ((uint64_t)op[1] * 17502397674372328731UL) + ((uint64_t)op[2] * 15165559855966190264UL) + ((uint64_t)op[3] * 11337178660009888593UL) + ((((uint64_t)op[4] * 15349282852981703588UL) + ((uint64_t)op[5] * 15039242982457981184UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 15039242982457981184UL) + ((uint64_t)op[1] * 7307750830529502813UL) + ((uint64_t)op[2] * 17502397674372328731UL) + ((uint64_t)op[3] * 15165559855966190264UL) + ((uint64_t)op[4] * 11337178660009888593UL) + ((uint64_t)op[5] * 12251821632253855560UL);
	tmp_q[5] = ((uint64_t)op[0] * 15349282852981703588UL) + ((uint64_t)op[1] * 15039242982457981184UL) + ((uint64_t)op[2] * 7307750830529502813UL) + ((uint64_t)op[3] * 17502397674372328731UL) + ((uint64_t)op[4] * 15165559855966190264UL) + ((uint64_t)op[5] * 11337178660009888593UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4385489177797L) + ((((int128)tmp_q[1] * 1313541179372L) + ((int128)tmp_q[2] * 1265463542125L) + ((int128)tmp_q[3] * 72185955151L) + ((int128)tmp_q[4] * 1627017180311L) + ((int128)tmp_q[5] * 4683584214342L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 4683584214342L) - ((int128)tmp_q[1] * 4385489177797L) + ((((int128)tmp_q[2] * 1313541179372L) + ((int128)tmp_q[3] * 1265463542125L) + ((int128)tmp_q[4] * 72185955151L) + ((int128)tmp_q[5] * 1627017180311L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 1627017180311L) + ((int128)tmp_q[1] * 4683584214342L) - ((int128)tmp_q[2] * 4385489177797L) + ((((int128)tmp_q[3] * 1313541179372L) + ((int128)tmp_q[4] * 1265463542125L) + ((int128)tmp_q[5] * 72185955151L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 72185955151L) + ((int128)tmp_q[1] * 1627017180311L) + ((int128)tmp_q[2] * 4683584214342L) - ((int128)tmp_q[3] * 4385489177797L) + ((((int128)tmp_q[4] * 1313541179372L) + ((int128)tmp_q[5] * 1265463542125L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 1265463542125L) + ((int128)tmp_q[1] * 72185955151L) + ((int128)tmp_q[2] * 1627017180311L) + ((int128)tmp_q[3] * 4683584214342L) - ((int128)tmp_q[4] * 4385489177797L) + ((int128)tmp_q[5] * 2627082358744L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1313541179372L) + ((int128)tmp_q[1] * 1265463542125L) + ((int128)tmp_q[2] * 72185955151L) + ((int128)tmp_q[3] * 1627017180311L) + ((int128)tmp_q[4] * 4683584214342L) - ((int128)tmp_q[5] * 4385489177797L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

