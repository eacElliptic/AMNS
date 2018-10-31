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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9625339970277592055UL) + ((((uint64_t)op[1] * 3154969694548689497UL) + ((uint64_t)op[2] * 11434481810808401973UL) + ((uint64_t)op[3] * 5177244435080169204UL) + ((uint64_t)op[4] * 17123497413820320169UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 17123497413820320169UL) + ((uint64_t)op[1] * 9625339970277592055UL) + ((((uint64_t)op[2] * 3154969694548689497UL) + ((uint64_t)op[3] * 11434481810808401973UL) + ((uint64_t)op[4] * 5177244435080169204UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 5177244435080169204UL) + ((uint64_t)op[1] * 17123497413820320169UL) + ((uint64_t)op[2] * 9625339970277592055UL) + ((((uint64_t)op[3] * 3154969694548689497UL) + ((uint64_t)op[4] * 11434481810808401973UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 11434481810808401973UL) + ((uint64_t)op[1] * 5177244435080169204UL) + ((uint64_t)op[2] * 17123497413820320169UL) + ((uint64_t)op[3] * 9625339970277592055UL) + ((uint64_t)op[4] * 12619878778194757988UL);
	tmp_q[4] = ((uint64_t)op[0] * 3154969694548689497UL) + ((uint64_t)op[1] * 11434481810808401973UL) + ((uint64_t)op[2] * 5177244435080169204UL) + ((uint64_t)op[3] * 17123497413820320169UL) + ((uint64_t)op[4] * 9625339970277592055UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 287472562495L) + ((((int128)tmp_q[1] * 59128153544L) - ((int128)tmp_q[2] * 123104500710L) - ((int128)tmp_q[3] * 22284903031L) - ((int128)tmp_q[4] * 182431291251L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 182431291251L) - ((int128)tmp_q[1] * 287472562495L) + ((((int128)tmp_q[2] * 59128153544L) - ((int128)tmp_q[3] * 123104500710L) - ((int128)tmp_q[4] * 22284903031L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 22284903031L) - ((int128)tmp_q[1] * 182431291251L) - ((int128)tmp_q[2] * 287472562495L) + ((((int128)tmp_q[3] * 59128153544L) - ((int128)tmp_q[4] * 123104500710L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 123104500710L) - ((int128)tmp_q[1] * 22284903031L) - ((int128)tmp_q[2] * 182431291251L) - ((int128)tmp_q[3] * 287472562495L) + ((int128)tmp_q[4] * 236512614176L);
	tmp_zero[4] = ((int128)tmp_q[0] * 59128153544L) - ((int128)tmp_q[1] * 123104500710L) - ((int128)tmp_q[2] * 22284903031L) - ((int128)tmp_q[3] * 182431291251L) - ((int128)tmp_q[4] * 287472562495L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

