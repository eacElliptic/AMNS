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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10174242020965357758UL) + ((((uint64_t)op[1] * 5825884965231322132UL) + ((uint64_t)op[2] * 7422492898454529783UL) + ((uint64_t)op[3] * 221425258653370777UL) + ((uint64_t)op[4] * 2333700255738950699UL) + ((uint64_t)op[5] * 10811191460477622003UL) + ((uint64_t)op[6] * 18069544045833391061UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 18069544045833391061UL) + ((uint64_t)op[1] * 10174242020965357758UL) + ((((uint64_t)op[2] * 5825884965231322132UL) + ((uint64_t)op[3] * 7422492898454529783UL) + ((uint64_t)op[4] * 221425258653370777UL) + ((uint64_t)op[5] * 2333700255738950699UL) + ((uint64_t)op[6] * 10811191460477622003UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 10811191460477622003UL) + ((uint64_t)op[1] * 18069544045833391061UL) + ((uint64_t)op[2] * 10174242020965357758UL) + ((((uint64_t)op[3] * 5825884965231322132UL) + ((uint64_t)op[4] * 7422492898454529783UL) + ((uint64_t)op[5] * 221425258653370777UL) + ((uint64_t)op[6] * 2333700255738950699UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 2333700255738950699UL) + ((uint64_t)op[1] * 10811191460477622003UL) + ((uint64_t)op[2] * 18069544045833391061UL) + ((uint64_t)op[3] * 10174242020965357758UL) + ((((uint64_t)op[4] * 5825884965231322132UL) + ((uint64_t)op[5] * 7422492898454529783UL) + ((uint64_t)op[6] * 221425258653370777UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 221425258653370777UL) + ((uint64_t)op[1] * 2333700255738950699UL) + ((uint64_t)op[2] * 10811191460477622003UL) + ((uint64_t)op[3] * 18069544045833391061UL) + ((uint64_t)op[4] * 10174242020965357758UL) + ((((uint64_t)op[5] * 5825884965231322132UL) + ((uint64_t)op[6] * 7422492898454529783UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 7422492898454529783UL) + ((uint64_t)op[1] * 221425258653370777UL) + ((uint64_t)op[2] * 2333700255738950699UL) + ((uint64_t)op[3] * 10811191460477622003UL) + ((uint64_t)op[4] * 18069544045833391061UL) + ((uint64_t)op[5] * 10174242020965357758UL) + ((uint64_t)op[6] * 10682680752447059044UL);
	tmp_q[6] = ((uint64_t)op[0] * 5825884965231322132UL) + ((uint64_t)op[1] * 7422492898454529783UL) + ((uint64_t)op[2] * 221425258653370777UL) + ((uint64_t)op[3] * 2333700255738950699UL) + ((uint64_t)op[4] * 10811191460477622003UL) + ((uint64_t)op[5] * 18069544045833391061UL) + ((uint64_t)op[6] * 10174242020965357758UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 19359666618L) + ((((int128)tmp_q[1] * 21753474761L) - ((int128)tmp_q[2] * 45010561532L) + ((int128)tmp_q[3] * 75455209293L) + ((int128)tmp_q[4] * 19275066314L) - ((int128)tmp_q[5] * 29298756077L) - ((int128)tmp_q[6] * 18941078262L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 18941078262L) - ((int128)tmp_q[1] * 19359666618L) + ((((int128)tmp_q[2] * 21753474761L) - ((int128)tmp_q[3] * 45010561532L) + ((int128)tmp_q[4] * 75455209293L) + ((int128)tmp_q[5] * 19275066314L) - ((int128)tmp_q[6] * 29298756077L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 29298756077L) - ((int128)tmp_q[1] * 18941078262L) - ((int128)tmp_q[2] * 19359666618L) + ((((int128)tmp_q[3] * 21753474761L) - ((int128)tmp_q[4] * 45010561532L) + ((int128)tmp_q[5] * 75455209293L) + ((int128)tmp_q[6] * 19275066314L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 19275066314L) - ((int128)tmp_q[1] * 29298756077L) - ((int128)tmp_q[2] * 18941078262L) - ((int128)tmp_q[3] * 19359666618L) + ((((int128)tmp_q[4] * 21753474761L) - ((int128)tmp_q[5] * 45010561532L) + ((int128)tmp_q[6] * 75455209293L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 75455209293L) + ((int128)tmp_q[1] * 19275066314L) - ((int128)tmp_q[2] * 29298756077L) - ((int128)tmp_q[3] * 18941078262L) - ((int128)tmp_q[4] * 19359666618L) + ((((int128)tmp_q[5] * 21753474761L) - ((int128)tmp_q[6] * 45010561532L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 45010561532L) + ((int128)tmp_q[1] * 75455209293L) + ((int128)tmp_q[2] * 19275066314L) - ((int128)tmp_q[3] * 29298756077L) - ((int128)tmp_q[4] * 18941078262L) - ((int128)tmp_q[5] * 19359666618L) + ((int128)tmp_q[6] * 108767373805L);
	tmp_zero[6] = ((int128)tmp_q[0] * 21753474761L) - ((int128)tmp_q[1] * 45010561532L) + ((int128)tmp_q[2] * 75455209293L) + ((int128)tmp_q[3] * 19275066314L) - ((int128)tmp_q[4] * 29298756077L) - ((int128)tmp_q[5] * 18941078262L) - ((int128)tmp_q[6] * 19359666618L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

