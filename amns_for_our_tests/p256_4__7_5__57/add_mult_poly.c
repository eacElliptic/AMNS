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
	tmp_q[0] = ((uint64_t)op[0] * 17208272469417034426UL) + ((((uint64_t)op[1] * 9183635126891911600UL) + ((uint64_t)op[2] * 5983315204410911834UL) + ((uint64_t)op[3] * 8490702294073102620UL) + ((uint64_t)op[4] * 15481157229426120147UL) + ((uint64_t)op[5] * 15099814262289494593UL) + ((uint64_t)op[6] * 965451799795375081UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 965451799795375081UL) + ((uint64_t)op[1] * 17208272469417034426UL) + ((((uint64_t)op[2] * 9183635126891911600UL) + ((uint64_t)op[3] * 5983315204410911834UL) + ((uint64_t)op[4] * 8490702294073102620UL) + ((uint64_t)op[5] * 15481157229426120147UL) + ((uint64_t)op[6] * 15099814262289494593UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 15099814262289494593UL) + ((uint64_t)op[1] * 965451799795375081UL) + ((uint64_t)op[2] * 17208272469417034426UL) + ((((uint64_t)op[3] * 9183635126891911600UL) + ((uint64_t)op[4] * 5983315204410911834UL) + ((uint64_t)op[5] * 8490702294073102620UL) + ((uint64_t)op[6] * 15481157229426120147UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 15481157229426120147UL) + ((uint64_t)op[1] * 15099814262289494593UL) + ((uint64_t)op[2] * 965451799795375081UL) + ((uint64_t)op[3] * 17208272469417034426UL) + ((((uint64_t)op[4] * 9183635126891911600UL) + ((uint64_t)op[5] * 5983315204410911834UL) + ((uint64_t)op[6] * 8490702294073102620UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 8490702294073102620UL) + ((uint64_t)op[1] * 15481157229426120147UL) + ((uint64_t)op[2] * 15099814262289494593UL) + ((uint64_t)op[3] * 965451799795375081UL) + ((uint64_t)op[4] * 17208272469417034426UL) + ((((uint64_t)op[5] * 9183635126891911600UL) + ((uint64_t)op[6] * 5983315204410911834UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 5983315204410911834UL) + ((uint64_t)op[1] * 8490702294073102620UL) + ((uint64_t)op[2] * 15481157229426120147UL) + ((uint64_t)op[3] * 15099814262289494593UL) + ((uint64_t)op[4] * 965451799795375081UL) + ((uint64_t)op[5] * 17208272469417034426UL) + ((uint64_t)op[6] * 9024687487040454768UL);
	tmp_q[6] = ((uint64_t)op[0] * 9183635126891911600UL) + ((uint64_t)op[1] * 5983315204410911834UL) + ((uint64_t)op[2] * 8490702294073102620UL) + ((uint64_t)op[3] * 15481157229426120147UL) + ((uint64_t)op[4] * 15099814262289494593UL) + ((uint64_t)op[5] * 965451799795375081UL) + ((uint64_t)op[6] * 17208272469417034426UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 59518262506L) + ((-((int128)tmp_q[1] * 9948010311L) + ((int128)tmp_q[2] * 5479368277L) - ((int128)tmp_q[3] * 17360264019L) - ((int128)tmp_q[4] * 447104450L) + ((int128)tmp_q[5] * 82423771201L) - ((int128)tmp_q[6] * 35317652761L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 35317652761L) + ((int128)tmp_q[1] * 59518262506L) + ((-((int128)tmp_q[2] * 9948010311L) + ((int128)tmp_q[3] * 5479368277L) - ((int128)tmp_q[4] * 17360264019L) - ((int128)tmp_q[5] * 447104450L) + ((int128)tmp_q[6] * 82423771201L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 82423771201L) - ((int128)tmp_q[1] * 35317652761L) + ((int128)tmp_q[2] * 59518262506L) + ((-((int128)tmp_q[3] * 9948010311L) + ((int128)tmp_q[4] * 5479368277L) - ((int128)tmp_q[5] * 17360264019L) - ((int128)tmp_q[6] * 447104450L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 447104450L) + ((int128)tmp_q[1] * 82423771201L) - ((int128)tmp_q[2] * 35317652761L) + ((int128)tmp_q[3] * 59518262506L) + ((-((int128)tmp_q[4] * 9948010311L) + ((int128)tmp_q[5] * 5479368277L) - ((int128)tmp_q[6] * 17360264019L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 17360264019L) - ((int128)tmp_q[1] * 447104450L) + ((int128)tmp_q[2] * 82423771201L) - ((int128)tmp_q[3] * 35317652761L) + ((int128)tmp_q[4] * 59518262506L) + ((-((int128)tmp_q[5] * 9948010311L) + ((int128)tmp_q[6] * 5479368277L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 5479368277L) - ((int128)tmp_q[1] * 17360264019L) - ((int128)tmp_q[2] * 447104450L) + ((int128)tmp_q[3] * 82423771201L) - ((int128)tmp_q[4] * 35317652761L) + ((int128)tmp_q[5] * 59518262506L) - ((int128)tmp_q[6] * 49740051555L);
	tmp_zero[6] = -((int128)tmp_q[0] * 9948010311L) + ((int128)tmp_q[1] * 5479368277L) - ((int128)tmp_q[2] * 17360264019L) - ((int128)tmp_q[3] * 447104450L) + ((int128)tmp_q[4] * 82423771201L) - ((int128)tmp_q[5] * 35317652761L) + ((int128)tmp_q[6] * 59518262506L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

