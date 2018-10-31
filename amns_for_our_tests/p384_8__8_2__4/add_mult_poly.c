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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6947435083814511235UL) + ((((uint64_t)op[1] * 9770372864476978738UL) + ((uint64_t)op[2] * 14209761983692068137UL) + ((uint64_t)op[3] * 2119473630284002257UL) + ((uint64_t)op[4] * 15347398701892804427UL) + ((uint64_t)op[5] * 11762724499419598552UL) + ((uint64_t)op[6] * 15472632839002102486UL) + ((uint64_t)op[7] * 6854342985082272333UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 6854342985082272333UL) + ((uint64_t)op[1] * 6947435083814511235UL) + ((((uint64_t)op[2] * 9770372864476978738UL) + ((uint64_t)op[3] * 14209761983692068137UL) + ((uint64_t)op[4] * 2119473630284002257UL) + ((uint64_t)op[5] * 15347398701892804427UL) + ((uint64_t)op[6] * 11762724499419598552UL) + ((uint64_t)op[7] * 15472632839002102486UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 15472632839002102486UL) + ((uint64_t)op[1] * 6854342985082272333UL) + ((uint64_t)op[2] * 6947435083814511235UL) + ((((uint64_t)op[3] * 9770372864476978738UL) + ((uint64_t)op[4] * 14209761983692068137UL) + ((uint64_t)op[5] * 2119473630284002257UL) + ((uint64_t)op[6] * 15347398701892804427UL) + ((uint64_t)op[7] * 11762724499419598552UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 11762724499419598552UL) + ((uint64_t)op[1] * 15472632839002102486UL) + ((uint64_t)op[2] * 6854342985082272333UL) + ((uint64_t)op[3] * 6947435083814511235UL) + ((((uint64_t)op[4] * 9770372864476978738UL) + ((uint64_t)op[5] * 14209761983692068137UL) + ((uint64_t)op[6] * 2119473630284002257UL) + ((uint64_t)op[7] * 15347398701892804427UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 15347398701892804427UL) + ((uint64_t)op[1] * 11762724499419598552UL) + ((uint64_t)op[2] * 15472632839002102486UL) + ((uint64_t)op[3] * 6854342985082272333UL) + ((uint64_t)op[4] * 6947435083814511235UL) + ((((uint64_t)op[5] * 9770372864476978738UL) + ((uint64_t)op[6] * 14209761983692068137UL) + ((uint64_t)op[7] * 2119473630284002257UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 2119473630284002257UL) + ((uint64_t)op[1] * 15347398701892804427UL) + ((uint64_t)op[2] * 11762724499419598552UL) + ((uint64_t)op[3] * 15472632839002102486UL) + ((uint64_t)op[4] * 6854342985082272333UL) + ((uint64_t)op[5] * 6947435083814511235UL) + ((((uint64_t)op[6] * 9770372864476978738UL) + ((uint64_t)op[7] * 14209761983692068137UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 14209761983692068137UL) + ((uint64_t)op[1] * 2119473630284002257UL) + ((uint64_t)op[2] * 15347398701892804427UL) + ((uint64_t)op[3] * 11762724499419598552UL) + ((uint64_t)op[4] * 15472632839002102486UL) + ((uint64_t)op[5] * 6854342985082272333UL) + ((uint64_t)op[6] * 6947435083814511235UL) + ((uint64_t)op[7] * 1094001655244405860UL);
	tmp_q[7] = ((uint64_t)op[0] * 9770372864476978738UL) + ((uint64_t)op[1] * 14209761983692068137UL) + ((uint64_t)op[2] * 2119473630284002257UL) + ((uint64_t)op[3] * 15347398701892804427UL) + ((uint64_t)op[4] * 11762724499419598552UL) + ((uint64_t)op[5] * 15472632839002102486UL) + ((uint64_t)op[6] * 6854342985082272333UL) + ((uint64_t)op[7] * 6947435083814511235UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 70530311598747L) + ((((int128)tmp_q[1] * 27718890850608L) + ((int128)tmp_q[2] * 21718546795343L) + ((int128)tmp_q[3] * 45839101977218L) - ((int128)tmp_q[4] * 56210694289854L) - ((int128)tmp_q[5] * 164661027344785L) - ((int128)tmp_q[6] * 108292518342957L) + ((int128)tmp_q[7] * 66897507096887L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 66897507096887L) - ((int128)tmp_q[1] * 70530311598747L) + ((((int128)tmp_q[2] * 27718890850608L) + ((int128)tmp_q[3] * 21718546795343L) + ((int128)tmp_q[4] * 45839101977218L) - ((int128)tmp_q[5] * 56210694289854L) - ((int128)tmp_q[6] * 164661027344785L) - ((int128)tmp_q[7] * 108292518342957L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 108292518342957L) + ((int128)tmp_q[1] * 66897507096887L) - ((int128)tmp_q[2] * 70530311598747L) + ((((int128)tmp_q[3] * 27718890850608L) + ((int128)tmp_q[4] * 21718546795343L) + ((int128)tmp_q[5] * 45839101977218L) - ((int128)tmp_q[6] * 56210694289854L) - ((int128)tmp_q[7] * 164661027344785L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 164661027344785L) - ((int128)tmp_q[1] * 108292518342957L) + ((int128)tmp_q[2] * 66897507096887L) - ((int128)tmp_q[3] * 70530311598747L) + ((((int128)tmp_q[4] * 27718890850608L) + ((int128)tmp_q[5] * 21718546795343L) + ((int128)tmp_q[6] * 45839101977218L) - ((int128)tmp_q[7] * 56210694289854L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 56210694289854L) - ((int128)tmp_q[1] * 164661027344785L) - ((int128)tmp_q[2] * 108292518342957L) + ((int128)tmp_q[3] * 66897507096887L) - ((int128)tmp_q[4] * 70530311598747L) + ((((int128)tmp_q[5] * 27718890850608L) + ((int128)tmp_q[6] * 21718546795343L) + ((int128)tmp_q[7] * 45839101977218L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 45839101977218L) - ((int128)tmp_q[1] * 56210694289854L) - ((int128)tmp_q[2] * 164661027344785L) - ((int128)tmp_q[3] * 108292518342957L) + ((int128)tmp_q[4] * 66897507096887L) - ((int128)tmp_q[5] * 70530311598747L) + ((((int128)tmp_q[6] * 27718890850608L) + ((int128)tmp_q[7] * 21718546795343L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 21718546795343L) + ((int128)tmp_q[1] * 45839101977218L) - ((int128)tmp_q[2] * 56210694289854L) - ((int128)tmp_q[3] * 164661027344785L) - ((int128)tmp_q[4] * 108292518342957L) + ((int128)tmp_q[5] * 66897507096887L) - ((int128)tmp_q[6] * 70530311598747L) + ((int128)tmp_q[7] * 55437781701216L);
	tmp_zero[7] = ((int128)tmp_q[0] * 27718890850608L) + ((int128)tmp_q[1] * 21718546795343L) + ((int128)tmp_q[2] * 45839101977218L) - ((int128)tmp_q[3] * 56210694289854L) - ((int128)tmp_q[4] * 164661027344785L) - ((int128)tmp_q[5] * 108292518342957L) + ((int128)tmp_q[6] * 66897507096887L) - ((int128)tmp_q[7] * 70530311598747L);

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

