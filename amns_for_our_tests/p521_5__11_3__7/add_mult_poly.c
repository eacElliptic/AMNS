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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 6);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12178686980926803898UL) + ((((uint64_t)op[1] * 3921892155269741073UL) + ((uint64_t)op[2] * 4165074587242253640UL) + ((uint64_t)op[3] * 2420979191559364243UL) + ((uint64_t)op[4] * 10380532057912079867UL) + ((uint64_t)op[5] * 6518110624721712258UL) + ((uint64_t)op[6] * 12255811099067596375UL) + ((uint64_t)op[7] * 753560797783495508UL) + ((uint64_t)op[8] * 14499649081365227296UL) + ((uint64_t)op[9] * 1118774341385841315UL) + ((uint64_t)op[10] * 12382886609617983744UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 12382886609617983744UL) + ((uint64_t)op[1] * 12178686980926803898UL) + ((((uint64_t)op[2] * 3921892155269741073UL) + ((uint64_t)op[3] * 4165074587242253640UL) + ((uint64_t)op[4] * 2420979191559364243UL) + ((uint64_t)op[5] * 10380532057912079867UL) + ((uint64_t)op[6] * 6518110624721712258UL) + ((uint64_t)op[7] * 12255811099067596375UL) + ((uint64_t)op[8] * 753560797783495508UL) + ((uint64_t)op[9] * 14499649081365227296UL) + ((uint64_t)op[10] * 1118774341385841315UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 1118774341385841315UL) + ((uint64_t)op[1] * 12382886609617983744UL) + ((uint64_t)op[2] * 12178686980926803898UL) + ((((uint64_t)op[3] * 3921892155269741073UL) + ((uint64_t)op[4] * 4165074587242253640UL) + ((uint64_t)op[5] * 2420979191559364243UL) + ((uint64_t)op[6] * 10380532057912079867UL) + ((uint64_t)op[7] * 6518110624721712258UL) + ((uint64_t)op[8] * 12255811099067596375UL) + ((uint64_t)op[9] * 753560797783495508UL) + ((uint64_t)op[10] * 14499649081365227296UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 14499649081365227296UL) + ((uint64_t)op[1] * 1118774341385841315UL) + ((uint64_t)op[2] * 12382886609617983744UL) + ((uint64_t)op[3] * 12178686980926803898UL) + ((((uint64_t)op[4] * 3921892155269741073UL) + ((uint64_t)op[5] * 4165074587242253640UL) + ((uint64_t)op[6] * 2420979191559364243UL) + ((uint64_t)op[7] * 10380532057912079867UL) + ((uint64_t)op[8] * 6518110624721712258UL) + ((uint64_t)op[9] * 12255811099067596375UL) + ((uint64_t)op[10] * 753560797783495508UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 753560797783495508UL) + ((uint64_t)op[1] * 14499649081365227296UL) + ((uint64_t)op[2] * 1118774341385841315UL) + ((uint64_t)op[3] * 12382886609617983744UL) + ((uint64_t)op[4] * 12178686980926803898UL) + ((((uint64_t)op[5] * 3921892155269741073UL) + ((uint64_t)op[6] * 4165074587242253640UL) + ((uint64_t)op[7] * 2420979191559364243UL) + ((uint64_t)op[8] * 10380532057912079867UL) + ((uint64_t)op[9] * 6518110624721712258UL) + ((uint64_t)op[10] * 12255811099067596375UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 12255811099067596375UL) + ((uint64_t)op[1] * 753560797783495508UL) + ((uint64_t)op[2] * 14499649081365227296UL) + ((uint64_t)op[3] * 1118774341385841315UL) + ((uint64_t)op[4] * 12382886609617983744UL) + ((uint64_t)op[5] * 12178686980926803898UL) + ((((uint64_t)op[6] * 3921892155269741073UL) + ((uint64_t)op[7] * 4165074587242253640UL) + ((uint64_t)op[8] * 2420979191559364243UL) + ((uint64_t)op[9] * 10380532057912079867UL) + ((uint64_t)op[10] * 6518110624721712258UL)) * 3);
	tmp_q[6] = ((uint64_t)op[0] * 6518110624721712258UL) + ((uint64_t)op[1] * 12255811099067596375UL) + ((uint64_t)op[2] * 753560797783495508UL) + ((uint64_t)op[3] * 14499649081365227296UL) + ((uint64_t)op[4] * 1118774341385841315UL) + ((uint64_t)op[5] * 12382886609617983744UL) + ((uint64_t)op[6] * 12178686980926803898UL) + ((((uint64_t)op[7] * 3921892155269741073UL) + ((uint64_t)op[8] * 4165074587242253640UL) + ((uint64_t)op[9] * 2420979191559364243UL) + ((uint64_t)op[10] * 10380532057912079867UL)) * 3);
	tmp_q[7] = ((uint64_t)op[0] * 10380532057912079867UL) + ((uint64_t)op[1] * 6518110624721712258UL) + ((uint64_t)op[2] * 12255811099067596375UL) + ((uint64_t)op[3] * 753560797783495508UL) + ((uint64_t)op[4] * 14499649081365227296UL) + ((uint64_t)op[5] * 1118774341385841315UL) + ((uint64_t)op[6] * 12382886609617983744UL) + ((uint64_t)op[7] * 12178686980926803898UL) + ((((uint64_t)op[8] * 3921892155269741073UL) + ((uint64_t)op[9] * 4165074587242253640UL) + ((uint64_t)op[10] * 2420979191559364243UL)) * 3);
	tmp_q[8] = ((uint64_t)op[0] * 2420979191559364243UL) + ((uint64_t)op[1] * 10380532057912079867UL) + ((uint64_t)op[2] * 6518110624721712258UL) + ((uint64_t)op[3] * 12255811099067596375UL) + ((uint64_t)op[4] * 753560797783495508UL) + ((uint64_t)op[5] * 14499649081365227296UL) + ((uint64_t)op[6] * 1118774341385841315UL) + ((uint64_t)op[7] * 12382886609617983744UL) + ((uint64_t)op[8] * 12178686980926803898UL) + ((((uint64_t)op[9] * 3921892155269741073UL) + ((uint64_t)op[10] * 4165074587242253640UL)) * 3);
	tmp_q[9] = ((uint64_t)op[0] * 4165074587242253640UL) + ((uint64_t)op[1] * 2420979191559364243UL) + ((uint64_t)op[2] * 10380532057912079867UL) + ((uint64_t)op[3] * 6518110624721712258UL) + ((uint64_t)op[4] * 12255811099067596375UL) + ((uint64_t)op[5] * 753560797783495508UL) + ((uint64_t)op[6] * 14499649081365227296UL) + ((uint64_t)op[7] * 1118774341385841315UL) + ((uint64_t)op[8] * 12382886609617983744UL) + ((uint64_t)op[9] * 12178686980926803898UL) + ((uint64_t)op[10] * 11765676465809223219UL);
	tmp_q[10] = ((uint64_t)op[0] * 3921892155269741073UL) + ((uint64_t)op[1] * 4165074587242253640UL) + ((uint64_t)op[2] * 2420979191559364243UL) + ((uint64_t)op[3] * 10380532057912079867UL) + ((uint64_t)op[4] * 6518110624721712258UL) + ((uint64_t)op[5] * 12255811099067596375UL) + ((uint64_t)op[6] * 753560797783495508UL) + ((uint64_t)op[7] * 14499649081365227296UL) + ((uint64_t)op[8] * 1118774341385841315UL) + ((uint64_t)op[9] * 12382886609617983744UL) + ((uint64_t)op[10] * 12178686980926803898UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 80500964039420L) + ((((int128)tmp_q[1] * 28276721435743L) - ((int128)tmp_q[2] * 50740828703245L) + ((int128)tmp_q[3] * 39568225254921L) + ((int128)tmp_q[4] * 110886638402004L) - ((int128)tmp_q[5] * 28790696800883L) - ((int128)tmp_q[6] * 76580064459937L) - ((int128)tmp_q[7] * 29293904003929L) + ((int128)tmp_q[8] * 36117511063631L) + ((int128)tmp_q[9] * 45798961747593L) + ((int128)tmp_q[10] * 2477982111881L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 2477982111881L) + ((int128)tmp_q[1] * 80500964039420L) + ((((int128)tmp_q[2] * 28276721435743L) - ((int128)tmp_q[3] * 50740828703245L) + ((int128)tmp_q[4] * 39568225254921L) + ((int128)tmp_q[5] * 110886638402004L) - ((int128)tmp_q[6] * 28790696800883L) - ((int128)tmp_q[7] * 76580064459937L) - ((int128)tmp_q[8] * 29293904003929L) + ((int128)tmp_q[9] * 36117511063631L) + ((int128)tmp_q[10] * 45798961747593L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 45798961747593L) + ((int128)tmp_q[1] * 2477982111881L) + ((int128)tmp_q[2] * 80500964039420L) + ((((int128)tmp_q[3] * 28276721435743L) - ((int128)tmp_q[4] * 50740828703245L) + ((int128)tmp_q[5] * 39568225254921L) + ((int128)tmp_q[6] * 110886638402004L) - ((int128)tmp_q[7] * 28790696800883L) - ((int128)tmp_q[8] * 76580064459937L) - ((int128)tmp_q[9] * 29293904003929L) + ((int128)tmp_q[10] * 36117511063631L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 36117511063631L) + ((int128)tmp_q[1] * 45798961747593L) + ((int128)tmp_q[2] * 2477982111881L) + ((int128)tmp_q[3] * 80500964039420L) + ((((int128)tmp_q[4] * 28276721435743L) - ((int128)tmp_q[5] * 50740828703245L) + ((int128)tmp_q[6] * 39568225254921L) + ((int128)tmp_q[7] * 110886638402004L) - ((int128)tmp_q[8] * 28790696800883L) - ((int128)tmp_q[9] * 76580064459937L) - ((int128)tmp_q[10] * 29293904003929L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 29293904003929L) + ((int128)tmp_q[1] * 36117511063631L) + ((int128)tmp_q[2] * 45798961747593L) + ((int128)tmp_q[3] * 2477982111881L) + ((int128)tmp_q[4] * 80500964039420L) + ((((int128)tmp_q[5] * 28276721435743L) - ((int128)tmp_q[6] * 50740828703245L) + ((int128)tmp_q[7] * 39568225254921L) + ((int128)tmp_q[8] * 110886638402004L) - ((int128)tmp_q[9] * 28790696800883L) - ((int128)tmp_q[10] * 76580064459937L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 76580064459937L) - ((int128)tmp_q[1] * 29293904003929L) + ((int128)tmp_q[2] * 36117511063631L) + ((int128)tmp_q[3] * 45798961747593L) + ((int128)tmp_q[4] * 2477982111881L) + ((int128)tmp_q[5] * 80500964039420L) + ((((int128)tmp_q[6] * 28276721435743L) - ((int128)tmp_q[7] * 50740828703245L) + ((int128)tmp_q[8] * 39568225254921L) + ((int128)tmp_q[9] * 110886638402004L) - ((int128)tmp_q[10] * 28790696800883L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 28790696800883L) - ((int128)tmp_q[1] * 76580064459937L) - ((int128)tmp_q[2] * 29293904003929L) + ((int128)tmp_q[3] * 36117511063631L) + ((int128)tmp_q[4] * 45798961747593L) + ((int128)tmp_q[5] * 2477982111881L) + ((int128)tmp_q[6] * 80500964039420L) + ((((int128)tmp_q[7] * 28276721435743L) - ((int128)tmp_q[8] * 50740828703245L) + ((int128)tmp_q[9] * 39568225254921L) + ((int128)tmp_q[10] * 110886638402004L)) * 3);
	tmp_zero[7] = ((int128)tmp_q[0] * 110886638402004L) - ((int128)tmp_q[1] * 28790696800883L) - ((int128)tmp_q[2] * 76580064459937L) - ((int128)tmp_q[3] * 29293904003929L) + ((int128)tmp_q[4] * 36117511063631L) + ((int128)tmp_q[5] * 45798961747593L) + ((int128)tmp_q[6] * 2477982111881L) + ((int128)tmp_q[7] * 80500964039420L) + ((((int128)tmp_q[8] * 28276721435743L) - ((int128)tmp_q[9] * 50740828703245L) + ((int128)tmp_q[10] * 39568225254921L)) * 3);
	tmp_zero[8] = ((int128)tmp_q[0] * 39568225254921L) + ((int128)tmp_q[1] * 110886638402004L) - ((int128)tmp_q[2] * 28790696800883L) - ((int128)tmp_q[3] * 76580064459937L) - ((int128)tmp_q[4] * 29293904003929L) + ((int128)tmp_q[5] * 36117511063631L) + ((int128)tmp_q[6] * 45798961747593L) + ((int128)tmp_q[7] * 2477982111881L) + ((int128)tmp_q[8] * 80500964039420L) + ((((int128)tmp_q[9] * 28276721435743L) - ((int128)tmp_q[10] * 50740828703245L)) * 3);
	tmp_zero[9] = -((int128)tmp_q[0] * 50740828703245L) + ((int128)tmp_q[1] * 39568225254921L) + ((int128)tmp_q[2] * 110886638402004L) - ((int128)tmp_q[3] * 28790696800883L) - ((int128)tmp_q[4] * 76580064459937L) - ((int128)tmp_q[5] * 29293904003929L) + ((int128)tmp_q[6] * 36117511063631L) + ((int128)tmp_q[7] * 45798961747593L) + ((int128)tmp_q[8] * 2477982111881L) + ((int128)tmp_q[9] * 80500964039420L) + ((int128)tmp_q[10] * 84830164307229L);
	tmp_zero[10] = ((int128)tmp_q[0] * 28276721435743L) - ((int128)tmp_q[1] * 50740828703245L) + ((int128)tmp_q[2] * 39568225254921L) + ((int128)tmp_q[3] * 110886638402004L) - ((int128)tmp_q[4] * 28790696800883L) - ((int128)tmp_q[5] * 76580064459937L) - ((int128)tmp_q[6] * 29293904003929L) + ((int128)tmp_q[7] * 36117511063631L) + ((int128)tmp_q[8] * 45798961747593L) + ((int128)tmp_q[9] * 2477982111881L) + ((int128)tmp_q[10] * 80500964039420L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
}

