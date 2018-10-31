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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 6);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 6);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 12);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11097655111185198091UL) + ((((uint64_t)op[1] * 5292568649962404132UL) + ((uint64_t)op[2] * 8918740459997593835UL) + ((uint64_t)op[3] * 11453451762639718088UL) + ((uint64_t)op[4] * 12441027704540564052UL) + ((uint64_t)op[5] * 7809816444702418696UL) + ((uint64_t)op[6] * 2558192755572601554UL) + ((uint64_t)op[7] * 9761403141522032877UL) + ((uint64_t)op[8] * 8970021237714896753UL) + ((uint64_t)op[9] * 16275354438909726196UL) + ((uint64_t)op[10] * 11734091413848547229UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 11734091413848547229UL) + ((uint64_t)op[1] * 11097655111185198091UL) + ((((uint64_t)op[2] * 5292568649962404132UL) + ((uint64_t)op[3] * 8918740459997593835UL) + ((uint64_t)op[4] * 11453451762639718088UL) + ((uint64_t)op[5] * 12441027704540564052UL) + ((uint64_t)op[6] * 7809816444702418696UL) + ((uint64_t)op[7] * 2558192755572601554UL) + ((uint64_t)op[8] * 9761403141522032877UL) + ((uint64_t)op[9] * 8970021237714896753UL) + ((uint64_t)op[10] * 16275354438909726196UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 16275354438909726196UL) + ((uint64_t)op[1] * 11734091413848547229UL) + ((uint64_t)op[2] * 11097655111185198091UL) + ((((uint64_t)op[3] * 5292568649962404132UL) + ((uint64_t)op[4] * 8918740459997593835UL) + ((uint64_t)op[5] * 11453451762639718088UL) + ((uint64_t)op[6] * 12441027704540564052UL) + ((uint64_t)op[7] * 7809816444702418696UL) + ((uint64_t)op[8] * 2558192755572601554UL) + ((uint64_t)op[9] * 9761403141522032877UL) + ((uint64_t)op[10] * 8970021237714896753UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 8970021237714896753UL) + ((uint64_t)op[1] * 16275354438909726196UL) + ((uint64_t)op[2] * 11734091413848547229UL) + ((uint64_t)op[3] * 11097655111185198091UL) + ((((uint64_t)op[4] * 5292568649962404132UL) + ((uint64_t)op[5] * 8918740459997593835UL) + ((uint64_t)op[6] * 11453451762639718088UL) + ((uint64_t)op[7] * 12441027704540564052UL) + ((uint64_t)op[8] * 7809816444702418696UL) + ((uint64_t)op[9] * 2558192755572601554UL) + ((uint64_t)op[10] * 9761403141522032877UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 9761403141522032877UL) + ((uint64_t)op[1] * 8970021237714896753UL) + ((uint64_t)op[2] * 16275354438909726196UL) + ((uint64_t)op[3] * 11734091413848547229UL) + ((uint64_t)op[4] * 11097655111185198091UL) + ((((uint64_t)op[5] * 5292568649962404132UL) + ((uint64_t)op[6] * 8918740459997593835UL) + ((uint64_t)op[7] * 11453451762639718088UL) + ((uint64_t)op[8] * 12441027704540564052UL) + ((uint64_t)op[9] * 7809816444702418696UL) + ((uint64_t)op[10] * 2558192755572601554UL)) * 6);
	tmp_q[5] = ((uint64_t)op[0] * 2558192755572601554UL) + ((uint64_t)op[1] * 9761403141522032877UL) + ((uint64_t)op[2] * 8970021237714896753UL) + ((uint64_t)op[3] * 16275354438909726196UL) + ((uint64_t)op[4] * 11734091413848547229UL) + ((uint64_t)op[5] * 11097655111185198091UL) + ((((uint64_t)op[6] * 5292568649962404132UL) + ((uint64_t)op[7] * 8918740459997593835UL) + ((uint64_t)op[8] * 11453451762639718088UL) + ((uint64_t)op[9] * 12441027704540564052UL) + ((uint64_t)op[10] * 7809816444702418696UL)) * 6);
	tmp_q[6] = ((uint64_t)op[0] * 7809816444702418696UL) + ((uint64_t)op[1] * 2558192755572601554UL) + ((uint64_t)op[2] * 9761403141522032877UL) + ((uint64_t)op[3] * 8970021237714896753UL) + ((uint64_t)op[4] * 16275354438909726196UL) + ((uint64_t)op[5] * 11734091413848547229UL) + ((uint64_t)op[6] * 11097655111185198091UL) + ((((uint64_t)op[7] * 5292568649962404132UL) + ((uint64_t)op[8] * 8918740459997593835UL) + ((uint64_t)op[9] * 11453451762639718088UL) + ((uint64_t)op[10] * 12441027704540564052UL)) * 6);
	tmp_q[7] = ((uint64_t)op[0] * 12441027704540564052UL) + ((uint64_t)op[1] * 7809816444702418696UL) + ((uint64_t)op[2] * 2558192755572601554UL) + ((uint64_t)op[3] * 9761403141522032877UL) + ((uint64_t)op[4] * 8970021237714896753UL) + ((uint64_t)op[5] * 16275354438909726196UL) + ((uint64_t)op[6] * 11734091413848547229UL) + ((uint64_t)op[7] * 11097655111185198091UL) + ((((uint64_t)op[8] * 5292568649962404132UL) + ((uint64_t)op[9] * 8918740459997593835UL) + ((uint64_t)op[10] * 11453451762639718088UL)) * 6);
	tmp_q[8] = ((uint64_t)op[0] * 11453451762639718088UL) + ((uint64_t)op[1] * 12441027704540564052UL) + ((uint64_t)op[2] * 7809816444702418696UL) + ((uint64_t)op[3] * 2558192755572601554UL) + ((uint64_t)op[4] * 9761403141522032877UL) + ((uint64_t)op[5] * 8970021237714896753UL) + ((uint64_t)op[6] * 16275354438909726196UL) + ((uint64_t)op[7] * 11734091413848547229UL) + ((uint64_t)op[8] * 11097655111185198091UL) + ((((uint64_t)op[9] * 5292568649962404132UL) + ((uint64_t)op[10] * 8918740459997593835UL)) * 6);
	tmp_q[9] = ((uint64_t)op[0] * 8918740459997593835UL) + ((uint64_t)op[1] * 11453451762639718088UL) + ((uint64_t)op[2] * 12441027704540564052UL) + ((uint64_t)op[3] * 7809816444702418696UL) + ((uint64_t)op[4] * 2558192755572601554UL) + ((uint64_t)op[5] * 9761403141522032877UL) + ((uint64_t)op[6] * 8970021237714896753UL) + ((uint64_t)op[7] * 16275354438909726196UL) + ((uint64_t)op[8] * 11734091413848547229UL) + ((uint64_t)op[9] * 11097655111185198091UL) + ((uint64_t)op[10] * 13308667826064873176UL);
	tmp_q[10] = ((uint64_t)op[0] * 5292568649962404132UL) + ((uint64_t)op[1] * 8918740459997593835UL) + ((uint64_t)op[2] * 11453451762639718088UL) + ((uint64_t)op[3] * 12441027704540564052UL) + ((uint64_t)op[4] * 7809816444702418696UL) + ((uint64_t)op[5] * 2558192755572601554UL) + ((uint64_t)op[6] * 9761403141522032877UL) + ((uint64_t)op[7] * 8970021237714896753UL) + ((uint64_t)op[8] * 16275354438909726196UL) + ((uint64_t)op[9] * 11734091413848547229UL) + ((uint64_t)op[10] * 11097655111185198091UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 45857299257587L) + ((-((int128)tmp_q[1] * 3067004662702L) + ((int128)tmp_q[2] * 24721429741363L) + ((int128)tmp_q[3] * 92547658633005L) + ((int128)tmp_q[4] * 26039476679081L) + ((int128)tmp_q[5] * 16834229309457L) - ((int128)tmp_q[6] * 36884383287444L) + ((int128)tmp_q[7] * 9378406501188L) - ((int128)tmp_q[8] * 55436987451312L) - ((int128)tmp_q[9] * 59629570328827L) - ((int128)tmp_q[10] * 42239518931065L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 42239518931065L) + ((int128)tmp_q[1] * 45857299257587L) + ((-((int128)tmp_q[2] * 3067004662702L) + ((int128)tmp_q[3] * 24721429741363L) + ((int128)tmp_q[4] * 92547658633005L) + ((int128)tmp_q[5] * 26039476679081L) + ((int128)tmp_q[6] * 16834229309457L) - ((int128)tmp_q[7] * 36884383287444L) + ((int128)tmp_q[8] * 9378406501188L) - ((int128)tmp_q[9] * 55436987451312L) - ((int128)tmp_q[10] * 59629570328827L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 59629570328827L) - ((int128)tmp_q[1] * 42239518931065L) + ((int128)tmp_q[2] * 45857299257587L) + ((-((int128)tmp_q[3] * 3067004662702L) + ((int128)tmp_q[4] * 24721429741363L) + ((int128)tmp_q[5] * 92547658633005L) + ((int128)tmp_q[6] * 26039476679081L) + ((int128)tmp_q[7] * 16834229309457L) - ((int128)tmp_q[8] * 36884383287444L) + ((int128)tmp_q[9] * 9378406501188L) - ((int128)tmp_q[10] * 55436987451312L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 55436987451312L) - ((int128)tmp_q[1] * 59629570328827L) - ((int128)tmp_q[2] * 42239518931065L) + ((int128)tmp_q[3] * 45857299257587L) + ((-((int128)tmp_q[4] * 3067004662702L) + ((int128)tmp_q[5] * 24721429741363L) + ((int128)tmp_q[6] * 92547658633005L) + ((int128)tmp_q[7] * 26039476679081L) + ((int128)tmp_q[8] * 16834229309457L) - ((int128)tmp_q[9] * 36884383287444L) + ((int128)tmp_q[10] * 9378406501188L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 9378406501188L) - ((int128)tmp_q[1] * 55436987451312L) - ((int128)tmp_q[2] * 59629570328827L) - ((int128)tmp_q[3] * 42239518931065L) + ((int128)tmp_q[4] * 45857299257587L) + ((-((int128)tmp_q[5] * 3067004662702L) + ((int128)tmp_q[6] * 24721429741363L) + ((int128)tmp_q[7] * 92547658633005L) + ((int128)tmp_q[8] * 26039476679081L) + ((int128)tmp_q[9] * 16834229309457L) - ((int128)tmp_q[10] * 36884383287444L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 36884383287444L) + ((int128)tmp_q[1] * 9378406501188L) - ((int128)tmp_q[2] * 55436987451312L) - ((int128)tmp_q[3] * 59629570328827L) - ((int128)tmp_q[4] * 42239518931065L) + ((int128)tmp_q[5] * 45857299257587L) + ((-((int128)tmp_q[6] * 3067004662702L) + ((int128)tmp_q[7] * 24721429741363L) + ((int128)tmp_q[8] * 92547658633005L) + ((int128)tmp_q[9] * 26039476679081L) + ((int128)tmp_q[10] * 16834229309457L)) * 6);
	tmp_zero[6] = ((int128)tmp_q[0] * 16834229309457L) - ((int128)tmp_q[1] * 36884383287444L) + ((int128)tmp_q[2] * 9378406501188L) - ((int128)tmp_q[3] * 55436987451312L) - ((int128)tmp_q[4] * 59629570328827L) - ((int128)tmp_q[5] * 42239518931065L) + ((int128)tmp_q[6] * 45857299257587L) + ((-((int128)tmp_q[7] * 3067004662702L) + ((int128)tmp_q[8] * 24721429741363L) + ((int128)tmp_q[9] * 92547658633005L) + ((int128)tmp_q[10] * 26039476679081L)) * 6);
	tmp_zero[7] = ((int128)tmp_q[0] * 26039476679081L) + ((int128)tmp_q[1] * 16834229309457L) - ((int128)tmp_q[2] * 36884383287444L) + ((int128)tmp_q[3] * 9378406501188L) - ((int128)tmp_q[4] * 55436987451312L) - ((int128)tmp_q[5] * 59629570328827L) - ((int128)tmp_q[6] * 42239518931065L) + ((int128)tmp_q[7] * 45857299257587L) + ((-((int128)tmp_q[8] * 3067004662702L) + ((int128)tmp_q[9] * 24721429741363L) + ((int128)tmp_q[10] * 92547658633005L)) * 6);
	tmp_zero[8] = ((int128)tmp_q[0] * 92547658633005L) + ((int128)tmp_q[1] * 26039476679081L) + ((int128)tmp_q[2] * 16834229309457L) - ((int128)tmp_q[3] * 36884383287444L) + ((int128)tmp_q[4] * 9378406501188L) - ((int128)tmp_q[5] * 55436987451312L) - ((int128)tmp_q[6] * 59629570328827L) - ((int128)tmp_q[7] * 42239518931065L) + ((int128)tmp_q[8] * 45857299257587L) + ((-((int128)tmp_q[9] * 3067004662702L) + ((int128)tmp_q[10] * 24721429741363L)) * 6);
	tmp_zero[9] = ((int128)tmp_q[0] * 24721429741363L) + ((int128)tmp_q[1] * 92547658633005L) + ((int128)tmp_q[2] * 26039476679081L) + ((int128)tmp_q[3] * 16834229309457L) - ((int128)tmp_q[4] * 36884383287444L) + ((int128)tmp_q[5] * 9378406501188L) - ((int128)tmp_q[6] * 55436987451312L) - ((int128)tmp_q[7] * 59629570328827L) - ((int128)tmp_q[8] * 42239518931065L) + ((int128)tmp_q[9] * 45857299257587L) - ((int128)tmp_q[10] * 18402027976212L);
	tmp_zero[10] = -((int128)tmp_q[0] * 3067004662702L) + ((int128)tmp_q[1] * 24721429741363L) + ((int128)tmp_q[2] * 92547658633005L) + ((int128)tmp_q[3] * 26039476679081L) + ((int128)tmp_q[4] * 16834229309457L) - ((int128)tmp_q[5] * 36884383287444L) + ((int128)tmp_q[6] * 9378406501188L) - ((int128)tmp_q[7] * 55436987451312L) - ((int128)tmp_q[8] * 59629570328827L) - ((int128)tmp_q[9] * 42239518931065L) + ((int128)tmp_q[10] * 45857299257587L);

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

