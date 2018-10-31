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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9651222573915829841UL) + ((((uint64_t)op[1] * 5753638542761860156UL) + ((uint64_t)op[2] * 10511923083017508152UL) + ((uint64_t)op[3] * 265428510690105081UL) + ((uint64_t)op[4] * 4253552797201822710UL) + ((uint64_t)op[5] * 9145369782445894570UL) + ((uint64_t)op[6] * 1731123949587645165UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 1731123949587645165UL) + ((uint64_t)op[1] * 9651222573915829841UL) + ((((uint64_t)op[2] * 5753638542761860156UL) + ((uint64_t)op[3] * 10511923083017508152UL) + ((uint64_t)op[4] * 265428510690105081UL) + ((uint64_t)op[5] * 4253552797201822710UL) + ((uint64_t)op[6] * 9145369782445894570UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 9145369782445894570UL) + ((uint64_t)op[1] * 1731123949587645165UL) + ((uint64_t)op[2] * 9651222573915829841UL) + ((((uint64_t)op[3] * 5753638542761860156UL) + ((uint64_t)op[4] * 10511923083017508152UL) + ((uint64_t)op[5] * 265428510690105081UL) + ((uint64_t)op[6] * 4253552797201822710UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 4253552797201822710UL) + ((uint64_t)op[1] * 9145369782445894570UL) + ((uint64_t)op[2] * 1731123949587645165UL) + ((uint64_t)op[3] * 9651222573915829841UL) + ((((uint64_t)op[4] * 5753638542761860156UL) + ((uint64_t)op[5] * 10511923083017508152UL) + ((uint64_t)op[6] * 265428510690105081UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 265428510690105081UL) + ((uint64_t)op[1] * 4253552797201822710UL) + ((uint64_t)op[2] * 9145369782445894570UL) + ((uint64_t)op[3] * 1731123949587645165UL) + ((uint64_t)op[4] * 9651222573915829841UL) + ((((uint64_t)op[5] * 5753638542761860156UL) + ((uint64_t)op[6] * 10511923083017508152UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 10511923083017508152UL) + ((uint64_t)op[1] * 265428510690105081UL) + ((uint64_t)op[2] * 4253552797201822710UL) + ((uint64_t)op[3] * 9145369782445894570UL) + ((uint64_t)op[4] * 1731123949587645165UL) + ((uint64_t)op[5] * 9651222573915829841UL) + ((uint64_t)op[6] * 17260915628285580468UL);
	tmp_q[6] = ((uint64_t)op[0] * 5753638542761860156UL) + ((uint64_t)op[1] * 10511923083017508152UL) + ((uint64_t)op[2] * 265428510690105081UL) + ((uint64_t)op[3] * 4253552797201822710UL) + ((uint64_t)op[4] * 9145369782445894570UL) + ((uint64_t)op[5] * 1731123949587645165UL) + ((uint64_t)op[6] * 9651222573915829841UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 44399671745L) + ((-((int128)tmp_q[1] * 28540065379L) - ((int128)tmp_q[2] * 46343196905L) + ((int128)tmp_q[3] * 36804653346L) - ((int128)tmp_q[4] * 12787622821L) + ((int128)tmp_q[5] * 36652568766L) + ((int128)tmp_q[6] * 11147535337L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 11147535337L) + ((int128)tmp_q[1] * 44399671745L) + ((-((int128)tmp_q[2] * 28540065379L) - ((int128)tmp_q[3] * 46343196905L) + ((int128)tmp_q[4] * 36804653346L) - ((int128)tmp_q[5] * 12787622821L) + ((int128)tmp_q[6] * 36652568766L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 36652568766L) + ((int128)tmp_q[1] * 11147535337L) + ((int128)tmp_q[2] * 44399671745L) + ((-((int128)tmp_q[3] * 28540065379L) - ((int128)tmp_q[4] * 46343196905L) + ((int128)tmp_q[5] * 36804653346L) - ((int128)tmp_q[6] * 12787622821L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 12787622821L) + ((int128)tmp_q[1] * 36652568766L) + ((int128)tmp_q[2] * 11147535337L) + ((int128)tmp_q[3] * 44399671745L) + ((-((int128)tmp_q[4] * 28540065379L) - ((int128)tmp_q[5] * 46343196905L) + ((int128)tmp_q[6] * 36804653346L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 36804653346L) - ((int128)tmp_q[1] * 12787622821L) + ((int128)tmp_q[2] * 36652568766L) + ((int128)tmp_q[3] * 11147535337L) + ((int128)tmp_q[4] * 44399671745L) + ((-((int128)tmp_q[5] * 28540065379L) - ((int128)tmp_q[6] * 46343196905L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 46343196905L) + ((int128)tmp_q[1] * 36804653346L) - ((int128)tmp_q[2] * 12787622821L) + ((int128)tmp_q[3] * 36652568766L) + ((int128)tmp_q[4] * 11147535337L) + ((int128)tmp_q[5] * 44399671745L) - ((int128)tmp_q[6] * 85620196137L);
	tmp_zero[6] = -((int128)tmp_q[0] * 28540065379L) - ((int128)tmp_q[1] * 46343196905L) + ((int128)tmp_q[2] * 36804653346L) - ((int128)tmp_q[3] * 12787622821L) + ((int128)tmp_q[4] * 36652568766L) + ((int128)tmp_q[5] * 11147535337L) + ((int128)tmp_q[6] * 44399671745L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

