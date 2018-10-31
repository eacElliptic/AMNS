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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 10);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15340048683174975360UL) + ((((uint64_t)op[1] * 89172972158469345UL) + ((uint64_t)op[2] * 15628113062811420988UL) + ((uint64_t)op[3] * 15463093778958002028UL) + ((uint64_t)op[4] * 12560055974760274288UL) + ((uint64_t)op[5] * 16574904321224587880UL) + ((uint64_t)op[6] * 11857916992991493751UL) + ((uint64_t)op[7] * 9751189577001604121UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 9751189577001604121UL) + ((uint64_t)op[1] * 15340048683174975360UL) + ((((uint64_t)op[2] * 89172972158469345UL) + ((uint64_t)op[3] * 15628113062811420988UL) + ((uint64_t)op[4] * 15463093778958002028UL) + ((uint64_t)op[5] * 12560055974760274288UL) + ((uint64_t)op[6] * 16574904321224587880UL) + ((uint64_t)op[7] * 11857916992991493751UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 11857916992991493751UL) + ((uint64_t)op[1] * 9751189577001604121UL) + ((uint64_t)op[2] * 15340048683174975360UL) + ((((uint64_t)op[3] * 89172972158469345UL) + ((uint64_t)op[4] * 15628113062811420988UL) + ((uint64_t)op[5] * 15463093778958002028UL) + ((uint64_t)op[6] * 12560055974760274288UL) + ((uint64_t)op[7] * 16574904321224587880UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 16574904321224587880UL) + ((uint64_t)op[1] * 11857916992991493751UL) + ((uint64_t)op[2] * 9751189577001604121UL) + ((uint64_t)op[3] * 15340048683174975360UL) + ((((uint64_t)op[4] * 89172972158469345UL) + ((uint64_t)op[5] * 15628113062811420988UL) + ((uint64_t)op[6] * 15463093778958002028UL) + ((uint64_t)op[7] * 12560055974760274288UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 12560055974760274288UL) + ((uint64_t)op[1] * 16574904321224587880UL) + ((uint64_t)op[2] * 11857916992991493751UL) + ((uint64_t)op[3] * 9751189577001604121UL) + ((uint64_t)op[4] * 15340048683174975360UL) + ((((uint64_t)op[5] * 89172972158469345UL) + ((uint64_t)op[6] * 15628113062811420988UL) + ((uint64_t)op[7] * 15463093778958002028UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 15463093778958002028UL) + ((uint64_t)op[1] * 12560055974760274288UL) + ((uint64_t)op[2] * 16574904321224587880UL) + ((uint64_t)op[3] * 11857916992991493751UL) + ((uint64_t)op[4] * 9751189577001604121UL) + ((uint64_t)op[5] * 15340048683174975360UL) + ((((uint64_t)op[6] * 89172972158469345UL) + ((uint64_t)op[7] * 15628113062811420988UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 15628113062811420988UL) + ((uint64_t)op[1] * 15463093778958002028UL) + ((uint64_t)op[2] * 12560055974760274288UL) + ((uint64_t)op[3] * 16574904321224587880UL) + ((uint64_t)op[4] * 11857916992991493751UL) + ((uint64_t)op[5] * 9751189577001604121UL) + ((uint64_t)op[6] * 15340048683174975360UL) + ((uint64_t)op[7] * 445864860792346725UL);
	tmp_q[7] = ((uint64_t)op[0] * 89172972158469345UL) + ((uint64_t)op[1] * 15628113062811420988UL) + ((uint64_t)op[2] * 15463093778958002028UL) + ((uint64_t)op[3] * 12560055974760274288UL) + ((uint64_t)op[4] * 16574904321224587880UL) + ((uint64_t)op[5] * 11857916992991493751UL) + ((uint64_t)op[6] * 9751189577001604121UL) + ((uint64_t)op[7] * 15340048683174975360UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 107479050683953L) + ((-((int128)tmp_q[1] * 118879013713447L) + ((int128)tmp_q[2] * 94407951889327L) - ((int128)tmp_q[3] * 151002856059284L) + ((int128)tmp_q[4] * 10409329727817L) - ((int128)tmp_q[5] * 18673325862860L) - ((int128)tmp_q[6] * 187198495377832L) + ((int128)tmp_q[7] * 142402613664745L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 142402613664745L) + ((int128)tmp_q[1] * 107479050683953L) + ((-((int128)tmp_q[2] * 118879013713447L) + ((int128)tmp_q[3] * 94407951889327L) - ((int128)tmp_q[4] * 151002856059284L) + ((int128)tmp_q[5] * 10409329727817L) - ((int128)tmp_q[6] * 18673325862860L) - ((int128)tmp_q[7] * 187198495377832L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 187198495377832L) + ((int128)tmp_q[1] * 142402613664745L) + ((int128)tmp_q[2] * 107479050683953L) + ((-((int128)tmp_q[3] * 118879013713447L) + ((int128)tmp_q[4] * 94407951889327L) - ((int128)tmp_q[5] * 151002856059284L) + ((int128)tmp_q[6] * 10409329727817L) - ((int128)tmp_q[7] * 18673325862860L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 18673325862860L) - ((int128)tmp_q[1] * 187198495377832L) + ((int128)tmp_q[2] * 142402613664745L) + ((int128)tmp_q[3] * 107479050683953L) + ((-((int128)tmp_q[4] * 118879013713447L) + ((int128)tmp_q[5] * 94407951889327L) - ((int128)tmp_q[6] * 151002856059284L) + ((int128)tmp_q[7] * 10409329727817L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 10409329727817L) - ((int128)tmp_q[1] * 18673325862860L) - ((int128)tmp_q[2] * 187198495377832L) + ((int128)tmp_q[3] * 142402613664745L) + ((int128)tmp_q[4] * 107479050683953L) + ((-((int128)tmp_q[5] * 118879013713447L) + ((int128)tmp_q[6] * 94407951889327L) - ((int128)tmp_q[7] * 151002856059284L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 151002856059284L) + ((int128)tmp_q[1] * 10409329727817L) - ((int128)tmp_q[2] * 18673325862860L) - ((int128)tmp_q[3] * 187198495377832L) + ((int128)tmp_q[4] * 142402613664745L) + ((int128)tmp_q[5] * 107479050683953L) + ((-((int128)tmp_q[6] * 118879013713447L) + ((int128)tmp_q[7] * 94407951889327L)) * 5);
	tmp_zero[6] = ((int128)tmp_q[0] * 94407951889327L) - ((int128)tmp_q[1] * 151002856059284L) + ((int128)tmp_q[2] * 10409329727817L) - ((int128)tmp_q[3] * 18673325862860L) - ((int128)tmp_q[4] * 187198495377832L) + ((int128)tmp_q[5] * 142402613664745L) + ((int128)tmp_q[6] * 107479050683953L) - ((int128)tmp_q[7] * 594395068567235L);
	tmp_zero[7] = -((int128)tmp_q[0] * 118879013713447L) + ((int128)tmp_q[1] * 94407951889327L) - ((int128)tmp_q[2] * 151002856059284L) + ((int128)tmp_q[3] * 10409329727817L) - ((int128)tmp_q[4] * 18673325862860L) - ((int128)tmp_q[5] * 187198495377832L) + ((int128)tmp_q[6] * 142402613664745L) + ((int128)tmp_q[7] * 107479050683953L);

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

