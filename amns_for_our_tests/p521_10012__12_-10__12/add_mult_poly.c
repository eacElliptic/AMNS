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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 10);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 10);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 10);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 10);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 10);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 10);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 10);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 10);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 10);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 10);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 10);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 20);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 20);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 20);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 20);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 10);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 20);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 10);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4571784853829693401UL) + ((((uint64_t)op[1] * 10018962020266657274UL) + ((uint64_t)op[2] * 5168529153368479498UL) + ((uint64_t)op[3] * 13444313339678805381UL) + ((uint64_t)op[4] * 16219044007440786299UL) + ((uint64_t)op[5] * 10789543072613395374UL) + ((uint64_t)op[6] * 4295185087568869185UL) + ((uint64_t)op[7] * 16565965787515404503UL) + ((uint64_t)op[8] * 5491159961322652503UL) + ((uint64_t)op[9] * 10016581184578952196UL) + ((uint64_t)op[10] * 14900748754645278934UL) + ((uint64_t)op[11] * 5785747452141223828UL)) * 18446744073709551606);
	tmp_q[1] = ((uint64_t)op[0] * 5785747452141223828UL) + ((uint64_t)op[1] * 4571784853829693401UL) + ((((uint64_t)op[2] * 10018962020266657274UL) + ((uint64_t)op[3] * 5168529153368479498UL) + ((uint64_t)op[4] * 13444313339678805381UL) + ((uint64_t)op[5] * 16219044007440786299UL) + ((uint64_t)op[6] * 10789543072613395374UL) + ((uint64_t)op[7] * 4295185087568869185UL) + ((uint64_t)op[8] * 16565965787515404503UL) + ((uint64_t)op[9] * 5491159961322652503UL) + ((uint64_t)op[10] * 10016581184578952196UL) + ((uint64_t)op[11] * 14900748754645278934UL)) * 18446744073709551606);
	tmp_q[2] = ((uint64_t)op[0] * 14900748754645278934UL) + ((uint64_t)op[1] * 5785747452141223828UL) + ((uint64_t)op[2] * 4571784853829693401UL) + ((((uint64_t)op[3] * 10018962020266657274UL) + ((uint64_t)op[4] * 5168529153368479498UL) + ((uint64_t)op[5] * 13444313339678805381UL) + ((uint64_t)op[6] * 16219044007440786299UL) + ((uint64_t)op[7] * 10789543072613395374UL) + ((uint64_t)op[8] * 4295185087568869185UL) + ((uint64_t)op[9] * 16565965787515404503UL) + ((uint64_t)op[10] * 5491159961322652503UL) + ((uint64_t)op[11] * 10016581184578952196UL)) * 18446744073709551606);
	tmp_q[3] = ((uint64_t)op[0] * 10016581184578952196UL) + ((uint64_t)op[1] * 14900748754645278934UL) + ((uint64_t)op[2] * 5785747452141223828UL) + ((uint64_t)op[3] * 4571784853829693401UL) + ((((uint64_t)op[4] * 10018962020266657274UL) + ((uint64_t)op[5] * 5168529153368479498UL) + ((uint64_t)op[6] * 13444313339678805381UL) + ((uint64_t)op[7] * 16219044007440786299UL) + ((uint64_t)op[8] * 10789543072613395374UL) + ((uint64_t)op[9] * 4295185087568869185UL) + ((uint64_t)op[10] * 16565965787515404503UL) + ((uint64_t)op[11] * 5491159961322652503UL)) * 18446744073709551606);
	tmp_q[4] = ((uint64_t)op[0] * 5491159961322652503UL) + ((uint64_t)op[1] * 10016581184578952196UL) + ((uint64_t)op[2] * 14900748754645278934UL) + ((uint64_t)op[3] * 5785747452141223828UL) + ((uint64_t)op[4] * 4571784853829693401UL) + ((((uint64_t)op[5] * 10018962020266657274UL) + ((uint64_t)op[6] * 5168529153368479498UL) + ((uint64_t)op[7] * 13444313339678805381UL) + ((uint64_t)op[8] * 16219044007440786299UL) + ((uint64_t)op[9] * 10789543072613395374UL) + ((uint64_t)op[10] * 4295185087568869185UL) + ((uint64_t)op[11] * 16565965787515404503UL)) * 18446744073709551606);
	tmp_q[5] = ((uint64_t)op[0] * 16565965787515404503UL) + ((uint64_t)op[1] * 5491159961322652503UL) + ((uint64_t)op[2] * 10016581184578952196UL) + ((uint64_t)op[3] * 14900748754645278934UL) + ((uint64_t)op[4] * 5785747452141223828UL) + ((uint64_t)op[5] * 4571784853829693401UL) + ((((uint64_t)op[6] * 10018962020266657274UL) + ((uint64_t)op[7] * 5168529153368479498UL) + ((uint64_t)op[8] * 13444313339678805381UL) + ((uint64_t)op[9] * 16219044007440786299UL) + ((uint64_t)op[10] * 10789543072613395374UL) + ((uint64_t)op[11] * 4295185087568869185UL)) * 18446744073709551606);
	tmp_q[6] = ((uint64_t)op[0] * 4295185087568869185UL) + ((uint64_t)op[1] * 16565965787515404503UL) + ((uint64_t)op[2] * 5491159961322652503UL) + ((uint64_t)op[3] * 10016581184578952196UL) + ((uint64_t)op[4] * 14900748754645278934UL) + ((uint64_t)op[5] * 5785747452141223828UL) + ((uint64_t)op[6] * 4571784853829693401UL) + ((((uint64_t)op[7] * 10018962020266657274UL) + ((uint64_t)op[8] * 5168529153368479498UL) + ((uint64_t)op[9] * 13444313339678805381UL) + ((uint64_t)op[10] * 16219044007440786299UL) + ((uint64_t)op[11] * 10789543072613395374UL)) * 18446744073709551606);
	tmp_q[7] = ((uint64_t)op[0] * 10789543072613395374UL) + ((uint64_t)op[1] * 4295185087568869185UL) + ((uint64_t)op[2] * 16565965787515404503UL) + ((uint64_t)op[3] * 5491159961322652503UL) + ((uint64_t)op[4] * 10016581184578952196UL) + ((uint64_t)op[5] * 14900748754645278934UL) + ((uint64_t)op[6] * 5785747452141223828UL) + ((uint64_t)op[7] * 4571784853829693401UL) + ((((uint64_t)op[8] * 10018962020266657274UL) + ((uint64_t)op[9] * 5168529153368479498UL) + ((uint64_t)op[10] * 13444313339678805381UL) + ((uint64_t)op[11] * 16219044007440786299UL)) * 18446744073709551606);
	tmp_q[8] = ((uint64_t)op[0] * 16219044007440786299UL) + ((uint64_t)op[1] * 10789543072613395374UL) + ((uint64_t)op[2] * 4295185087568869185UL) + ((uint64_t)op[3] * 16565965787515404503UL) + ((uint64_t)op[4] * 5491159961322652503UL) + ((uint64_t)op[5] * 10016581184578952196UL) + ((uint64_t)op[6] * 14900748754645278934UL) + ((uint64_t)op[7] * 5785747452141223828UL) + ((uint64_t)op[8] * 4571784853829693401UL) + ((((uint64_t)op[9] * 10018962020266657274UL) + ((uint64_t)op[10] * 5168529153368479498UL) + ((uint64_t)op[11] * 13444313339678805381UL)) * 18446744073709551606);
	tmp_q[9] = ((uint64_t)op[0] * 13444313339678805381UL) + ((uint64_t)op[1] * 16219044007440786299UL) + ((uint64_t)op[2] * 10789543072613395374UL) + ((uint64_t)op[3] * 4295185087568869185UL) + ((uint64_t)op[4] * 16565965787515404503UL) + ((uint64_t)op[5] * 5491159961322652503UL) + ((uint64_t)op[6] * 10016581184578952196UL) + ((uint64_t)op[7] * 14900748754645278934UL) + ((uint64_t)op[8] * 5785747452141223828UL) + ((uint64_t)op[9] * 4571784853829693401UL) + ((((uint64_t)op[10] * 10018962020266657274UL) + ((uint64_t)op[11] * 5168529153368479498UL)) * 18446744073709551606);
	tmp_q[10] = ((uint64_t)op[0] * 5168529153368479498UL) + ((uint64_t)op[1] * 13444313339678805381UL) + ((uint64_t)op[2] * 16219044007440786299UL) + ((uint64_t)op[3] * 10789543072613395374UL) + ((uint64_t)op[4] * 4295185087568869185UL) + ((uint64_t)op[5] * 16565965787515404503UL) + ((uint64_t)op[6] * 5491159961322652503UL) + ((uint64_t)op[7] * 10016581184578952196UL) + ((uint64_t)op[8] * 14900748754645278934UL) + ((uint64_t)op[9] * 5785747452141223828UL) + ((uint64_t)op[10] * 4571784853829693401UL) + ((uint64_t)op[11] * 10490844239590736956UL);
	tmp_q[11] = ((uint64_t)op[0] * 10018962020266657274UL) + ((uint64_t)op[1] * 5168529153368479498UL) + ((uint64_t)op[2] * 13444313339678805381UL) + ((uint64_t)op[3] * 16219044007440786299UL) + ((uint64_t)op[4] * 10789543072613395374UL) + ((uint64_t)op[5] * 4295185087568869185UL) + ((uint64_t)op[6] * 16565965787515404503UL) + ((uint64_t)op[7] * 5491159961322652503UL) + ((uint64_t)op[8] * 10016581184578952196UL) + ((uint64_t)op[9] * 14900748754645278934UL) + ((uint64_t)op[10] * 5785747452141223828UL) + ((uint64_t)op[11] * 4571784853829693401UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 5868169780971L) - ((-((int128)tmp_q[1] * 3354125399230L) - ((int128)tmp_q[2] * 1162084541777L) - ((int128)tmp_q[3] * 3369654567083L) + ((int128)tmp_q[4] * 4141609118370L) - ((int128)tmp_q[5] * 5492226004432L) + ((int128)tmp_q[6] * 136330842539L) - ((int128)tmp_q[7] * 4139106621025L) + ((int128)tmp_q[8] * 40831642937L) - ((int128)tmp_q[9] * 34241695394L) - ((int128)tmp_q[10] * 15298643674L) - ((int128)tmp_q[11] * 7341466365674L)) * 10);
	tmp_zero[1] = -((int128)tmp_q[0] * 7341466365674L) + ((int128)tmp_q[1] * 5868169780971L) - ((-((int128)tmp_q[2] * 3354125399230L) - ((int128)tmp_q[3] * 1162084541777L) - ((int128)tmp_q[4] * 3369654567083L) + ((int128)tmp_q[5] * 4141609118370L) - ((int128)tmp_q[6] * 5492226004432L) + ((int128)tmp_q[7] * 136330842539L) - ((int128)tmp_q[8] * 4139106621025L) + ((int128)tmp_q[9] * 40831642937L) - ((int128)tmp_q[10] * 34241695394L) - ((int128)tmp_q[11] * 15298643674L)) * 10);
	tmp_zero[2] = -((int128)tmp_q[0] * 15298643674L) - ((int128)tmp_q[1] * 7341466365674L) + ((int128)tmp_q[2] * 5868169780971L) - ((-((int128)tmp_q[3] * 3354125399230L) - ((int128)tmp_q[4] * 1162084541777L) - ((int128)tmp_q[5] * 3369654567083L) + ((int128)tmp_q[6] * 4141609118370L) - ((int128)tmp_q[7] * 5492226004432L) + ((int128)tmp_q[8] * 136330842539L) - ((int128)tmp_q[9] * 4139106621025L) + ((int128)tmp_q[10] * 40831642937L) - ((int128)tmp_q[11] * 34241695394L)) * 10);
	tmp_zero[3] = -((int128)tmp_q[0] * 34241695394L) - ((int128)tmp_q[1] * 15298643674L) - ((int128)tmp_q[2] * 7341466365674L) + ((int128)tmp_q[3] * 5868169780971L) - ((-((int128)tmp_q[4] * 3354125399230L) - ((int128)tmp_q[5] * 1162084541777L) - ((int128)tmp_q[6] * 3369654567083L) + ((int128)tmp_q[7] * 4141609118370L) - ((int128)tmp_q[8] * 5492226004432L) + ((int128)tmp_q[9] * 136330842539L) - ((int128)tmp_q[10] * 4139106621025L) + ((int128)tmp_q[11] * 40831642937L)) * 10);
	tmp_zero[4] = ((int128)tmp_q[0] * 40831642937L) - ((int128)tmp_q[1] * 34241695394L) - ((int128)tmp_q[2] * 15298643674L) - ((int128)tmp_q[3] * 7341466365674L) + ((int128)tmp_q[4] * 5868169780971L) - ((-((int128)tmp_q[5] * 3354125399230L) - ((int128)tmp_q[6] * 1162084541777L) - ((int128)tmp_q[7] * 3369654567083L) + ((int128)tmp_q[8] * 4141609118370L) - ((int128)tmp_q[9] * 5492226004432L) + ((int128)tmp_q[10] * 136330842539L) - ((int128)tmp_q[11] * 4139106621025L)) * 10);
	tmp_zero[5] = -((int128)tmp_q[0] * 4139106621025L) + ((int128)tmp_q[1] * 40831642937L) - ((int128)tmp_q[2] * 34241695394L) - ((int128)tmp_q[3] * 15298643674L) - ((int128)tmp_q[4] * 7341466365674L) + ((int128)tmp_q[5] * 5868169780971L) - ((-((int128)tmp_q[6] * 3354125399230L) - ((int128)tmp_q[7] * 1162084541777L) - ((int128)tmp_q[8] * 3369654567083L) + ((int128)tmp_q[9] * 4141609118370L) - ((int128)tmp_q[10] * 5492226004432L) + ((int128)tmp_q[11] * 136330842539L)) * 10);
	tmp_zero[6] = ((int128)tmp_q[0] * 136330842539L) - ((int128)tmp_q[1] * 4139106621025L) + ((int128)tmp_q[2] * 40831642937L) - ((int128)tmp_q[3] * 34241695394L) - ((int128)tmp_q[4] * 15298643674L) - ((int128)tmp_q[5] * 7341466365674L) + ((int128)tmp_q[6] * 5868169780971L) - ((-((int128)tmp_q[7] * 3354125399230L) - ((int128)tmp_q[8] * 1162084541777L) - ((int128)tmp_q[9] * 3369654567083L) + ((int128)tmp_q[10] * 4141609118370L) - ((int128)tmp_q[11] * 5492226004432L)) * 10);
	tmp_zero[7] = -((int128)tmp_q[0] * 5492226004432L) + ((int128)tmp_q[1] * 136330842539L) - ((int128)tmp_q[2] * 4139106621025L) + ((int128)tmp_q[3] * 40831642937L) - ((int128)tmp_q[4] * 34241695394L) - ((int128)tmp_q[5] * 15298643674L) - ((int128)tmp_q[6] * 7341466365674L) + ((int128)tmp_q[7] * 5868169780971L) - ((-((int128)tmp_q[8] * 3354125399230L) - ((int128)tmp_q[9] * 1162084541777L) - ((int128)tmp_q[10] * 3369654567083L) + ((int128)tmp_q[11] * 4141609118370L)) * 10);
	tmp_zero[8] = ((int128)tmp_q[0] * 4141609118370L) - ((int128)tmp_q[1] * 5492226004432L) + ((int128)tmp_q[2] * 136330842539L) - ((int128)tmp_q[3] * 4139106621025L) + ((int128)tmp_q[4] * 40831642937L) - ((int128)tmp_q[5] * 34241695394L) - ((int128)tmp_q[6] * 15298643674L) - ((int128)tmp_q[7] * 7341466365674L) + ((int128)tmp_q[8] * 5868169780971L) - ((-((int128)tmp_q[9] * 3354125399230L) - ((int128)tmp_q[10] * 1162084541777L) - ((int128)tmp_q[11] * 3369654567083L)) * 10);
	tmp_zero[9] = -((int128)tmp_q[0] * 3369654567083L) + ((int128)tmp_q[1] * 4141609118370L) - ((int128)tmp_q[2] * 5492226004432L) + ((int128)tmp_q[3] * 136330842539L) - ((int128)tmp_q[4] * 4139106621025L) + ((int128)tmp_q[5] * 40831642937L) - ((int128)tmp_q[6] * 34241695394L) - ((int128)tmp_q[7] * 15298643674L) - ((int128)tmp_q[8] * 7341466365674L) + ((int128)tmp_q[9] * 5868169780971L) - ((-((int128)tmp_q[10] * 3354125399230L) - ((int128)tmp_q[11] * 1162084541777L)) * 10);
	tmp_zero[10] = -((int128)tmp_q[0] * 1162084541777L) - ((int128)tmp_q[1] * 3369654567083L) + ((int128)tmp_q[2] * 4141609118370L) - ((int128)tmp_q[3] * 5492226004432L) + ((int128)tmp_q[4] * 136330842539L) - ((int128)tmp_q[5] * 4139106621025L) + ((int128)tmp_q[6] * 40831642937L) - ((int128)tmp_q[7] * 34241695394L) - ((int128)tmp_q[8] * 15298643674L) - ((int128)tmp_q[9] * 7341466365674L) + ((int128)tmp_q[10] * 5868169780971L) + ((int128)tmp_q[11] * 33541253992300L);
	tmp_zero[11] = -((int128)tmp_q[0] * 3354125399230L) - ((int128)tmp_q[1] * 1162084541777L) - ((int128)tmp_q[2] * 3369654567083L) + ((int128)tmp_q[3] * 4141609118370L) - ((int128)tmp_q[4] * 5492226004432L) + ((int128)tmp_q[5] * 136330842539L) - ((int128)tmp_q[6] * 4139106621025L) + ((int128)tmp_q[7] * 40831642937L) - ((int128)tmp_q[8] * 34241695394L) - ((int128)tmp_q[9] * 15298643674L) - ((int128)tmp_q[10] * 7341466365674L) + ((int128)tmp_q[11] * 5868169780971L);

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
	rop[11] = (op[11] + tmp_zero[11]) >> WORD_SIZE;
}

