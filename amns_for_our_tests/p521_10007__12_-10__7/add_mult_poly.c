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
	tmp_q[0] = ((uint64_t)op[0] * 15431124270551605419UL) + ((((uint64_t)op[1] * 6744729342025392123UL) + ((uint64_t)op[2] * 16189961347825308534UL) + ((uint64_t)op[3] * 17844050137228712259UL) + ((uint64_t)op[4] * 11446284856851719951UL) + ((uint64_t)op[5] * 908747623513467604UL) + ((uint64_t)op[6] * 3674476697559281267UL) + ((uint64_t)op[7] * 13759344822164504210UL) + ((uint64_t)op[8] * 1952328131235802638UL) + ((uint64_t)op[9] * 8131561547005457202UL) + ((uint64_t)op[10] * 14323211572726630239UL) + ((uint64_t)op[11] * 5086946768604074256UL)) * 18446744073709551606);
	tmp_q[1] = ((uint64_t)op[0] * 5086946768604074256UL) + ((uint64_t)op[1] * 15431124270551605419UL) + ((((uint64_t)op[2] * 6744729342025392123UL) + ((uint64_t)op[3] * 16189961347825308534UL) + ((uint64_t)op[4] * 17844050137228712259UL) + ((uint64_t)op[5] * 11446284856851719951UL) + ((uint64_t)op[6] * 908747623513467604UL) + ((uint64_t)op[7] * 3674476697559281267UL) + ((uint64_t)op[8] * 13759344822164504210UL) + ((uint64_t)op[9] * 1952328131235802638UL) + ((uint64_t)op[10] * 8131561547005457202UL) + ((uint64_t)op[11] * 14323211572726630239UL)) * 18446744073709551606);
	tmp_q[2] = ((uint64_t)op[0] * 14323211572726630239UL) + ((uint64_t)op[1] * 5086946768604074256UL) + ((uint64_t)op[2] * 15431124270551605419UL) + ((((uint64_t)op[3] * 6744729342025392123UL) + ((uint64_t)op[4] * 16189961347825308534UL) + ((uint64_t)op[5] * 17844050137228712259UL) + ((uint64_t)op[6] * 11446284856851719951UL) + ((uint64_t)op[7] * 908747623513467604UL) + ((uint64_t)op[8] * 3674476697559281267UL) + ((uint64_t)op[9] * 13759344822164504210UL) + ((uint64_t)op[10] * 1952328131235802638UL) + ((uint64_t)op[11] * 8131561547005457202UL)) * 18446744073709551606);
	tmp_q[3] = ((uint64_t)op[0] * 8131561547005457202UL) + ((uint64_t)op[1] * 14323211572726630239UL) + ((uint64_t)op[2] * 5086946768604074256UL) + ((uint64_t)op[3] * 15431124270551605419UL) + ((((uint64_t)op[4] * 6744729342025392123UL) + ((uint64_t)op[5] * 16189961347825308534UL) + ((uint64_t)op[6] * 17844050137228712259UL) + ((uint64_t)op[7] * 11446284856851719951UL) + ((uint64_t)op[8] * 908747623513467604UL) + ((uint64_t)op[9] * 3674476697559281267UL) + ((uint64_t)op[10] * 13759344822164504210UL) + ((uint64_t)op[11] * 1952328131235802638UL)) * 18446744073709551606);
	tmp_q[4] = ((uint64_t)op[0] * 1952328131235802638UL) + ((uint64_t)op[1] * 8131561547005457202UL) + ((uint64_t)op[2] * 14323211572726630239UL) + ((uint64_t)op[3] * 5086946768604074256UL) + ((uint64_t)op[4] * 15431124270551605419UL) + ((((uint64_t)op[5] * 6744729342025392123UL) + ((uint64_t)op[6] * 16189961347825308534UL) + ((uint64_t)op[7] * 17844050137228712259UL) + ((uint64_t)op[8] * 11446284856851719951UL) + ((uint64_t)op[9] * 908747623513467604UL) + ((uint64_t)op[10] * 3674476697559281267UL) + ((uint64_t)op[11] * 13759344822164504210UL)) * 18446744073709551606);
	tmp_q[5] = ((uint64_t)op[0] * 13759344822164504210UL) + ((uint64_t)op[1] * 1952328131235802638UL) + ((uint64_t)op[2] * 8131561547005457202UL) + ((uint64_t)op[3] * 14323211572726630239UL) + ((uint64_t)op[4] * 5086946768604074256UL) + ((uint64_t)op[5] * 15431124270551605419UL) + ((((uint64_t)op[6] * 6744729342025392123UL) + ((uint64_t)op[7] * 16189961347825308534UL) + ((uint64_t)op[8] * 17844050137228712259UL) + ((uint64_t)op[9] * 11446284856851719951UL) + ((uint64_t)op[10] * 908747623513467604UL) + ((uint64_t)op[11] * 3674476697559281267UL)) * 18446744073709551606);
	tmp_q[6] = ((uint64_t)op[0] * 3674476697559281267UL) + ((uint64_t)op[1] * 13759344822164504210UL) + ((uint64_t)op[2] * 1952328131235802638UL) + ((uint64_t)op[3] * 8131561547005457202UL) + ((uint64_t)op[4] * 14323211572726630239UL) + ((uint64_t)op[5] * 5086946768604074256UL) + ((uint64_t)op[6] * 15431124270551605419UL) + ((((uint64_t)op[7] * 6744729342025392123UL) + ((uint64_t)op[8] * 16189961347825308534UL) + ((uint64_t)op[9] * 17844050137228712259UL) + ((uint64_t)op[10] * 11446284856851719951UL) + ((uint64_t)op[11] * 908747623513467604UL)) * 18446744073709551606);
	tmp_q[7] = ((uint64_t)op[0] * 908747623513467604UL) + ((uint64_t)op[1] * 3674476697559281267UL) + ((uint64_t)op[2] * 13759344822164504210UL) + ((uint64_t)op[3] * 1952328131235802638UL) + ((uint64_t)op[4] * 8131561547005457202UL) + ((uint64_t)op[5] * 14323211572726630239UL) + ((uint64_t)op[6] * 5086946768604074256UL) + ((uint64_t)op[7] * 15431124270551605419UL) + ((((uint64_t)op[8] * 6744729342025392123UL) + ((uint64_t)op[9] * 16189961347825308534UL) + ((uint64_t)op[10] * 17844050137228712259UL) + ((uint64_t)op[11] * 11446284856851719951UL)) * 18446744073709551606);
	tmp_q[8] = ((uint64_t)op[0] * 11446284856851719951UL) + ((uint64_t)op[1] * 908747623513467604UL) + ((uint64_t)op[2] * 3674476697559281267UL) + ((uint64_t)op[3] * 13759344822164504210UL) + ((uint64_t)op[4] * 1952328131235802638UL) + ((uint64_t)op[5] * 8131561547005457202UL) + ((uint64_t)op[6] * 14323211572726630239UL) + ((uint64_t)op[7] * 5086946768604074256UL) + ((uint64_t)op[8] * 15431124270551605419UL) + ((((uint64_t)op[9] * 6744729342025392123UL) + ((uint64_t)op[10] * 16189961347825308534UL) + ((uint64_t)op[11] * 17844050137228712259UL)) * 18446744073709551606);
	tmp_q[9] = ((uint64_t)op[0] * 17844050137228712259UL) + ((uint64_t)op[1] * 11446284856851719951UL) + ((uint64_t)op[2] * 908747623513467604UL) + ((uint64_t)op[3] * 3674476697559281267UL) + ((uint64_t)op[4] * 13759344822164504210UL) + ((uint64_t)op[5] * 1952328131235802638UL) + ((uint64_t)op[6] * 8131561547005457202UL) + ((uint64_t)op[7] * 14323211572726630239UL) + ((uint64_t)op[8] * 5086946768604074256UL) + ((uint64_t)op[9] * 15431124270551605419UL) + ((((uint64_t)op[10] * 6744729342025392123UL) + ((uint64_t)op[11] * 16189961347825308534UL)) * 18446744073709551606);
	tmp_q[10] = ((uint64_t)op[0] * 16189961347825308534UL) + ((uint64_t)op[1] * 17844050137228712259UL) + ((uint64_t)op[2] * 11446284856851719951UL) + ((uint64_t)op[3] * 908747623513467604UL) + ((uint64_t)op[4] * 3674476697559281267UL) + ((uint64_t)op[5] * 13759344822164504210UL) + ((uint64_t)op[6] * 1952328131235802638UL) + ((uint64_t)op[7] * 8131561547005457202UL) + ((uint64_t)op[8] * 14323211572726630239UL) + ((uint64_t)op[9] * 5086946768604074256UL) + ((uint64_t)op[10] * 15431124270551605419UL) + ((uint64_t)op[11] * 6339682874584285234UL);
	tmp_q[11] = ((uint64_t)op[0] * 6744729342025392123UL) + ((uint64_t)op[1] * 16189961347825308534UL) + ((uint64_t)op[2] * 17844050137228712259UL) + ((uint64_t)op[3] * 11446284856851719951UL) + ((uint64_t)op[4] * 908747623513467604UL) + ((uint64_t)op[5] * 3674476697559281267UL) + ((uint64_t)op[6] * 13759344822164504210UL) + ((uint64_t)op[7] * 1952328131235802638UL) + ((uint64_t)op[8] * 8131561547005457202UL) + ((uint64_t)op[9] * 14323211572726630239UL) + ((uint64_t)op[10] * 5086946768604074256UL) + ((uint64_t)op[11] * 15431124270551605419UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 803596199881L) - ((((int128)tmp_q[1] * 580598060439L) - ((int128)tmp_q[2] * 4264163915826L) - ((int128)tmp_q[3] * 5597658947451L) - ((int128)tmp_q[4] * 4379795745878L) + ((int128)tmp_q[5] * 2924129035188L) - ((int128)tmp_q[6] * 2845378044860L) - ((int128)tmp_q[7] * 3420111681664L) - ((int128)tmp_q[8] * 7107432341503L) - ((int128)tmp_q[9] * 2812157261752L) - ((int128)tmp_q[10] * 1545932451075L) + ((int128)tmp_q[11] * 1370812617974L)) * 10);
	tmp_zero[1] = ((int128)tmp_q[0] * 1370812617974L) - ((int128)tmp_q[1] * 803596199881L) - ((((int128)tmp_q[2] * 580598060439L) - ((int128)tmp_q[3] * 4264163915826L) - ((int128)tmp_q[4] * 5597658947451L) - ((int128)tmp_q[5] * 4379795745878L) + ((int128)tmp_q[6] * 2924129035188L) - ((int128)tmp_q[7] * 2845378044860L) - ((int128)tmp_q[8] * 3420111681664L) - ((int128)tmp_q[9] * 7107432341503L) - ((int128)tmp_q[10] * 2812157261752L) - ((int128)tmp_q[11] * 1545932451075L)) * 10);
	tmp_zero[2] = -((int128)tmp_q[0] * 1545932451075L) + ((int128)tmp_q[1] * 1370812617974L) - ((int128)tmp_q[2] * 803596199881L) - ((((int128)tmp_q[3] * 580598060439L) - ((int128)tmp_q[4] * 4264163915826L) - ((int128)tmp_q[5] * 5597658947451L) - ((int128)tmp_q[6] * 4379795745878L) + ((int128)tmp_q[7] * 2924129035188L) - ((int128)tmp_q[8] * 2845378044860L) - ((int128)tmp_q[9] * 3420111681664L) - ((int128)tmp_q[10] * 7107432341503L) - ((int128)tmp_q[11] * 2812157261752L)) * 10);
	tmp_zero[3] = -((int128)tmp_q[0] * 2812157261752L) - ((int128)tmp_q[1] * 1545932451075L) + ((int128)tmp_q[2] * 1370812617974L) - ((int128)tmp_q[3] * 803596199881L) - ((((int128)tmp_q[4] * 580598060439L) - ((int128)tmp_q[5] * 4264163915826L) - ((int128)tmp_q[6] * 5597658947451L) - ((int128)tmp_q[7] * 4379795745878L) + ((int128)tmp_q[8] * 2924129035188L) - ((int128)tmp_q[9] * 2845378044860L) - ((int128)tmp_q[10] * 3420111681664L) - ((int128)tmp_q[11] * 7107432341503L)) * 10);
	tmp_zero[4] = -((int128)tmp_q[0] * 7107432341503L) - ((int128)tmp_q[1] * 2812157261752L) - ((int128)tmp_q[2] * 1545932451075L) + ((int128)tmp_q[3] * 1370812617974L) - ((int128)tmp_q[4] * 803596199881L) - ((((int128)tmp_q[5] * 580598060439L) - ((int128)tmp_q[6] * 4264163915826L) - ((int128)tmp_q[7] * 5597658947451L) - ((int128)tmp_q[8] * 4379795745878L) + ((int128)tmp_q[9] * 2924129035188L) - ((int128)tmp_q[10] * 2845378044860L) - ((int128)tmp_q[11] * 3420111681664L)) * 10);
	tmp_zero[5] = -((int128)tmp_q[0] * 3420111681664L) - ((int128)tmp_q[1] * 7107432341503L) - ((int128)tmp_q[2] * 2812157261752L) - ((int128)tmp_q[3] * 1545932451075L) + ((int128)tmp_q[4] * 1370812617974L) - ((int128)tmp_q[5] * 803596199881L) - ((((int128)tmp_q[6] * 580598060439L) - ((int128)tmp_q[7] * 4264163915826L) - ((int128)tmp_q[8] * 5597658947451L) - ((int128)tmp_q[9] * 4379795745878L) + ((int128)tmp_q[10] * 2924129035188L) - ((int128)tmp_q[11] * 2845378044860L)) * 10);
	tmp_zero[6] = -((int128)tmp_q[0] * 2845378044860L) - ((int128)tmp_q[1] * 3420111681664L) - ((int128)tmp_q[2] * 7107432341503L) - ((int128)tmp_q[3] * 2812157261752L) - ((int128)tmp_q[4] * 1545932451075L) + ((int128)tmp_q[5] * 1370812617974L) - ((int128)tmp_q[6] * 803596199881L) - ((((int128)tmp_q[7] * 580598060439L) - ((int128)tmp_q[8] * 4264163915826L) - ((int128)tmp_q[9] * 5597658947451L) - ((int128)tmp_q[10] * 4379795745878L) + ((int128)tmp_q[11] * 2924129035188L)) * 10);
	tmp_zero[7] = ((int128)tmp_q[0] * 2924129035188L) - ((int128)tmp_q[1] * 2845378044860L) - ((int128)tmp_q[2] * 3420111681664L) - ((int128)tmp_q[3] * 7107432341503L) - ((int128)tmp_q[4] * 2812157261752L) - ((int128)tmp_q[5] * 1545932451075L) + ((int128)tmp_q[6] * 1370812617974L) - ((int128)tmp_q[7] * 803596199881L) - ((((int128)tmp_q[8] * 580598060439L) - ((int128)tmp_q[9] * 4264163915826L) - ((int128)tmp_q[10] * 5597658947451L) - ((int128)tmp_q[11] * 4379795745878L)) * 10);
	tmp_zero[8] = -((int128)tmp_q[0] * 4379795745878L) + ((int128)tmp_q[1] * 2924129035188L) - ((int128)tmp_q[2] * 2845378044860L) - ((int128)tmp_q[3] * 3420111681664L) - ((int128)tmp_q[4] * 7107432341503L) - ((int128)tmp_q[5] * 2812157261752L) - ((int128)tmp_q[6] * 1545932451075L) + ((int128)tmp_q[7] * 1370812617974L) - ((int128)tmp_q[8] * 803596199881L) - ((((int128)tmp_q[9] * 580598060439L) - ((int128)tmp_q[10] * 4264163915826L) - ((int128)tmp_q[11] * 5597658947451L)) * 10);
	tmp_zero[9] = -((int128)tmp_q[0] * 5597658947451L) - ((int128)tmp_q[1] * 4379795745878L) + ((int128)tmp_q[2] * 2924129035188L) - ((int128)tmp_q[3] * 2845378044860L) - ((int128)tmp_q[4] * 3420111681664L) - ((int128)tmp_q[5] * 7107432341503L) - ((int128)tmp_q[6] * 2812157261752L) - ((int128)tmp_q[7] * 1545932451075L) + ((int128)tmp_q[8] * 1370812617974L) - ((int128)tmp_q[9] * 803596199881L) - ((((int128)tmp_q[10] * 580598060439L) - ((int128)tmp_q[11] * 4264163915826L)) * 10);
	tmp_zero[10] = -((int128)tmp_q[0] * 4264163915826L) - ((int128)tmp_q[1] * 5597658947451L) - ((int128)tmp_q[2] * 4379795745878L) + ((int128)tmp_q[3] * 2924129035188L) - ((int128)tmp_q[4] * 2845378044860L) - ((int128)tmp_q[5] * 3420111681664L) - ((int128)tmp_q[6] * 7107432341503L) - ((int128)tmp_q[7] * 2812157261752L) - ((int128)tmp_q[8] * 1545932451075L) + ((int128)tmp_q[9] * 1370812617974L) - ((int128)tmp_q[10] * 803596199881L) - ((int128)tmp_q[11] * 5805980604390L);
	tmp_zero[11] = ((int128)tmp_q[0] * 580598060439L) - ((int128)tmp_q[1] * 4264163915826L) - ((int128)tmp_q[2] * 5597658947451L) - ((int128)tmp_q[3] * 4379795745878L) + ((int128)tmp_q[4] * 2924129035188L) - ((int128)tmp_q[5] * 2845378044860L) - ((int128)tmp_q[6] * 3420111681664L) - ((int128)tmp_q[7] * 7107432341503L) - ((int128)tmp_q[8] * 2812157261752L) - ((int128)tmp_q[9] * 1545932451075L) + ((int128)tmp_q[10] * 1370812617974L) - ((int128)tmp_q[11] * 803596199881L);

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

