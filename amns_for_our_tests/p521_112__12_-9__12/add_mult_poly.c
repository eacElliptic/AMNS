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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 9);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 9);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 9);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 9);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 9);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 9);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 9);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 9);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 9);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 9);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 9);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 9);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 18);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 9);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 18);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 9);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 18);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 9);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 18);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 9);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 18);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 9);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 41458622030470971UL) + ((((uint64_t)op[1] * 15265390351117160521UL) + ((uint64_t)op[2] * 326700523315903953UL) + ((uint64_t)op[3] * 12825064626874971187UL) + ((uint64_t)op[4] * 4814587162290322127UL) + ((uint64_t)op[5] * 7514595851264824095UL) + ((uint64_t)op[6] * 5435229718218315520UL) + ((uint64_t)op[7] * 4909456695852457600UL) + ((uint64_t)op[8] * 3153486730549308732UL) + ((uint64_t)op[9] * 7281846818456319898UL) + ((uint64_t)op[10] * 11712592278866228615UL) + ((uint64_t)op[11] * 3344981151613424656UL)) * 18446744073709551607);
	tmp_q[1] = ((uint64_t)op[0] * 3344981151613424656UL) + ((uint64_t)op[1] * 41458622030470971UL) + ((((uint64_t)op[2] * 15265390351117160521UL) + ((uint64_t)op[3] * 326700523315903953UL) + ((uint64_t)op[4] * 12825064626874971187UL) + ((uint64_t)op[5] * 4814587162290322127UL) + ((uint64_t)op[6] * 7514595851264824095UL) + ((uint64_t)op[7] * 5435229718218315520UL) + ((uint64_t)op[8] * 4909456695852457600UL) + ((uint64_t)op[9] * 3153486730549308732UL) + ((uint64_t)op[10] * 7281846818456319898UL) + ((uint64_t)op[11] * 11712592278866228615UL)) * 18446744073709551607);
	tmp_q[2] = ((uint64_t)op[0] * 11712592278866228615UL) + ((uint64_t)op[1] * 3344981151613424656UL) + ((uint64_t)op[2] * 41458622030470971UL) + ((((uint64_t)op[3] * 15265390351117160521UL) + ((uint64_t)op[4] * 326700523315903953UL) + ((uint64_t)op[5] * 12825064626874971187UL) + ((uint64_t)op[6] * 4814587162290322127UL) + ((uint64_t)op[7] * 7514595851264824095UL) + ((uint64_t)op[8] * 5435229718218315520UL) + ((uint64_t)op[9] * 4909456695852457600UL) + ((uint64_t)op[10] * 3153486730549308732UL) + ((uint64_t)op[11] * 7281846818456319898UL)) * 18446744073709551607);
	tmp_q[3] = ((uint64_t)op[0] * 7281846818456319898UL) + ((uint64_t)op[1] * 11712592278866228615UL) + ((uint64_t)op[2] * 3344981151613424656UL) + ((uint64_t)op[3] * 41458622030470971UL) + ((((uint64_t)op[4] * 15265390351117160521UL) + ((uint64_t)op[5] * 326700523315903953UL) + ((uint64_t)op[6] * 12825064626874971187UL) + ((uint64_t)op[7] * 4814587162290322127UL) + ((uint64_t)op[8] * 7514595851264824095UL) + ((uint64_t)op[9] * 5435229718218315520UL) + ((uint64_t)op[10] * 4909456695852457600UL) + ((uint64_t)op[11] * 3153486730549308732UL)) * 18446744073709551607);
	tmp_q[4] = ((uint64_t)op[0] * 3153486730549308732UL) + ((uint64_t)op[1] * 7281846818456319898UL) + ((uint64_t)op[2] * 11712592278866228615UL) + ((uint64_t)op[3] * 3344981151613424656UL) + ((uint64_t)op[4] * 41458622030470971UL) + ((((uint64_t)op[5] * 15265390351117160521UL) + ((uint64_t)op[6] * 326700523315903953UL) + ((uint64_t)op[7] * 12825064626874971187UL) + ((uint64_t)op[8] * 4814587162290322127UL) + ((uint64_t)op[9] * 7514595851264824095UL) + ((uint64_t)op[10] * 5435229718218315520UL) + ((uint64_t)op[11] * 4909456695852457600UL)) * 18446744073709551607);
	tmp_q[5] = ((uint64_t)op[0] * 4909456695852457600UL) + ((uint64_t)op[1] * 3153486730549308732UL) + ((uint64_t)op[2] * 7281846818456319898UL) + ((uint64_t)op[3] * 11712592278866228615UL) + ((uint64_t)op[4] * 3344981151613424656UL) + ((uint64_t)op[5] * 41458622030470971UL) + ((((uint64_t)op[6] * 15265390351117160521UL) + ((uint64_t)op[7] * 326700523315903953UL) + ((uint64_t)op[8] * 12825064626874971187UL) + ((uint64_t)op[9] * 4814587162290322127UL) + ((uint64_t)op[10] * 7514595851264824095UL) + ((uint64_t)op[11] * 5435229718218315520UL)) * 18446744073709551607);
	tmp_q[6] = ((uint64_t)op[0] * 5435229718218315520UL) + ((uint64_t)op[1] * 4909456695852457600UL) + ((uint64_t)op[2] * 3153486730549308732UL) + ((uint64_t)op[3] * 7281846818456319898UL) + ((uint64_t)op[4] * 11712592278866228615UL) + ((uint64_t)op[5] * 3344981151613424656UL) + ((uint64_t)op[6] * 41458622030470971UL) + ((((uint64_t)op[7] * 15265390351117160521UL) + ((uint64_t)op[8] * 326700523315903953UL) + ((uint64_t)op[9] * 12825064626874971187UL) + ((uint64_t)op[10] * 4814587162290322127UL) + ((uint64_t)op[11] * 7514595851264824095UL)) * 18446744073709551607);
	tmp_q[7] = ((uint64_t)op[0] * 7514595851264824095UL) + ((uint64_t)op[1] * 5435229718218315520UL) + ((uint64_t)op[2] * 4909456695852457600UL) + ((uint64_t)op[3] * 3153486730549308732UL) + ((uint64_t)op[4] * 7281846818456319898UL) + ((uint64_t)op[5] * 11712592278866228615UL) + ((uint64_t)op[6] * 3344981151613424656UL) + ((uint64_t)op[7] * 41458622030470971UL) + ((((uint64_t)op[8] * 15265390351117160521UL) + ((uint64_t)op[9] * 326700523315903953UL) + ((uint64_t)op[10] * 12825064626874971187UL) + ((uint64_t)op[11] * 4814587162290322127UL)) * 18446744073709551607);
	tmp_q[8] = ((uint64_t)op[0] * 4814587162290322127UL) + ((uint64_t)op[1] * 7514595851264824095UL) + ((uint64_t)op[2] * 5435229718218315520UL) + ((uint64_t)op[3] * 4909456695852457600UL) + ((uint64_t)op[4] * 3153486730549308732UL) + ((uint64_t)op[5] * 7281846818456319898UL) + ((uint64_t)op[6] * 11712592278866228615UL) + ((uint64_t)op[7] * 3344981151613424656UL) + ((uint64_t)op[8] * 41458622030470971UL) + ((((uint64_t)op[9] * 15265390351117160521UL) + ((uint64_t)op[10] * 326700523315903953UL) + ((uint64_t)op[11] * 12825064626874971187UL)) * 18446744073709551607);
	tmp_q[9] = ((uint64_t)op[0] * 12825064626874971187UL) + ((uint64_t)op[1] * 4814587162290322127UL) + ((uint64_t)op[2] * 7514595851264824095UL) + ((uint64_t)op[3] * 5435229718218315520UL) + ((uint64_t)op[4] * 4909456695852457600UL) + ((uint64_t)op[5] * 3153486730549308732UL) + ((uint64_t)op[6] * 7281846818456319898UL) + ((uint64_t)op[7] * 11712592278866228615UL) + ((uint64_t)op[8] * 3344981151613424656UL) + ((uint64_t)op[9] * 41458622030470971UL) + ((((uint64_t)op[10] * 15265390351117160521UL) + ((uint64_t)op[11] * 326700523315903953UL)) * 18446744073709551607);
	tmp_q[10] = ((uint64_t)op[0] * 326700523315903953UL) + ((uint64_t)op[1] * 12825064626874971187UL) + ((uint64_t)op[2] * 4814587162290322127UL) + ((uint64_t)op[3] * 7514595851264824095UL) + ((uint64_t)op[4] * 5435229718218315520UL) + ((uint64_t)op[5] * 4909456695852457600UL) + ((uint64_t)op[6] * 3153486730549308732UL) + ((uint64_t)op[7] * 7281846818456319898UL) + ((uint64_t)op[8] * 11712592278866228615UL) + ((uint64_t)op[9] * 3344981151613424656UL) + ((uint64_t)op[10] * 41458622030470971UL) + ((uint64_t)op[11] * 10185439429621968239UL);
	tmp_q[11] = ((uint64_t)op[0] * 15265390351117160521UL) + ((uint64_t)op[1] * 326700523315903953UL) + ((uint64_t)op[2] * 12825064626874971187UL) + ((uint64_t)op[3] * 4814587162290322127UL) + ((uint64_t)op[4] * 7514595851264824095UL) + ((uint64_t)op[5] * 5435229718218315520UL) + ((uint64_t)op[6] * 4909456695852457600UL) + ((uint64_t)op[7] * 3153486730549308732UL) + ((uint64_t)op[8] * 7281846818456319898UL) + ((uint64_t)op[9] * 11712592278866228615UL) + ((uint64_t)op[10] * 3344981151613424656UL) + ((uint64_t)op[11] * 41458622030470971UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1998894887168L) - ((-((int128)tmp_q[1] * 1001515882109L) + ((int128)tmp_q[2] * 1654455341667L) + ((int128)tmp_q[3] * 4437593908487L) - ((int128)tmp_q[4] * 2852566735783L) - ((int128)tmp_q[5] * 5592943170870L) - ((int128)tmp_q[6] * 2233511772555L) + ((int128)tmp_q[7] * 7568347825490L) + ((int128)tmp_q[8] * 728633710545L) - ((int128)tmp_q[9] * 7829488419734L) + ((int128)tmp_q[10] * 5922369628380L) + ((int128)tmp_q[11] * 4245546247235L)) * 9);
	tmp_zero[1] = ((int128)tmp_q[0] * 4245546247235L) + ((int128)tmp_q[1] * 1998894887168L) - ((-((int128)tmp_q[2] * 1001515882109L) + ((int128)tmp_q[3] * 1654455341667L) + ((int128)tmp_q[4] * 4437593908487L) - ((int128)tmp_q[5] * 2852566735783L) - ((int128)tmp_q[6] * 5592943170870L) - ((int128)tmp_q[7] * 2233511772555L) + ((int128)tmp_q[8] * 7568347825490L) + ((int128)tmp_q[9] * 728633710545L) - ((int128)tmp_q[10] * 7829488419734L) + ((int128)tmp_q[11] * 5922369628380L)) * 9);
	tmp_zero[2] = ((int128)tmp_q[0] * 5922369628380L) + ((int128)tmp_q[1] * 4245546247235L) + ((int128)tmp_q[2] * 1998894887168L) - ((-((int128)tmp_q[3] * 1001515882109L) + ((int128)tmp_q[4] * 1654455341667L) + ((int128)tmp_q[5] * 4437593908487L) - ((int128)tmp_q[6] * 2852566735783L) - ((int128)tmp_q[7] * 5592943170870L) - ((int128)tmp_q[8] * 2233511772555L) + ((int128)tmp_q[9] * 7568347825490L) + ((int128)tmp_q[10] * 728633710545L) - ((int128)tmp_q[11] * 7829488419734L)) * 9);
	tmp_zero[3] = -((int128)tmp_q[0] * 7829488419734L) + ((int128)tmp_q[1] * 5922369628380L) + ((int128)tmp_q[2] * 4245546247235L) + ((int128)tmp_q[3] * 1998894887168L) - ((-((int128)tmp_q[4] * 1001515882109L) + ((int128)tmp_q[5] * 1654455341667L) + ((int128)tmp_q[6] * 4437593908487L) - ((int128)tmp_q[7] * 2852566735783L) - ((int128)tmp_q[8] * 5592943170870L) - ((int128)tmp_q[9] * 2233511772555L) + ((int128)tmp_q[10] * 7568347825490L) + ((int128)tmp_q[11] * 728633710545L)) * 9);
	tmp_zero[4] = ((int128)tmp_q[0] * 728633710545L) - ((int128)tmp_q[1] * 7829488419734L) + ((int128)tmp_q[2] * 5922369628380L) + ((int128)tmp_q[3] * 4245546247235L) + ((int128)tmp_q[4] * 1998894887168L) - ((-((int128)tmp_q[5] * 1001515882109L) + ((int128)tmp_q[6] * 1654455341667L) + ((int128)tmp_q[7] * 4437593908487L) - ((int128)tmp_q[8] * 2852566735783L) - ((int128)tmp_q[9] * 5592943170870L) - ((int128)tmp_q[10] * 2233511772555L) + ((int128)tmp_q[11] * 7568347825490L)) * 9);
	tmp_zero[5] = ((int128)tmp_q[0] * 7568347825490L) + ((int128)tmp_q[1] * 728633710545L) - ((int128)tmp_q[2] * 7829488419734L) + ((int128)tmp_q[3] * 5922369628380L) + ((int128)tmp_q[4] * 4245546247235L) + ((int128)tmp_q[5] * 1998894887168L) - ((-((int128)tmp_q[6] * 1001515882109L) + ((int128)tmp_q[7] * 1654455341667L) + ((int128)tmp_q[8] * 4437593908487L) - ((int128)tmp_q[9] * 2852566735783L) - ((int128)tmp_q[10] * 5592943170870L) - ((int128)tmp_q[11] * 2233511772555L)) * 9);
	tmp_zero[6] = -((int128)tmp_q[0] * 2233511772555L) + ((int128)tmp_q[1] * 7568347825490L) + ((int128)tmp_q[2] * 728633710545L) - ((int128)tmp_q[3] * 7829488419734L) + ((int128)tmp_q[4] * 5922369628380L) + ((int128)tmp_q[5] * 4245546247235L) + ((int128)tmp_q[6] * 1998894887168L) - ((-((int128)tmp_q[7] * 1001515882109L) + ((int128)tmp_q[8] * 1654455341667L) + ((int128)tmp_q[9] * 4437593908487L) - ((int128)tmp_q[10] * 2852566735783L) - ((int128)tmp_q[11] * 5592943170870L)) * 9);
	tmp_zero[7] = -((int128)tmp_q[0] * 5592943170870L) - ((int128)tmp_q[1] * 2233511772555L) + ((int128)tmp_q[2] * 7568347825490L) + ((int128)tmp_q[3] * 728633710545L) - ((int128)tmp_q[4] * 7829488419734L) + ((int128)tmp_q[5] * 5922369628380L) + ((int128)tmp_q[6] * 4245546247235L) + ((int128)tmp_q[7] * 1998894887168L) - ((-((int128)tmp_q[8] * 1001515882109L) + ((int128)tmp_q[9] * 1654455341667L) + ((int128)tmp_q[10] * 4437593908487L) - ((int128)tmp_q[11] * 2852566735783L)) * 9);
	tmp_zero[8] = -((int128)tmp_q[0] * 2852566735783L) - ((int128)tmp_q[1] * 5592943170870L) - ((int128)tmp_q[2] * 2233511772555L) + ((int128)tmp_q[3] * 7568347825490L) + ((int128)tmp_q[4] * 728633710545L) - ((int128)tmp_q[5] * 7829488419734L) + ((int128)tmp_q[6] * 5922369628380L) + ((int128)tmp_q[7] * 4245546247235L) + ((int128)tmp_q[8] * 1998894887168L) - ((-((int128)tmp_q[9] * 1001515882109L) + ((int128)tmp_q[10] * 1654455341667L) + ((int128)tmp_q[11] * 4437593908487L)) * 9);
	tmp_zero[9] = ((int128)tmp_q[0] * 4437593908487L) - ((int128)tmp_q[1] * 2852566735783L) - ((int128)tmp_q[2] * 5592943170870L) - ((int128)tmp_q[3] * 2233511772555L) + ((int128)tmp_q[4] * 7568347825490L) + ((int128)tmp_q[5] * 728633710545L) - ((int128)tmp_q[6] * 7829488419734L) + ((int128)tmp_q[7] * 5922369628380L) + ((int128)tmp_q[8] * 4245546247235L) + ((int128)tmp_q[9] * 1998894887168L) - ((-((int128)tmp_q[10] * 1001515882109L) + ((int128)tmp_q[11] * 1654455341667L)) * 9);
	tmp_zero[10] = ((int128)tmp_q[0] * 1654455341667L) + ((int128)tmp_q[1] * 4437593908487L) - ((int128)tmp_q[2] * 2852566735783L) - ((int128)tmp_q[3] * 5592943170870L) - ((int128)tmp_q[4] * 2233511772555L) + ((int128)tmp_q[5] * 7568347825490L) + ((int128)tmp_q[6] * 728633710545L) - ((int128)tmp_q[7] * 7829488419734L) + ((int128)tmp_q[8] * 5922369628380L) + ((int128)tmp_q[9] * 4245546247235L) + ((int128)tmp_q[10] * 1998894887168L) + ((int128)tmp_q[11] * 9013642938981L);
	tmp_zero[11] = -((int128)tmp_q[0] * 1001515882109L) + ((int128)tmp_q[1] * 1654455341667L) + ((int128)tmp_q[2] * 4437593908487L) - ((int128)tmp_q[3] * 2852566735783L) - ((int128)tmp_q[4] * 5592943170870L) - ((int128)tmp_q[5] * 2233511772555L) + ((int128)tmp_q[6] * 7568347825490L) + ((int128)tmp_q[7] * 728633710545L) - ((int128)tmp_q[8] * 7829488419734L) + ((int128)tmp_q[9] * 5922369628380L) + ((int128)tmp_q[10] * 4245546247235L) + ((int128)tmp_q[11] * 1998894887168L);

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

