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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) << 4);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) << 4);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) << 4);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) << 4);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) << 4);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) << 4);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) << 4);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) << 4);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) << 4);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) << 4);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) << 4);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) << 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) << 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) << 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) << 4);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) << 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) << 4);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) << 5);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) << 4);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11006936700485570137UL) + ((((uint64_t)op[1] * 11426515914107012300UL) + ((uint64_t)op[2] * 5619274444686996823UL) + ((uint64_t)op[3] * 16241926334950006912UL) + ((uint64_t)op[4] * 1388654943435253604UL) + ((uint64_t)op[5] * 1453109792944704396UL) + ((uint64_t)op[6] * 12710492834476651738UL) + ((uint64_t)op[7] * 989886695836914904UL) + ((uint64_t)op[8] * 6851573905507941151UL) + ((uint64_t)op[9] * 4360416026152177361UL) + ((uint64_t)op[10] * 276319644861830145UL) + ((uint64_t)op[11] * 7773059531132371638UL)) * 18446744073709551600);
	tmp_q[1] = ((uint64_t)op[0] * 7773059531132371638UL) + ((uint64_t)op[1] * 11006936700485570137UL) + ((((uint64_t)op[2] * 11426515914107012300UL) + ((uint64_t)op[3] * 5619274444686996823UL) + ((uint64_t)op[4] * 16241926334950006912UL) + ((uint64_t)op[5] * 1388654943435253604UL) + ((uint64_t)op[6] * 1453109792944704396UL) + ((uint64_t)op[7] * 12710492834476651738UL) + ((uint64_t)op[8] * 989886695836914904UL) + ((uint64_t)op[9] * 6851573905507941151UL) + ((uint64_t)op[10] * 4360416026152177361UL) + ((uint64_t)op[11] * 276319644861830145UL)) * 18446744073709551600);
	tmp_q[2] = ((uint64_t)op[0] * 276319644861830145UL) + ((uint64_t)op[1] * 7773059531132371638UL) + ((uint64_t)op[2] * 11006936700485570137UL) + ((((uint64_t)op[3] * 11426515914107012300UL) + ((uint64_t)op[4] * 5619274444686996823UL) + ((uint64_t)op[5] * 16241926334950006912UL) + ((uint64_t)op[6] * 1388654943435253604UL) + ((uint64_t)op[7] * 1453109792944704396UL) + ((uint64_t)op[8] * 12710492834476651738UL) + ((uint64_t)op[9] * 989886695836914904UL) + ((uint64_t)op[10] * 6851573905507941151UL) + ((uint64_t)op[11] * 4360416026152177361UL)) * 18446744073709551600);
	tmp_q[3] = ((uint64_t)op[0] * 4360416026152177361UL) + ((uint64_t)op[1] * 276319644861830145UL) + ((uint64_t)op[2] * 7773059531132371638UL) + ((uint64_t)op[3] * 11006936700485570137UL) + ((((uint64_t)op[4] * 11426515914107012300UL) + ((uint64_t)op[5] * 5619274444686996823UL) + ((uint64_t)op[6] * 16241926334950006912UL) + ((uint64_t)op[7] * 1388654943435253604UL) + ((uint64_t)op[8] * 1453109792944704396UL) + ((uint64_t)op[9] * 12710492834476651738UL) + ((uint64_t)op[10] * 989886695836914904UL) + ((uint64_t)op[11] * 6851573905507941151UL)) * 18446744073709551600);
	tmp_q[4] = ((uint64_t)op[0] * 6851573905507941151UL) + ((uint64_t)op[1] * 4360416026152177361UL) + ((uint64_t)op[2] * 276319644861830145UL) + ((uint64_t)op[3] * 7773059531132371638UL) + ((uint64_t)op[4] * 11006936700485570137UL) + ((((uint64_t)op[5] * 11426515914107012300UL) + ((uint64_t)op[6] * 5619274444686996823UL) + ((uint64_t)op[7] * 16241926334950006912UL) + ((uint64_t)op[8] * 1388654943435253604UL) + ((uint64_t)op[9] * 1453109792944704396UL) + ((uint64_t)op[10] * 12710492834476651738UL) + ((uint64_t)op[11] * 989886695836914904UL)) * 18446744073709551600);
	tmp_q[5] = ((uint64_t)op[0] * 989886695836914904UL) + ((uint64_t)op[1] * 6851573905507941151UL) + ((uint64_t)op[2] * 4360416026152177361UL) + ((uint64_t)op[3] * 276319644861830145UL) + ((uint64_t)op[4] * 7773059531132371638UL) + ((uint64_t)op[5] * 11006936700485570137UL) + ((((uint64_t)op[6] * 11426515914107012300UL) + ((uint64_t)op[7] * 5619274444686996823UL) + ((uint64_t)op[8] * 16241926334950006912UL) + ((uint64_t)op[9] * 1388654943435253604UL) + ((uint64_t)op[10] * 1453109792944704396UL) + ((uint64_t)op[11] * 12710492834476651738UL)) * 18446744073709551600);
	tmp_q[6] = ((uint64_t)op[0] * 12710492834476651738UL) + ((uint64_t)op[1] * 989886695836914904UL) + ((uint64_t)op[2] * 6851573905507941151UL) + ((uint64_t)op[3] * 4360416026152177361UL) + ((uint64_t)op[4] * 276319644861830145UL) + ((uint64_t)op[5] * 7773059531132371638UL) + ((uint64_t)op[6] * 11006936700485570137UL) + ((((uint64_t)op[7] * 11426515914107012300UL) + ((uint64_t)op[8] * 5619274444686996823UL) + ((uint64_t)op[9] * 16241926334950006912UL) + ((uint64_t)op[10] * 1388654943435253604UL) + ((uint64_t)op[11] * 1453109792944704396UL)) * 18446744073709551600);
	tmp_q[7] = ((uint64_t)op[0] * 1453109792944704396UL) + ((uint64_t)op[1] * 12710492834476651738UL) + ((uint64_t)op[2] * 989886695836914904UL) + ((uint64_t)op[3] * 6851573905507941151UL) + ((uint64_t)op[4] * 4360416026152177361UL) + ((uint64_t)op[5] * 276319644861830145UL) + ((uint64_t)op[6] * 7773059531132371638UL) + ((uint64_t)op[7] * 11006936700485570137UL) + ((((uint64_t)op[8] * 11426515914107012300UL) + ((uint64_t)op[9] * 5619274444686996823UL) + ((uint64_t)op[10] * 16241926334950006912UL) + ((uint64_t)op[11] * 1388654943435253604UL)) * 18446744073709551600);
	tmp_q[8] = ((uint64_t)op[0] * 1388654943435253604UL) + ((uint64_t)op[1] * 1453109792944704396UL) + ((uint64_t)op[2] * 12710492834476651738UL) + ((uint64_t)op[3] * 989886695836914904UL) + ((uint64_t)op[4] * 6851573905507941151UL) + ((uint64_t)op[5] * 4360416026152177361UL) + ((uint64_t)op[6] * 276319644861830145UL) + ((uint64_t)op[7] * 7773059531132371638UL) + ((uint64_t)op[8] * 11006936700485570137UL) + ((((uint64_t)op[9] * 11426515914107012300UL) + ((uint64_t)op[10] * 5619274444686996823UL) + ((uint64_t)op[11] * 16241926334950006912UL)) * 18446744073709551600);
	tmp_q[9] = ((uint64_t)op[0] * 16241926334950006912UL) + ((uint64_t)op[1] * 1388654943435253604UL) + ((uint64_t)op[2] * 1453109792944704396UL) + ((uint64_t)op[3] * 12710492834476651738UL) + ((uint64_t)op[4] * 989886695836914904UL) + ((uint64_t)op[5] * 6851573905507941151UL) + ((uint64_t)op[6] * 4360416026152177361UL) + ((uint64_t)op[7] * 276319644861830145UL) + ((uint64_t)op[8] * 7773059531132371638UL) + ((uint64_t)op[9] * 11006936700485570137UL) + ((((uint64_t)op[10] * 11426515914107012300UL) + ((uint64_t)op[11] * 5619274444686996823UL)) * 18446744073709551600);
	tmp_q[10] = ((uint64_t)op[0] * 5619274444686996823UL) + ((uint64_t)op[1] * 16241926334950006912UL) + ((uint64_t)op[2] * 1388654943435253604UL) + ((uint64_t)op[3] * 1453109792944704396UL) + ((uint64_t)op[4] * 12710492834476651738UL) + ((uint64_t)op[5] * 989886695836914904UL) + ((uint64_t)op[6] * 6851573905507941151UL) + ((uint64_t)op[7] * 4360416026152177361UL) + ((uint64_t)op[8] * 276319644861830145UL) + ((uint64_t)op[9] * 7773059531132371638UL) + ((uint64_t)op[10] * 11006936700485570137UL) + ((uint64_t)op[11] * 1643186111383319360UL);
	tmp_q[11] = ((uint64_t)op[0] * 11426515914107012300UL) + ((uint64_t)op[1] * 5619274444686996823UL) + ((uint64_t)op[2] * 16241926334950006912UL) + ((uint64_t)op[3] * 1388654943435253604UL) + ((uint64_t)op[4] * 1453109792944704396UL) + ((uint64_t)op[5] * 12710492834476651738UL) + ((uint64_t)op[6] * 989886695836914904UL) + ((uint64_t)op[7] * 6851573905507941151UL) + ((uint64_t)op[8] * 4360416026152177361UL) + ((uint64_t)op[9] * 276319644861830145UL) + ((uint64_t)op[10] * 7773059531132371638UL) + ((uint64_t)op[11] * 11006936700485570137UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 5847793091049L) - ((-((int128)tmp_q[1] * 4265963647408L) + ((int128)tmp_q[2] * 144388755532L) - ((int128)tmp_q[3] * 1053680603469L) - ((int128)tmp_q[4] * 1360444911798L) + ((int128)tmp_q[5] * 4942555107807L) + ((int128)tmp_q[6] * 6323514805348L) - ((int128)tmp_q[7] * 3003693672080L) - ((int128)tmp_q[8] * 1479492219354L) + ((int128)tmp_q[9] * 4738686363117L) - ((int128)tmp_q[10] * 1622921654739L) + ((int128)tmp_q[11] * 1125950440566L)) * 16);
	tmp_zero[1] = ((int128)tmp_q[0] * 1125950440566L) - ((int128)tmp_q[1] * 5847793091049L) - ((-((int128)tmp_q[2] * 4265963647408L) + ((int128)tmp_q[3] * 144388755532L) - ((int128)tmp_q[4] * 1053680603469L) - ((int128)tmp_q[5] * 1360444911798L) + ((int128)tmp_q[6] * 4942555107807L) + ((int128)tmp_q[7] * 6323514805348L) - ((int128)tmp_q[8] * 3003693672080L) - ((int128)tmp_q[9] * 1479492219354L) + ((int128)tmp_q[10] * 4738686363117L) - ((int128)tmp_q[11] * 1622921654739L)) * 16);
	tmp_zero[2] = -((int128)tmp_q[0] * 1622921654739L) + ((int128)tmp_q[1] * 1125950440566L) - ((int128)tmp_q[2] * 5847793091049L) - ((-((int128)tmp_q[3] * 4265963647408L) + ((int128)tmp_q[4] * 144388755532L) - ((int128)tmp_q[5] * 1053680603469L) - ((int128)tmp_q[6] * 1360444911798L) + ((int128)tmp_q[7] * 4942555107807L) + ((int128)tmp_q[8] * 6323514805348L) - ((int128)tmp_q[9] * 3003693672080L) - ((int128)tmp_q[10] * 1479492219354L) + ((int128)tmp_q[11] * 4738686363117L)) * 16);
	tmp_zero[3] = ((int128)tmp_q[0] * 4738686363117L) - ((int128)tmp_q[1] * 1622921654739L) + ((int128)tmp_q[2] * 1125950440566L) - ((int128)tmp_q[3] * 5847793091049L) - ((-((int128)tmp_q[4] * 4265963647408L) + ((int128)tmp_q[5] * 144388755532L) - ((int128)tmp_q[6] * 1053680603469L) - ((int128)tmp_q[7] * 1360444911798L) + ((int128)tmp_q[8] * 4942555107807L) + ((int128)tmp_q[9] * 6323514805348L) - ((int128)tmp_q[10] * 3003693672080L) - ((int128)tmp_q[11] * 1479492219354L)) * 16);
	tmp_zero[4] = -((int128)tmp_q[0] * 1479492219354L) + ((int128)tmp_q[1] * 4738686363117L) - ((int128)tmp_q[2] * 1622921654739L) + ((int128)tmp_q[3] * 1125950440566L) - ((int128)tmp_q[4] * 5847793091049L) - ((-((int128)tmp_q[5] * 4265963647408L) + ((int128)tmp_q[6] * 144388755532L) - ((int128)tmp_q[7] * 1053680603469L) - ((int128)tmp_q[8] * 1360444911798L) + ((int128)tmp_q[9] * 4942555107807L) + ((int128)tmp_q[10] * 6323514805348L) - ((int128)tmp_q[11] * 3003693672080L)) * 16);
	tmp_zero[5] = -((int128)tmp_q[0] * 3003693672080L) - ((int128)tmp_q[1] * 1479492219354L) + ((int128)tmp_q[2] * 4738686363117L) - ((int128)tmp_q[3] * 1622921654739L) + ((int128)tmp_q[4] * 1125950440566L) - ((int128)tmp_q[5] * 5847793091049L) - ((-((int128)tmp_q[6] * 4265963647408L) + ((int128)tmp_q[7] * 144388755532L) - ((int128)tmp_q[8] * 1053680603469L) - ((int128)tmp_q[9] * 1360444911798L) + ((int128)tmp_q[10] * 4942555107807L) + ((int128)tmp_q[11] * 6323514805348L)) * 16);
	tmp_zero[6] = ((int128)tmp_q[0] * 6323514805348L) - ((int128)tmp_q[1] * 3003693672080L) - ((int128)tmp_q[2] * 1479492219354L) + ((int128)tmp_q[3] * 4738686363117L) - ((int128)tmp_q[4] * 1622921654739L) + ((int128)tmp_q[5] * 1125950440566L) - ((int128)tmp_q[6] * 5847793091049L) - ((-((int128)tmp_q[7] * 4265963647408L) + ((int128)tmp_q[8] * 144388755532L) - ((int128)tmp_q[9] * 1053680603469L) - ((int128)tmp_q[10] * 1360444911798L) + ((int128)tmp_q[11] * 4942555107807L)) * 16);
	tmp_zero[7] = ((int128)tmp_q[0] * 4942555107807L) + ((int128)tmp_q[1] * 6323514805348L) - ((int128)tmp_q[2] * 3003693672080L) - ((int128)tmp_q[3] * 1479492219354L) + ((int128)tmp_q[4] * 4738686363117L) - ((int128)tmp_q[5] * 1622921654739L) + ((int128)tmp_q[6] * 1125950440566L) - ((int128)tmp_q[7] * 5847793091049L) - ((-((int128)tmp_q[8] * 4265963647408L) + ((int128)tmp_q[9] * 144388755532L) - ((int128)tmp_q[10] * 1053680603469L) - ((int128)tmp_q[11] * 1360444911798L)) * 16);
	tmp_zero[8] = -((int128)tmp_q[0] * 1360444911798L) + ((int128)tmp_q[1] * 4942555107807L) + ((int128)tmp_q[2] * 6323514805348L) - ((int128)tmp_q[3] * 3003693672080L) - ((int128)tmp_q[4] * 1479492219354L) + ((int128)tmp_q[5] * 4738686363117L) - ((int128)tmp_q[6] * 1622921654739L) + ((int128)tmp_q[7] * 1125950440566L) - ((int128)tmp_q[8] * 5847793091049L) - ((-((int128)tmp_q[9] * 4265963647408L) + ((int128)tmp_q[10] * 144388755532L) - ((int128)tmp_q[11] * 1053680603469L)) * 16);
	tmp_zero[9] = -((int128)tmp_q[0] * 1053680603469L) - ((int128)tmp_q[1] * 1360444911798L) + ((int128)tmp_q[2] * 4942555107807L) + ((int128)tmp_q[3] * 6323514805348L) - ((int128)tmp_q[4] * 3003693672080L) - ((int128)tmp_q[5] * 1479492219354L) + ((int128)tmp_q[6] * 4738686363117L) - ((int128)tmp_q[7] * 1622921654739L) + ((int128)tmp_q[8] * 1125950440566L) - ((int128)tmp_q[9] * 5847793091049L) - ((-((int128)tmp_q[10] * 4265963647408L) + ((int128)tmp_q[11] * 144388755532L)) * 16);
	tmp_zero[10] = ((int128)tmp_q[0] * 144388755532L) - ((int128)tmp_q[1] * 1053680603469L) - ((int128)tmp_q[2] * 1360444911798L) + ((int128)tmp_q[3] * 4942555107807L) + ((int128)tmp_q[4] * 6323514805348L) - ((int128)tmp_q[5] * 3003693672080L) - ((int128)tmp_q[6] * 1479492219354L) + ((int128)tmp_q[7] * 4738686363117L) - ((int128)tmp_q[8] * 1622921654739L) + ((int128)tmp_q[9] * 1125950440566L) - ((int128)tmp_q[10] * 5847793091049L) + ((int128)tmp_q[11] * 68255418358528L);
	tmp_zero[11] = -((int128)tmp_q[0] * 4265963647408L) + ((int128)tmp_q[1] * 144388755532L) - ((int128)tmp_q[2] * 1053680603469L) - ((int128)tmp_q[3] * 1360444911798L) + ((int128)tmp_q[4] * 4942555107807L) + ((int128)tmp_q[5] * 6323514805348L) - ((int128)tmp_q[6] * 3003693672080L) - ((int128)tmp_q[7] * 1479492219354L) + ((int128)tmp_q[8] * 4738686363117L) - ((int128)tmp_q[9] * 1622921654739L) + ((int128)tmp_q[10] * 1125950440566L) - ((int128)tmp_q[11] * 5847793091049L);

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

