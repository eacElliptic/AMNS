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
	tmp_q[0] = ((uint64_t)op[0] * 8203490476046241450UL) + ((((uint64_t)op[1] * 17160412308995579215UL) + ((uint64_t)op[2] * 11882851899042216395UL) + ((uint64_t)op[3] * 11484644786867148233UL) + ((uint64_t)op[4] * 9822345591995209869UL) + ((uint64_t)op[5] * 16048266049312751605UL) + ((uint64_t)op[6] * 7368009420086348896UL) + ((uint64_t)op[7] * 18174958489028537175UL) + ((uint64_t)op[8] * 6766332347828823146UL) + ((uint64_t)op[9] * 17310486234607353687UL) + ((uint64_t)op[10] * 560719363052726780UL) + ((uint64_t)op[11] * 5328883357746293756UL)) * 18446744073709551607);
	tmp_q[1] = ((uint64_t)op[0] * 5328883357746293756UL) + ((uint64_t)op[1] * 8203490476046241450UL) + ((((uint64_t)op[2] * 17160412308995579215UL) + ((uint64_t)op[3] * 11882851899042216395UL) + ((uint64_t)op[4] * 11484644786867148233UL) + ((uint64_t)op[5] * 9822345591995209869UL) + ((uint64_t)op[6] * 16048266049312751605UL) + ((uint64_t)op[7] * 7368009420086348896UL) + ((uint64_t)op[8] * 18174958489028537175UL) + ((uint64_t)op[9] * 6766332347828823146UL) + ((uint64_t)op[10] * 17310486234607353687UL) + ((uint64_t)op[11] * 560719363052726780UL)) * 18446744073709551607);
	tmp_q[2] = ((uint64_t)op[0] * 560719363052726780UL) + ((uint64_t)op[1] * 5328883357746293756UL) + ((uint64_t)op[2] * 8203490476046241450UL) + ((((uint64_t)op[3] * 17160412308995579215UL) + ((uint64_t)op[4] * 11882851899042216395UL) + ((uint64_t)op[5] * 11484644786867148233UL) + ((uint64_t)op[6] * 9822345591995209869UL) + ((uint64_t)op[7] * 16048266049312751605UL) + ((uint64_t)op[8] * 7368009420086348896UL) + ((uint64_t)op[9] * 18174958489028537175UL) + ((uint64_t)op[10] * 6766332347828823146UL) + ((uint64_t)op[11] * 17310486234607353687UL)) * 18446744073709551607);
	tmp_q[3] = ((uint64_t)op[0] * 17310486234607353687UL) + ((uint64_t)op[1] * 560719363052726780UL) + ((uint64_t)op[2] * 5328883357746293756UL) + ((uint64_t)op[3] * 8203490476046241450UL) + ((((uint64_t)op[4] * 17160412308995579215UL) + ((uint64_t)op[5] * 11882851899042216395UL) + ((uint64_t)op[6] * 11484644786867148233UL) + ((uint64_t)op[7] * 9822345591995209869UL) + ((uint64_t)op[8] * 16048266049312751605UL) + ((uint64_t)op[9] * 7368009420086348896UL) + ((uint64_t)op[10] * 18174958489028537175UL) + ((uint64_t)op[11] * 6766332347828823146UL)) * 18446744073709551607);
	tmp_q[4] = ((uint64_t)op[0] * 6766332347828823146UL) + ((uint64_t)op[1] * 17310486234607353687UL) + ((uint64_t)op[2] * 560719363052726780UL) + ((uint64_t)op[3] * 5328883357746293756UL) + ((uint64_t)op[4] * 8203490476046241450UL) + ((((uint64_t)op[5] * 17160412308995579215UL) + ((uint64_t)op[6] * 11882851899042216395UL) + ((uint64_t)op[7] * 11484644786867148233UL) + ((uint64_t)op[8] * 9822345591995209869UL) + ((uint64_t)op[9] * 16048266049312751605UL) + ((uint64_t)op[10] * 7368009420086348896UL) + ((uint64_t)op[11] * 18174958489028537175UL)) * 18446744073709551607);
	tmp_q[5] = ((uint64_t)op[0] * 18174958489028537175UL) + ((uint64_t)op[1] * 6766332347828823146UL) + ((uint64_t)op[2] * 17310486234607353687UL) + ((uint64_t)op[3] * 560719363052726780UL) + ((uint64_t)op[4] * 5328883357746293756UL) + ((uint64_t)op[5] * 8203490476046241450UL) + ((((uint64_t)op[6] * 17160412308995579215UL) + ((uint64_t)op[7] * 11882851899042216395UL) + ((uint64_t)op[8] * 11484644786867148233UL) + ((uint64_t)op[9] * 9822345591995209869UL) + ((uint64_t)op[10] * 16048266049312751605UL) + ((uint64_t)op[11] * 7368009420086348896UL)) * 18446744073709551607);
	tmp_q[6] = ((uint64_t)op[0] * 7368009420086348896UL) + ((uint64_t)op[1] * 18174958489028537175UL) + ((uint64_t)op[2] * 6766332347828823146UL) + ((uint64_t)op[3] * 17310486234607353687UL) + ((uint64_t)op[4] * 560719363052726780UL) + ((uint64_t)op[5] * 5328883357746293756UL) + ((uint64_t)op[6] * 8203490476046241450UL) + ((((uint64_t)op[7] * 17160412308995579215UL) + ((uint64_t)op[8] * 11882851899042216395UL) + ((uint64_t)op[9] * 11484644786867148233UL) + ((uint64_t)op[10] * 9822345591995209869UL) + ((uint64_t)op[11] * 16048266049312751605UL)) * 18446744073709551607);
	tmp_q[7] = ((uint64_t)op[0] * 16048266049312751605UL) + ((uint64_t)op[1] * 7368009420086348896UL) + ((uint64_t)op[2] * 18174958489028537175UL) + ((uint64_t)op[3] * 6766332347828823146UL) + ((uint64_t)op[4] * 17310486234607353687UL) + ((uint64_t)op[5] * 560719363052726780UL) + ((uint64_t)op[6] * 5328883357746293756UL) + ((uint64_t)op[7] * 8203490476046241450UL) + ((((uint64_t)op[8] * 17160412308995579215UL) + ((uint64_t)op[9] * 11882851899042216395UL) + ((uint64_t)op[10] * 11484644786867148233UL) + ((uint64_t)op[11] * 9822345591995209869UL)) * 18446744073709551607);
	tmp_q[8] = ((uint64_t)op[0] * 9822345591995209869UL) + ((uint64_t)op[1] * 16048266049312751605UL) + ((uint64_t)op[2] * 7368009420086348896UL) + ((uint64_t)op[3] * 18174958489028537175UL) + ((uint64_t)op[4] * 6766332347828823146UL) + ((uint64_t)op[5] * 17310486234607353687UL) + ((uint64_t)op[6] * 560719363052726780UL) + ((uint64_t)op[7] * 5328883357746293756UL) + ((uint64_t)op[8] * 8203490476046241450UL) + ((((uint64_t)op[9] * 17160412308995579215UL) + ((uint64_t)op[10] * 11882851899042216395UL) + ((uint64_t)op[11] * 11484644786867148233UL)) * 18446744073709551607);
	tmp_q[9] = ((uint64_t)op[0] * 11484644786867148233UL) + ((uint64_t)op[1] * 9822345591995209869UL) + ((uint64_t)op[2] * 16048266049312751605UL) + ((uint64_t)op[3] * 7368009420086348896UL) + ((uint64_t)op[4] * 18174958489028537175UL) + ((uint64_t)op[5] * 6766332347828823146UL) + ((uint64_t)op[6] * 17310486234607353687UL) + ((uint64_t)op[7] * 560719363052726780UL) + ((uint64_t)op[8] * 5328883357746293756UL) + ((uint64_t)op[9] * 8203490476046241450UL) + ((((uint64_t)op[10] * 17160412308995579215UL) + ((uint64_t)op[11] * 11882851899042216395UL)) * 18446744073709551607);
	tmp_q[10] = ((uint64_t)op[0] * 11882851899042216395UL) + ((uint64_t)op[1] * 11484644786867148233UL) + ((uint64_t)op[2] * 9822345591995209869UL) + ((uint64_t)op[3] * 16048266049312751605UL) + ((uint64_t)op[4] * 7368009420086348896UL) + ((uint64_t)op[5] * 18174958489028537175UL) + ((uint64_t)op[6] * 6766332347828823146UL) + ((uint64_t)op[7] * 17310486234607353687UL) + ((uint64_t)op[8] * 560719363052726780UL) + ((uint64_t)op[9] * 5328883357746293756UL) + ((uint64_t)op[10] * 8203490476046241450UL) + ((uint64_t)op[11] * 11576985882425751609UL);
	tmp_q[11] = ((uint64_t)op[0] * 17160412308995579215UL) + ((uint64_t)op[1] * 11882851899042216395UL) + ((uint64_t)op[2] * 11484644786867148233UL) + ((uint64_t)op[3] * 9822345591995209869UL) + ((uint64_t)op[4] * 16048266049312751605UL) + ((uint64_t)op[5] * 7368009420086348896UL) + ((uint64_t)op[6] * 18174958489028537175UL) + ((uint64_t)op[7] * 6766332347828823146UL) + ((uint64_t)op[8] * 17310486234607353687UL) + ((uint64_t)op[9] * 560719363052726780UL) + ((uint64_t)op[10] * 5328883357746293756UL) + ((uint64_t)op[11] * 8203490476046241450UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1449025555520L) - ((((int128)tmp_q[1] * 2905279822727L) - ((int128)tmp_q[2] * 6394742504133L) + ((int128)tmp_q[3] * 5259267357046L) - ((int128)tmp_q[4] * 6547114088201L) - ((int128)tmp_q[5] * 674926074022L) + ((int128)tmp_q[6] * 1329761433245L) - ((int128)tmp_q[7] * 1324343798645L) + ((int128)tmp_q[8] * 4026951854344L) - ((int128)tmp_q[9] * 423288691315L) + ((int128)tmp_q[10] * 2847689439733L) + ((int128)tmp_q[11] * 5411350517120L)) * 9);
	tmp_zero[1] = ((int128)tmp_q[0] * 5411350517120L) - ((int128)tmp_q[1] * 1449025555520L) - ((((int128)tmp_q[2] * 2905279822727L) - ((int128)tmp_q[3] * 6394742504133L) + ((int128)tmp_q[4] * 5259267357046L) - ((int128)tmp_q[5] * 6547114088201L) - ((int128)tmp_q[6] * 674926074022L) + ((int128)tmp_q[7] * 1329761433245L) - ((int128)tmp_q[8] * 1324343798645L) + ((int128)tmp_q[9] * 4026951854344L) - ((int128)tmp_q[10] * 423288691315L) + ((int128)tmp_q[11] * 2847689439733L)) * 9);
	tmp_zero[2] = ((int128)tmp_q[0] * 2847689439733L) + ((int128)tmp_q[1] * 5411350517120L) - ((int128)tmp_q[2] * 1449025555520L) - ((((int128)tmp_q[3] * 2905279822727L) - ((int128)tmp_q[4] * 6394742504133L) + ((int128)tmp_q[5] * 5259267357046L) - ((int128)tmp_q[6] * 6547114088201L) - ((int128)tmp_q[7] * 674926074022L) + ((int128)tmp_q[8] * 1329761433245L) - ((int128)tmp_q[9] * 1324343798645L) + ((int128)tmp_q[10] * 4026951854344L) - ((int128)tmp_q[11] * 423288691315L)) * 9);
	tmp_zero[3] = -((int128)tmp_q[0] * 423288691315L) + ((int128)tmp_q[1] * 2847689439733L) + ((int128)tmp_q[2] * 5411350517120L) - ((int128)tmp_q[3] * 1449025555520L) - ((((int128)tmp_q[4] * 2905279822727L) - ((int128)tmp_q[5] * 6394742504133L) + ((int128)tmp_q[6] * 5259267357046L) - ((int128)tmp_q[7] * 6547114088201L) - ((int128)tmp_q[8] * 674926074022L) + ((int128)tmp_q[9] * 1329761433245L) - ((int128)tmp_q[10] * 1324343798645L) + ((int128)tmp_q[11] * 4026951854344L)) * 9);
	tmp_zero[4] = ((int128)tmp_q[0] * 4026951854344L) - ((int128)tmp_q[1] * 423288691315L) + ((int128)tmp_q[2] * 2847689439733L) + ((int128)tmp_q[3] * 5411350517120L) - ((int128)tmp_q[4] * 1449025555520L) - ((((int128)tmp_q[5] * 2905279822727L) - ((int128)tmp_q[6] * 6394742504133L) + ((int128)tmp_q[7] * 5259267357046L) - ((int128)tmp_q[8] * 6547114088201L) - ((int128)tmp_q[9] * 674926074022L) + ((int128)tmp_q[10] * 1329761433245L) - ((int128)tmp_q[11] * 1324343798645L)) * 9);
	tmp_zero[5] = -((int128)tmp_q[0] * 1324343798645L) + ((int128)tmp_q[1] * 4026951854344L) - ((int128)tmp_q[2] * 423288691315L) + ((int128)tmp_q[3] * 2847689439733L) + ((int128)tmp_q[4] * 5411350517120L) - ((int128)tmp_q[5] * 1449025555520L) - ((((int128)tmp_q[6] * 2905279822727L) - ((int128)tmp_q[7] * 6394742504133L) + ((int128)tmp_q[8] * 5259267357046L) - ((int128)tmp_q[9] * 6547114088201L) - ((int128)tmp_q[10] * 674926074022L) + ((int128)tmp_q[11] * 1329761433245L)) * 9);
	tmp_zero[6] = ((int128)tmp_q[0] * 1329761433245L) - ((int128)tmp_q[1] * 1324343798645L) + ((int128)tmp_q[2] * 4026951854344L) - ((int128)tmp_q[3] * 423288691315L) + ((int128)tmp_q[4] * 2847689439733L) + ((int128)tmp_q[5] * 5411350517120L) - ((int128)tmp_q[6] * 1449025555520L) - ((((int128)tmp_q[7] * 2905279822727L) - ((int128)tmp_q[8] * 6394742504133L) + ((int128)tmp_q[9] * 5259267357046L) - ((int128)tmp_q[10] * 6547114088201L) - ((int128)tmp_q[11] * 674926074022L)) * 9);
	tmp_zero[7] = -((int128)tmp_q[0] * 674926074022L) + ((int128)tmp_q[1] * 1329761433245L) - ((int128)tmp_q[2] * 1324343798645L) + ((int128)tmp_q[3] * 4026951854344L) - ((int128)tmp_q[4] * 423288691315L) + ((int128)tmp_q[5] * 2847689439733L) + ((int128)tmp_q[6] * 5411350517120L) - ((int128)tmp_q[7] * 1449025555520L) - ((((int128)tmp_q[8] * 2905279822727L) - ((int128)tmp_q[9] * 6394742504133L) + ((int128)tmp_q[10] * 5259267357046L) - ((int128)tmp_q[11] * 6547114088201L)) * 9);
	tmp_zero[8] = -((int128)tmp_q[0] * 6547114088201L) - ((int128)tmp_q[1] * 674926074022L) + ((int128)tmp_q[2] * 1329761433245L) - ((int128)tmp_q[3] * 1324343798645L) + ((int128)tmp_q[4] * 4026951854344L) - ((int128)tmp_q[5] * 423288691315L) + ((int128)tmp_q[6] * 2847689439733L) + ((int128)tmp_q[7] * 5411350517120L) - ((int128)tmp_q[8] * 1449025555520L) - ((((int128)tmp_q[9] * 2905279822727L) - ((int128)tmp_q[10] * 6394742504133L) + ((int128)tmp_q[11] * 5259267357046L)) * 9);
	tmp_zero[9] = ((int128)tmp_q[0] * 5259267357046L) - ((int128)tmp_q[1] * 6547114088201L) - ((int128)tmp_q[2] * 674926074022L) + ((int128)tmp_q[3] * 1329761433245L) - ((int128)tmp_q[4] * 1324343798645L) + ((int128)tmp_q[5] * 4026951854344L) - ((int128)tmp_q[6] * 423288691315L) + ((int128)tmp_q[7] * 2847689439733L) + ((int128)tmp_q[8] * 5411350517120L) - ((int128)tmp_q[9] * 1449025555520L) - ((((int128)tmp_q[10] * 2905279822727L) - ((int128)tmp_q[11] * 6394742504133L)) * 9);
	tmp_zero[10] = -((int128)tmp_q[0] * 6394742504133L) + ((int128)tmp_q[1] * 5259267357046L) - ((int128)tmp_q[2] * 6547114088201L) - ((int128)tmp_q[3] * 674926074022L) + ((int128)tmp_q[4] * 1329761433245L) - ((int128)tmp_q[5] * 1324343798645L) + ((int128)tmp_q[6] * 4026951854344L) - ((int128)tmp_q[7] * 423288691315L) + ((int128)tmp_q[8] * 2847689439733L) + ((int128)tmp_q[9] * 5411350517120L) - ((int128)tmp_q[10] * 1449025555520L) - ((int128)tmp_q[11] * 26147518404543L);
	tmp_zero[11] = ((int128)tmp_q[0] * 2905279822727L) - ((int128)tmp_q[1] * 6394742504133L) + ((int128)tmp_q[2] * 5259267357046L) - ((int128)tmp_q[3] * 6547114088201L) - ((int128)tmp_q[4] * 674926074022L) + ((int128)tmp_q[5] * 1329761433245L) - ((int128)tmp_q[6] * 1324343798645L) + ((int128)tmp_q[7] * 4026951854344L) - ((int128)tmp_q[8] * 423288691315L) + ((int128)tmp_q[9] * 2847689439733L) + ((int128)tmp_q[10] * 5411350517120L) - ((int128)tmp_q[11] * 1449025555520L);

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

