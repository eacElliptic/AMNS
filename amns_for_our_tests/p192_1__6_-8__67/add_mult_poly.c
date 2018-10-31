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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13817073792408655875UL) + ((((uint64_t)op[1] * 13660838411095690563UL) + ((uint64_t)op[2] * 983635324110093331UL) + ((uint64_t)op[3] * 17493802813262985780UL) + ((uint64_t)op[4] * 9352330536088158828UL) + ((uint64_t)op[5] * 12844681507364880714UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 12844681507364880714UL) + ((uint64_t)op[1] * 13817073792408655875UL) + ((((uint64_t)op[2] * 13660838411095690563UL) + ((uint64_t)op[3] * 983635324110093331UL) + ((uint64_t)op[4] * 17493802813262985780UL) + ((uint64_t)op[5] * 9352330536088158828UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 9352330536088158828UL) + ((uint64_t)op[1] * 12844681507364880714UL) + ((uint64_t)op[2] * 13817073792408655875UL) + ((((uint64_t)op[3] * 13660838411095690563UL) + ((uint64_t)op[4] * 983635324110093331UL) + ((uint64_t)op[5] * 17493802813262985780UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 17493802813262985780UL) + ((uint64_t)op[1] * 9352330536088158828UL) + ((uint64_t)op[2] * 12844681507364880714UL) + ((uint64_t)op[3] * 13817073792408655875UL) + ((((uint64_t)op[4] * 13660838411095690563UL) + ((uint64_t)op[5] * 983635324110093331UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 983635324110093331UL) + ((uint64_t)op[1] * 17493802813262985780UL) + ((uint64_t)op[2] * 9352330536088158828UL) + ((uint64_t)op[3] * 12844681507364880714UL) + ((uint64_t)op[4] * 13817073792408655875UL) + ((uint64_t)op[5] * 1393757153491785192UL);
	tmp_q[5] = ((uint64_t)op[0] * 13660838411095690563UL) + ((uint64_t)op[1] * 983635324110093331UL) + ((uint64_t)op[2] * 17493802813262985780UL) + ((uint64_t)op[3] * 9352330536088158828UL) + ((uint64_t)op[4] * 12844681507364880714UL) + ((uint64_t)op[5] * 13817073792408655875UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3007266304725L) - ((-((int128)tmp_q[1] * 43776884343305L) - ((int128)tmp_q[2] * 41557792305773L) + ((int128)tmp_q[3] * 13437978417452L) - ((int128)tmp_q[4] * 14764363543432L) + ((int128)tmp_q[5] * 17663210294554L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 17663210294554L) + ((int128)tmp_q[1] * 3007266304725L) - ((-((int128)tmp_q[2] * 43776884343305L) - ((int128)tmp_q[3] * 41557792305773L) + ((int128)tmp_q[4] * 13437978417452L) - ((int128)tmp_q[5] * 14764363543432L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 14764363543432L) + ((int128)tmp_q[1] * 17663210294554L) + ((int128)tmp_q[2] * 3007266304725L) - ((-((int128)tmp_q[3] * 43776884343305L) - ((int128)tmp_q[4] * 41557792305773L) + ((int128)tmp_q[5] * 13437978417452L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 13437978417452L) - ((int128)tmp_q[1] * 14764363543432L) + ((int128)tmp_q[2] * 17663210294554L) + ((int128)tmp_q[3] * 3007266304725L) - ((-((int128)tmp_q[4] * 43776884343305L) - ((int128)tmp_q[5] * 41557792305773L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 41557792305773L) + ((int128)tmp_q[1] * 13437978417452L) - ((int128)tmp_q[2] * 14764363543432L) + ((int128)tmp_q[3] * 17663210294554L) + ((int128)tmp_q[4] * 3007266304725L) + ((int128)tmp_q[5] * 350215074746440L);
	tmp_zero[5] = -((int128)tmp_q[0] * 43776884343305L) - ((int128)tmp_q[1] * 41557792305773L) + ((int128)tmp_q[2] * 13437978417452L) - ((int128)tmp_q[3] * 14764363543432L) + ((int128)tmp_q[4] * 17663210294554L) + ((int128)tmp_q[5] * 3007266304725L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

