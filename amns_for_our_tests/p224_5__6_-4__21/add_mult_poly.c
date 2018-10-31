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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 642205319251612793UL) + ((((uint64_t)op[1] * 12897050057890916975UL) + ((uint64_t)op[2] * 4597241507308786362UL) + ((uint64_t)op[3] * 862399274729255717UL) + ((uint64_t)op[4] * 13605469519729058258UL) + ((uint64_t)op[5] * 15281358231284911200UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 15281358231284911200UL) + ((uint64_t)op[1] * 642205319251612793UL) + ((((uint64_t)op[2] * 12897050057890916975UL) + ((uint64_t)op[3] * 4597241507308786362UL) + ((uint64_t)op[4] * 862399274729255717UL) + ((uint64_t)op[5] * 13605469519729058258UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 13605469519729058258UL) + ((uint64_t)op[1] * 15281358231284911200UL) + ((uint64_t)op[2] * 642205319251612793UL) + ((((uint64_t)op[3] * 12897050057890916975UL) + ((uint64_t)op[4] * 4597241507308786362UL) + ((uint64_t)op[5] * 862399274729255717UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 862399274729255717UL) + ((uint64_t)op[1] * 13605469519729058258UL) + ((uint64_t)op[2] * 15281358231284911200UL) + ((uint64_t)op[3] * 642205319251612793UL) + ((((uint64_t)op[4] * 12897050057890916975UL) + ((uint64_t)op[5] * 4597241507308786362UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 4597241507308786362UL) + ((uint64_t)op[1] * 862399274729255717UL) + ((uint64_t)op[2] * 13605469519729058258UL) + ((uint64_t)op[3] * 15281358231284911200UL) + ((uint64_t)op[4] * 642205319251612793UL) + ((uint64_t)op[5] * 3752031989564986948UL);
	tmp_q[5] = ((uint64_t)op[0] * 12897050057890916975UL) + ((uint64_t)op[1] * 4597241507308786362UL) + ((uint64_t)op[2] * 862399274729255717UL) + ((uint64_t)op[3] * 13605469519729058258UL) + ((uint64_t)op[4] * 15281358231284911200UL) + ((uint64_t)op[5] * 642205319251612793UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 56990056523L) - ((((int128)tmp_q[1] * 35901275447L) - ((int128)tmp_q[2] * 64723076974L) - ((int128)tmp_q[3] * 14332937535L) + ((int128)tmp_q[4] * 91311449186L) + ((int128)tmp_q[5] * 62130236576L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 62130236576L) + ((int128)tmp_q[1] * 56990056523L) - ((((int128)tmp_q[2] * 35901275447L) - ((int128)tmp_q[3] * 64723076974L) - ((int128)tmp_q[4] * 14332937535L) + ((int128)tmp_q[5] * 91311449186L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 91311449186L) + ((int128)tmp_q[1] * 62130236576L) + ((int128)tmp_q[2] * 56990056523L) - ((((int128)tmp_q[3] * 35901275447L) - ((int128)tmp_q[4] * 64723076974L) - ((int128)tmp_q[5] * 14332937535L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 14332937535L) + ((int128)tmp_q[1] * 91311449186L) + ((int128)tmp_q[2] * 62130236576L) + ((int128)tmp_q[3] * 56990056523L) - ((((int128)tmp_q[4] * 35901275447L) - ((int128)tmp_q[5] * 64723076974L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 64723076974L) - ((int128)tmp_q[1] * 14332937535L) + ((int128)tmp_q[2] * 91311449186L) + ((int128)tmp_q[3] * 62130236576L) + ((int128)tmp_q[4] * 56990056523L) - ((int128)tmp_q[5] * 143605101788L);
	tmp_zero[5] = ((int128)tmp_q[0] * 35901275447L) - ((int128)tmp_q[1] * 64723076974L) - ((int128)tmp_q[2] * 14332937535L) + ((int128)tmp_q[3] * 91311449186L) + ((int128)tmp_q[4] * 62130236576L) + ((int128)tmp_q[5] * 56990056523L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

