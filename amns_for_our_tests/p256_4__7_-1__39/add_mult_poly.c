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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - ((int128)pa[6] * pb[6]);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - ((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - ((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - ((int128)pa[6] * pa[6]);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9097175789421631260UL) + ((((uint64_t)op[1] * 15163691219178405610UL) + ((uint64_t)op[2] * 7781265775120356656UL) + ((uint64_t)op[3] * 2653727700720135763UL) + ((uint64_t)op[4] * 1163650877949278016UL) + ((uint64_t)op[5] * 16592940551822490080UL) + ((uint64_t)op[6] * 16115874534363476424UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 16115874534363476424UL) + ((uint64_t)op[1] * 9097175789421631260UL) + ((((uint64_t)op[2] * 15163691219178405610UL) + ((uint64_t)op[3] * 7781265775120356656UL) + ((uint64_t)op[4] * 2653727700720135763UL) + ((uint64_t)op[5] * 1163650877949278016UL) + ((uint64_t)op[6] * 16592940551822490080UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 16592940551822490080UL) + ((uint64_t)op[1] * 16115874534363476424UL) + ((uint64_t)op[2] * 9097175789421631260UL) + ((((uint64_t)op[3] * 15163691219178405610UL) + ((uint64_t)op[4] * 7781265775120356656UL) + ((uint64_t)op[5] * 2653727700720135763UL) + ((uint64_t)op[6] * 1163650877949278016UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 1163650877949278016UL) + ((uint64_t)op[1] * 16592940551822490080UL) + ((uint64_t)op[2] * 16115874534363476424UL) + ((uint64_t)op[3] * 9097175789421631260UL) + ((((uint64_t)op[4] * 15163691219178405610UL) + ((uint64_t)op[5] * 7781265775120356656UL) + ((uint64_t)op[6] * 2653727700720135763UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 2653727700720135763UL) + ((uint64_t)op[1] * 1163650877949278016UL) + ((uint64_t)op[2] * 16592940551822490080UL) + ((uint64_t)op[3] * 16115874534363476424UL) + ((uint64_t)op[4] * 9097175789421631260UL) + ((((uint64_t)op[5] * 15163691219178405610UL) + ((uint64_t)op[6] * 7781265775120356656UL)) * 18446744073709551615);
	tmp_q[5] = ((uint64_t)op[0] * 7781265775120356656UL) + ((uint64_t)op[1] * 2653727700720135763UL) + ((uint64_t)op[2] * 1163650877949278016UL) + ((uint64_t)op[3] * 16592940551822490080UL) + ((uint64_t)op[4] * 16115874534363476424UL) + ((uint64_t)op[5] * 9097175789421631260UL) + ((uint64_t)op[6] * 3283052854531146006UL);
	tmp_q[6] = ((uint64_t)op[0] * 15163691219178405610UL) + ((uint64_t)op[1] * 7781265775120356656UL) + ((uint64_t)op[2] * 2653727700720135763UL) + ((uint64_t)op[3] * 1163650877949278016UL) + ((uint64_t)op[4] * 16592940551822490080UL) + ((uint64_t)op[5] * 16115874534363476424UL) + ((uint64_t)op[6] * 9097175789421631260UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3098690330764L) - (((int128)tmp_q[1] * 620763150620L) - ((int128)tmp_q[2] * 2553712912282L) - ((int128)tmp_q[3] * 728708749728L) - ((int128)tmp_q[4] * 1027952052021L) - ((int128)tmp_q[5] * 4377094995368L) + ((int128)tmp_q[6] * 2195314700592L));
	tmp_zero[1] = ((int128)tmp_q[0] * 2195314700592L) + ((int128)tmp_q[1] * 3098690330764L) - (((int128)tmp_q[2] * 620763150620L) - ((int128)tmp_q[3] * 2553712912282L) - ((int128)tmp_q[4] * 728708749728L) - ((int128)tmp_q[5] * 1027952052021L) - ((int128)tmp_q[6] * 4377094995368L));
	tmp_zero[2] = -((int128)tmp_q[0] * 4377094995368L) + ((int128)tmp_q[1] * 2195314700592L) + ((int128)tmp_q[2] * 3098690330764L) - (((int128)tmp_q[3] * 620763150620L) - ((int128)tmp_q[4] * 2553712912282L) - ((int128)tmp_q[5] * 728708749728L) - ((int128)tmp_q[6] * 1027952052021L));
	tmp_zero[3] = -((int128)tmp_q[0] * 1027952052021L) - ((int128)tmp_q[1] * 4377094995368L) + ((int128)tmp_q[2] * 2195314700592L) + ((int128)tmp_q[3] * 3098690330764L) - (((int128)tmp_q[4] * 620763150620L) - ((int128)tmp_q[5] * 2553712912282L) - ((int128)tmp_q[6] * 728708749728L));
	tmp_zero[4] = -((int128)tmp_q[0] * 728708749728L) - ((int128)tmp_q[1] * 1027952052021L) - ((int128)tmp_q[2] * 4377094995368L) + ((int128)tmp_q[3] * 2195314700592L) + ((int128)tmp_q[4] * 3098690330764L) - (((int128)tmp_q[5] * 620763150620L) - ((int128)tmp_q[6] * 2553712912282L));
	tmp_zero[5] = -((int128)tmp_q[0] * 2553712912282L) - ((int128)tmp_q[1] * 728708749728L) - ((int128)tmp_q[2] * 1027952052021L) - ((int128)tmp_q[3] * 4377094995368L) + ((int128)tmp_q[4] * 2195314700592L) + ((int128)tmp_q[5] * 3098690330764L) - ((int128)tmp_q[6] * 620763150620L);
	tmp_zero[6] = ((int128)tmp_q[0] * 620763150620L) - ((int128)tmp_q[1] * 2553712912282L) - ((int128)tmp_q[2] * 728708749728L) - ((int128)tmp_q[3] * 1027952052021L) - ((int128)tmp_q[4] * 4377094995368L) + ((int128)tmp_q[5] * 2195314700592L) + ((int128)tmp_q[6] * 3098690330764L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

