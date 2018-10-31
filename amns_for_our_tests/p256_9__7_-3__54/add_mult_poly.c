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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8718002830747094624UL) + ((((uint64_t)op[1] * 7490886961034760495UL) + ((uint64_t)op[2] * 4924958866171939771UL) + ((uint64_t)op[3] * 18240298273481830261UL) + ((uint64_t)op[4] * 6273530143475699304UL) + ((uint64_t)op[5] * 12982892758676393654UL) + ((uint64_t)op[6] * 12866136014748737450UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 12866136014748737450UL) + ((uint64_t)op[1] * 8718002830747094624UL) + ((((uint64_t)op[2] * 7490886961034760495UL) + ((uint64_t)op[3] * 4924958866171939771UL) + ((uint64_t)op[4] * 18240298273481830261UL) + ((uint64_t)op[5] * 6273530143475699304UL) + ((uint64_t)op[6] * 12982892758676393654UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 12982892758676393654UL) + ((uint64_t)op[1] * 12866136014748737450UL) + ((uint64_t)op[2] * 8718002830747094624UL) + ((((uint64_t)op[3] * 7490886961034760495UL) + ((uint64_t)op[4] * 4924958866171939771UL) + ((uint64_t)op[5] * 18240298273481830261UL) + ((uint64_t)op[6] * 6273530143475699304UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 6273530143475699304UL) + ((uint64_t)op[1] * 12982892758676393654UL) + ((uint64_t)op[2] * 12866136014748737450UL) + ((uint64_t)op[3] * 8718002830747094624UL) + ((((uint64_t)op[4] * 7490886961034760495UL) + ((uint64_t)op[5] * 4924958866171939771UL) + ((uint64_t)op[6] * 18240298273481830261UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 18240298273481830261UL) + ((uint64_t)op[1] * 6273530143475699304UL) + ((uint64_t)op[2] * 12982892758676393654UL) + ((uint64_t)op[3] * 12866136014748737450UL) + ((uint64_t)op[4] * 8718002830747094624UL) + ((((uint64_t)op[5] * 7490886961034760495UL) + ((uint64_t)op[6] * 4924958866171939771UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 4924958866171939771UL) + ((uint64_t)op[1] * 18240298273481830261UL) + ((uint64_t)op[2] * 6273530143475699304UL) + ((uint64_t)op[3] * 12982892758676393654UL) + ((uint64_t)op[4] * 12866136014748737450UL) + ((uint64_t)op[5] * 8718002830747094624UL) + ((uint64_t)op[6] * 14420827264314821747UL);
	tmp_q[6] = ((uint64_t)op[0] * 7490886961034760495UL) + ((uint64_t)op[1] * 4924958866171939771UL) + ((uint64_t)op[2] * 18240298273481830261UL) + ((uint64_t)op[3] * 6273530143475699304UL) + ((uint64_t)op[4] * 12982892758676393654UL) + ((uint64_t)op[5] * 12866136014748737450UL) + ((uint64_t)op[6] * 8718002830747094624UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 34138945060L) - ((((int128)tmp_q[1] * 7806061941L) - ((int128)tmp_q[2] * 27296873441L) + ((int128)tmp_q[3] * 31852772946L) + ((int128)tmp_q[4] * 8607793079L) + ((int128)tmp_q[5] * 33100591733L) - ((int128)tmp_q[6] * 25835466485L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 25835466485L) - ((int128)tmp_q[1] * 34138945060L) - ((((int128)tmp_q[2] * 7806061941L) - ((int128)tmp_q[3] * 27296873441L) + ((int128)tmp_q[4] * 31852772946L) + ((int128)tmp_q[5] * 8607793079L) + ((int128)tmp_q[6] * 33100591733L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 33100591733L) - ((int128)tmp_q[1] * 25835466485L) - ((int128)tmp_q[2] * 34138945060L) - ((((int128)tmp_q[3] * 7806061941L) - ((int128)tmp_q[4] * 27296873441L) + ((int128)tmp_q[5] * 31852772946L) + ((int128)tmp_q[6] * 8607793079L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 8607793079L) + ((int128)tmp_q[1] * 33100591733L) - ((int128)tmp_q[2] * 25835466485L) - ((int128)tmp_q[3] * 34138945060L) - ((((int128)tmp_q[4] * 7806061941L) - ((int128)tmp_q[5] * 27296873441L) + ((int128)tmp_q[6] * 31852772946L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 31852772946L) + ((int128)tmp_q[1] * 8607793079L) + ((int128)tmp_q[2] * 33100591733L) - ((int128)tmp_q[3] * 25835466485L) - ((int128)tmp_q[4] * 34138945060L) - ((((int128)tmp_q[5] * 7806061941L) - ((int128)tmp_q[6] * 27296873441L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 27296873441L) + ((int128)tmp_q[1] * 31852772946L) + ((int128)tmp_q[2] * 8607793079L) + ((int128)tmp_q[3] * 33100591733L) - ((int128)tmp_q[4] * 25835466485L) - ((int128)tmp_q[5] * 34138945060L) - ((int128)tmp_q[6] * 23418185823L);
	tmp_zero[6] = ((int128)tmp_q[0] * 7806061941L) - ((int128)tmp_q[1] * 27296873441L) + ((int128)tmp_q[2] * 31852772946L) + ((int128)tmp_q[3] * 8607793079L) + ((int128)tmp_q[4] * 33100591733L) - ((int128)tmp_q[5] * 25835466485L) - ((int128)tmp_q[6] * 34138945060L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

