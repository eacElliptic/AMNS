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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 6);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 6);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 12);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14271637060425691457UL) + ((((uint64_t)op[1] * 2591291328803453785UL) + ((uint64_t)op[2] * 15249416208318797670UL) + ((uint64_t)op[3] * 14628224601508414051UL) + ((uint64_t)op[4] * 15411742774886212473UL) + ((uint64_t)op[5] * 9458343025437417664UL) + ((uint64_t)op[6] * 15246660488433576457UL) + ((uint64_t)op[7] * 10391159287510716697UL) + ((uint64_t)op[8] * 1183786549986879857UL) + ((uint64_t)op[9] * 16245961435328838072UL) + ((uint64_t)op[10] * 7256247525221878254UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 7256247525221878254UL) + ((uint64_t)op[1] * 14271637060425691457UL) + ((((uint64_t)op[2] * 2591291328803453785UL) + ((uint64_t)op[3] * 15249416208318797670UL) + ((uint64_t)op[4] * 14628224601508414051UL) + ((uint64_t)op[5] * 15411742774886212473UL) + ((uint64_t)op[6] * 9458343025437417664UL) + ((uint64_t)op[7] * 15246660488433576457UL) + ((uint64_t)op[8] * 10391159287510716697UL) + ((uint64_t)op[9] * 1183786549986879857UL) + ((uint64_t)op[10] * 16245961435328838072UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 16245961435328838072UL) + ((uint64_t)op[1] * 7256247525221878254UL) + ((uint64_t)op[2] * 14271637060425691457UL) + ((((uint64_t)op[3] * 2591291328803453785UL) + ((uint64_t)op[4] * 15249416208318797670UL) + ((uint64_t)op[5] * 14628224601508414051UL) + ((uint64_t)op[6] * 15411742774886212473UL) + ((uint64_t)op[7] * 9458343025437417664UL) + ((uint64_t)op[8] * 15246660488433576457UL) + ((uint64_t)op[9] * 10391159287510716697UL) + ((uint64_t)op[10] * 1183786549986879857UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 1183786549986879857UL) + ((uint64_t)op[1] * 16245961435328838072UL) + ((uint64_t)op[2] * 7256247525221878254UL) + ((uint64_t)op[3] * 14271637060425691457UL) + ((((uint64_t)op[4] * 2591291328803453785UL) + ((uint64_t)op[5] * 15249416208318797670UL) + ((uint64_t)op[6] * 14628224601508414051UL) + ((uint64_t)op[7] * 15411742774886212473UL) + ((uint64_t)op[8] * 9458343025437417664UL) + ((uint64_t)op[9] * 15246660488433576457UL) + ((uint64_t)op[10] * 10391159287510716697UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 10391159287510716697UL) + ((uint64_t)op[1] * 1183786549986879857UL) + ((uint64_t)op[2] * 16245961435328838072UL) + ((uint64_t)op[3] * 7256247525221878254UL) + ((uint64_t)op[4] * 14271637060425691457UL) + ((((uint64_t)op[5] * 2591291328803453785UL) + ((uint64_t)op[6] * 15249416208318797670UL) + ((uint64_t)op[7] * 14628224601508414051UL) + ((uint64_t)op[8] * 15411742774886212473UL) + ((uint64_t)op[9] * 9458343025437417664UL) + ((uint64_t)op[10] * 15246660488433576457UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 15246660488433576457UL) + ((uint64_t)op[1] * 10391159287510716697UL) + ((uint64_t)op[2] * 1183786549986879857UL) + ((uint64_t)op[3] * 16245961435328838072UL) + ((uint64_t)op[4] * 7256247525221878254UL) + ((uint64_t)op[5] * 14271637060425691457UL) + ((((uint64_t)op[6] * 2591291328803453785UL) + ((uint64_t)op[7] * 15249416208318797670UL) + ((uint64_t)op[8] * 14628224601508414051UL) + ((uint64_t)op[9] * 15411742774886212473UL) + ((uint64_t)op[10] * 9458343025437417664UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 9458343025437417664UL) + ((uint64_t)op[1] * 15246660488433576457UL) + ((uint64_t)op[2] * 10391159287510716697UL) + ((uint64_t)op[3] * 1183786549986879857UL) + ((uint64_t)op[4] * 16245961435328838072UL) + ((uint64_t)op[5] * 7256247525221878254UL) + ((uint64_t)op[6] * 14271637060425691457UL) + ((((uint64_t)op[7] * 2591291328803453785UL) + ((uint64_t)op[8] * 15249416208318797670UL) + ((uint64_t)op[9] * 14628224601508414051UL) + ((uint64_t)op[10] * 15411742774886212473UL)) * 18446744073709551610);
	tmp_q[7] = ((uint64_t)op[0] * 15411742774886212473UL) + ((uint64_t)op[1] * 9458343025437417664UL) + ((uint64_t)op[2] * 15246660488433576457UL) + ((uint64_t)op[3] * 10391159287510716697UL) + ((uint64_t)op[4] * 1183786549986879857UL) + ((uint64_t)op[5] * 16245961435328838072UL) + ((uint64_t)op[6] * 7256247525221878254UL) + ((uint64_t)op[7] * 14271637060425691457UL) + ((((uint64_t)op[8] * 2591291328803453785UL) + ((uint64_t)op[9] * 15249416208318797670UL) + ((uint64_t)op[10] * 14628224601508414051UL)) * 18446744073709551610);
	tmp_q[8] = ((uint64_t)op[0] * 14628224601508414051UL) + ((uint64_t)op[1] * 15411742774886212473UL) + ((uint64_t)op[2] * 9458343025437417664UL) + ((uint64_t)op[3] * 15246660488433576457UL) + ((uint64_t)op[4] * 10391159287510716697UL) + ((uint64_t)op[5] * 1183786549986879857UL) + ((uint64_t)op[6] * 16245961435328838072UL) + ((uint64_t)op[7] * 7256247525221878254UL) + ((uint64_t)op[8] * 14271637060425691457UL) + ((((uint64_t)op[9] * 2591291328803453785UL) + ((uint64_t)op[10] * 15249416208318797670UL)) * 18446744073709551610);
	tmp_q[9] = ((uint64_t)op[0] * 15249416208318797670UL) + ((uint64_t)op[1] * 14628224601508414051UL) + ((uint64_t)op[2] * 15411742774886212473UL) + ((uint64_t)op[3] * 9458343025437417664UL) + ((uint64_t)op[4] * 15246660488433576457UL) + ((uint64_t)op[5] * 10391159287510716697UL) + ((uint64_t)op[6] * 1183786549986879857UL) + ((uint64_t)op[7] * 16245961435328838072UL) + ((uint64_t)op[8] * 7256247525221878254UL) + ((uint64_t)op[9] * 14271637060425691457UL) + ((uint64_t)op[10] * 2898996100888828906UL);
	tmp_q[10] = ((uint64_t)op[0] * 2591291328803453785UL) + ((uint64_t)op[1] * 15249416208318797670UL) + ((uint64_t)op[2] * 14628224601508414051UL) + ((uint64_t)op[3] * 15411742774886212473UL) + ((uint64_t)op[4] * 9458343025437417664UL) + ((uint64_t)op[5] * 15246660488433576457UL) + ((uint64_t)op[6] * 10391159287510716697UL) + ((uint64_t)op[7] * 1183786549986879857UL) + ((uint64_t)op[8] * 16245961435328838072UL) + ((uint64_t)op[9] * 7256247525221878254UL) + ((uint64_t)op[10] * 14271637060425691457UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 32724806318805L) - ((((int128)tmp_q[1] * 55451179918261L) - ((int128)tmp_q[2] * 65821444144727L) - ((int128)tmp_q[3] * 34390061316752L) - ((int128)tmp_q[4] * 95718679835719L) + ((int128)tmp_q[5] * 77145923305277L) - ((int128)tmp_q[6] * 44214330099425L) + ((int128)tmp_q[7] * 8502397298063L) - ((int128)tmp_q[8] * 107272776163117L) + ((int128)tmp_q[9] * 48205789289162L) - ((int128)tmp_q[10] * 68498386294958L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 68498386294958L) - ((int128)tmp_q[1] * 32724806318805L) - ((((int128)tmp_q[2] * 55451179918261L) - ((int128)tmp_q[3] * 65821444144727L) - ((int128)tmp_q[4] * 34390061316752L) - ((int128)tmp_q[5] * 95718679835719L) + ((int128)tmp_q[6] * 77145923305277L) - ((int128)tmp_q[7] * 44214330099425L) + ((int128)tmp_q[8] * 8502397298063L) - ((int128)tmp_q[9] * 107272776163117L) + ((int128)tmp_q[10] * 48205789289162L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 48205789289162L) - ((int128)tmp_q[1] * 68498386294958L) - ((int128)tmp_q[2] * 32724806318805L) - ((((int128)tmp_q[3] * 55451179918261L) - ((int128)tmp_q[4] * 65821444144727L) - ((int128)tmp_q[5] * 34390061316752L) - ((int128)tmp_q[6] * 95718679835719L) + ((int128)tmp_q[7] * 77145923305277L) - ((int128)tmp_q[8] * 44214330099425L) + ((int128)tmp_q[9] * 8502397298063L) - ((int128)tmp_q[10] * 107272776163117L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 107272776163117L) + ((int128)tmp_q[1] * 48205789289162L) - ((int128)tmp_q[2] * 68498386294958L) - ((int128)tmp_q[3] * 32724806318805L) - ((((int128)tmp_q[4] * 55451179918261L) - ((int128)tmp_q[5] * 65821444144727L) - ((int128)tmp_q[6] * 34390061316752L) - ((int128)tmp_q[7] * 95718679835719L) + ((int128)tmp_q[8] * 77145923305277L) - ((int128)tmp_q[9] * 44214330099425L) + ((int128)tmp_q[10] * 8502397298063L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 8502397298063L) - ((int128)tmp_q[1] * 107272776163117L) + ((int128)tmp_q[2] * 48205789289162L) - ((int128)tmp_q[3] * 68498386294958L) - ((int128)tmp_q[4] * 32724806318805L) - ((((int128)tmp_q[5] * 55451179918261L) - ((int128)tmp_q[6] * 65821444144727L) - ((int128)tmp_q[7] * 34390061316752L) - ((int128)tmp_q[8] * 95718679835719L) + ((int128)tmp_q[9] * 77145923305277L) - ((int128)tmp_q[10] * 44214330099425L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 44214330099425L) + ((int128)tmp_q[1] * 8502397298063L) - ((int128)tmp_q[2] * 107272776163117L) + ((int128)tmp_q[3] * 48205789289162L) - ((int128)tmp_q[4] * 68498386294958L) - ((int128)tmp_q[5] * 32724806318805L) - ((((int128)tmp_q[6] * 55451179918261L) - ((int128)tmp_q[7] * 65821444144727L) - ((int128)tmp_q[8] * 34390061316752L) - ((int128)tmp_q[9] * 95718679835719L) + ((int128)tmp_q[10] * 77145923305277L)) * 6);
	tmp_zero[6] = ((int128)tmp_q[0] * 77145923305277L) - ((int128)tmp_q[1] * 44214330099425L) + ((int128)tmp_q[2] * 8502397298063L) - ((int128)tmp_q[3] * 107272776163117L) + ((int128)tmp_q[4] * 48205789289162L) - ((int128)tmp_q[5] * 68498386294958L) - ((int128)tmp_q[6] * 32724806318805L) - ((((int128)tmp_q[7] * 55451179918261L) - ((int128)tmp_q[8] * 65821444144727L) - ((int128)tmp_q[9] * 34390061316752L) - ((int128)tmp_q[10] * 95718679835719L)) * 6);
	tmp_zero[7] = -((int128)tmp_q[0] * 95718679835719L) + ((int128)tmp_q[1] * 77145923305277L) - ((int128)tmp_q[2] * 44214330099425L) + ((int128)tmp_q[3] * 8502397298063L) - ((int128)tmp_q[4] * 107272776163117L) + ((int128)tmp_q[5] * 48205789289162L) - ((int128)tmp_q[6] * 68498386294958L) - ((int128)tmp_q[7] * 32724806318805L) - ((((int128)tmp_q[8] * 55451179918261L) - ((int128)tmp_q[9] * 65821444144727L) - ((int128)tmp_q[10] * 34390061316752L)) * 6);
	tmp_zero[8] = -((int128)tmp_q[0] * 34390061316752L) - ((int128)tmp_q[1] * 95718679835719L) + ((int128)tmp_q[2] * 77145923305277L) - ((int128)tmp_q[3] * 44214330099425L) + ((int128)tmp_q[4] * 8502397298063L) - ((int128)tmp_q[5] * 107272776163117L) + ((int128)tmp_q[6] * 48205789289162L) - ((int128)tmp_q[7] * 68498386294958L) - ((int128)tmp_q[8] * 32724806318805L) - ((((int128)tmp_q[9] * 55451179918261L) - ((int128)tmp_q[10] * 65821444144727L)) * 6);
	tmp_zero[9] = -((int128)tmp_q[0] * 65821444144727L) - ((int128)tmp_q[1] * 34390061316752L) - ((int128)tmp_q[2] * 95718679835719L) + ((int128)tmp_q[3] * 77145923305277L) - ((int128)tmp_q[4] * 44214330099425L) + ((int128)tmp_q[5] * 8502397298063L) - ((int128)tmp_q[6] * 107272776163117L) + ((int128)tmp_q[7] * 48205789289162L) - ((int128)tmp_q[8] * 68498386294958L) - ((int128)tmp_q[9] * 32724806318805L) - ((int128)tmp_q[10] * 332707079509566L);
	tmp_zero[10] = ((int128)tmp_q[0] * 55451179918261L) - ((int128)tmp_q[1] * 65821444144727L) - ((int128)tmp_q[2] * 34390061316752L) - ((int128)tmp_q[3] * 95718679835719L) + ((int128)tmp_q[4] * 77145923305277L) - ((int128)tmp_q[5] * 44214330099425L) + ((int128)tmp_q[6] * 8502397298063L) - ((int128)tmp_q[7] * 107272776163117L) + ((int128)tmp_q[8] * 48205789289162L) - ((int128)tmp_q[9] * 68498386294958L) - ((int128)tmp_q[10] * 32724806318805L);

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
}

