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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2537303214894462607UL) + ((((uint64_t)op[1] * 7723399312630549152UL) + ((uint64_t)op[2] * 18126978537761904455UL) + ((uint64_t)op[3] * 10869296744367155141UL) + ((uint64_t)op[4] * 10318861558953384607UL) + ((uint64_t)op[5] * 10083863870967937326UL) + ((uint64_t)op[6] * 12662466292265012026UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 12662466292265012026UL) + ((uint64_t)op[1] * 2537303214894462607UL) + ((((uint64_t)op[2] * 7723399312630549152UL) + ((uint64_t)op[3] * 18126978537761904455UL) + ((uint64_t)op[4] * 10869296744367155141UL) + ((uint64_t)op[5] * 10318861558953384607UL) + ((uint64_t)op[6] * 10083863870967937326UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 10083863870967937326UL) + ((uint64_t)op[1] * 12662466292265012026UL) + ((uint64_t)op[2] * 2537303214894462607UL) + ((((uint64_t)op[3] * 7723399312630549152UL) + ((uint64_t)op[4] * 18126978537761904455UL) + ((uint64_t)op[5] * 10869296744367155141UL) + ((uint64_t)op[6] * 10318861558953384607UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 10318861558953384607UL) + ((uint64_t)op[1] * 10083863870967937326UL) + ((uint64_t)op[2] * 12662466292265012026UL) + ((uint64_t)op[3] * 2537303214894462607UL) + ((((uint64_t)op[4] * 7723399312630549152UL) + ((uint64_t)op[5] * 18126978537761904455UL) + ((uint64_t)op[6] * 10869296744367155141UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 10869296744367155141UL) + ((uint64_t)op[1] * 10318861558953384607UL) + ((uint64_t)op[2] * 10083863870967937326UL) + ((uint64_t)op[3] * 12662466292265012026UL) + ((uint64_t)op[4] * 2537303214894462607UL) + ((((uint64_t)op[5] * 7723399312630549152UL) + ((uint64_t)op[6] * 18126978537761904455UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 18126978537761904455UL) + ((uint64_t)op[1] * 10869296744367155141UL) + ((uint64_t)op[2] * 10318861558953384607UL) + ((uint64_t)op[3] * 10083863870967937326UL) + ((uint64_t)op[4] * 12662466292265012026UL) + ((uint64_t)op[5] * 2537303214894462607UL) + ((uint64_t)op[6] * 15446798625261098304UL);
	tmp_q[6] = ((uint64_t)op[0] * 7723399312630549152UL) + ((uint64_t)op[1] * 18126978537761904455UL) + ((uint64_t)op[2] * 10869296744367155141UL) + ((uint64_t)op[3] * 10318861558953384607UL) + ((uint64_t)op[4] * 10083863870967937326UL) + ((uint64_t)op[5] * 12662466292265012026UL) + ((uint64_t)op[6] * 2537303214894462607UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 8093988352370851L) + ((((int128)tmp_q[1] * 1203168232124289L) + ((int128)tmp_q[2] * 15761391193697267L) - ((int128)tmp_q[3] * 6258541789542303L) + ((int128)tmp_q[4] * 1570674034259535L) - ((int128)tmp_q[5] * 14809568966510552L) + ((int128)tmp_q[6] * 14911090468445604L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 14911090468445604L) - ((int128)tmp_q[1] * 8093988352370851L) + ((((int128)tmp_q[2] * 1203168232124289L) + ((int128)tmp_q[3] * 15761391193697267L) - ((int128)tmp_q[4] * 6258541789542303L) + ((int128)tmp_q[5] * 1570674034259535L) - ((int128)tmp_q[6] * 14809568966510552L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 14809568966510552L) + ((int128)tmp_q[1] * 14911090468445604L) - ((int128)tmp_q[2] * 8093988352370851L) + ((((int128)tmp_q[3] * 1203168232124289L) + ((int128)tmp_q[4] * 15761391193697267L) - ((int128)tmp_q[5] * 6258541789542303L) + ((int128)tmp_q[6] * 1570674034259535L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 1570674034259535L) - ((int128)tmp_q[1] * 14809568966510552L) + ((int128)tmp_q[2] * 14911090468445604L) - ((int128)tmp_q[3] * 8093988352370851L) + ((((int128)tmp_q[4] * 1203168232124289L) + ((int128)tmp_q[5] * 15761391193697267L) - ((int128)tmp_q[6] * 6258541789542303L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 6258541789542303L) + ((int128)tmp_q[1] * 1570674034259535L) - ((int128)tmp_q[2] * 14809568966510552L) + ((int128)tmp_q[3] * 14911090468445604L) - ((int128)tmp_q[4] * 8093988352370851L) + ((((int128)tmp_q[5] * 1203168232124289L) + ((int128)tmp_q[6] * 15761391193697267L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 15761391193697267L) - ((int128)tmp_q[1] * 6258541789542303L) + ((int128)tmp_q[2] * 1570674034259535L) - ((int128)tmp_q[3] * 14809568966510552L) + ((int128)tmp_q[4] * 14911090468445604L) - ((int128)tmp_q[5] * 8093988352370851L) + ((int128)tmp_q[6] * 2406336464248578L);
	tmp_zero[6] = ((int128)tmp_q[0] * 1203168232124289L) + ((int128)tmp_q[1] * 15761391193697267L) - ((int128)tmp_q[2] * 6258541789542303L) + ((int128)tmp_q[3] * 1570674034259535L) - ((int128)tmp_q[4] * 14809568966510552L) + ((int128)tmp_q[5] * 14911090468445604L) - ((int128)tmp_q[6] * 8093988352370851L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

