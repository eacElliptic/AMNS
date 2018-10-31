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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17916611061853057888UL) + ((((uint64_t)op[1] * 10231162015678562956UL) + ((uint64_t)op[2] * 6643269424425549475UL) + ((uint64_t)op[3] * 3717852105350205698UL) + ((uint64_t)op[4] * 2577202517601382547UL) + ((uint64_t)op[5] * 7740295898300413269UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 7740295898300413269UL) + ((uint64_t)op[1] * 17916611061853057888UL) + ((((uint64_t)op[2] * 10231162015678562956UL) + ((uint64_t)op[3] * 6643269424425549475UL) + ((uint64_t)op[4] * 3717852105350205698UL) + ((uint64_t)op[5] * 2577202517601382547UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 2577202517601382547UL) + ((uint64_t)op[1] * 7740295898300413269UL) + ((uint64_t)op[2] * 17916611061853057888UL) + ((((uint64_t)op[3] * 10231162015678562956UL) + ((uint64_t)op[4] * 6643269424425549475UL) + ((uint64_t)op[5] * 3717852105350205698UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 3717852105350205698UL) + ((uint64_t)op[1] * 2577202517601382547UL) + ((uint64_t)op[2] * 7740295898300413269UL) + ((uint64_t)op[3] * 17916611061853057888UL) + ((((uint64_t)op[4] * 10231162015678562956UL) + ((uint64_t)op[5] * 6643269424425549475UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 6643269424425549475UL) + ((uint64_t)op[1] * 3717852105350205698UL) + ((uint64_t)op[2] * 2577202517601382547UL) + ((uint64_t)op[3] * 7740295898300413269UL) + ((uint64_t)op[4] * 17916611061853057888UL) + ((uint64_t)op[5] * 16277901888621285844UL);
	tmp_q[5] = ((uint64_t)op[0] * 10231162015678562956UL) + ((uint64_t)op[1] * 6643269424425549475UL) + ((uint64_t)op[2] * 3717852105350205698UL) + ((uint64_t)op[3] * 2577202517601382547UL) + ((uint64_t)op[4] * 7740295898300413269UL) + ((uint64_t)op[5] * 17916611061853057888UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 42222613735L) + ((((int128)tmp_q[1] * 22983674472L) + ((int128)tmp_q[2] * 69671335025L) + ((int128)tmp_q[3] * 12836637659L) + ((int128)tmp_q[4] * 71735232144L) - ((int128)tmp_q[5] * 55411179774L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 55411179774L) + ((int128)tmp_q[1] * 42222613735L) + ((((int128)tmp_q[2] * 22983674472L) + ((int128)tmp_q[3] * 69671335025L) + ((int128)tmp_q[4] * 12836637659L) + ((int128)tmp_q[5] * 71735232144L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 71735232144L) - ((int128)tmp_q[1] * 55411179774L) + ((int128)tmp_q[2] * 42222613735L) + ((((int128)tmp_q[3] * 22983674472L) + ((int128)tmp_q[4] * 69671335025L) + ((int128)tmp_q[5] * 12836637659L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 12836637659L) + ((int128)tmp_q[1] * 71735232144L) - ((int128)tmp_q[2] * 55411179774L) + ((int128)tmp_q[3] * 42222613735L) + ((((int128)tmp_q[4] * 22983674472L) + ((int128)tmp_q[5] * 69671335025L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 69671335025L) + ((int128)tmp_q[1] * 12836637659L) + ((int128)tmp_q[2] * 71735232144L) - ((int128)tmp_q[3] * 55411179774L) + ((int128)tmp_q[4] * 42222613735L) + ((int128)tmp_q[5] * 160885721304L);
	tmp_zero[5] = ((int128)tmp_q[0] * 22983674472L) + ((int128)tmp_q[1] * 69671335025L) + ((int128)tmp_q[2] * 12836637659L) + ((int128)tmp_q[3] * 71735232144L) - ((int128)tmp_q[4] * 55411179774L) + ((int128)tmp_q[5] * 42222613735L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

