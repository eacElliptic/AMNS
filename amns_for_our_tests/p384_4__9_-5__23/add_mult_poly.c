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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[8] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7348449846292026884UL) + ((((uint64_t)op[1] * 15536830961903983339UL) + ((uint64_t)op[2] * 10168235157789149452UL) + ((uint64_t)op[3] * 11116162907409091367UL) + ((uint64_t)op[4] * 4572241908450385724UL) + ((uint64_t)op[5] * 8476344114449901833UL) + ((uint64_t)op[6] * 17819273278027349621UL) + ((uint64_t)op[7] * 7480397238287676342UL) + ((uint64_t)op[8] * 2684487481200131803UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 2684487481200131803UL) + ((uint64_t)op[1] * 7348449846292026884UL) + ((((uint64_t)op[2] * 15536830961903983339UL) + ((uint64_t)op[3] * 10168235157789149452UL) + ((uint64_t)op[4] * 11116162907409091367UL) + ((uint64_t)op[5] * 4572241908450385724UL) + ((uint64_t)op[6] * 8476344114449901833UL) + ((uint64_t)op[7] * 17819273278027349621UL) + ((uint64_t)op[8] * 7480397238287676342UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 7480397238287676342UL) + ((uint64_t)op[1] * 2684487481200131803UL) + ((uint64_t)op[2] * 7348449846292026884UL) + ((((uint64_t)op[3] * 15536830961903983339UL) + ((uint64_t)op[4] * 10168235157789149452UL) + ((uint64_t)op[5] * 11116162907409091367UL) + ((uint64_t)op[6] * 4572241908450385724UL) + ((uint64_t)op[7] * 8476344114449901833UL) + ((uint64_t)op[8] * 17819273278027349621UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 17819273278027349621UL) + ((uint64_t)op[1] * 7480397238287676342UL) + ((uint64_t)op[2] * 2684487481200131803UL) + ((uint64_t)op[3] * 7348449846292026884UL) + ((((uint64_t)op[4] * 15536830961903983339UL) + ((uint64_t)op[5] * 10168235157789149452UL) + ((uint64_t)op[6] * 11116162907409091367UL) + ((uint64_t)op[7] * 4572241908450385724UL) + ((uint64_t)op[8] * 8476344114449901833UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 8476344114449901833UL) + ((uint64_t)op[1] * 17819273278027349621UL) + ((uint64_t)op[2] * 7480397238287676342UL) + ((uint64_t)op[3] * 2684487481200131803UL) + ((uint64_t)op[4] * 7348449846292026884UL) + ((((uint64_t)op[5] * 15536830961903983339UL) + ((uint64_t)op[6] * 10168235157789149452UL) + ((uint64_t)op[7] * 11116162907409091367UL) + ((uint64_t)op[8] * 4572241908450385724UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 4572241908450385724UL) + ((uint64_t)op[1] * 8476344114449901833UL) + ((uint64_t)op[2] * 17819273278027349621UL) + ((uint64_t)op[3] * 7480397238287676342UL) + ((uint64_t)op[4] * 2684487481200131803UL) + ((uint64_t)op[5] * 7348449846292026884UL) + ((((uint64_t)op[6] * 15536830961903983339UL) + ((uint64_t)op[7] * 10168235157789149452UL) + ((uint64_t)op[8] * 11116162907409091367UL)) * 18446744073709551611);
	tmp_q[6] = ((uint64_t)op[0] * 11116162907409091367UL) + ((uint64_t)op[1] * 4572241908450385724UL) + ((uint64_t)op[2] * 8476344114449901833UL) + ((uint64_t)op[3] * 17819273278027349621UL) + ((uint64_t)op[4] * 7480397238287676342UL) + ((uint64_t)op[5] * 2684487481200131803UL) + ((uint64_t)op[6] * 7348449846292026884UL) + ((((uint64_t)op[7] * 15536830961903983339UL) + ((uint64_t)op[8] * 10168235157789149452UL)) * 18446744073709551611);
	tmp_q[7] = ((uint64_t)op[0] * 10168235157789149452UL) + ((uint64_t)op[1] * 11116162907409091367UL) + ((uint64_t)op[2] * 4572241908450385724UL) + ((uint64_t)op[3] * 8476344114449901833UL) + ((uint64_t)op[4] * 17819273278027349621UL) + ((uint64_t)op[5] * 7480397238287676342UL) + ((uint64_t)op[6] * 2684487481200131803UL) + ((uint64_t)op[7] * 7348449846292026884UL) + ((uint64_t)op[8] * 14549565559027841385UL);
	tmp_q[8] = ((uint64_t)op[0] * 15536830961903983339UL) + ((uint64_t)op[1] * 10168235157789149452UL) + ((uint64_t)op[2] * 11116162907409091367UL) + ((uint64_t)op[3] * 4572241908450385724UL) + ((uint64_t)op[4] * 8476344114449901833UL) + ((uint64_t)op[5] * 17819273278027349621UL) + ((uint64_t)op[6] * 7480397238287676342UL) + ((uint64_t)op[7] * 2684487481200131803UL) + ((uint64_t)op[8] * 7348449846292026884UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1145512888063L) - ((-((int128)tmp_q[1] * 350407859959L) + ((int128)tmp_q[2] * 640096580557L) - ((int128)tmp_q[3] * 2085376770554L) - ((int128)tmp_q[4] * 3740925448206L) + ((int128)tmp_q[5] * 3296587593667L) - ((int128)tmp_q[6] * 3100790827871L) + ((int128)tmp_q[7] * 543109219935L) + ((int128)tmp_q[8] * 3725615184291L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 3725615184291L) - ((int128)tmp_q[1] * 1145512888063L) - ((-((int128)tmp_q[2] * 350407859959L) + ((int128)tmp_q[3] * 640096580557L) - ((int128)tmp_q[4] * 2085376770554L) - ((int128)tmp_q[5] * 3740925448206L) + ((int128)tmp_q[6] * 3296587593667L) - ((int128)tmp_q[7] * 3100790827871L) + ((int128)tmp_q[8] * 543109219935L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 543109219935L) + ((int128)tmp_q[1] * 3725615184291L) - ((int128)tmp_q[2] * 1145512888063L) - ((-((int128)tmp_q[3] * 350407859959L) + ((int128)tmp_q[4] * 640096580557L) - ((int128)tmp_q[5] * 2085376770554L) - ((int128)tmp_q[6] * 3740925448206L) + ((int128)tmp_q[7] * 3296587593667L) - ((int128)tmp_q[8] * 3100790827871L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 3100790827871L) + ((int128)tmp_q[1] * 543109219935L) + ((int128)tmp_q[2] * 3725615184291L) - ((int128)tmp_q[3] * 1145512888063L) - ((-((int128)tmp_q[4] * 350407859959L) + ((int128)tmp_q[5] * 640096580557L) - ((int128)tmp_q[6] * 2085376770554L) - ((int128)tmp_q[7] * 3740925448206L) + ((int128)tmp_q[8] * 3296587593667L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 3296587593667L) - ((int128)tmp_q[1] * 3100790827871L) + ((int128)tmp_q[2] * 543109219935L) + ((int128)tmp_q[3] * 3725615184291L) - ((int128)tmp_q[4] * 1145512888063L) - ((-((int128)tmp_q[5] * 350407859959L) + ((int128)tmp_q[6] * 640096580557L) - ((int128)tmp_q[7] * 2085376770554L) - ((int128)tmp_q[8] * 3740925448206L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 3740925448206L) + ((int128)tmp_q[1] * 3296587593667L) - ((int128)tmp_q[2] * 3100790827871L) + ((int128)tmp_q[3] * 543109219935L) + ((int128)tmp_q[4] * 3725615184291L) - ((int128)tmp_q[5] * 1145512888063L) - ((-((int128)tmp_q[6] * 350407859959L) + ((int128)tmp_q[7] * 640096580557L) - ((int128)tmp_q[8] * 2085376770554L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 2085376770554L) - ((int128)tmp_q[1] * 3740925448206L) + ((int128)tmp_q[2] * 3296587593667L) - ((int128)tmp_q[3] * 3100790827871L) + ((int128)tmp_q[4] * 543109219935L) + ((int128)tmp_q[5] * 3725615184291L) - ((int128)tmp_q[6] * 1145512888063L) - ((-((int128)tmp_q[7] * 350407859959L) + ((int128)tmp_q[8] * 640096580557L)) * 5);
	tmp_zero[7] = ((int128)tmp_q[0] * 640096580557L) - ((int128)tmp_q[1] * 2085376770554L) - ((int128)tmp_q[2] * 3740925448206L) + ((int128)tmp_q[3] * 3296587593667L) - ((int128)tmp_q[4] * 3100790827871L) + ((int128)tmp_q[5] * 543109219935L) + ((int128)tmp_q[6] * 3725615184291L) - ((int128)tmp_q[7] * 1145512888063L) + ((int128)tmp_q[8] * 1752039299795L);
	tmp_zero[8] = -((int128)tmp_q[0] * 350407859959L) + ((int128)tmp_q[1] * 640096580557L) - ((int128)tmp_q[2] * 2085376770554L) - ((int128)tmp_q[3] * 3740925448206L) + ((int128)tmp_q[4] * 3296587593667L) - ((int128)tmp_q[5] * 3100790827871L) + ((int128)tmp_q[6] * 543109219935L) + ((int128)tmp_q[7] * 3725615184291L) - ((int128)tmp_q[8] * 1145512888063L);

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
}

