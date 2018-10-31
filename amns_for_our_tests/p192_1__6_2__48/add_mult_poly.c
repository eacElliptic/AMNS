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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11515849461469159803UL) + ((((uint64_t)op[1] * 7139790015543895368UL) + ((uint64_t)op[2] * 7735855574573029118UL) + ((uint64_t)op[3] * 12745499591907786096UL) + ((uint64_t)op[4] * 5780565205597168959UL) + ((uint64_t)op[5] * 12582968618133923338UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 12582968618133923338UL) + ((uint64_t)op[1] * 11515849461469159803UL) + ((((uint64_t)op[2] * 7139790015543895368UL) + ((uint64_t)op[3] * 7735855574573029118UL) + ((uint64_t)op[4] * 12745499591907786096UL) + ((uint64_t)op[5] * 5780565205597168959UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 5780565205597168959UL) + ((uint64_t)op[1] * 12582968618133923338UL) + ((uint64_t)op[2] * 11515849461469159803UL) + ((((uint64_t)op[3] * 7139790015543895368UL) + ((uint64_t)op[4] * 7735855574573029118UL) + ((uint64_t)op[5] * 12745499591907786096UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 12745499591907786096UL) + ((uint64_t)op[1] * 5780565205597168959UL) + ((uint64_t)op[2] * 12582968618133923338UL) + ((uint64_t)op[3] * 11515849461469159803UL) + ((((uint64_t)op[4] * 7139790015543895368UL) + ((uint64_t)op[5] * 7735855574573029118UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 7735855574573029118UL) + ((uint64_t)op[1] * 12745499591907786096UL) + ((uint64_t)op[2] * 5780565205597168959UL) + ((uint64_t)op[3] * 12582968618133923338UL) + ((uint64_t)op[4] * 11515849461469159803UL) + ((uint64_t)op[5] * 14279580031087790736UL);
	tmp_q[5] = ((uint64_t)op[0] * 7139790015543895368UL) + ((uint64_t)op[1] * 7735855574573029118UL) + ((uint64_t)op[2] * 12745499591907786096UL) + ((uint64_t)op[3] * 5780565205597168959UL) + ((uint64_t)op[4] * 12582968618133923338UL) + ((uint64_t)op[5] * 11515849461469159803UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4041134633L) + ((-((int128)tmp_q[1] * 106741122L) - ((int128)tmp_q[2] * 430505703L) + ((int128)tmp_q[3] * 250321352L) + ((int128)tmp_q[4] * 1501879333L) - ((int128)tmp_q[5] * 530496462L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 530496462L) - ((int128)tmp_q[1] * 4041134633L) + ((-((int128)tmp_q[2] * 106741122L) - ((int128)tmp_q[3] * 430505703L) + ((int128)tmp_q[4] * 250321352L) + ((int128)tmp_q[5] * 1501879333L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 1501879333L) - ((int128)tmp_q[1] * 530496462L) - ((int128)tmp_q[2] * 4041134633L) + ((-((int128)tmp_q[3] * 106741122L) - ((int128)tmp_q[4] * 430505703L) + ((int128)tmp_q[5] * 250321352L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 250321352L) + ((int128)tmp_q[1] * 1501879333L) - ((int128)tmp_q[2] * 530496462L) - ((int128)tmp_q[3] * 4041134633L) + ((-((int128)tmp_q[4] * 106741122L) - ((int128)tmp_q[5] * 430505703L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 430505703L) + ((int128)tmp_q[1] * 250321352L) + ((int128)tmp_q[2] * 1501879333L) - ((int128)tmp_q[3] * 530496462L) - ((int128)tmp_q[4] * 4041134633L) - ((int128)tmp_q[5] * 213482244L);
	tmp_zero[5] = -((int128)tmp_q[0] * 106741122L) - ((int128)tmp_q[1] * 430505703L) + ((int128)tmp_q[2] * 250321352L) + ((int128)tmp_q[3] * 1501879333L) - ((int128)tmp_q[4] * 530496462L) - ((int128)tmp_q[5] * 4041134633L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

