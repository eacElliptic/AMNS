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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - ((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - ((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - ((int128)pa[5] * pa[5]);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18132246096884682486UL) + ((((uint64_t)op[1] * 5440569427598520072UL) + ((uint64_t)op[2] * 875876523328253781UL) + ((uint64_t)op[3] * 17867723223757732254UL) + ((uint64_t)op[4] * 15318773805471025944UL) + ((uint64_t)op[5] * 5048456166675391536UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 5048456166675391536UL) + ((uint64_t)op[1] * 18132246096884682486UL) + ((((uint64_t)op[2] * 5440569427598520072UL) + ((uint64_t)op[3] * 875876523328253781UL) + ((uint64_t)op[4] * 17867723223757732254UL) + ((uint64_t)op[5] * 15318773805471025944UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 15318773805471025944UL) + ((uint64_t)op[1] * 5048456166675391536UL) + ((uint64_t)op[2] * 18132246096884682486UL) + ((((uint64_t)op[3] * 5440569427598520072UL) + ((uint64_t)op[4] * 875876523328253781UL) + ((uint64_t)op[5] * 17867723223757732254UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 17867723223757732254UL) + ((uint64_t)op[1] * 15318773805471025944UL) + ((uint64_t)op[2] * 5048456166675391536UL) + ((uint64_t)op[3] * 18132246096884682486UL) + ((((uint64_t)op[4] * 5440569427598520072UL) + ((uint64_t)op[5] * 875876523328253781UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 875876523328253781UL) + ((uint64_t)op[1] * 17867723223757732254UL) + ((uint64_t)op[2] * 15318773805471025944UL) + ((uint64_t)op[3] * 5048456166675391536UL) + ((uint64_t)op[4] * 18132246096884682486UL) + ((uint64_t)op[5] * 13006174646111031544UL);
	tmp_q[5] = ((uint64_t)op[0] * 5440569427598520072UL) + ((uint64_t)op[1] * 875876523328253781UL) + ((uint64_t)op[2] * 17867723223757732254UL) + ((uint64_t)op[3] * 15318773805471025944UL) + ((uint64_t)op[4] * 5048456166675391536UL) + ((uint64_t)op[5] * 18132246096884682486UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1265765539839688L) - (-((int128)tmp_q[1] * 38015694078157824L) + ((int128)tmp_q[2] * 8522449700161638L) - ((int128)tmp_q[3] * 33997051074781200L) + ((int128)tmp_q[4] * 9788215240001325L) + ((int128)tmp_q[5] * 4018643003376626L));
	tmp_zero[1] = ((int128)tmp_q[0] * 4018643003376626L) + ((int128)tmp_q[1] * 1265765539839688L) - (-((int128)tmp_q[2] * 38015694078157824L) + ((int128)tmp_q[3] * 8522449700161638L) - ((int128)tmp_q[4] * 33997051074781200L) + ((int128)tmp_q[5] * 9788215240001325L));
	tmp_zero[2] = ((int128)tmp_q[0] * 9788215240001325L) + ((int128)tmp_q[1] * 4018643003376626L) + ((int128)tmp_q[2] * 1265765539839688L) - (-((int128)tmp_q[3] * 38015694078157824L) + ((int128)tmp_q[4] * 8522449700161638L) - ((int128)tmp_q[5] * 33997051074781200L));
	tmp_zero[3] = -((int128)tmp_q[0] * 33997051074781200L) + ((int128)tmp_q[1] * 9788215240001325L) + ((int128)tmp_q[2] * 4018643003376626L) + ((int128)tmp_q[3] * 1265765539839688L) - (-((int128)tmp_q[4] * 38015694078157824L) + ((int128)tmp_q[5] * 8522449700161638L));
	tmp_zero[4] = ((int128)tmp_q[0] * 8522449700161638L) - ((int128)tmp_q[1] * 33997051074781200L) + ((int128)tmp_q[2] * 9788215240001325L) + ((int128)tmp_q[3] * 4018643003376626L) + ((int128)tmp_q[4] * 1265765539839688L) + ((int128)tmp_q[5] * 38015694078157824L);
	tmp_zero[5] = -((int128)tmp_q[0] * 38015694078157824L) + ((int128)tmp_q[1] * 8522449700161638L) - ((int128)tmp_q[2] * 33997051074781200L) + ((int128)tmp_q[3] * 9788215240001325L) + ((int128)tmp_q[4] * 4018643003376626L) + ((int128)tmp_q[5] * 1265765539839688L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

