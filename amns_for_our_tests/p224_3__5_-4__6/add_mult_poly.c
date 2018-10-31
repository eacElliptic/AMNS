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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7548658831482033119UL) + ((((uint64_t)op[1] * 7384592812649969774UL) + ((uint64_t)op[2] * 15323491107099124023UL) + ((uint64_t)op[3] * 904864896135321782UL) + ((uint64_t)op[4] * 10751242838372843300UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 10751242838372843300UL) + ((uint64_t)op[1] * 7548658831482033119UL) + ((((uint64_t)op[2] * 7384592812649969774UL) + ((uint64_t)op[3] * 15323491107099124023UL) + ((uint64_t)op[4] * 904864896135321782UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 904864896135321782UL) + ((uint64_t)op[1] * 10751242838372843300UL) + ((uint64_t)op[2] * 7548658831482033119UL) + ((((uint64_t)op[3] * 7384592812649969774UL) + ((uint64_t)op[4] * 15323491107099124023UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 15323491107099124023UL) + ((uint64_t)op[1] * 904864896135321782UL) + ((uint64_t)op[2] * 10751242838372843300UL) + ((uint64_t)op[3] * 7548658831482033119UL) + ((uint64_t)op[4] * 7355116896819224136UL);
	tmp_q[4] = ((uint64_t)op[0] * 7384592812649969774UL) + ((uint64_t)op[1] * 15323491107099124023UL) + ((uint64_t)op[2] * 904864896135321782UL) + ((uint64_t)op[3] * 10751242838372843300UL) + ((uint64_t)op[4] * 7548658831482033119UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1209541147953L) - ((((int128)tmp_q[1] * 7178908109134L) + ((int128)tmp_q[2] * 21664755790511L) - ((int128)tmp_q[3] * 5276591573562L) + ((int128)tmp_q[4] * 1561279631712L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 1561279631712L) + ((int128)tmp_q[1] * 1209541147953L) - ((((int128)tmp_q[2] * 7178908109134L) + ((int128)tmp_q[3] * 21664755790511L) - ((int128)tmp_q[4] * 5276591573562L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 5276591573562L) + ((int128)tmp_q[1] * 1561279631712L) + ((int128)tmp_q[2] * 1209541147953L) - ((((int128)tmp_q[3] * 7178908109134L) + ((int128)tmp_q[4] * 21664755790511L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 21664755790511L) - ((int128)tmp_q[1] * 5276591573562L) + ((int128)tmp_q[2] * 1561279631712L) + ((int128)tmp_q[3] * 1209541147953L) - ((int128)tmp_q[4] * 28715632436536L);
	tmp_zero[4] = ((int128)tmp_q[0] * 7178908109134L) + ((int128)tmp_q[1] * 21664755790511L) - ((int128)tmp_q[2] * 5276591573562L) + ((int128)tmp_q[3] * 1561279631712L) + ((int128)tmp_q[4] * 1209541147953L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

