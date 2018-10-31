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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5717764958283850773UL) + ((((uint64_t)op[1] * 8392589489802897368UL) + ((uint64_t)op[2] * 16870712854057143662UL) + ((uint64_t)op[3] * 16776434306912751365UL) + ((uint64_t)op[4] * 13879467374483062934UL) + ((uint64_t)op[5] * 12628833036188415361UL) + ((uint64_t)op[6] * 13572423582393799122UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 13572423582393799122UL) + ((uint64_t)op[1] * 5717764958283850773UL) + ((((uint64_t)op[2] * 8392589489802897368UL) + ((uint64_t)op[3] * 16870712854057143662UL) + ((uint64_t)op[4] * 16776434306912751365UL) + ((uint64_t)op[5] * 13879467374483062934UL) + ((uint64_t)op[6] * 12628833036188415361UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 12628833036188415361UL) + ((uint64_t)op[1] * 13572423582393799122UL) + ((uint64_t)op[2] * 5717764958283850773UL) + ((((uint64_t)op[3] * 8392589489802897368UL) + ((uint64_t)op[4] * 16870712854057143662UL) + ((uint64_t)op[5] * 16776434306912751365UL) + ((uint64_t)op[6] * 13879467374483062934UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 13879467374483062934UL) + ((uint64_t)op[1] * 12628833036188415361UL) + ((uint64_t)op[2] * 13572423582393799122UL) + ((uint64_t)op[3] * 5717764958283850773UL) + ((((uint64_t)op[4] * 8392589489802897368UL) + ((uint64_t)op[5] * 16870712854057143662UL) + ((uint64_t)op[6] * 16776434306912751365UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 16776434306912751365UL) + ((uint64_t)op[1] * 13879467374483062934UL) + ((uint64_t)op[2] * 12628833036188415361UL) + ((uint64_t)op[3] * 13572423582393799122UL) + ((uint64_t)op[4] * 5717764958283850773UL) + ((((uint64_t)op[5] * 8392589489802897368UL) + ((uint64_t)op[6] * 16870712854057143662UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 16870712854057143662UL) + ((uint64_t)op[1] * 16776434306912751365UL) + ((uint64_t)op[2] * 13879467374483062934UL) + ((uint64_t)op[3] * 12628833036188415361UL) + ((uint64_t)op[4] * 13572423582393799122UL) + ((uint64_t)op[5] * 5717764958283850773UL) + ((uint64_t)op[6] * 5069459301595383608UL);
	tmp_q[6] = ((uint64_t)op[0] * 8392589489802897368UL) + ((uint64_t)op[1] * 16870712854057143662UL) + ((uint64_t)op[2] * 16776434306912751365UL) + ((uint64_t)op[3] * 13879467374483062934UL) + ((uint64_t)op[4] * 12628833036188415361UL) + ((uint64_t)op[5] * 13572423582393799122UL) + ((uint64_t)op[6] * 5717764958283850773UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 13476037697L) + ((-((int128)tmp_q[1] * 2848753469L) - ((int128)tmp_q[2] * 9302692833L) - ((int128)tmp_q[3] * 23996886561L) + ((int128)tmp_q[4] * 12509537517L) + ((int128)tmp_q[5] * 52160073774L) - ((int128)tmp_q[6] * 59726180996L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 59726180996L) - ((int128)tmp_q[1] * 13476037697L) + ((-((int128)tmp_q[2] * 2848753469L) - ((int128)tmp_q[3] * 9302692833L) - ((int128)tmp_q[4] * 23996886561L) + ((int128)tmp_q[5] * 12509537517L) + ((int128)tmp_q[6] * 52160073774L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 52160073774L) - ((int128)tmp_q[1] * 59726180996L) - ((int128)tmp_q[2] * 13476037697L) + ((-((int128)tmp_q[3] * 2848753469L) - ((int128)tmp_q[4] * 9302692833L) - ((int128)tmp_q[5] * 23996886561L) + ((int128)tmp_q[6] * 12509537517L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 12509537517L) + ((int128)tmp_q[1] * 52160073774L) - ((int128)tmp_q[2] * 59726180996L) - ((int128)tmp_q[3] * 13476037697L) + ((-((int128)tmp_q[4] * 2848753469L) - ((int128)tmp_q[5] * 9302692833L) - ((int128)tmp_q[6] * 23996886561L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 23996886561L) + ((int128)tmp_q[1] * 12509537517L) + ((int128)tmp_q[2] * 52160073774L) - ((int128)tmp_q[3] * 59726180996L) - ((int128)tmp_q[4] * 13476037697L) + ((-((int128)tmp_q[5] * 2848753469L) - ((int128)tmp_q[6] * 9302692833L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 9302692833L) - ((int128)tmp_q[1] * 23996886561L) + ((int128)tmp_q[2] * 12509537517L) + ((int128)tmp_q[3] * 52160073774L) - ((int128)tmp_q[4] * 59726180996L) - ((int128)tmp_q[5] * 13476037697L) - ((int128)tmp_q[6] * 14243767345L);
	tmp_zero[6] = -((int128)tmp_q[0] * 2848753469L) - ((int128)tmp_q[1] * 9302692833L) - ((int128)tmp_q[2] * 23996886561L) + ((int128)tmp_q[3] * 12509537517L) + ((int128)tmp_q[4] * 52160073774L) - ((int128)tmp_q[5] * 59726180996L) - ((int128)tmp_q[6] * 13476037697L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

