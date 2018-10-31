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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[3] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1107060251482376027UL) + ((((uint64_t)op[1] * 1739419374733238386UL) + ((uint64_t)op[2] * 3077772582076886366UL) + ((uint64_t)op[3] * 8335203623384588822UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 8335203623384588822UL) + ((uint64_t)op[1] * 1107060251482376027UL) + ((((uint64_t)op[2] * 1739419374733238386UL) + ((uint64_t)op[3] * 3077772582076886366UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 3077772582076886366UL) + ((uint64_t)op[1] * 8335203623384588822UL) + ((uint64_t)op[2] * 1107060251482376027UL) + ((uint64_t)op[3] * 5218258124199715158UL);
	tmp_q[3] = ((uint64_t)op[0] * 1739419374733238386UL) + ((uint64_t)op[1] * 3077772582076886366UL) + ((uint64_t)op[2] * 8335203623384588822UL) + ((uint64_t)op[3] * 1107060251482376027UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 225246656385601L) + ((-((int128)tmp_q[1] * 116792325329502L) - ((int128)tmp_q[2] * 45453137138538L) + ((int128)tmp_q[3] * 150596746167126L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 150596746167126L) + ((int128)tmp_q[1] * 225246656385601L) + ((-((int128)tmp_q[2] * 116792325329502L) - ((int128)tmp_q[3] * 45453137138538L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 45453137138538L) + ((int128)tmp_q[1] * 150596746167126L) + ((int128)tmp_q[2] * 225246656385601L) - ((int128)tmp_q[3] * 350376975988506L);
	tmp_zero[3] = -((int128)tmp_q[0] * 116792325329502L) - ((int128)tmp_q[1] * 45453137138538L) + ((int128)tmp_q[2] * 150596746167126L) + ((int128)tmp_q[3] * 225246656385601L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

