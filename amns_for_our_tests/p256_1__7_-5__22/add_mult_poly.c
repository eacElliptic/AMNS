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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 741676892932634579UL) + ((((uint64_t)op[1] * 842799006679111839UL) + ((uint64_t)op[2] * 2768506922109124417UL) + ((uint64_t)op[3] * 3567942826394631927UL) + ((uint64_t)op[4] * 14907655503118622182UL) + ((uint64_t)op[5] * 17677432238978783618UL) + ((uint64_t)op[6] * 6521272494927823297UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 6521272494927823297UL) + ((uint64_t)op[1] * 741676892932634579UL) + ((((uint64_t)op[2] * 842799006679111839UL) + ((uint64_t)op[3] * 2768506922109124417UL) + ((uint64_t)op[4] * 3567942826394631927UL) + ((uint64_t)op[5] * 14907655503118622182UL) + ((uint64_t)op[6] * 17677432238978783618UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 17677432238978783618UL) + ((uint64_t)op[1] * 6521272494927823297UL) + ((uint64_t)op[2] * 741676892932634579UL) + ((((uint64_t)op[3] * 842799006679111839UL) + ((uint64_t)op[4] * 2768506922109124417UL) + ((uint64_t)op[5] * 3567942826394631927UL) + ((uint64_t)op[6] * 14907655503118622182UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 14907655503118622182UL) + ((uint64_t)op[1] * 17677432238978783618UL) + ((uint64_t)op[2] * 6521272494927823297UL) + ((uint64_t)op[3] * 741676892932634579UL) + ((((uint64_t)op[4] * 842799006679111839UL) + ((uint64_t)op[5] * 2768506922109124417UL) + ((uint64_t)op[6] * 3567942826394631927UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 3567942826394631927UL) + ((uint64_t)op[1] * 14907655503118622182UL) + ((uint64_t)op[2] * 17677432238978783618UL) + ((uint64_t)op[3] * 6521272494927823297UL) + ((uint64_t)op[4] * 741676892932634579UL) + ((((uint64_t)op[5] * 842799006679111839UL) + ((uint64_t)op[6] * 2768506922109124417UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 2768506922109124417UL) + ((uint64_t)op[1] * 3567942826394631927UL) + ((uint64_t)op[2] * 14907655503118622182UL) + ((uint64_t)op[3] * 17677432238978783618UL) + ((uint64_t)op[4] * 6521272494927823297UL) + ((uint64_t)op[5] * 741676892932634579UL) + ((uint64_t)op[6] * 14232749040313992421UL);
	tmp_q[6] = ((uint64_t)op[0] * 842799006679111839UL) + ((uint64_t)op[1] * 2768506922109124417UL) + ((uint64_t)op[2] * 3567942826394631927UL) + ((uint64_t)op[3] * 14907655503118622182UL) + ((uint64_t)op[4] * 17677432238978783618UL) + ((uint64_t)op[5] * 6521272494927823297UL) + ((uint64_t)op[6] * 741676892932634579UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 29327398008L) - ((((int128)tmp_q[1] * 16148604483L) - ((int128)tmp_q[2] * 17730104348L) - ((int128)tmp_q[3] * 57752124940L) + ((int128)tmp_q[4] * 26942071753L) - ((int128)tmp_q[5] * 15215581060L) + ((int128)tmp_q[6] * 1067654649L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 1067654649L) + ((int128)tmp_q[1] * 29327398008L) - ((((int128)tmp_q[2] * 16148604483L) - ((int128)tmp_q[3] * 17730104348L) - ((int128)tmp_q[4] * 57752124940L) + ((int128)tmp_q[5] * 26942071753L) - ((int128)tmp_q[6] * 15215581060L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 15215581060L) + ((int128)tmp_q[1] * 1067654649L) + ((int128)tmp_q[2] * 29327398008L) - ((((int128)tmp_q[3] * 16148604483L) - ((int128)tmp_q[4] * 17730104348L) - ((int128)tmp_q[5] * 57752124940L) + ((int128)tmp_q[6] * 26942071753L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 26942071753L) - ((int128)tmp_q[1] * 15215581060L) + ((int128)tmp_q[2] * 1067654649L) + ((int128)tmp_q[3] * 29327398008L) - ((((int128)tmp_q[4] * 16148604483L) - ((int128)tmp_q[5] * 17730104348L) - ((int128)tmp_q[6] * 57752124940L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 57752124940L) + ((int128)tmp_q[1] * 26942071753L) - ((int128)tmp_q[2] * 15215581060L) + ((int128)tmp_q[3] * 1067654649L) + ((int128)tmp_q[4] * 29327398008L) - ((((int128)tmp_q[5] * 16148604483L) - ((int128)tmp_q[6] * 17730104348L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 17730104348L) - ((int128)tmp_q[1] * 57752124940L) + ((int128)tmp_q[2] * 26942071753L) - ((int128)tmp_q[3] * 15215581060L) + ((int128)tmp_q[4] * 1067654649L) + ((int128)tmp_q[5] * 29327398008L) - ((int128)tmp_q[6] * 80743022415L);
	tmp_zero[6] = ((int128)tmp_q[0] * 16148604483L) - ((int128)tmp_q[1] * 17730104348L) - ((int128)tmp_q[2] * 57752124940L) + ((int128)tmp_q[3] * 26942071753L) - ((int128)tmp_q[4] * 15215581060L) + ((int128)tmp_q[5] * 1067654649L) + ((int128)tmp_q[6] * 29327398008L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

