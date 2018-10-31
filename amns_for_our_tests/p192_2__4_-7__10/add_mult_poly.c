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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[3] + (int128)pa[3] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[3] * pa[1]) << 1) + (int128)pa[2] * pa[2]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[3] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8476779484070440184UL) + ((((uint64_t)op[1] * 4517039474427438618UL) + ((uint64_t)op[2] * 17546553678431452673UL) + ((uint64_t)op[3] * 7773482846847573814UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 7773482846847573814UL) + ((uint64_t)op[1] * 8476779484070440184UL) + ((((uint64_t)op[2] * 4517039474427438618UL) + ((uint64_t)op[3] * 17546553678431452673UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 17546553678431452673UL) + ((uint64_t)op[1] * 7773482846847573814UL) + ((uint64_t)op[2] * 8476779484070440184UL) + ((uint64_t)op[3] * 5274211826427032906UL);
	tmp_q[3] = ((uint64_t)op[0] * 4517039474427438618UL) + ((uint64_t)op[1] * 17546553678431452673UL) + ((uint64_t)op[2] * 7773482846847573814UL) + ((uint64_t)op[3] * 8476779484070440184UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 50833363719216L) - ((-((int128)tmp_q[1] * 101940984143126L) - ((int128)tmp_q[2] * 70617337956225L) - ((int128)tmp_q[3] * 102023011163546L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 102023011163546L) - ((int128)tmp_q[1] * 50833363719216L) - ((-((int128)tmp_q[2] * 101940984143126L) - ((int128)tmp_q[3] * 70617337956225L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 70617337956225L) - ((int128)tmp_q[1] * 102023011163546L) - ((int128)tmp_q[2] * 50833363719216L) + ((int128)tmp_q[3] * 713586889001882L);
	tmp_zero[3] = -((int128)tmp_q[0] * 101940984143126L) - ((int128)tmp_q[1] * 70617337956225L) - ((int128)tmp_q[2] * 102023011163546L) - ((int128)tmp_q[3] * 50833363719216L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

