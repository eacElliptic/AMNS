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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17132633658236386835UL) + ((((uint64_t)op[1] * 16428680986350263173UL) + ((uint64_t)op[2] * 16920461383436411875UL) + ((uint64_t)op[3] * 4636016568459497345UL) + ((uint64_t)op[4] * 1558935545621611554UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 1558935545621611554UL) + ((uint64_t)op[1] * 17132633658236386835UL) + ((((uint64_t)op[2] * 16428680986350263173UL) + ((uint64_t)op[3] * 16920461383436411875UL) + ((uint64_t)op[4] * 4636016568459497345UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 4636016568459497345UL) + ((uint64_t)op[1] * 1558935545621611554UL) + ((uint64_t)op[2] * 17132633658236386835UL) + ((((uint64_t)op[3] * 16428680986350263173UL) + ((uint64_t)op[4] * 16920461383436411875UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 16920461383436411875UL) + ((uint64_t)op[1] * 4636016568459497345UL) + ((uint64_t)op[2] * 1558935545621611554UL) + ((uint64_t)op[3] * 17132633658236386835UL) + ((uint64_t)op[4] * 6338365549553820958UL);
	tmp_q[4] = ((uint64_t)op[0] * 16428680986350263173UL) + ((uint64_t)op[1] * 16920461383436411875UL) + ((uint64_t)op[2] * 4636016568459497345UL) + ((uint64_t)op[3] * 1558935545621611554UL) + ((uint64_t)op[4] * 17132633658236386835UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 133502423745L) + ((-((int128)tmp_q[1] * 138038349524L) + ((int128)tmp_q[2] * 33003052679L) - ((int128)tmp_q[3] * 28370688325L) + ((int128)tmp_q[4] * 212679779034L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 212679779034L) + ((int128)tmp_q[1] * 133502423745L) + ((-((int128)tmp_q[2] * 138038349524L) + ((int128)tmp_q[3] * 33003052679L) - ((int128)tmp_q[4] * 28370688325L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 28370688325L) + ((int128)tmp_q[1] * 212679779034L) + ((int128)tmp_q[2] * 133502423745L) + ((-((int128)tmp_q[3] * 138038349524L) + ((int128)tmp_q[4] * 33003052679L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 33003052679L) - ((int128)tmp_q[1] * 28370688325L) + ((int128)tmp_q[2] * 212679779034L) + ((int128)tmp_q[3] * 133502423745L) - ((int128)tmp_q[4] * 828230097144L);
	tmp_zero[4] = -((int128)tmp_q[0] * 138038349524L) + ((int128)tmp_q[1] * 33003052679L) - ((int128)tmp_q[2] * 28370688325L) + ((int128)tmp_q[3] * 212679779034L) + ((int128)tmp_q[4] * 133502423745L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

