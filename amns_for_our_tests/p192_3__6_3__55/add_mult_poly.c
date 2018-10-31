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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4842687503373353822UL) + ((((uint64_t)op[1] * 16685439728553560634UL) + ((uint64_t)op[2] * 6049841194951080633UL) + ((uint64_t)op[3] * 9354784725626734126UL) + ((uint64_t)op[4] * 18283019325419418300UL) + ((uint64_t)op[5] * 3979020378399871780UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 3979020378399871780UL) + ((uint64_t)op[1] * 4842687503373353822UL) + ((((uint64_t)op[2] * 16685439728553560634UL) + ((uint64_t)op[3] * 6049841194951080633UL) + ((uint64_t)op[4] * 9354784725626734126UL) + ((uint64_t)op[5] * 18283019325419418300UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 18283019325419418300UL) + ((uint64_t)op[1] * 3979020378399871780UL) + ((uint64_t)op[2] * 4842687503373353822UL) + ((((uint64_t)op[3] * 16685439728553560634UL) + ((uint64_t)op[4] * 6049841194951080633UL) + ((uint64_t)op[5] * 9354784725626734126UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 9354784725626734126UL) + ((uint64_t)op[1] * 18283019325419418300UL) + ((uint64_t)op[2] * 3979020378399871780UL) + ((uint64_t)op[3] * 4842687503373353822UL) + ((((uint64_t)op[4] * 16685439728553560634UL) + ((uint64_t)op[5] * 6049841194951080633UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 6049841194951080633UL) + ((uint64_t)op[1] * 9354784725626734126UL) + ((uint64_t)op[2] * 18283019325419418300UL) + ((uint64_t)op[3] * 3979020378399871780UL) + ((uint64_t)op[4] * 4842687503373353822UL) + ((uint64_t)op[5] * 13162831038241578670UL);
	tmp_q[5] = ((uint64_t)op[0] * 16685439728553560634UL) + ((uint64_t)op[1] * 6049841194951080633UL) + ((uint64_t)op[2] * 9354784725626734126UL) + ((uint64_t)op[3] * 18283019325419418300UL) + ((uint64_t)op[4] * 3979020378399871780UL) + ((uint64_t)op[5] * 4842687503373353822UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 409272548L) + ((((int128)tmp_q[1] * 845985188L) + ((int128)tmp_q[2] * 1047339698L) + ((int128)tmp_q[3] * 1212841934L) + ((int128)tmp_q[4] * 3210186453L) + ((int128)tmp_q[5] * 904767658L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 904767658L) - ((int128)tmp_q[1] * 409272548L) + ((((int128)tmp_q[2] * 845985188L) + ((int128)tmp_q[3] * 1047339698L) + ((int128)tmp_q[4] * 1212841934L) + ((int128)tmp_q[5] * 3210186453L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 3210186453L) + ((int128)tmp_q[1] * 904767658L) - ((int128)tmp_q[2] * 409272548L) + ((((int128)tmp_q[3] * 845985188L) + ((int128)tmp_q[4] * 1047339698L) + ((int128)tmp_q[5] * 1212841934L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 1212841934L) + ((int128)tmp_q[1] * 3210186453L) + ((int128)tmp_q[2] * 904767658L) - ((int128)tmp_q[3] * 409272548L) + ((((int128)tmp_q[4] * 845985188L) + ((int128)tmp_q[5] * 1047339698L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 1047339698L) + ((int128)tmp_q[1] * 1212841934L) + ((int128)tmp_q[2] * 3210186453L) + ((int128)tmp_q[3] * 904767658L) - ((int128)tmp_q[4] * 409272548L) + ((int128)tmp_q[5] * 2537955564L);
	tmp_zero[5] = ((int128)tmp_q[0] * 845985188L) + ((int128)tmp_q[1] * 1047339698L) + ((int128)tmp_q[2] * 1212841934L) + ((int128)tmp_q[3] * 3210186453L) + ((int128)tmp_q[4] * 904767658L) - ((int128)tmp_q[5] * 409272548L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

