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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15085822155764637614UL) + ((((uint64_t)op[1] * 8191657366720884696UL) + ((uint64_t)op[2] * 485059667976147754UL) + ((uint64_t)op[3] * 9545240971946732668UL) + ((uint64_t)op[4] * 7216823613677579664UL) + ((uint64_t)op[5] * 16059369633361222573UL) + ((uint64_t)op[6] * 17424652497843315302UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 17424652497843315302UL) + ((uint64_t)op[1] * 15085822155764637614UL) + ((((uint64_t)op[2] * 8191657366720884696UL) + ((uint64_t)op[3] * 485059667976147754UL) + ((uint64_t)op[4] * 9545240971946732668UL) + ((uint64_t)op[5] * 7216823613677579664UL) + ((uint64_t)op[6] * 16059369633361222573UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 16059369633361222573UL) + ((uint64_t)op[1] * 17424652497843315302UL) + ((uint64_t)op[2] * 15085822155764637614UL) + ((((uint64_t)op[3] * 8191657366720884696UL) + ((uint64_t)op[4] * 485059667976147754UL) + ((uint64_t)op[5] * 9545240971946732668UL) + ((uint64_t)op[6] * 7216823613677579664UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 7216823613677579664UL) + ((uint64_t)op[1] * 16059369633361222573UL) + ((uint64_t)op[2] * 17424652497843315302UL) + ((uint64_t)op[3] * 15085822155764637614UL) + ((((uint64_t)op[4] * 8191657366720884696UL) + ((uint64_t)op[5] * 485059667976147754UL) + ((uint64_t)op[6] * 9545240971946732668UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 9545240971946732668UL) + ((uint64_t)op[1] * 7216823613677579664UL) + ((uint64_t)op[2] * 16059369633361222573UL) + ((uint64_t)op[3] * 17424652497843315302UL) + ((uint64_t)op[4] * 15085822155764637614UL) + ((((uint64_t)op[5] * 8191657366720884696UL) + ((uint64_t)op[6] * 485059667976147754UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 485059667976147754UL) + ((uint64_t)op[1] * 9545240971946732668UL) + ((uint64_t)op[2] * 7216823613677579664UL) + ((uint64_t)op[3] * 16059369633361222573UL) + ((uint64_t)op[4] * 17424652497843315302UL) + ((uint64_t)op[5] * 15085822155764637614UL) + ((uint64_t)op[6] * 6128228026453102472UL);
	tmp_q[6] = ((uint64_t)op[0] * 8191657366720884696UL) + ((uint64_t)op[1] * 485059667976147754UL) + ((uint64_t)op[2] * 9545240971946732668UL) + ((uint64_t)op[3] * 7216823613677579664UL) + ((uint64_t)op[4] * 16059369633361222573UL) + ((uint64_t)op[5] * 17424652497843315302UL) + ((uint64_t)op[6] * 15085822155764637614UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 35589714412L) + ((-((int128)tmp_q[1] * 43020431256L) + ((int128)tmp_q[2] * 14063396985L) + ((int128)tmp_q[3] * 41846719638L) + ((int128)tmp_q[4] * 68053576070L) + ((int128)tmp_q[5] * 71004160928L) - ((int128)tmp_q[6] * 32483546490L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 32483546490L) - ((int128)tmp_q[1] * 35589714412L) + ((-((int128)tmp_q[2] * 43020431256L) + ((int128)tmp_q[3] * 14063396985L) + ((int128)tmp_q[4] * 41846719638L) + ((int128)tmp_q[5] * 68053576070L) + ((int128)tmp_q[6] * 71004160928L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 71004160928L) - ((int128)tmp_q[1] * 32483546490L) - ((int128)tmp_q[2] * 35589714412L) + ((-((int128)tmp_q[3] * 43020431256L) + ((int128)tmp_q[4] * 14063396985L) + ((int128)tmp_q[5] * 41846719638L) + ((int128)tmp_q[6] * 68053576070L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 68053576070L) + ((int128)tmp_q[1] * 71004160928L) - ((int128)tmp_q[2] * 32483546490L) - ((int128)tmp_q[3] * 35589714412L) + ((-((int128)tmp_q[4] * 43020431256L) + ((int128)tmp_q[5] * 14063396985L) + ((int128)tmp_q[6] * 41846719638L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 41846719638L) + ((int128)tmp_q[1] * 68053576070L) + ((int128)tmp_q[2] * 71004160928L) - ((int128)tmp_q[3] * 32483546490L) - ((int128)tmp_q[4] * 35589714412L) + ((-((int128)tmp_q[5] * 43020431256L) + ((int128)tmp_q[6] * 14063396985L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 14063396985L) + ((int128)tmp_q[1] * 41846719638L) + ((int128)tmp_q[2] * 68053576070L) + ((int128)tmp_q[3] * 71004160928L) - ((int128)tmp_q[4] * 32483546490L) - ((int128)tmp_q[5] * 35589714412L) - ((int128)tmp_q[6] * 129061293768L);
	tmp_zero[6] = -((int128)tmp_q[0] * 43020431256L) + ((int128)tmp_q[1] * 14063396985L) + ((int128)tmp_q[2] * 41846719638L) + ((int128)tmp_q[3] * 68053576070L) + ((int128)tmp_q[4] * 71004160928L) - ((int128)tmp_q[5] * 32483546490L) - ((int128)tmp_q[6] * 35589714412L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

