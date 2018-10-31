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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14038851869899999509UL) + ((((uint64_t)op[1] * 5464171001065146238UL) + ((uint64_t)op[2] * 5414994660439508063UL) + ((uint64_t)op[3] * 14573318395261585762UL) + ((uint64_t)op[4] * 1087406978345860525UL) + ((uint64_t)op[5] * 1878814335062012664UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 1878814335062012664UL) + ((uint64_t)op[1] * 14038851869899999509UL) + ((((uint64_t)op[2] * 5464171001065146238UL) + ((uint64_t)op[3] * 5414994660439508063UL) + ((uint64_t)op[4] * 14573318395261585762UL) + ((uint64_t)op[5] * 1087406978345860525UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 1087406978345860525UL) + ((uint64_t)op[1] * 1878814335062012664UL) + ((uint64_t)op[2] * 14038851869899999509UL) + ((((uint64_t)op[3] * 5464171001065146238UL) + ((uint64_t)op[4] * 5414994660439508063UL) + ((uint64_t)op[5] * 14573318395261585762UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 14573318395261585762UL) + ((uint64_t)op[1] * 1087406978345860525UL) + ((uint64_t)op[2] * 1878814335062012664UL) + ((uint64_t)op[3] * 14038851869899999509UL) + ((((uint64_t)op[4] * 5464171001065146238UL) + ((uint64_t)op[5] * 5414994660439508063UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 5414994660439508063UL) + ((uint64_t)op[1] * 14573318395261585762UL) + ((uint64_t)op[2] * 1087406978345860525UL) + ((uint64_t)op[3] * 1878814335062012664UL) + ((uint64_t)op[4] * 14038851869899999509UL) + ((uint64_t)op[5] * 14338281932681325812UL);
	tmp_q[5] = ((uint64_t)op[0] * 5464171001065146238UL) + ((uint64_t)op[1] * 5414994660439508063UL) + ((uint64_t)op[2] * 14573318395261585762UL) + ((uint64_t)op[3] * 1087406978345860525UL) + ((uint64_t)op[4] * 1878814335062012664UL) + ((uint64_t)op[5] * 14038851869899999509UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2446776023047L) + ((-((int128)tmp_q[1] * 1669226341358L) - ((int128)tmp_q[2] * 3992287762686L) + ((int128)tmp_q[3] * 1505334292214L) - ((int128)tmp_q[4] * 1008665686981L) - ((int128)tmp_q[5] * 210536488628L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 210536488628L) - ((int128)tmp_q[1] * 2446776023047L) + ((-((int128)tmp_q[2] * 1669226341358L) - ((int128)tmp_q[3] * 3992287762686L) + ((int128)tmp_q[4] * 1505334292214L) - ((int128)tmp_q[5] * 1008665686981L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 1008665686981L) - ((int128)tmp_q[1] * 210536488628L) - ((int128)tmp_q[2] * 2446776023047L) + ((-((int128)tmp_q[3] * 1669226341358L) - ((int128)tmp_q[4] * 3992287762686L) + ((int128)tmp_q[5] * 1505334292214L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 1505334292214L) - ((int128)tmp_q[1] * 1008665686981L) - ((int128)tmp_q[2] * 210536488628L) - ((int128)tmp_q[3] * 2446776023047L) + ((-((int128)tmp_q[4] * 1669226341358L) - ((int128)tmp_q[5] * 3992287762686L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 3992287762686L) + ((int128)tmp_q[1] * 1505334292214L) - ((int128)tmp_q[2] * 1008665686981L) - ((int128)tmp_q[3] * 210536488628L) - ((int128)tmp_q[4] * 2446776023047L) - ((int128)tmp_q[5] * 10015358048148L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1669226341358L) - ((int128)tmp_q[1] * 3992287762686L) + ((int128)tmp_q[2] * 1505334292214L) - ((int128)tmp_q[3] * 1008665686981L) - ((int128)tmp_q[4] * 210536488628L) - ((int128)tmp_q[5] * 2446776023047L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

