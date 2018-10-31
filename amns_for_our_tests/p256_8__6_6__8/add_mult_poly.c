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
	tmp_q[0] = ((uint64_t)op[0] * 6827173719340415423UL) + ((((uint64_t)op[1] * 17112576043856387913UL) + ((uint64_t)op[2] * 10086105280530914197UL) + ((uint64_t)op[3] * 3665043501370130409UL) + ((uint64_t)op[4] * 7911276826571450412UL) + ((uint64_t)op[5] * 2754478810157455281UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 2754478810157455281UL) + ((uint64_t)op[1] * 6827173719340415423UL) + ((((uint64_t)op[2] * 17112576043856387913UL) + ((uint64_t)op[3] * 10086105280530914197UL) + ((uint64_t)op[4] * 3665043501370130409UL) + ((uint64_t)op[5] * 7911276826571450412UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 7911276826571450412UL) + ((uint64_t)op[1] * 2754478810157455281UL) + ((uint64_t)op[2] * 6827173719340415423UL) + ((((uint64_t)op[3] * 17112576043856387913UL) + ((uint64_t)op[4] * 10086105280530914197UL) + ((uint64_t)op[5] * 3665043501370130409UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 3665043501370130409UL) + ((uint64_t)op[1] * 7911276826571450412UL) + ((uint64_t)op[2] * 2754478810157455281UL) + ((uint64_t)op[3] * 6827173719340415423UL) + ((((uint64_t)op[4] * 17112576043856387913UL) + ((uint64_t)op[5] * 10086105280530914197UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 10086105280530914197UL) + ((uint64_t)op[1] * 3665043501370130409UL) + ((uint64_t)op[2] * 7911276826571450412UL) + ((uint64_t)op[3] * 2754478810157455281UL) + ((uint64_t)op[4] * 6827173719340415423UL) + ((uint64_t)op[5] * 10441735894590569398UL);
	tmp_q[5] = ((uint64_t)op[0] * 17112576043856387913UL) + ((uint64_t)op[1] * 10086105280530914197UL) + ((uint64_t)op[2] * 3665043501370130409UL) + ((uint64_t)op[3] * 7911276826571450412UL) + ((uint64_t)op[4] * 2754478810157455281UL) + ((uint64_t)op[5] * 6827173719340415423UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1591028306937L) + ((-((int128)tmp_q[1] * 4955395539291L) - ((int128)tmp_q[2] * 340059026126L) + ((int128)tmp_q[3] * 1607399273592L) + ((int128)tmp_q[4] * 465200038667L) - ((int128)tmp_q[5] * 1816461565011L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 1816461565011L) - ((int128)tmp_q[1] * 1591028306937L) + ((-((int128)tmp_q[2] * 4955395539291L) - ((int128)tmp_q[3] * 340059026126L) + ((int128)tmp_q[4] * 1607399273592L) + ((int128)tmp_q[5] * 465200038667L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 465200038667L) - ((int128)tmp_q[1] * 1816461565011L) - ((int128)tmp_q[2] * 1591028306937L) + ((-((int128)tmp_q[3] * 4955395539291L) - ((int128)tmp_q[4] * 340059026126L) + ((int128)tmp_q[5] * 1607399273592L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 1607399273592L) + ((int128)tmp_q[1] * 465200038667L) - ((int128)tmp_q[2] * 1816461565011L) - ((int128)tmp_q[3] * 1591028306937L) + ((-((int128)tmp_q[4] * 4955395539291L) - ((int128)tmp_q[5] * 340059026126L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 340059026126L) + ((int128)tmp_q[1] * 1607399273592L) + ((int128)tmp_q[2] * 465200038667L) - ((int128)tmp_q[3] * 1816461565011L) - ((int128)tmp_q[4] * 1591028306937L) - ((int128)tmp_q[5] * 29732373235746L);
	tmp_zero[5] = -((int128)tmp_q[0] * 4955395539291L) - ((int128)tmp_q[1] * 340059026126L) + ((int128)tmp_q[2] * 1607399273592L) + ((int128)tmp_q[3] * 465200038667L) - ((int128)tmp_q[4] * 1816461565011L) - ((int128)tmp_q[5] * 1591028306937L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

