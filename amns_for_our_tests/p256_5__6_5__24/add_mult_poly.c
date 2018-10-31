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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8561154301888035710UL) + ((((uint64_t)op[1] * 11044729357414625580UL) + ((uint64_t)op[2] * 16699390260457613622UL) + ((uint64_t)op[3] * 12610541151632757591UL) + ((uint64_t)op[4] * 7654094079010521176UL) + ((uint64_t)op[5] * 1753696230593894872UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 1753696230593894872UL) + ((uint64_t)op[1] * 8561154301888035710UL) + ((((uint64_t)op[2] * 11044729357414625580UL) + ((uint64_t)op[3] * 16699390260457613622UL) + ((uint64_t)op[4] * 12610541151632757591UL) + ((uint64_t)op[5] * 7654094079010521176UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 7654094079010521176UL) + ((uint64_t)op[1] * 1753696230593894872UL) + ((uint64_t)op[2] * 8561154301888035710UL) + ((((uint64_t)op[3] * 11044729357414625580UL) + ((uint64_t)op[4] * 16699390260457613622UL) + ((uint64_t)op[5] * 12610541151632757591UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 12610541151632757591UL) + ((uint64_t)op[1] * 7654094079010521176UL) + ((uint64_t)op[2] * 1753696230593894872UL) + ((uint64_t)op[3] * 8561154301888035710UL) + ((((uint64_t)op[4] * 11044729357414625580UL) + ((uint64_t)op[5] * 16699390260457613622UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 16699390260457613622UL) + ((uint64_t)op[1] * 12610541151632757591UL) + ((uint64_t)op[2] * 7654094079010521176UL) + ((uint64_t)op[3] * 1753696230593894872UL) + ((uint64_t)op[4] * 8561154301888035710UL) + ((uint64_t)op[5] * 18330158639654024668UL);
	tmp_q[5] = ((uint64_t)op[0] * 11044729357414625580UL) + ((uint64_t)op[1] * 16699390260457613622UL) + ((uint64_t)op[2] * 12610541151632757591UL) + ((uint64_t)op[3] * 7654094079010521176UL) + ((uint64_t)op[4] * 1753696230593894872UL) + ((uint64_t)op[5] * 8561154301888035710UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 5333179162438L) + ((-((int128)tmp_q[1] * 1366389371200L) + ((int128)tmp_q[2] * 4306934376966L) - ((int128)tmp_q[3] * 2140496208263L) - ((int128)tmp_q[4] * 3386928323952L) - ((int128)tmp_q[5] * 3976456106800L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 3976456106800L) + ((int128)tmp_q[1] * 5333179162438L) + ((-((int128)tmp_q[2] * 1366389371200L) + ((int128)tmp_q[3] * 4306934376966L) - ((int128)tmp_q[4] * 2140496208263L) - ((int128)tmp_q[5] * 3386928323952L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 3386928323952L) - ((int128)tmp_q[1] * 3976456106800L) + ((int128)tmp_q[2] * 5333179162438L) + ((-((int128)tmp_q[3] * 1366389371200L) + ((int128)tmp_q[4] * 4306934376966L) - ((int128)tmp_q[5] * 2140496208263L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 2140496208263L) - ((int128)tmp_q[1] * 3386928323952L) - ((int128)tmp_q[2] * 3976456106800L) + ((int128)tmp_q[3] * 5333179162438L) + ((-((int128)tmp_q[4] * 1366389371200L) + ((int128)tmp_q[5] * 4306934376966L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 4306934376966L) - ((int128)tmp_q[1] * 2140496208263L) - ((int128)tmp_q[2] * 3386928323952L) - ((int128)tmp_q[3] * 3976456106800L) + ((int128)tmp_q[4] * 5333179162438L) - ((int128)tmp_q[5] * 6831946856000L);
	tmp_zero[5] = -((int128)tmp_q[0] * 1366389371200L) + ((int128)tmp_q[1] * 4306934376966L) - ((int128)tmp_q[2] * 2140496208263L) - ((int128)tmp_q[3] * 3386928323952L) - ((int128)tmp_q[4] * 3976456106800L) + ((int128)tmp_q[5] * 5333179162438L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

