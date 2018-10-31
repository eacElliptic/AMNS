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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + ((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + ((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + ((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + ((int128)pa[4] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + ((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + ((int128)pa[4] * pa[4]);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 71996709448111527UL) + ((uint64_t)op[1] * 16942726803087307363UL) + ((uint64_t)op[2] * 18172008861773825924UL) + ((uint64_t)op[3] * 8657452982393209059UL) + ((uint64_t)op[4] * 11496046864426200974UL);
	tmp_q[1] = ((uint64_t)op[0] * 11496046864426200974UL) + ((uint64_t)op[1] * 71996709448111527UL) + ((uint64_t)op[2] * 16942726803087307363UL) + ((uint64_t)op[3] * 18172008861773825924UL) + ((uint64_t)op[4] * 8657452982393209059UL);
	tmp_q[2] = ((uint64_t)op[0] * 8657452982393209059UL) + ((uint64_t)op[1] * 11496046864426200974UL) + ((uint64_t)op[2] * 71996709448111527UL) + ((uint64_t)op[3] * 16942726803087307363UL) + ((uint64_t)op[4] * 18172008861773825924UL);
	tmp_q[3] = ((uint64_t)op[0] * 18172008861773825924UL) + ((uint64_t)op[1] * 8657452982393209059UL) + ((uint64_t)op[2] * 11496046864426200974UL) + ((uint64_t)op[3] * 71996709448111527UL) + ((uint64_t)op[4] * 16942726803087307363UL);
	tmp_q[4] = ((uint64_t)op[0] * 16942726803087307363UL) + ((uint64_t)op[1] * 18172008861773825924UL) + ((uint64_t)op[2] * 8657452982393209059UL) + ((uint64_t)op[3] * 11496046864426200974UL) + ((uint64_t)op[4] * 71996709448111527UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 97035963086700L) + (-((int128)tmp_q[1] * 34846675209L) - ((int128)tmp_q[2] * 124967671044055L) + ((int128)tmp_q[3] * 184083504189551L) + ((int128)tmp_q[4] * 37954976616414L));
	tmp_zero[1] = ((int128)tmp_q[0] * 37954976616414L) - ((int128)tmp_q[1] * 97035963086700L) + (-((int128)tmp_q[2] * 34846675209L) - ((int128)tmp_q[3] * 124967671044055L) + ((int128)tmp_q[4] * 184083504189551L));
	tmp_zero[2] = ((int128)tmp_q[0] * 184083504189551L) + ((int128)tmp_q[1] * 37954976616414L) - ((int128)tmp_q[2] * 97035963086700L) + (-((int128)tmp_q[3] * 34846675209L) - ((int128)tmp_q[4] * 124967671044055L));
	tmp_zero[3] = -((int128)tmp_q[0] * 124967671044055L) + ((int128)tmp_q[1] * 184083504189551L) + ((int128)tmp_q[2] * 37954976616414L) - ((int128)tmp_q[3] * 97035963086700L) - ((int128)tmp_q[4] * 34846675209L);
	tmp_zero[4] = -((int128)tmp_q[0] * 34846675209L) - ((int128)tmp_q[1] * 124967671044055L) + ((int128)tmp_q[2] * 184083504189551L) + ((int128)tmp_q[3] * 37954976616414L) - ((int128)tmp_q[4] * 97035963086700L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

