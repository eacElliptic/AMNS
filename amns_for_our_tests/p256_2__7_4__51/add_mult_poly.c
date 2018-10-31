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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8735231032896795599UL) + ((((uint64_t)op[1] * 6442758020575755363UL) + ((uint64_t)op[2] * 9524335271146308589UL) + ((uint64_t)op[3] * 393795604101892667UL) + ((uint64_t)op[4] * 1823028408051073808UL) + ((uint64_t)op[5] * 8283541846087647528UL) + ((uint64_t)op[6] * 10846822680092352945UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 10846822680092352945UL) + ((uint64_t)op[1] * 8735231032896795599UL) + ((((uint64_t)op[2] * 6442758020575755363UL) + ((uint64_t)op[3] * 9524335271146308589UL) + ((uint64_t)op[4] * 393795604101892667UL) + ((uint64_t)op[5] * 1823028408051073808UL) + ((uint64_t)op[6] * 8283541846087647528UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 8283541846087647528UL) + ((uint64_t)op[1] * 10846822680092352945UL) + ((uint64_t)op[2] * 8735231032896795599UL) + ((((uint64_t)op[3] * 6442758020575755363UL) + ((uint64_t)op[4] * 9524335271146308589UL) + ((uint64_t)op[5] * 393795604101892667UL) + ((uint64_t)op[6] * 1823028408051073808UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 1823028408051073808UL) + ((uint64_t)op[1] * 8283541846087647528UL) + ((uint64_t)op[2] * 10846822680092352945UL) + ((uint64_t)op[3] * 8735231032896795599UL) + ((((uint64_t)op[4] * 6442758020575755363UL) + ((uint64_t)op[5] * 9524335271146308589UL) + ((uint64_t)op[6] * 393795604101892667UL)) * 4);
	tmp_q[4] = ((uint64_t)op[0] * 393795604101892667UL) + ((uint64_t)op[1] * 1823028408051073808UL) + ((uint64_t)op[2] * 8283541846087647528UL) + ((uint64_t)op[3] * 10846822680092352945UL) + ((uint64_t)op[4] * 8735231032896795599UL) + ((((uint64_t)op[5] * 6442758020575755363UL) + ((uint64_t)op[6] * 9524335271146308589UL)) * 4);
	tmp_q[5] = ((uint64_t)op[0] * 9524335271146308589UL) + ((uint64_t)op[1] * 393795604101892667UL) + ((uint64_t)op[2] * 1823028408051073808UL) + ((uint64_t)op[3] * 8283541846087647528UL) + ((uint64_t)op[4] * 10846822680092352945UL) + ((uint64_t)op[5] * 8735231032896795599UL) + ((uint64_t)op[6] * 7324288008593469836UL);
	tmp_q[6] = ((uint64_t)op[0] * 6442758020575755363UL) + ((uint64_t)op[1] * 9524335271146308589UL) + ((uint64_t)op[2] * 393795604101892667UL) + ((uint64_t)op[3] * 1823028408051073808UL) + ((uint64_t)op[4] * 8283541846087647528UL) + ((uint64_t)op[5] * 10846822680092352945UL) + ((uint64_t)op[6] * 8735231032896795599UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 7516240569L) + ((((int128)tmp_q[1] * 790749923L) - ((int128)tmp_q[2] * 56149188924L) + ((int128)tmp_q[3] * 31358456528L) + ((int128)tmp_q[4] * 18270374553L) - ((int128)tmp_q[5] * 1831598195L) - ((int128)tmp_q[6] * 16871487351L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 16871487351L) + ((int128)tmp_q[1] * 7516240569L) + ((((int128)tmp_q[2] * 790749923L) - ((int128)tmp_q[3] * 56149188924L) + ((int128)tmp_q[4] * 31358456528L) + ((int128)tmp_q[5] * 18270374553L) - ((int128)tmp_q[6] * 1831598195L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 1831598195L) - ((int128)tmp_q[1] * 16871487351L) + ((int128)tmp_q[2] * 7516240569L) + ((((int128)tmp_q[3] * 790749923L) - ((int128)tmp_q[4] * 56149188924L) + ((int128)tmp_q[5] * 31358456528L) + ((int128)tmp_q[6] * 18270374553L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 18270374553L) - ((int128)tmp_q[1] * 1831598195L) - ((int128)tmp_q[2] * 16871487351L) + ((int128)tmp_q[3] * 7516240569L) + ((((int128)tmp_q[4] * 790749923L) - ((int128)tmp_q[5] * 56149188924L) + ((int128)tmp_q[6] * 31358456528L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 31358456528L) + ((int128)tmp_q[1] * 18270374553L) - ((int128)tmp_q[2] * 1831598195L) - ((int128)tmp_q[3] * 16871487351L) + ((int128)tmp_q[4] * 7516240569L) + ((((int128)tmp_q[5] * 790749923L) - ((int128)tmp_q[6] * 56149188924L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 56149188924L) + ((int128)tmp_q[1] * 31358456528L) + ((int128)tmp_q[2] * 18270374553L) - ((int128)tmp_q[3] * 1831598195L) - ((int128)tmp_q[4] * 16871487351L) + ((int128)tmp_q[5] * 7516240569L) + ((int128)tmp_q[6] * 3162999692L);
	tmp_zero[6] = ((int128)tmp_q[0] * 790749923L) - ((int128)tmp_q[1] * 56149188924L) + ((int128)tmp_q[2] * 31358456528L) + ((int128)tmp_q[3] * 18270374553L) - ((int128)tmp_q[4] * 1831598195L) - ((int128)tmp_q[5] * 16871487351L) + ((int128)tmp_q[6] * 7516240569L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

