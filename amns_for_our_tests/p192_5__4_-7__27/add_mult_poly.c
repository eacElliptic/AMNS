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
	tmp_q[0] = ((uint64_t)op[0] * 7683737246485837200UL) + ((((uint64_t)op[1] * 10039402637776582496UL) + ((uint64_t)op[2] * 12348502763410823565UL) + ((uint64_t)op[3] * 9196857665890010426UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 9196857665890010426UL) + ((uint64_t)op[1] * 7683737246485837200UL) + ((((uint64_t)op[2] * 10039402637776582496UL) + ((uint64_t)op[3] * 12348502763410823565UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 12348502763410823565UL) + ((uint64_t)op[1] * 9196857665890010426UL) + ((uint64_t)op[2] * 7683737246485837200UL) + ((uint64_t)op[3] * 3511157830402128992UL);
	tmp_q[3] = ((uint64_t)op[0] * 10039402637776582496UL) + ((uint64_t)op[1] * 12348502763410823565UL) + ((uint64_t)op[2] * 9196857665890010426UL) + ((uint64_t)op[3] * 7683737246485837200UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 71050938225636L) - ((((int128)tmp_q[1] * 34999585411016L) + ((int128)tmp_q[2] * 9117208816963L) + ((int128)tmp_q[3] * 111715201213082L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 111715201213082L) - ((int128)tmp_q[1] * 71050938225636L) - ((((int128)tmp_q[2] * 34999585411016L) + ((int128)tmp_q[3] * 9117208816963L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 9117208816963L) + ((int128)tmp_q[1] * 111715201213082L) - ((int128)tmp_q[2] * 71050938225636L) - ((int128)tmp_q[3] * 244997097877112L);
	tmp_zero[3] = ((int128)tmp_q[0] * 34999585411016L) + ((int128)tmp_q[1] * 9117208816963L) + ((int128)tmp_q[2] * 111715201213082L) - ((int128)tmp_q[3] * 71050938225636L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
}

