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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9680381270770780307UL) + ((((uint64_t)op[1] * 14510413457574190928UL) + ((uint64_t)op[2] * 15104215421466475500UL) + ((uint64_t)op[3] * 16146217190789674036UL) + ((uint64_t)op[4] * 4808881821881701328UL) + ((uint64_t)op[5] * 17808453575077831874UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 17808453575077831874UL) + ((uint64_t)op[1] * 9680381270770780307UL) + ((((uint64_t)op[2] * 14510413457574190928UL) + ((uint64_t)op[3] * 15104215421466475500UL) + ((uint64_t)op[4] * 16146217190789674036UL) + ((uint64_t)op[5] * 4808881821881701328UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 4808881821881701328UL) + ((uint64_t)op[1] * 17808453575077831874UL) + ((uint64_t)op[2] * 9680381270770780307UL) + ((((uint64_t)op[3] * 14510413457574190928UL) + ((uint64_t)op[4] * 15104215421466475500UL) + ((uint64_t)op[5] * 16146217190789674036UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 16146217190789674036UL) + ((uint64_t)op[1] * 4808881821881701328UL) + ((uint64_t)op[2] * 17808453575077831874UL) + ((uint64_t)op[3] * 9680381270770780307UL) + ((((uint64_t)op[4] * 14510413457574190928UL) + ((uint64_t)op[5] * 15104215421466475500UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 15104215421466475500UL) + ((uint64_t)op[1] * 16146217190789674036UL) + ((uint64_t)op[2] * 4808881821881701328UL) + ((uint64_t)op[3] * 17808453575077831874UL) + ((uint64_t)op[4] * 9680381270770780307UL) + ((uint64_t)op[5] * 1234909006967251824UL);
	tmp_q[5] = ((uint64_t)op[0] * 14510413457574190928UL) + ((uint64_t)op[1] * 15104215421466475500UL) + ((uint64_t)op[2] * 16146217190789674036UL) + ((uint64_t)op[3] * 4808881821881701328UL) + ((uint64_t)op[4] * 17808453575077831874UL) + ((uint64_t)op[5] * 9680381270770780307UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 584946074053L) - ((((int128)tmp_q[1] * 3712266381328L) + ((int128)tmp_q[2] * 148221166572L) + ((int128)tmp_q[3] * 3638753454652L) + ((int128)tmp_q[4] * 6448336975924L) - ((int128)tmp_q[5] * 22777805902L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 22777805902L) + ((int128)tmp_q[1] * 584946074053L) - ((((int128)tmp_q[2] * 3712266381328L) + ((int128)tmp_q[3] * 148221166572L) + ((int128)tmp_q[4] * 3638753454652L) + ((int128)tmp_q[5] * 6448336975924L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 6448336975924L) - ((int128)tmp_q[1] * 22777805902L) + ((int128)tmp_q[2] * 584946074053L) - ((((int128)tmp_q[3] * 3712266381328L) + ((int128)tmp_q[4] * 148221166572L) + ((int128)tmp_q[5] * 3638753454652L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 3638753454652L) + ((int128)tmp_q[1] * 6448336975924L) - ((int128)tmp_q[2] * 22777805902L) + ((int128)tmp_q[3] * 584946074053L) - ((((int128)tmp_q[4] * 3712266381328L) + ((int128)tmp_q[5] * 148221166572L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 148221166572L) + ((int128)tmp_q[1] * 3638753454652L) + ((int128)tmp_q[2] * 6448336975924L) - ((int128)tmp_q[3] * 22777805902L) + ((int128)tmp_q[4] * 584946074053L) - ((int128)tmp_q[5] * 18561331906640L);
	tmp_zero[5] = ((int128)tmp_q[0] * 3712266381328L) + ((int128)tmp_q[1] * 148221166572L) + ((int128)tmp_q[2] * 3638753454652L) + ((int128)tmp_q[3] * 6448336975924L) - ((int128)tmp_q[4] * 22777805902L) + ((int128)tmp_q[5] * 584946074053L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

