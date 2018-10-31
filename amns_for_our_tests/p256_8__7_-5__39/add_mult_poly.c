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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1806799135077294067UL) + ((((uint64_t)op[1] * 6025982956261339846UL) + ((uint64_t)op[2] * 3826797987146243823UL) + ((uint64_t)op[3] * 3402412088961204283UL) + ((uint64_t)op[4] * 3693034966571312981UL) + ((uint64_t)op[5] * 7933916339959028336UL) + ((uint64_t)op[6] * 10957193811830801065UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 10957193811830801065UL) + ((uint64_t)op[1] * 1806799135077294067UL) + ((((uint64_t)op[2] * 6025982956261339846UL) + ((uint64_t)op[3] * 3826797987146243823UL) + ((uint64_t)op[4] * 3402412088961204283UL) + ((uint64_t)op[5] * 3693034966571312981UL) + ((uint64_t)op[6] * 7933916339959028336UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 7933916339959028336UL) + ((uint64_t)op[1] * 10957193811830801065UL) + ((uint64_t)op[2] * 1806799135077294067UL) + ((((uint64_t)op[3] * 6025982956261339846UL) + ((uint64_t)op[4] * 3826797987146243823UL) + ((uint64_t)op[5] * 3402412088961204283UL) + ((uint64_t)op[6] * 3693034966571312981UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 3693034966571312981UL) + ((uint64_t)op[1] * 7933916339959028336UL) + ((uint64_t)op[2] * 10957193811830801065UL) + ((uint64_t)op[3] * 1806799135077294067UL) + ((((uint64_t)op[4] * 6025982956261339846UL) + ((uint64_t)op[5] * 3826797987146243823UL) + ((uint64_t)op[6] * 3402412088961204283UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 3402412088961204283UL) + ((uint64_t)op[1] * 3693034966571312981UL) + ((uint64_t)op[2] * 7933916339959028336UL) + ((uint64_t)op[3] * 10957193811830801065UL) + ((uint64_t)op[4] * 1806799135077294067UL) + ((((uint64_t)op[5] * 6025982956261339846UL) + ((uint64_t)op[6] * 3826797987146243823UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 3826797987146243823UL) + ((uint64_t)op[1] * 3402412088961204283UL) + ((uint64_t)op[2] * 3693034966571312981UL) + ((uint64_t)op[3] * 7933916339959028336UL) + ((uint64_t)op[4] * 10957193811830801065UL) + ((uint64_t)op[5] * 1806799135077294067UL) + ((uint64_t)op[6] * 6763573366112404002UL);
	tmp_q[6] = ((uint64_t)op[0] * 6025982956261339846UL) + ((uint64_t)op[1] * 3826797987146243823UL) + ((uint64_t)op[2] * 3402412088961204283UL) + ((uint64_t)op[3] * 3693034966571312981UL) + ((uint64_t)op[4] * 7933916339959028336UL) + ((uint64_t)op[5] * 10957193811830801065UL) + ((uint64_t)op[6] * 1806799135077294067UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 30784616384L) - ((-((int128)tmp_q[1] * 1399666676L) + ((int128)tmp_q[2] * 43119343908L) - ((int128)tmp_q[3] * 6983708235L) + ((int128)tmp_q[4] * 30308902191L) - ((int128)tmp_q[5] * 12495349839L) - ((int128)tmp_q[6] * 39726739874L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 39726739874L) - ((int128)tmp_q[1] * 30784616384L) - ((-((int128)tmp_q[2] * 1399666676L) + ((int128)tmp_q[3] * 43119343908L) - ((int128)tmp_q[4] * 6983708235L) + ((int128)tmp_q[5] * 30308902191L) - ((int128)tmp_q[6] * 12495349839L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 12495349839L) - ((int128)tmp_q[1] * 39726739874L) - ((int128)tmp_q[2] * 30784616384L) - ((-((int128)tmp_q[3] * 1399666676L) + ((int128)tmp_q[4] * 43119343908L) - ((int128)tmp_q[5] * 6983708235L) + ((int128)tmp_q[6] * 30308902191L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 30308902191L) - ((int128)tmp_q[1] * 12495349839L) - ((int128)tmp_q[2] * 39726739874L) - ((int128)tmp_q[3] * 30784616384L) - ((-((int128)tmp_q[4] * 1399666676L) + ((int128)tmp_q[5] * 43119343908L) - ((int128)tmp_q[6] * 6983708235L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 6983708235L) + ((int128)tmp_q[1] * 30308902191L) - ((int128)tmp_q[2] * 12495349839L) - ((int128)tmp_q[3] * 39726739874L) - ((int128)tmp_q[4] * 30784616384L) - ((-((int128)tmp_q[5] * 1399666676L) + ((int128)tmp_q[6] * 43119343908L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 43119343908L) - ((int128)tmp_q[1] * 6983708235L) + ((int128)tmp_q[2] * 30308902191L) - ((int128)tmp_q[3] * 12495349839L) - ((int128)tmp_q[4] * 39726739874L) - ((int128)tmp_q[5] * 30784616384L) + ((int128)tmp_q[6] * 6998333380L);
	tmp_zero[6] = -((int128)tmp_q[0] * 1399666676L) + ((int128)tmp_q[1] * 43119343908L) - ((int128)tmp_q[2] * 6983708235L) + ((int128)tmp_q[3] * 30308902191L) - ((int128)tmp_q[4] * 12495349839L) - ((int128)tmp_q[5] * 39726739874L) - ((int128)tmp_q[6] * 30784616384L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

