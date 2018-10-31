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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4189864291938278059UL) + ((((uint64_t)op[1] * 15543465530966380043UL) + ((uint64_t)op[2] * 11636010811982156669UL) + ((uint64_t)op[3] * 8793593957425316821UL) + ((uint64_t)op[4] * 3521235415284290228UL) + ((uint64_t)op[5] * 18278209585865011994UL) + ((uint64_t)op[6] * 5036860976926730749UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 5036860976926730749UL) + ((uint64_t)op[1] * 4189864291938278059UL) + ((((uint64_t)op[2] * 15543465530966380043UL) + ((uint64_t)op[3] * 11636010811982156669UL) + ((uint64_t)op[4] * 8793593957425316821UL) + ((uint64_t)op[5] * 3521235415284290228UL) + ((uint64_t)op[6] * 18278209585865011994UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 18278209585865011994UL) + ((uint64_t)op[1] * 5036860976926730749UL) + ((uint64_t)op[2] * 4189864291938278059UL) + ((((uint64_t)op[3] * 15543465530966380043UL) + ((uint64_t)op[4] * 11636010811982156669UL) + ((uint64_t)op[5] * 8793593957425316821UL) + ((uint64_t)op[6] * 3521235415284290228UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 3521235415284290228UL) + ((uint64_t)op[1] * 18278209585865011994UL) + ((uint64_t)op[2] * 5036860976926730749UL) + ((uint64_t)op[3] * 4189864291938278059UL) + ((((uint64_t)op[4] * 15543465530966380043UL) + ((uint64_t)op[5] * 11636010811982156669UL) + ((uint64_t)op[6] * 8793593957425316821UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 8793593957425316821UL) + ((uint64_t)op[1] * 3521235415284290228UL) + ((uint64_t)op[2] * 18278209585865011994UL) + ((uint64_t)op[3] * 5036860976926730749UL) + ((uint64_t)op[4] * 4189864291938278059UL) + ((((uint64_t)op[5] * 15543465530966380043UL) + ((uint64_t)op[6] * 11636010811982156669UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 11636010811982156669UL) + ((uint64_t)op[1] * 8793593957425316821UL) + ((uint64_t)op[2] * 3521235415284290228UL) + ((uint64_t)op[3] * 18278209585865011994UL) + ((uint64_t)op[4] * 5036860976926730749UL) + ((uint64_t)op[5] * 4189864291938278059UL) + ((uint64_t)op[6] * 3930351359993693751UL);
	tmp_q[6] = ((uint64_t)op[0] * 15543465530966380043UL) + ((uint64_t)op[1] * 11636010811982156669UL) + ((uint64_t)op[2] * 8793593957425316821UL) + ((uint64_t)op[3] * 3521235415284290228UL) + ((uint64_t)op[4] * 18278209585865011994UL) + ((uint64_t)op[5] * 5036860976926730749UL) + ((uint64_t)op[6] * 4189864291938278059UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 336963428L) + ((-((int128)tmp_q[1] * 6442875129L) - ((int128)tmp_q[2] * 29290342442L) - ((int128)tmp_q[3] * 53200200230L) - ((int128)tmp_q[4] * 810802379L) - ((int128)tmp_q[5] * 21301664748L) - ((int128)tmp_q[6] * 33158350591L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 33158350591L) + ((int128)tmp_q[1] * 336963428L) + ((-((int128)tmp_q[2] * 6442875129L) - ((int128)tmp_q[3] * 29290342442L) - ((int128)tmp_q[4] * 53200200230L) - ((int128)tmp_q[5] * 810802379L) - ((int128)tmp_q[6] * 21301664748L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 21301664748L) - ((int128)tmp_q[1] * 33158350591L) + ((int128)tmp_q[2] * 336963428L) + ((-((int128)tmp_q[3] * 6442875129L) - ((int128)tmp_q[4] * 29290342442L) - ((int128)tmp_q[5] * 53200200230L) - ((int128)tmp_q[6] * 810802379L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 810802379L) - ((int128)tmp_q[1] * 21301664748L) - ((int128)tmp_q[2] * 33158350591L) + ((int128)tmp_q[3] * 336963428L) + ((-((int128)tmp_q[4] * 6442875129L) - ((int128)tmp_q[5] * 29290342442L) - ((int128)tmp_q[6] * 53200200230L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 53200200230L) - ((int128)tmp_q[1] * 810802379L) - ((int128)tmp_q[2] * 21301664748L) - ((int128)tmp_q[3] * 33158350591L) + ((int128)tmp_q[4] * 336963428L) + ((-((int128)tmp_q[5] * 6442875129L) - ((int128)tmp_q[6] * 29290342442L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 29290342442L) - ((int128)tmp_q[1] * 53200200230L) - ((int128)tmp_q[2] * 810802379L) - ((int128)tmp_q[3] * 21301664748L) - ((int128)tmp_q[4] * 33158350591L) + ((int128)tmp_q[5] * 336963428L) - ((int128)tmp_q[6] * 32214375645L);
	tmp_zero[6] = -((int128)tmp_q[0] * 6442875129L) - ((int128)tmp_q[1] * 29290342442L) - ((int128)tmp_q[2] * 53200200230L) - ((int128)tmp_q[3] * 810802379L) - ((int128)tmp_q[4] * 21301664748L) - ((int128)tmp_q[5] * 33158350591L) + ((int128)tmp_q[6] * 336963428L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

