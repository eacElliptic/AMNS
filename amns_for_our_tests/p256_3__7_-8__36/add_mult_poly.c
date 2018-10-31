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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17741252984987301723UL) + ((((uint64_t)op[1] * 16789867694554058763UL) + ((uint64_t)op[2] * 12942357375628860052UL) + ((uint64_t)op[3] * 4619672982504592877UL) + ((uint64_t)op[4] * 16206664422495578326UL) + ((uint64_t)op[5] * 2704444378802727317UL) + ((uint64_t)op[6] * 14366363203566896811UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 14366363203566896811UL) + ((uint64_t)op[1] * 17741252984987301723UL) + ((((uint64_t)op[2] * 16789867694554058763UL) + ((uint64_t)op[3] * 12942357375628860052UL) + ((uint64_t)op[4] * 4619672982504592877UL) + ((uint64_t)op[5] * 16206664422495578326UL) + ((uint64_t)op[6] * 2704444378802727317UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 2704444378802727317UL) + ((uint64_t)op[1] * 14366363203566896811UL) + ((uint64_t)op[2] * 17741252984987301723UL) + ((((uint64_t)op[3] * 16789867694554058763UL) + ((uint64_t)op[4] * 12942357375628860052UL) + ((uint64_t)op[5] * 4619672982504592877UL) + ((uint64_t)op[6] * 16206664422495578326UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 16206664422495578326UL) + ((uint64_t)op[1] * 2704444378802727317UL) + ((uint64_t)op[2] * 14366363203566896811UL) + ((uint64_t)op[3] * 17741252984987301723UL) + ((((uint64_t)op[4] * 16789867694554058763UL) + ((uint64_t)op[5] * 12942357375628860052UL) + ((uint64_t)op[6] * 4619672982504592877UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 4619672982504592877UL) + ((uint64_t)op[1] * 16206664422495578326UL) + ((uint64_t)op[2] * 2704444378802727317UL) + ((uint64_t)op[3] * 14366363203566896811UL) + ((uint64_t)op[4] * 17741252984987301723UL) + ((((uint64_t)op[5] * 16789867694554058763UL) + ((uint64_t)op[6] * 12942357375628860052UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 12942357375628860052UL) + ((uint64_t)op[1] * 4619672982504592877UL) + ((uint64_t)op[2] * 16206664422495578326UL) + ((uint64_t)op[3] * 2704444378802727317UL) + ((uint64_t)op[4] * 14366363203566896811UL) + ((uint64_t)op[5] * 17741252984987301723UL) + ((uint64_t)op[6] * 13255011033243942824UL);
	tmp_q[6] = ((uint64_t)op[0] * 16789867694554058763UL) + ((uint64_t)op[1] * 12942357375628860052UL) + ((uint64_t)op[2] * 4619672982504592877UL) + ((uint64_t)op[3] * 16206664422495578326UL) + ((uint64_t)op[4] * 2704444378802727317UL) + ((uint64_t)op[5] * 14366363203566896811UL) + ((uint64_t)op[6] * 17741252984987301723UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 35531810069L) - ((-((int128)tmp_q[1] * 57370752371L) + ((int128)tmp_q[2] * 23473434272L) - ((int128)tmp_q[3] * 2349939494L) + ((int128)tmp_q[4] * 12729493047L) + ((int128)tmp_q[5] * 39901823882L) - ((int128)tmp_q[6] * 68221227405L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 68221227405L) + ((int128)tmp_q[1] * 35531810069L) - ((-((int128)tmp_q[2] * 57370752371L) + ((int128)tmp_q[3] * 23473434272L) - ((int128)tmp_q[4] * 2349939494L) + ((int128)tmp_q[5] * 12729493047L) + ((int128)tmp_q[6] * 39901823882L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 39901823882L) - ((int128)tmp_q[1] * 68221227405L) + ((int128)tmp_q[2] * 35531810069L) - ((-((int128)tmp_q[3] * 57370752371L) + ((int128)tmp_q[4] * 23473434272L) - ((int128)tmp_q[5] * 2349939494L) + ((int128)tmp_q[6] * 12729493047L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 12729493047L) + ((int128)tmp_q[1] * 39901823882L) - ((int128)tmp_q[2] * 68221227405L) + ((int128)tmp_q[3] * 35531810069L) - ((-((int128)tmp_q[4] * 57370752371L) + ((int128)tmp_q[5] * 23473434272L) - ((int128)tmp_q[6] * 2349939494L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 2349939494L) + ((int128)tmp_q[1] * 12729493047L) + ((int128)tmp_q[2] * 39901823882L) - ((int128)tmp_q[3] * 68221227405L) + ((int128)tmp_q[4] * 35531810069L) - ((-((int128)tmp_q[5] * 57370752371L) + ((int128)tmp_q[6] * 23473434272L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 23473434272L) - ((int128)tmp_q[1] * 2349939494L) + ((int128)tmp_q[2] * 12729493047L) + ((int128)tmp_q[3] * 39901823882L) - ((int128)tmp_q[4] * 68221227405L) + ((int128)tmp_q[5] * 35531810069L) + ((int128)tmp_q[6] * 458966018968L);
	tmp_zero[6] = -((int128)tmp_q[0] * 57370752371L) + ((int128)tmp_q[1] * 23473434272L) - ((int128)tmp_q[2] * 2349939494L) + ((int128)tmp_q[3] * 12729493047L) + ((int128)tmp_q[4] * 39901823882L) - ((int128)tmp_q[5] * 68221227405L) + ((int128)tmp_q[6] * 35531810069L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

