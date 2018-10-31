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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2044177889323124265UL) + ((((uint64_t)op[1] * 9434368364614778902UL) + ((uint64_t)op[2] * 1581434713025789093UL) + ((uint64_t)op[3] * 17018685431747127981UL) + ((uint64_t)op[4] * 8441264421428323450UL) + ((uint64_t)op[5] * 15384500294950694730UL) + ((uint64_t)op[6] * 10500609018878516616UL) + ((uint64_t)op[7] * 1586173893041583994UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 1586173893041583994UL) + ((uint64_t)op[1] * 2044177889323124265UL) + ((((uint64_t)op[2] * 9434368364614778902UL) + ((uint64_t)op[3] * 1581434713025789093UL) + ((uint64_t)op[4] * 17018685431747127981UL) + ((uint64_t)op[5] * 8441264421428323450UL) + ((uint64_t)op[6] * 15384500294950694730UL) + ((uint64_t)op[7] * 10500609018878516616UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 10500609018878516616UL) + ((uint64_t)op[1] * 1586173893041583994UL) + ((uint64_t)op[2] * 2044177889323124265UL) + ((((uint64_t)op[3] * 9434368364614778902UL) + ((uint64_t)op[4] * 1581434713025789093UL) + ((uint64_t)op[5] * 17018685431747127981UL) + ((uint64_t)op[6] * 8441264421428323450UL) + ((uint64_t)op[7] * 15384500294950694730UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 15384500294950694730UL) + ((uint64_t)op[1] * 10500609018878516616UL) + ((uint64_t)op[2] * 1586173893041583994UL) + ((uint64_t)op[3] * 2044177889323124265UL) + ((((uint64_t)op[4] * 9434368364614778902UL) + ((uint64_t)op[5] * 1581434713025789093UL) + ((uint64_t)op[6] * 17018685431747127981UL) + ((uint64_t)op[7] * 8441264421428323450UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 8441264421428323450UL) + ((uint64_t)op[1] * 15384500294950694730UL) + ((uint64_t)op[2] * 10500609018878516616UL) + ((uint64_t)op[3] * 1586173893041583994UL) + ((uint64_t)op[4] * 2044177889323124265UL) + ((((uint64_t)op[5] * 9434368364614778902UL) + ((uint64_t)op[6] * 1581434713025789093UL) + ((uint64_t)op[7] * 17018685431747127981UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 17018685431747127981UL) + ((uint64_t)op[1] * 8441264421428323450UL) + ((uint64_t)op[2] * 15384500294950694730UL) + ((uint64_t)op[3] * 10500609018878516616UL) + ((uint64_t)op[4] * 1586173893041583994UL) + ((uint64_t)op[5] * 2044177889323124265UL) + ((((uint64_t)op[6] * 9434368364614778902UL) + ((uint64_t)op[7] * 1581434713025789093UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 1581434713025789093UL) + ((uint64_t)op[1] * 17018685431747127981UL) + ((uint64_t)op[2] * 8441264421428323450UL) + ((uint64_t)op[3] * 15384500294950694730UL) + ((uint64_t)op[4] * 10500609018878516616UL) + ((uint64_t)op[5] * 1586173893041583994UL) + ((uint64_t)op[6] * 2044177889323124265UL) + ((uint64_t)op[7] * 8590383053574766526UL);
	tmp_q[7] = ((uint64_t)op[0] * 9434368364614778902UL) + ((uint64_t)op[1] * 1581434713025789093UL) + ((uint64_t)op[2] * 17018685431747127981UL) + ((uint64_t)op[3] * 8441264421428323450UL) + ((uint64_t)op[4] * 15384500294950694730UL) + ((uint64_t)op[5] * 10500609018878516616UL) + ((uint64_t)op[6] * 1586173893041583994UL) + ((uint64_t)op[7] * 2044177889323124265UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 62768984328265L) - ((-((int128)tmp_q[1] * 100913590384774L) - ((int128)tmp_q[2] * 101166206996268L) - ((int128)tmp_q[3] * 164503040053101L) - ((int128)tmp_q[4] * 134866422782950L) - ((int128)tmp_q[5] * 42161002891739L) - ((int128)tmp_q[6] * 150277051644309L) + ((int128)tmp_q[7] * 27919231025627L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 27919231025627L) - ((int128)tmp_q[1] * 62768984328265L) - ((-((int128)tmp_q[2] * 100913590384774L) - ((int128)tmp_q[3] * 101166206996268L) - ((int128)tmp_q[4] * 164503040053101L) - ((int128)tmp_q[5] * 134866422782950L) - ((int128)tmp_q[6] * 42161002891739L) - ((int128)tmp_q[7] * 150277051644309L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 150277051644309L) + ((int128)tmp_q[1] * 27919231025627L) - ((int128)tmp_q[2] * 62768984328265L) - ((-((int128)tmp_q[3] * 100913590384774L) - ((int128)tmp_q[4] * 101166206996268L) - ((int128)tmp_q[5] * 164503040053101L) - ((int128)tmp_q[6] * 134866422782950L) - ((int128)tmp_q[7] * 42161002891739L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 42161002891739L) - ((int128)tmp_q[1] * 150277051644309L) + ((int128)tmp_q[2] * 27919231025627L) - ((int128)tmp_q[3] * 62768984328265L) - ((-((int128)tmp_q[4] * 100913590384774L) - ((int128)tmp_q[5] * 101166206996268L) - ((int128)tmp_q[6] * 164503040053101L) - ((int128)tmp_q[7] * 134866422782950L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 134866422782950L) - ((int128)tmp_q[1] * 42161002891739L) - ((int128)tmp_q[2] * 150277051644309L) + ((int128)tmp_q[3] * 27919231025627L) - ((int128)tmp_q[4] * 62768984328265L) - ((-((int128)tmp_q[5] * 100913590384774L) - ((int128)tmp_q[6] * 101166206996268L) - ((int128)tmp_q[7] * 164503040053101L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 164503040053101L) - ((int128)tmp_q[1] * 134866422782950L) - ((int128)tmp_q[2] * 42161002891739L) - ((int128)tmp_q[3] * 150277051644309L) + ((int128)tmp_q[4] * 27919231025627L) - ((int128)tmp_q[5] * 62768984328265L) - ((-((int128)tmp_q[6] * 100913590384774L) - ((int128)tmp_q[7] * 101166206996268L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 101166206996268L) - ((int128)tmp_q[1] * 164503040053101L) - ((int128)tmp_q[2] * 134866422782950L) - ((int128)tmp_q[3] * 42161002891739L) - ((int128)tmp_q[4] * 150277051644309L) + ((int128)tmp_q[5] * 27919231025627L) - ((int128)tmp_q[6] * 62768984328265L) + ((int128)tmp_q[7] * 302740771154322L);
	tmp_zero[7] = -((int128)tmp_q[0] * 100913590384774L) - ((int128)tmp_q[1] * 101166206996268L) - ((int128)tmp_q[2] * 164503040053101L) - ((int128)tmp_q[3] * 134866422782950L) - ((int128)tmp_q[4] * 42161002891739L) - ((int128)tmp_q[5] * 150277051644309L) + ((int128)tmp_q[6] * 27919231025627L) - ((int128)tmp_q[7] * 62768984328265L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

