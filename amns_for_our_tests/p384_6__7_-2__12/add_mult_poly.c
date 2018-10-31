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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 387652656273109645UL) + ((((uint64_t)op[1] * 14014798024086378517UL) + ((uint64_t)op[2] * 2127302728775053117UL) + ((uint64_t)op[3] * 631972674352485791UL) + ((uint64_t)op[4] * 10280807893939712498UL) + ((uint64_t)op[5] * 10376353615306438720UL) + ((uint64_t)op[6] * 1609334715803879078UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 1609334715803879078UL) + ((uint64_t)op[1] * 387652656273109645UL) + ((((uint64_t)op[2] * 14014798024086378517UL) + ((uint64_t)op[3] * 2127302728775053117UL) + ((uint64_t)op[4] * 631972674352485791UL) + ((uint64_t)op[5] * 10280807893939712498UL) + ((uint64_t)op[6] * 10376353615306438720UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 10376353615306438720UL) + ((uint64_t)op[1] * 1609334715803879078UL) + ((uint64_t)op[2] * 387652656273109645UL) + ((((uint64_t)op[3] * 14014798024086378517UL) + ((uint64_t)op[4] * 2127302728775053117UL) + ((uint64_t)op[5] * 631972674352485791UL) + ((uint64_t)op[6] * 10280807893939712498UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 10280807893939712498UL) + ((uint64_t)op[1] * 10376353615306438720UL) + ((uint64_t)op[2] * 1609334715803879078UL) + ((uint64_t)op[3] * 387652656273109645UL) + ((((uint64_t)op[4] * 14014798024086378517UL) + ((uint64_t)op[5] * 2127302728775053117UL) + ((uint64_t)op[6] * 631972674352485791UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 631972674352485791UL) + ((uint64_t)op[1] * 10280807893939712498UL) + ((uint64_t)op[2] * 10376353615306438720UL) + ((uint64_t)op[3] * 1609334715803879078UL) + ((uint64_t)op[4] * 387652656273109645UL) + ((((uint64_t)op[5] * 14014798024086378517UL) + ((uint64_t)op[6] * 2127302728775053117UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 2127302728775053117UL) + ((uint64_t)op[1] * 631972674352485791UL) + ((uint64_t)op[2] * 10280807893939712498UL) + ((uint64_t)op[3] * 10376353615306438720UL) + ((uint64_t)op[4] * 1609334715803879078UL) + ((uint64_t)op[5] * 387652656273109645UL) + ((uint64_t)op[6] * 8863892099246346198UL);
	tmp_q[6] = ((uint64_t)op[0] * 14014798024086378517UL) + ((uint64_t)op[1] * 2127302728775053117UL) + ((uint64_t)op[2] * 631972674352485791UL) + ((uint64_t)op[3] * 10280807893939712498UL) + ((uint64_t)op[4] * 10376353615306438720UL) + ((uint64_t)op[5] * 1609334715803879078UL) + ((uint64_t)op[6] * 387652656273109645UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 12567420566498443L) - ((-((int128)tmp_q[1] * 7718016092140613L) + ((int128)tmp_q[2] * 7488029599318645L) + ((int128)tmp_q[3] * 12669185237472679L) - ((int128)tmp_q[4] * 19645593250229748L) + ((int128)tmp_q[5] * 18215046445075816L) - ((int128)tmp_q[6] * 1082099879588100L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 1082099879588100L) + ((int128)tmp_q[1] * 12567420566498443L) - ((-((int128)tmp_q[2] * 7718016092140613L) + ((int128)tmp_q[3] * 7488029599318645L) + ((int128)tmp_q[4] * 12669185237472679L) - ((int128)tmp_q[5] * 19645593250229748L) + ((int128)tmp_q[6] * 18215046445075816L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 18215046445075816L) - ((int128)tmp_q[1] * 1082099879588100L) + ((int128)tmp_q[2] * 12567420566498443L) - ((-((int128)tmp_q[3] * 7718016092140613L) + ((int128)tmp_q[4] * 7488029599318645L) + ((int128)tmp_q[5] * 12669185237472679L) - ((int128)tmp_q[6] * 19645593250229748L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 19645593250229748L) + ((int128)tmp_q[1] * 18215046445075816L) - ((int128)tmp_q[2] * 1082099879588100L) + ((int128)tmp_q[3] * 12567420566498443L) - ((-((int128)tmp_q[4] * 7718016092140613L) + ((int128)tmp_q[5] * 7488029599318645L) + ((int128)tmp_q[6] * 12669185237472679L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 12669185237472679L) - ((int128)tmp_q[1] * 19645593250229748L) + ((int128)tmp_q[2] * 18215046445075816L) - ((int128)tmp_q[3] * 1082099879588100L) + ((int128)tmp_q[4] * 12567420566498443L) - ((-((int128)tmp_q[5] * 7718016092140613L) + ((int128)tmp_q[6] * 7488029599318645L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 7488029599318645L) + ((int128)tmp_q[1] * 12669185237472679L) - ((int128)tmp_q[2] * 19645593250229748L) + ((int128)tmp_q[3] * 18215046445075816L) - ((int128)tmp_q[4] * 1082099879588100L) + ((int128)tmp_q[5] * 12567420566498443L) + ((int128)tmp_q[6] * 15436032184281226L);
	tmp_zero[6] = -((int128)tmp_q[0] * 7718016092140613L) + ((int128)tmp_q[1] * 7488029599318645L) + ((int128)tmp_q[2] * 12669185237472679L) - ((int128)tmp_q[3] * 19645593250229748L) + ((int128)tmp_q[4] * 18215046445075816L) - ((int128)tmp_q[5] * 1082099879588100L) + ((int128)tmp_q[6] * 12567420566498443L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

