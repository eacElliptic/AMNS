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
	tmp_q[0] = ((uint64_t)op[0] * 16702812791252635676UL) + ((((uint64_t)op[1] * 13543199783177230578UL) + ((uint64_t)op[2] * 11343071444723754658UL) + ((uint64_t)op[3] * 11211854268670426922UL) + ((uint64_t)op[4] * 18327740045477195921UL) + ((uint64_t)op[5] * 9205789609956259801UL) + ((uint64_t)op[6] * 1871856922438044031UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 1871856922438044031UL) + ((uint64_t)op[1] * 16702812791252635676UL) + ((((uint64_t)op[2] * 13543199783177230578UL) + ((uint64_t)op[3] * 11343071444723754658UL) + ((uint64_t)op[4] * 11211854268670426922UL) + ((uint64_t)op[5] * 18327740045477195921UL) + ((uint64_t)op[6] * 9205789609956259801UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 9205789609956259801UL) + ((uint64_t)op[1] * 1871856922438044031UL) + ((uint64_t)op[2] * 16702812791252635676UL) + ((((uint64_t)op[3] * 13543199783177230578UL) + ((uint64_t)op[4] * 11343071444723754658UL) + ((uint64_t)op[5] * 11211854268670426922UL) + ((uint64_t)op[6] * 18327740045477195921UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 18327740045477195921UL) + ((uint64_t)op[1] * 9205789609956259801UL) + ((uint64_t)op[2] * 1871856922438044031UL) + ((uint64_t)op[3] * 16702812791252635676UL) + ((((uint64_t)op[4] * 13543199783177230578UL) + ((uint64_t)op[5] * 11343071444723754658UL) + ((uint64_t)op[6] * 11211854268670426922UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 11211854268670426922UL) + ((uint64_t)op[1] * 18327740045477195921UL) + ((uint64_t)op[2] * 9205789609956259801UL) + ((uint64_t)op[3] * 1871856922438044031UL) + ((uint64_t)op[4] * 16702812791252635676UL) + ((((uint64_t)op[5] * 13543199783177230578UL) + ((uint64_t)op[6] * 11343071444723754658UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 11343071444723754658UL) + ((uint64_t)op[1] * 11211854268670426922UL) + ((uint64_t)op[2] * 18327740045477195921UL) + ((uint64_t)op[3] * 9205789609956259801UL) + ((uint64_t)op[4] * 1871856922438044031UL) + ((uint64_t)op[5] * 16702812791252635676UL) + ((uint64_t)op[6] * 6070977378952053574UL);
	tmp_q[6] = ((uint64_t)op[0] * 13543199783177230578UL) + ((uint64_t)op[1] * 11343071444723754658UL) + ((uint64_t)op[2] * 11211854268670426922UL) + ((uint64_t)op[3] * 18327740045477195921UL) + ((uint64_t)op[4] * 9205789609956259801UL) + ((uint64_t)op[5] * 1871856922438044031UL) + ((uint64_t)op[6] * 16702812791252635676UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 11015746878L) - ((-((int128)tmp_q[1] * 46107464125L) - ((int128)tmp_q[2] * 12874411673L) - ((int128)tmp_q[3] * 39231469271L) + ((int128)tmp_q[4] * 22146362258L) + ((int128)tmp_q[5] * 33270474431L) - ((int128)tmp_q[6] * 16264130073L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 16264130073L) + ((int128)tmp_q[1] * 11015746878L) - ((-((int128)tmp_q[2] * 46107464125L) - ((int128)tmp_q[3] * 12874411673L) - ((int128)tmp_q[4] * 39231469271L) + ((int128)tmp_q[5] * 22146362258L) + ((int128)tmp_q[6] * 33270474431L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 33270474431L) - ((int128)tmp_q[1] * 16264130073L) + ((int128)tmp_q[2] * 11015746878L) - ((-((int128)tmp_q[3] * 46107464125L) - ((int128)tmp_q[4] * 12874411673L) - ((int128)tmp_q[5] * 39231469271L) + ((int128)tmp_q[6] * 22146362258L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 22146362258L) + ((int128)tmp_q[1] * 33270474431L) - ((int128)tmp_q[2] * 16264130073L) + ((int128)tmp_q[3] * 11015746878L) - ((-((int128)tmp_q[4] * 46107464125L) - ((int128)tmp_q[5] * 12874411673L) - ((int128)tmp_q[6] * 39231469271L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 39231469271L) + ((int128)tmp_q[1] * 22146362258L) + ((int128)tmp_q[2] * 33270474431L) - ((int128)tmp_q[3] * 16264130073L) + ((int128)tmp_q[4] * 11015746878L) - ((-((int128)tmp_q[5] * 46107464125L) - ((int128)tmp_q[6] * 12874411673L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 12874411673L) - ((int128)tmp_q[1] * 39231469271L) + ((int128)tmp_q[2] * 22146362258L) + ((int128)tmp_q[3] * 33270474431L) - ((int128)tmp_q[4] * 16264130073L) + ((int128)tmp_q[5] * 11015746878L) + ((int128)tmp_q[6] * 230537320625L);
	tmp_zero[6] = -((int128)tmp_q[0] * 46107464125L) - ((int128)tmp_q[1] * 12874411673L) - ((int128)tmp_q[2] * 39231469271L) + ((int128)tmp_q[3] * 22146362258L) + ((int128)tmp_q[4] * 33270474431L) - ((int128)tmp_q[5] * 16264130073L) + ((int128)tmp_q[6] * 11015746878L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

