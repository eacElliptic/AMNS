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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2239479898716951448UL) + ((((uint64_t)op[1] * 8970591127573584764UL) + ((uint64_t)op[2] * 12849825779777081443UL) + ((uint64_t)op[3] * 1819535322181279505UL) + ((uint64_t)op[4] * 14587184225417764878UL) + ((uint64_t)op[5] * 11298008975319909137UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 11298008975319909137UL) + ((uint64_t)op[1] * 2239479898716951448UL) + ((((uint64_t)op[2] * 8970591127573584764UL) + ((uint64_t)op[3] * 12849825779777081443UL) + ((uint64_t)op[4] * 1819535322181279505UL) + ((uint64_t)op[5] * 14587184225417764878UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 14587184225417764878UL) + ((uint64_t)op[1] * 11298008975319909137UL) + ((uint64_t)op[2] * 2239479898716951448UL) + ((((uint64_t)op[3] * 8970591127573584764UL) + ((uint64_t)op[4] * 12849825779777081443UL) + ((uint64_t)op[5] * 1819535322181279505UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 1819535322181279505UL) + ((uint64_t)op[1] * 14587184225417764878UL) + ((uint64_t)op[2] * 11298008975319909137UL) + ((uint64_t)op[3] * 2239479898716951448UL) + ((((uint64_t)op[4] * 8970591127573584764UL) + ((uint64_t)op[5] * 12849825779777081443UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 12849825779777081443UL) + ((uint64_t)op[1] * 1819535322181279505UL) + ((uint64_t)op[2] * 14587184225417764878UL) + ((uint64_t)op[3] * 11298008975319909137UL) + ((uint64_t)op[4] * 2239479898716951448UL) + ((uint64_t)op[5] * 10992838401823113116UL);
	tmp_q[5] = ((uint64_t)op[0] * 8970591127573584764UL) + ((uint64_t)op[1] * 12849825779777081443UL) + ((uint64_t)op[2] * 1819535322181279505UL) + ((uint64_t)op[3] * 14587184225417764878UL) + ((uint64_t)op[4] * 11298008975319909137UL) + ((uint64_t)op[5] * 2239479898716951448UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3494066825160L) - ((((int128)tmp_q[1] * 191080375670L) + ((int128)tmp_q[2] * 2344811361903L) - ((int128)tmp_q[3] * 3096993040993L) + ((int128)tmp_q[4] * 1060340560620L) + ((int128)tmp_q[5] * 2317601535825L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 2317601535825L) + ((int128)tmp_q[1] * 3494066825160L) - ((((int128)tmp_q[2] * 191080375670L) + ((int128)tmp_q[3] * 2344811361903L) - ((int128)tmp_q[4] * 3096993040993L) + ((int128)tmp_q[5] * 1060340560620L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 1060340560620L) + ((int128)tmp_q[1] * 2317601535825L) + ((int128)tmp_q[2] * 3494066825160L) - ((((int128)tmp_q[3] * 191080375670L) + ((int128)tmp_q[4] * 2344811361903L) - ((int128)tmp_q[5] * 3096993040993L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 3096993040993L) + ((int128)tmp_q[1] * 1060340560620L) + ((int128)tmp_q[2] * 2317601535825L) + ((int128)tmp_q[3] * 3494066825160L) - ((((int128)tmp_q[4] * 191080375670L) + ((int128)tmp_q[5] * 2344811361903L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 2344811361903L) - ((int128)tmp_q[1] * 3096993040993L) + ((int128)tmp_q[2] * 1060340560620L) + ((int128)tmp_q[3] * 2317601535825L) + ((int128)tmp_q[4] * 3494066825160L) - ((int128)tmp_q[5] * 1337562629690L);
	tmp_zero[5] = ((int128)tmp_q[0] * 191080375670L) + ((int128)tmp_q[1] * 2344811361903L) - ((int128)tmp_q[2] * 3096993040993L) + ((int128)tmp_q[3] * 1060340560620L) + ((int128)tmp_q[4] * 2317601535825L) + ((int128)tmp_q[5] * 3494066825160L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

