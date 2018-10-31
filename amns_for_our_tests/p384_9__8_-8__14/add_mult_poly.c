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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) << 4);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6976739322967624631UL) + ((((uint64_t)op[1] * 15161254208378211635UL) + ((uint64_t)op[2] * 16851541736958624525UL) + ((uint64_t)op[3] * 2456163355700635475UL) + ((uint64_t)op[4] * 8317837961600810275UL) + ((uint64_t)op[5] * 16776140340805393080UL) + ((uint64_t)op[6] * 8696446885227025071UL) + ((uint64_t)op[7] * 3142332157316163300UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 3142332157316163300UL) + ((uint64_t)op[1] * 6976739322967624631UL) + ((((uint64_t)op[2] * 15161254208378211635UL) + ((uint64_t)op[3] * 16851541736958624525UL) + ((uint64_t)op[4] * 2456163355700635475UL) + ((uint64_t)op[5] * 8317837961600810275UL) + ((uint64_t)op[6] * 16776140340805393080UL) + ((uint64_t)op[7] * 8696446885227025071UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 8696446885227025071UL) + ((uint64_t)op[1] * 3142332157316163300UL) + ((uint64_t)op[2] * 6976739322967624631UL) + ((((uint64_t)op[3] * 15161254208378211635UL) + ((uint64_t)op[4] * 16851541736958624525UL) + ((uint64_t)op[5] * 2456163355700635475UL) + ((uint64_t)op[6] * 8317837961600810275UL) + ((uint64_t)op[7] * 16776140340805393080UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 16776140340805393080UL) + ((uint64_t)op[1] * 8696446885227025071UL) + ((uint64_t)op[2] * 3142332157316163300UL) + ((uint64_t)op[3] * 6976739322967624631UL) + ((((uint64_t)op[4] * 15161254208378211635UL) + ((uint64_t)op[5] * 16851541736958624525UL) + ((uint64_t)op[6] * 2456163355700635475UL) + ((uint64_t)op[7] * 8317837961600810275UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 8317837961600810275UL) + ((uint64_t)op[1] * 16776140340805393080UL) + ((uint64_t)op[2] * 8696446885227025071UL) + ((uint64_t)op[3] * 3142332157316163300UL) + ((uint64_t)op[4] * 6976739322967624631UL) + ((((uint64_t)op[5] * 15161254208378211635UL) + ((uint64_t)op[6] * 16851541736958624525UL) + ((uint64_t)op[7] * 2456163355700635475UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 2456163355700635475UL) + ((uint64_t)op[1] * 8317837961600810275UL) + ((uint64_t)op[2] * 16776140340805393080UL) + ((uint64_t)op[3] * 8696446885227025071UL) + ((uint64_t)op[4] * 3142332157316163300UL) + ((uint64_t)op[5] * 6976739322967624631UL) + ((((uint64_t)op[6] * 15161254208378211635UL) + ((uint64_t)op[7] * 16851541736958624525UL)) * 18446744073709551608);
	tmp_q[6] = ((uint64_t)op[0] * 16851541736958624525UL) + ((uint64_t)op[1] * 2456163355700635475UL) + ((uint64_t)op[2] * 8317837961600810275UL) + ((uint64_t)op[3] * 16776140340805393080UL) + ((uint64_t)op[4] * 8696446885227025071UL) + ((uint64_t)op[5] * 3142332157316163300UL) + ((uint64_t)op[6] * 6976739322967624631UL) + ((uint64_t)op[7] * 7837174848941168232UL);
	tmp_q[7] = ((uint64_t)op[0] * 15161254208378211635UL) + ((uint64_t)op[1] * 16851541736958624525UL) + ((uint64_t)op[2] * 2456163355700635475UL) + ((uint64_t)op[3] * 8317837961600810275UL) + ((uint64_t)op[4] * 16776140340805393080UL) + ((uint64_t)op[5] * 8696446885227025071UL) + ((uint64_t)op[6] * 3142332157316163300UL) + ((uint64_t)op[7] * 6976739322967624631UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 114726571827343L) - ((((int128)tmp_q[1] * 112665561744861L) + ((int128)tmp_q[2] * 67070980661838L) - ((int128)tmp_q[3] * 7159080553849L) - ((int128)tmp_q[4] * 33869703731660L) - ((int128)tmp_q[5] * 46692318679144L) + ((int128)tmp_q[6] * 75940280883295L) + ((int128)tmp_q[7] * 110653332342524L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 110653332342524L) - ((int128)tmp_q[1] * 114726571827343L) - ((((int128)tmp_q[2] * 112665561744861L) + ((int128)tmp_q[3] * 67070980661838L) - ((int128)tmp_q[4] * 7159080553849L) - ((int128)tmp_q[5] * 33869703731660L) - ((int128)tmp_q[6] * 46692318679144L) + ((int128)tmp_q[7] * 75940280883295L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 75940280883295L) + ((int128)tmp_q[1] * 110653332342524L) - ((int128)tmp_q[2] * 114726571827343L) - ((((int128)tmp_q[3] * 112665561744861L) + ((int128)tmp_q[4] * 67070980661838L) - ((int128)tmp_q[5] * 7159080553849L) - ((int128)tmp_q[6] * 33869703731660L) - ((int128)tmp_q[7] * 46692318679144L)) * 8);
	tmp_zero[3] = -((int128)tmp_q[0] * 46692318679144L) + ((int128)tmp_q[1] * 75940280883295L) + ((int128)tmp_q[2] * 110653332342524L) - ((int128)tmp_q[3] * 114726571827343L) - ((((int128)tmp_q[4] * 112665561744861L) + ((int128)tmp_q[5] * 67070980661838L) - ((int128)tmp_q[6] * 7159080553849L) - ((int128)tmp_q[7] * 33869703731660L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 33869703731660L) - ((int128)tmp_q[1] * 46692318679144L) + ((int128)tmp_q[2] * 75940280883295L) + ((int128)tmp_q[3] * 110653332342524L) - ((int128)tmp_q[4] * 114726571827343L) - ((((int128)tmp_q[5] * 112665561744861L) + ((int128)tmp_q[6] * 67070980661838L) - ((int128)tmp_q[7] * 7159080553849L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 7159080553849L) - ((int128)tmp_q[1] * 33869703731660L) - ((int128)tmp_q[2] * 46692318679144L) + ((int128)tmp_q[3] * 75940280883295L) + ((int128)tmp_q[4] * 110653332342524L) - ((int128)tmp_q[5] * 114726571827343L) - ((((int128)tmp_q[6] * 112665561744861L) + ((int128)tmp_q[7] * 67070980661838L)) * 8);
	tmp_zero[6] = ((int128)tmp_q[0] * 67070980661838L) - ((int128)tmp_q[1] * 7159080553849L) - ((int128)tmp_q[2] * 33869703731660L) - ((int128)tmp_q[3] * 46692318679144L) + ((int128)tmp_q[4] * 75940280883295L) + ((int128)tmp_q[5] * 110653332342524L) - ((int128)tmp_q[6] * 114726571827343L) - ((int128)tmp_q[7] * 901324493958888L);
	tmp_zero[7] = ((int128)tmp_q[0] * 112665561744861L) + ((int128)tmp_q[1] * 67070980661838L) - ((int128)tmp_q[2] * 7159080553849L) - ((int128)tmp_q[3] * 33869703731660L) - ((int128)tmp_q[4] * 46692318679144L) + ((int128)tmp_q[5] * 75940280883295L) + ((int128)tmp_q[6] * 110653332342524L) - ((int128)tmp_q[7] * 114726571827343L);

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

