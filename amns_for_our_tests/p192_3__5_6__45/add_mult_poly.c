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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17142365947852960197UL) + ((((uint64_t)op[1] * 15457954402246966051UL) + ((uint64_t)op[2] * 2929536103038208709UL) + ((uint64_t)op[3] * 3157762988740752906UL) + ((uint64_t)op[4] * 11428244525180117033UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 11428244525180117033UL) + ((uint64_t)op[1] * 17142365947852960197UL) + ((((uint64_t)op[2] * 15457954402246966051UL) + ((uint64_t)op[3] * 2929536103038208709UL) + ((uint64_t)op[4] * 3157762988740752906UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 3157762988740752906UL) + ((uint64_t)op[1] * 11428244525180117033UL) + ((uint64_t)op[2] * 17142365947852960197UL) + ((((uint64_t)op[3] * 15457954402246966051UL) + ((uint64_t)op[4] * 2929536103038208709UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 2929536103038208709UL) + ((uint64_t)op[1] * 3157762988740752906UL) + ((uint64_t)op[2] * 11428244525180117033UL) + ((uint64_t)op[3] * 17142365947852960197UL) + ((uint64_t)op[4] * 514006044934038226UL);
	tmp_q[4] = ((uint64_t)op[0] * 15457954402246966051UL) + ((uint64_t)op[1] * 2929536103038208709UL) + ((uint64_t)op[2] * 3157762988740752906UL) + ((uint64_t)op[3] * 11428244525180117033UL) + ((uint64_t)op[4] * 17142365947852960197UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 29949319535L) + ((((int128)tmp_q[1] * 214113540890L) - ((int128)tmp_q[2] * 27436505512L) - ((int128)tmp_q[3] * 163768312837L) + ((int128)tmp_q[4] * 83970433263L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 83970433263L) + ((int128)tmp_q[1] * 29949319535L) + ((((int128)tmp_q[2] * 214113540890L) - ((int128)tmp_q[3] * 27436505512L) - ((int128)tmp_q[4] * 163768312837L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 163768312837L) + ((int128)tmp_q[1] * 83970433263L) + ((int128)tmp_q[2] * 29949319535L) + ((((int128)tmp_q[3] * 214113540890L) - ((int128)tmp_q[4] * 27436505512L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 27436505512L) - ((int128)tmp_q[1] * 163768312837L) + ((int128)tmp_q[2] * 83970433263L) + ((int128)tmp_q[3] * 29949319535L) + ((int128)tmp_q[4] * 1284681245340L);
	tmp_zero[4] = ((int128)tmp_q[0] * 214113540890L) - ((int128)tmp_q[1] * 27436505512L) - ((int128)tmp_q[2] * 163768312837L) + ((int128)tmp_q[3] * 83970433263L) + ((int128)tmp_q[4] * 29949319535L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

