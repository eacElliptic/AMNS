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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 12);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13951214725866049749UL) + ((((uint64_t)op[1] * 16283130561821061865UL) + ((uint64_t)op[2] * 16704168137851855070UL) + ((uint64_t)op[3] * 16891489960959778439UL) + ((uint64_t)op[4] * 6374952760598911734UL) + ((uint64_t)op[5] * 9388314692590074405UL) + ((uint64_t)op[6] * 3014950665452689330UL) + ((uint64_t)op[7] * 9781462088054834467UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 9781462088054834467UL) + ((uint64_t)op[1] * 13951214725866049749UL) + ((((uint64_t)op[2] * 16283130561821061865UL) + ((uint64_t)op[3] * 16704168137851855070UL) + ((uint64_t)op[4] * 16891489960959778439UL) + ((uint64_t)op[5] * 6374952760598911734UL) + ((uint64_t)op[6] * 9388314692590074405UL) + ((uint64_t)op[7] * 3014950665452689330UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 3014950665452689330UL) + ((uint64_t)op[1] * 9781462088054834467UL) + ((uint64_t)op[2] * 13951214725866049749UL) + ((((uint64_t)op[3] * 16283130561821061865UL) + ((uint64_t)op[4] * 16704168137851855070UL) + ((uint64_t)op[5] * 16891489960959778439UL) + ((uint64_t)op[6] * 6374952760598911734UL) + ((uint64_t)op[7] * 9388314692590074405UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 9388314692590074405UL) + ((uint64_t)op[1] * 3014950665452689330UL) + ((uint64_t)op[2] * 9781462088054834467UL) + ((uint64_t)op[3] * 13951214725866049749UL) + ((((uint64_t)op[4] * 16283130561821061865UL) + ((uint64_t)op[5] * 16704168137851855070UL) + ((uint64_t)op[6] * 16891489960959778439UL) + ((uint64_t)op[7] * 6374952760598911734UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 6374952760598911734UL) + ((uint64_t)op[1] * 9388314692590074405UL) + ((uint64_t)op[2] * 3014950665452689330UL) + ((uint64_t)op[3] * 9781462088054834467UL) + ((uint64_t)op[4] * 13951214725866049749UL) + ((((uint64_t)op[5] * 16283130561821061865UL) + ((uint64_t)op[6] * 16704168137851855070UL) + ((uint64_t)op[7] * 16891489960959778439UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 16891489960959778439UL) + ((uint64_t)op[1] * 6374952760598911734UL) + ((uint64_t)op[2] * 9388314692590074405UL) + ((uint64_t)op[3] * 3014950665452689330UL) + ((uint64_t)op[4] * 9781462088054834467UL) + ((uint64_t)op[5] * 13951214725866049749UL) + ((((uint64_t)op[6] * 16283130561821061865UL) + ((uint64_t)op[7] * 16704168137851855070UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 16704168137851855070UL) + ((uint64_t)op[1] * 16891489960959778439UL) + ((uint64_t)op[2] * 6374952760598911734UL) + ((uint64_t)op[3] * 9388314692590074405UL) + ((uint64_t)op[4] * 3014950665452689330UL) + ((uint64_t)op[5] * 9781462088054834467UL) + ((uint64_t)op[6] * 13951214725866049749UL) + ((uint64_t)op[7] * 12981681071330938506UL);
	tmp_q[7] = ((uint64_t)op[0] * 16283130561821061865UL) + ((uint64_t)op[1] * 16704168137851855070UL) + ((uint64_t)op[2] * 16891489960959778439UL) + ((uint64_t)op[3] * 6374952760598911734UL) + ((uint64_t)op[4] * 9388314692590074405UL) + ((uint64_t)op[5] * 3014950665452689330UL) + ((uint64_t)op[6] * 9781462088054834467UL) + ((uint64_t)op[7] * 13951214725866049749UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 91120044100035L) - ((((int128)tmp_q[1] * 27791613020141L) - ((int128)tmp_q[2] * 146966674447752L) + ((int128)tmp_q[3] * 139138105595901L) + ((int128)tmp_q[4] * 44867626039633L) - ((int128)tmp_q[5] * 97450480005170L) + ((int128)tmp_q[6] * 54224355166951L) + ((int128)tmp_q[7] * 4066024616241L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 4066024616241L) - ((int128)tmp_q[1] * 91120044100035L) - ((((int128)tmp_q[2] * 27791613020141L) - ((int128)tmp_q[3] * 146966674447752L) + ((int128)tmp_q[4] * 139138105595901L) + ((int128)tmp_q[5] * 44867626039633L) - ((int128)tmp_q[6] * 97450480005170L) + ((int128)tmp_q[7] * 54224355166951L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 54224355166951L) + ((int128)tmp_q[1] * 4066024616241L) - ((int128)tmp_q[2] * 91120044100035L) - ((((int128)tmp_q[3] * 27791613020141L) - ((int128)tmp_q[4] * 146966674447752L) + ((int128)tmp_q[5] * 139138105595901L) + ((int128)tmp_q[6] * 44867626039633L) - ((int128)tmp_q[7] * 97450480005170L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 97450480005170L) + ((int128)tmp_q[1] * 54224355166951L) + ((int128)tmp_q[2] * 4066024616241L) - ((int128)tmp_q[3] * 91120044100035L) - ((((int128)tmp_q[4] * 27791613020141L) - ((int128)tmp_q[5] * 146966674447752L) + ((int128)tmp_q[6] * 139138105595901L) + ((int128)tmp_q[7] * 44867626039633L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 44867626039633L) - ((int128)tmp_q[1] * 97450480005170L) + ((int128)tmp_q[2] * 54224355166951L) + ((int128)tmp_q[3] * 4066024616241L) - ((int128)tmp_q[4] * 91120044100035L) - ((((int128)tmp_q[5] * 27791613020141L) - ((int128)tmp_q[6] * 146966674447752L) + ((int128)tmp_q[7] * 139138105595901L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 139138105595901L) + ((int128)tmp_q[1] * 44867626039633L) - ((int128)tmp_q[2] * 97450480005170L) + ((int128)tmp_q[3] * 54224355166951L) + ((int128)tmp_q[4] * 4066024616241L) - ((int128)tmp_q[5] * 91120044100035L) - ((((int128)tmp_q[6] * 27791613020141L) - ((int128)tmp_q[7] * 146966674447752L)) * 6);
	tmp_zero[6] = -((int128)tmp_q[0] * 146966674447752L) + ((int128)tmp_q[1] * 139138105595901L) + ((int128)tmp_q[2] * 44867626039633L) - ((int128)tmp_q[3] * 97450480005170L) + ((int128)tmp_q[4] * 54224355166951L) + ((int128)tmp_q[5] * 4066024616241L) - ((int128)tmp_q[6] * 91120044100035L) - ((int128)tmp_q[7] * 166749678120846L);
	tmp_zero[7] = ((int128)tmp_q[0] * 27791613020141L) - ((int128)tmp_q[1] * 146966674447752L) + ((int128)tmp_q[2] * 139138105595901L) + ((int128)tmp_q[3] * 44867626039633L) - ((int128)tmp_q[4] * 97450480005170L) + ((int128)tmp_q[5] * 54224355166951L) + ((int128)tmp_q[6] * 4066024616241L) - ((int128)tmp_q[7] * 91120044100035L);

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

