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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9522687483690698887UL) + ((((uint64_t)op[1] * 2456413463247746969UL) + ((uint64_t)op[2] * 4277247053648686754UL) + ((uint64_t)op[3] * 7909343745345571009UL) + ((uint64_t)op[4] * 15873790942262494968UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 15873790942262494968UL) + ((uint64_t)op[1] * 9522687483690698887UL) + ((((uint64_t)op[2] * 2456413463247746969UL) + ((uint64_t)op[3] * 4277247053648686754UL) + ((uint64_t)op[4] * 7909343745345571009UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 7909343745345571009UL) + ((uint64_t)op[1] * 15873790942262494968UL) + ((uint64_t)op[2] * 9522687483690698887UL) + ((((uint64_t)op[3] * 2456413463247746969UL) + ((uint64_t)op[4] * 4277247053648686754UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 4277247053648686754UL) + ((uint64_t)op[1] * 7909343745345571009UL) + ((uint64_t)op[2] * 15873790942262494968UL) + ((uint64_t)op[3] * 9522687483690698887UL) + ((uint64_t)op[4] * 3708263294223069802UL);
	tmp_q[4] = ((uint64_t)op[0] * 2456413463247746969UL) + ((uint64_t)op[1] * 4277247053648686754UL) + ((uint64_t)op[2] * 7909343745345571009UL) + ((uint64_t)op[3] * 15873790942262494968UL) + ((uint64_t)op[4] * 9522687483690698887UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 475393770985745L) - ((-((int128)tmp_q[1] * 424382979195938L) + ((int128)tmp_q[2] * 1747444043942532L) + ((int128)tmp_q[3] * 374704265627497L) - ((int128)tmp_q[4] * 1694860462403202L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 1694860462403202L) + ((int128)tmp_q[1] * 475393770985745L) - ((-((int128)tmp_q[2] * 424382979195938L) + ((int128)tmp_q[3] * 1747444043942532L) + ((int128)tmp_q[4] * 374704265627497L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 374704265627497L) - ((int128)tmp_q[1] * 1694860462403202L) + ((int128)tmp_q[2] * 475393770985745L) - ((-((int128)tmp_q[3] * 424382979195938L) + ((int128)tmp_q[4] * 1747444043942532L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 1747444043942532L) + ((int128)tmp_q[1] * 374704265627497L) - ((int128)tmp_q[2] * 1694860462403202L) + ((int128)tmp_q[3] * 475393770985745L) + ((int128)tmp_q[4] * 2546297875175628L);
	tmp_zero[4] = -((int128)tmp_q[0] * 424382979195938L) + ((int128)tmp_q[1] * 1747444043942532L) + ((int128)tmp_q[2] * 374704265627497L) - ((int128)tmp_q[3] * 1694860462403202L) + ((int128)tmp_q[4] * 475393770985745L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

