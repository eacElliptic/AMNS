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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - ((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - ((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - ((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - ((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - ((int128)pa[5] * pb[5]);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - ((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - ((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - ((int128)pa[5] * pa[5]);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7515043656170482257UL) + ((((uint64_t)op[1] * 2041833301930030575UL) + ((uint64_t)op[2] * 6557903897941956082UL) + ((uint64_t)op[3] * 15913525380690250377UL) + ((uint64_t)op[4] * 14072947554112438338UL) + ((uint64_t)op[5] * 13871692078760219802UL)) * 18446744073709551615);
	tmp_q[1] = ((uint64_t)op[0] * 13871692078760219802UL) + ((uint64_t)op[1] * 7515043656170482257UL) + ((((uint64_t)op[2] * 2041833301930030575UL) + ((uint64_t)op[3] * 6557903897941956082UL) + ((uint64_t)op[4] * 15913525380690250377UL) + ((uint64_t)op[5] * 14072947554112438338UL)) * 18446744073709551615);
	tmp_q[2] = ((uint64_t)op[0] * 14072947554112438338UL) + ((uint64_t)op[1] * 13871692078760219802UL) + ((uint64_t)op[2] * 7515043656170482257UL) + ((((uint64_t)op[3] * 2041833301930030575UL) + ((uint64_t)op[4] * 6557903897941956082UL) + ((uint64_t)op[5] * 15913525380690250377UL)) * 18446744073709551615);
	tmp_q[3] = ((uint64_t)op[0] * 15913525380690250377UL) + ((uint64_t)op[1] * 14072947554112438338UL) + ((uint64_t)op[2] * 13871692078760219802UL) + ((uint64_t)op[3] * 7515043656170482257UL) + ((((uint64_t)op[4] * 2041833301930030575UL) + ((uint64_t)op[5] * 6557903897941956082UL)) * 18446744073709551615);
	tmp_q[4] = ((uint64_t)op[0] * 6557903897941956082UL) + ((uint64_t)op[1] * 15913525380690250377UL) + ((uint64_t)op[2] * 14072947554112438338UL) + ((uint64_t)op[3] * 13871692078760219802UL) + ((uint64_t)op[4] * 7515043656170482257UL) + ((uint64_t)op[5] * 16404910771779521041UL);
	tmp_q[5] = ((uint64_t)op[0] * 2041833301930030575UL) + ((uint64_t)op[1] * 6557903897941956082UL) + ((uint64_t)op[2] * 15913525380690250377UL) + ((uint64_t)op[3] * 14072947554112438338UL) + ((uint64_t)op[4] * 13871692078760219802UL) + ((uint64_t)op[5] * 7515043656170482257UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 87714942664860L) - (-((int128)tmp_q[1] * 2986300999319L) - ((int128)tmp_q[2] * 173175875636344L) + ((int128)tmp_q[3] * 60263868712698L) - ((int128)tmp_q[4] * 85460932971483L) + ((int128)tmp_q[5] * 63250169712017L));
	tmp_zero[1] = ((int128)tmp_q[0] * 63250169712017L) + ((int128)tmp_q[1] * 87714942664860L) - (-((int128)tmp_q[2] * 2986300999319L) - ((int128)tmp_q[3] * 173175875636344L) + ((int128)tmp_q[4] * 60263868712698L) - ((int128)tmp_q[5] * 85460932971483L));
	tmp_zero[2] = -((int128)tmp_q[0] * 85460932971483L) + ((int128)tmp_q[1] * 63250169712017L) + ((int128)tmp_q[2] * 87714942664860L) - (-((int128)tmp_q[3] * 2986300999319L) - ((int128)tmp_q[4] * 173175875636344L) + ((int128)tmp_q[5] * 60263868712698L));
	tmp_zero[3] = ((int128)tmp_q[0] * 60263868712698L) - ((int128)tmp_q[1] * 85460932971483L) + ((int128)tmp_q[2] * 63250169712017L) + ((int128)tmp_q[3] * 87714942664860L) - (-((int128)tmp_q[4] * 2986300999319L) - ((int128)tmp_q[5] * 173175875636344L));
	tmp_zero[4] = -((int128)tmp_q[0] * 173175875636344L) + ((int128)tmp_q[1] * 60263868712698L) - ((int128)tmp_q[2] * 85460932971483L) + ((int128)tmp_q[3] * 63250169712017L) + ((int128)tmp_q[4] * 87714942664860L) + ((int128)tmp_q[5] * 2986300999319L);
	tmp_zero[5] = -((int128)tmp_q[0] * 2986300999319L) - ((int128)tmp_q[1] * 173175875636344L) + ((int128)tmp_q[2] * 60263868712698L) - ((int128)tmp_q[3] * 85460932971483L) + ((int128)tmp_q[4] * 63250169712017L) + ((int128)tmp_q[5] * 87714942664860L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

