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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15028937207716359414UL) + ((((uint64_t)op[1] * 14042614774181344254UL) + ((uint64_t)op[2] * 16618728702623826186UL) + ((uint64_t)op[3] * 6235889632638747850UL) + ((uint64_t)op[4] * 18003848408318722713UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 18003848408318722713UL) + ((uint64_t)op[1] * 15028937207716359414UL) + ((((uint64_t)op[2] * 14042614774181344254UL) + ((uint64_t)op[3] * 16618728702623826186UL) + ((uint64_t)op[4] * 6235889632638747850UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 6235889632638747850UL) + ((uint64_t)op[1] * 18003848408318722713UL) + ((uint64_t)op[2] * 15028937207716359414UL) + ((((uint64_t)op[3] * 14042614774181344254UL) + ((uint64_t)op[4] * 16618728702623826186UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 16618728702623826186UL) + ((uint64_t)op[1] * 6235889632638747850UL) + ((uint64_t)op[2] * 18003848408318722713UL) + ((uint64_t)op[3] * 15028937207716359414UL) + ((uint64_t)op[4] * 6064583050721651698UL);
	tmp_q[4] = ((uint64_t)op[0] * 14042614774181344254UL) + ((uint64_t)op[1] * 16618728702623826186UL) + ((uint64_t)op[2] * 6235889632638747850UL) + ((uint64_t)op[3] * 18003848408318722713UL) + ((uint64_t)op[4] * 15028937207716359414UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 232884015234L) + ((-((int128)tmp_q[1] * 38285647423L) - ((int128)tmp_q[2] * 122947458106L) - ((int128)tmp_q[3] * 175997399054L) - ((int128)tmp_q[4] * 43858044026L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 43858044026L) - ((int128)tmp_q[1] * 232884015234L) + ((-((int128)tmp_q[2] * 38285647423L) - ((int128)tmp_q[3] * 122947458106L) - ((int128)tmp_q[4] * 175997399054L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 175997399054L) - ((int128)tmp_q[1] * 43858044026L) - ((int128)tmp_q[2] * 232884015234L) + ((-((int128)tmp_q[3] * 38285647423L) - ((int128)tmp_q[4] * 122947458106L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 122947458106L) - ((int128)tmp_q[1] * 175997399054L) - ((int128)tmp_q[2] * 43858044026L) - ((int128)tmp_q[3] * 232884015234L) - ((int128)tmp_q[4] * 267999531961L);
	tmp_zero[4] = -((int128)tmp_q[0] * 38285647423L) - ((int128)tmp_q[1] * 122947458106L) - ((int128)tmp_q[2] * 175997399054L) - ((int128)tmp_q[3] * 43858044026L) - ((int128)tmp_q[4] * 232884015234L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

