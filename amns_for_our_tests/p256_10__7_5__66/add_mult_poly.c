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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3729025387245575658UL) + ((((uint64_t)op[1] * 1766858154768292620UL) + ((uint64_t)op[2] * 10842432913510289513UL) + ((uint64_t)op[3] * 2312767378934400501UL) + ((uint64_t)op[4] * 2576345106490322972UL) + ((uint64_t)op[5] * 7293484046771382890UL) + ((uint64_t)op[6] * 7590950340415729187UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 7590950340415729187UL) + ((uint64_t)op[1] * 3729025387245575658UL) + ((((uint64_t)op[2] * 1766858154768292620UL) + ((uint64_t)op[3] * 10842432913510289513UL) + ((uint64_t)op[4] * 2312767378934400501UL) + ((uint64_t)op[5] * 2576345106490322972UL) + ((uint64_t)op[6] * 7293484046771382890UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 7293484046771382890UL) + ((uint64_t)op[1] * 7590950340415729187UL) + ((uint64_t)op[2] * 3729025387245575658UL) + ((((uint64_t)op[3] * 1766858154768292620UL) + ((uint64_t)op[4] * 10842432913510289513UL) + ((uint64_t)op[5] * 2312767378934400501UL) + ((uint64_t)op[6] * 2576345106490322972UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 2576345106490322972UL) + ((uint64_t)op[1] * 7293484046771382890UL) + ((uint64_t)op[2] * 7590950340415729187UL) + ((uint64_t)op[3] * 3729025387245575658UL) + ((((uint64_t)op[4] * 1766858154768292620UL) + ((uint64_t)op[5] * 10842432913510289513UL) + ((uint64_t)op[6] * 2312767378934400501UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 2312767378934400501UL) + ((uint64_t)op[1] * 2576345106490322972UL) + ((uint64_t)op[2] * 7293484046771382890UL) + ((uint64_t)op[3] * 7590950340415729187UL) + ((uint64_t)op[4] * 3729025387245575658UL) + ((((uint64_t)op[5] * 1766858154768292620UL) + ((uint64_t)op[6] * 10842432913510289513UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 10842432913510289513UL) + ((uint64_t)op[1] * 2312767378934400501UL) + ((uint64_t)op[2] * 2576345106490322972UL) + ((uint64_t)op[3] * 7293484046771382890UL) + ((uint64_t)op[4] * 7590950340415729187UL) + ((uint64_t)op[5] * 3729025387245575658UL) + ((uint64_t)op[6] * 8834290773841463100UL);
	tmp_q[6] = ((uint64_t)op[0] * 1766858154768292620UL) + ((uint64_t)op[1] * 10842432913510289513UL) + ((uint64_t)op[2] * 2312767378934400501UL) + ((uint64_t)op[3] * 2576345106490322972UL) + ((uint64_t)op[4] * 7293484046771382890UL) + ((uint64_t)op[5] * 7590950340415729187UL) + ((uint64_t)op[6] * 3729025387245575658UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3255954898L) + ((((int128)tmp_q[1] * 22050127123L) - ((int128)tmp_q[2] * 53853915402L) + ((int128)tmp_q[3] * 45543026389L) - ((int128)tmp_q[4] * 12724374103L) + ((int128)tmp_q[5] * 56410180325L) + ((int128)tmp_q[6] * 53607134937L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 53607134937L) - ((int128)tmp_q[1] * 3255954898L) + ((((int128)tmp_q[2] * 22050127123L) - ((int128)tmp_q[3] * 53853915402L) + ((int128)tmp_q[4] * 45543026389L) - ((int128)tmp_q[5] * 12724374103L) + ((int128)tmp_q[6] * 56410180325L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 56410180325L) + ((int128)tmp_q[1] * 53607134937L) - ((int128)tmp_q[2] * 3255954898L) + ((((int128)tmp_q[3] * 22050127123L) - ((int128)tmp_q[4] * 53853915402L) + ((int128)tmp_q[5] * 45543026389L) - ((int128)tmp_q[6] * 12724374103L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 12724374103L) + ((int128)tmp_q[1] * 56410180325L) + ((int128)tmp_q[2] * 53607134937L) - ((int128)tmp_q[3] * 3255954898L) + ((((int128)tmp_q[4] * 22050127123L) - ((int128)tmp_q[5] * 53853915402L) + ((int128)tmp_q[6] * 45543026389L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 45543026389L) - ((int128)tmp_q[1] * 12724374103L) + ((int128)tmp_q[2] * 56410180325L) + ((int128)tmp_q[3] * 53607134937L) - ((int128)tmp_q[4] * 3255954898L) + ((((int128)tmp_q[5] * 22050127123L) - ((int128)tmp_q[6] * 53853915402L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 53853915402L) + ((int128)tmp_q[1] * 45543026389L) - ((int128)tmp_q[2] * 12724374103L) + ((int128)tmp_q[3] * 56410180325L) + ((int128)tmp_q[4] * 53607134937L) - ((int128)tmp_q[5] * 3255954898L) + ((int128)tmp_q[6] * 110250635615L);
	tmp_zero[6] = ((int128)tmp_q[0] * 22050127123L) - ((int128)tmp_q[1] * 53853915402L) + ((int128)tmp_q[2] * 45543026389L) - ((int128)tmp_q[3] * 12724374103L) + ((int128)tmp_q[4] * 56410180325L) + ((int128)tmp_q[5] * 53607134937L) - ((int128)tmp_q[6] * 3255954898L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

