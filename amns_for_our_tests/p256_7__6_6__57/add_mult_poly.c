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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1192638223139869571UL) + ((((uint64_t)op[1] * 15885369977370716255UL) + ((uint64_t)op[2] * 398427257845542083UL) + ((uint64_t)op[3] * 12507817008275041893UL) + ((uint64_t)op[4] * 10108919349319166829UL) + ((uint64_t)op[5] * 13058271593386889916UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 13058271593386889916UL) + ((uint64_t)op[1] * 1192638223139869571UL) + ((((uint64_t)op[2] * 15885369977370716255UL) + ((uint64_t)op[3] * 398427257845542083UL) + ((uint64_t)op[4] * 12507817008275041893UL) + ((uint64_t)op[5] * 10108919349319166829UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 10108919349319166829UL) + ((uint64_t)op[1] * 13058271593386889916UL) + ((uint64_t)op[2] * 1192638223139869571UL) + ((((uint64_t)op[3] * 15885369977370716255UL) + ((uint64_t)op[4] * 398427257845542083UL) + ((uint64_t)op[5] * 12507817008275041893UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 12507817008275041893UL) + ((uint64_t)op[1] * 10108919349319166829UL) + ((uint64_t)op[2] * 13058271593386889916UL) + ((uint64_t)op[3] * 1192638223139869571UL) + ((((uint64_t)op[4] * 15885369977370716255UL) + ((uint64_t)op[5] * 398427257845542083UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 398427257845542083UL) + ((uint64_t)op[1] * 12507817008275041893UL) + ((uint64_t)op[2] * 10108919349319166829UL) + ((uint64_t)op[3] * 13058271593386889916UL) + ((uint64_t)op[4] * 1192638223139869571UL) + ((uint64_t)op[5] * 3078499495676539450UL);
	tmp_q[5] = ((uint64_t)op[0] * 15885369977370716255UL) + ((uint64_t)op[1] * 398427257845542083UL) + ((uint64_t)op[2] * 12507817008275041893UL) + ((uint64_t)op[3] * 10108919349319166829UL) + ((uint64_t)op[4] * 13058271593386889916UL) + ((uint64_t)op[5] * 1192638223139869571UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3658868272607L) + ((((int128)tmp_q[1] * 984749301979L) + ((int128)tmp_q[2] * 4240759455920L) + ((int128)tmp_q[3] * 79968219149L) - ((int128)tmp_q[4] * 276873815051L) - ((int128)tmp_q[5] * 1385028827670L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 1385028827670L) - ((int128)tmp_q[1] * 3658868272607L) + ((((int128)tmp_q[2] * 984749301979L) + ((int128)tmp_q[3] * 4240759455920L) + ((int128)tmp_q[4] * 79968219149L) - ((int128)tmp_q[5] * 276873815051L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 276873815051L) - ((int128)tmp_q[1] * 1385028827670L) - ((int128)tmp_q[2] * 3658868272607L) + ((((int128)tmp_q[3] * 984749301979L) + ((int128)tmp_q[4] * 4240759455920L) + ((int128)tmp_q[5] * 79968219149L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 79968219149L) - ((int128)tmp_q[1] * 276873815051L) - ((int128)tmp_q[2] * 1385028827670L) - ((int128)tmp_q[3] * 3658868272607L) + ((((int128)tmp_q[4] * 984749301979L) + ((int128)tmp_q[5] * 4240759455920L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 4240759455920L) + ((int128)tmp_q[1] * 79968219149L) - ((int128)tmp_q[2] * 276873815051L) - ((int128)tmp_q[3] * 1385028827670L) - ((int128)tmp_q[4] * 3658868272607L) + ((int128)tmp_q[5] * 5908495811874L);
	tmp_zero[5] = ((int128)tmp_q[0] * 984749301979L) + ((int128)tmp_q[1] * 4240759455920L) + ((int128)tmp_q[2] * 79968219149L) - ((int128)tmp_q[3] * 276873815051L) - ((int128)tmp_q[4] * 1385028827670L) - ((int128)tmp_q[5] * 3658868272607L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

