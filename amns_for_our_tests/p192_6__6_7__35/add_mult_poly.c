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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5134979562656136537UL) + ((((uint64_t)op[1] * 10402870928188549072UL) + ((uint64_t)op[2] * 12114008333730045240UL) + ((uint64_t)op[3] * 14216536388602988983UL) + ((uint64_t)op[4] * 11124723767791504586UL) + ((uint64_t)op[5] * 9946415518944176931UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 9946415518944176931UL) + ((uint64_t)op[1] * 5134979562656136537UL) + ((((uint64_t)op[2] * 10402870928188549072UL) + ((uint64_t)op[3] * 12114008333730045240UL) + ((uint64_t)op[4] * 14216536388602988983UL) + ((uint64_t)op[5] * 11124723767791504586UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 11124723767791504586UL) + ((uint64_t)op[1] * 9946415518944176931UL) + ((uint64_t)op[2] * 5134979562656136537UL) + ((((uint64_t)op[3] * 10402870928188549072UL) + ((uint64_t)op[4] * 12114008333730045240UL) + ((uint64_t)op[5] * 14216536388602988983UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 14216536388602988983UL) + ((uint64_t)op[1] * 11124723767791504586UL) + ((uint64_t)op[2] * 9946415518944176931UL) + ((uint64_t)op[3] * 5134979562656136537UL) + ((((uint64_t)op[4] * 10402870928188549072UL) + ((uint64_t)op[5] * 12114008333730045240UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 12114008333730045240UL) + ((uint64_t)op[1] * 14216536388602988983UL) + ((uint64_t)op[2] * 11124723767791504586UL) + ((uint64_t)op[3] * 9946415518944176931UL) + ((uint64_t)op[4] * 5134979562656136537UL) + ((uint64_t)op[5] * 17479864276191188656UL);
	tmp_q[5] = ((uint64_t)op[0] * 10402870928188549072UL) + ((uint64_t)op[1] * 12114008333730045240UL) + ((uint64_t)op[2] * 14216536388602988983UL) + ((uint64_t)op[3] * 11124723767791504586UL) + ((uint64_t)op[4] * 9946415518944176931UL) + ((uint64_t)op[5] * 5134979562656136537UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1701588160L) + ((-((int128)tmp_q[1] * 2512889063L) - ((int128)tmp_q[2] * 260790129L) - ((int128)tmp_q[3] * 676853320L) - ((int128)tmp_q[4] * 645161498L) - ((int128)tmp_q[5] * 671937229L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 671937229L) - ((int128)tmp_q[1] * 1701588160L) + ((-((int128)tmp_q[2] * 2512889063L) - ((int128)tmp_q[3] * 260790129L) - ((int128)tmp_q[4] * 676853320L) - ((int128)tmp_q[5] * 645161498L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 645161498L) - ((int128)tmp_q[1] * 671937229L) - ((int128)tmp_q[2] * 1701588160L) + ((-((int128)tmp_q[3] * 2512889063L) - ((int128)tmp_q[4] * 260790129L) - ((int128)tmp_q[5] * 676853320L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 676853320L) - ((int128)tmp_q[1] * 645161498L) - ((int128)tmp_q[2] * 671937229L) - ((int128)tmp_q[3] * 1701588160L) + ((-((int128)tmp_q[4] * 2512889063L) - ((int128)tmp_q[5] * 260790129L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 260790129L) - ((int128)tmp_q[1] * 676853320L) - ((int128)tmp_q[2] * 645161498L) - ((int128)tmp_q[3] * 671937229L) - ((int128)tmp_q[4] * 1701588160L) - ((int128)tmp_q[5] * 17590223441L);
	tmp_zero[5] = -((int128)tmp_q[0] * 2512889063L) - ((int128)tmp_q[1] * 260790129L) - ((int128)tmp_q[2] * 676853320L) - ((int128)tmp_q[3] * 645161498L) - ((int128)tmp_q[4] * 671937229L) - ((int128)tmp_q[5] * 1701588160L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

