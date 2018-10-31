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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2765297424945027721UL) + ((((uint64_t)op[1] * 3555067122660423139UL) + ((uint64_t)op[2] * 10109466307999210938UL) + ((uint64_t)op[3] * 2188661220420879889UL) + ((uint64_t)op[4] * 4401684544701201412UL) + ((uint64_t)op[5] * 1725218981391103106UL) + ((uint64_t)op[6] * 4164452271635163311UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 4164452271635163311UL) + ((uint64_t)op[1] * 2765297424945027721UL) + ((((uint64_t)op[2] * 3555067122660423139UL) + ((uint64_t)op[3] * 10109466307999210938UL) + ((uint64_t)op[4] * 2188661220420879889UL) + ((uint64_t)op[5] * 4401684544701201412UL) + ((uint64_t)op[6] * 1725218981391103106UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 1725218981391103106UL) + ((uint64_t)op[1] * 4164452271635163311UL) + ((uint64_t)op[2] * 2765297424945027721UL) + ((((uint64_t)op[3] * 3555067122660423139UL) + ((uint64_t)op[4] * 10109466307999210938UL) + ((uint64_t)op[5] * 2188661220420879889UL) + ((uint64_t)op[6] * 4401684544701201412UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 4401684544701201412UL) + ((uint64_t)op[1] * 1725218981391103106UL) + ((uint64_t)op[2] * 4164452271635163311UL) + ((uint64_t)op[3] * 2765297424945027721UL) + ((((uint64_t)op[4] * 3555067122660423139UL) + ((uint64_t)op[5] * 10109466307999210938UL) + ((uint64_t)op[6] * 2188661220420879889UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 2188661220420879889UL) + ((uint64_t)op[1] * 4401684544701201412UL) + ((uint64_t)op[2] * 1725218981391103106UL) + ((uint64_t)op[3] * 4164452271635163311UL) + ((uint64_t)op[4] * 2765297424945027721UL) + ((((uint64_t)op[5] * 3555067122660423139UL) + ((uint64_t)op[6] * 10109466307999210938UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 10109466307999210938UL) + ((uint64_t)op[1] * 2188661220420879889UL) + ((uint64_t)op[2] * 4401684544701201412UL) + ((uint64_t)op[3] * 1725218981391103106UL) + ((uint64_t)op[4] * 4164452271635163311UL) + ((uint64_t)op[5] * 2765297424945027721UL) + ((uint64_t)op[6] * 7110134245320846278UL);
	tmp_q[6] = ((uint64_t)op[0] * 3555067122660423139UL) + ((uint64_t)op[1] * 10109466307999210938UL) + ((uint64_t)op[2] * 2188661220420879889UL) + ((uint64_t)op[3] * 4401684544701201412UL) + ((uint64_t)op[4] * 1725218981391103106UL) + ((uint64_t)op[5] * 4164452271635163311UL) + ((uint64_t)op[6] * 2765297424945027721UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 9411102677882381L) + ((-((int128)tmp_q[1] * 14256015479250103L) - ((int128)tmp_q[2] * 13414328371894485L) - ((int128)tmp_q[3] * 5233419491280608L) + ((int128)tmp_q[4] * 11237856844201169L) + ((int128)tmp_q[5] * 3803499742111145L) + ((int128)tmp_q[6] * 14796589549258879L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 14796589549258879L) + ((int128)tmp_q[1] * 9411102677882381L) + ((-((int128)tmp_q[2] * 14256015479250103L) - ((int128)tmp_q[3] * 13414328371894485L) - ((int128)tmp_q[4] * 5233419491280608L) + ((int128)tmp_q[5] * 11237856844201169L) + ((int128)tmp_q[6] * 3803499742111145L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 3803499742111145L) + ((int128)tmp_q[1] * 14796589549258879L) + ((int128)tmp_q[2] * 9411102677882381L) + ((-((int128)tmp_q[3] * 14256015479250103L) - ((int128)tmp_q[4] * 13414328371894485L) - ((int128)tmp_q[5] * 5233419491280608L) + ((int128)tmp_q[6] * 11237856844201169L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 11237856844201169L) + ((int128)tmp_q[1] * 3803499742111145L) + ((int128)tmp_q[2] * 14796589549258879L) + ((int128)tmp_q[3] * 9411102677882381L) + ((-((int128)tmp_q[4] * 14256015479250103L) - ((int128)tmp_q[5] * 13414328371894485L) - ((int128)tmp_q[6] * 5233419491280608L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 5233419491280608L) + ((int128)tmp_q[1] * 11237856844201169L) + ((int128)tmp_q[2] * 3803499742111145L) + ((int128)tmp_q[3] * 14796589549258879L) + ((int128)tmp_q[4] * 9411102677882381L) + ((-((int128)tmp_q[5] * 14256015479250103L) - ((int128)tmp_q[6] * 13414328371894485L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 13414328371894485L) - ((int128)tmp_q[1] * 5233419491280608L) + ((int128)tmp_q[2] * 11237856844201169L) + ((int128)tmp_q[3] * 3803499742111145L) + ((int128)tmp_q[4] * 14796589549258879L) + ((int128)tmp_q[5] * 9411102677882381L) - ((int128)tmp_q[6] * 28512030958500206L);
	tmp_zero[6] = -((int128)tmp_q[0] * 14256015479250103L) - ((int128)tmp_q[1] * 13414328371894485L) - ((int128)tmp_q[2] * 5233419491280608L) + ((int128)tmp_q[3] * 11237856844201169L) + ((int128)tmp_q[4] * 3803499742111145L) + ((int128)tmp_q[5] * 14796589549258879L) + ((int128)tmp_q[6] * 9411102677882381L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

