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
	tmp_q[0] = ((uint64_t)op[0] * 15017881964218166852UL) + ((((uint64_t)op[1] * 17015217530661285082UL) + ((uint64_t)op[2] * 13852594961215715751UL) + ((uint64_t)op[3] * 5254129526572140651UL) + ((uint64_t)op[4] * 9839873785181111329UL) + ((uint64_t)op[5] * 463160036307106100UL) + ((uint64_t)op[6] * 6153551921839161738UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 6153551921839161738UL) + ((uint64_t)op[1] * 15017881964218166852UL) + ((((uint64_t)op[2] * 17015217530661285082UL) + ((uint64_t)op[3] * 13852594961215715751UL) + ((uint64_t)op[4] * 5254129526572140651UL) + ((uint64_t)op[5] * 9839873785181111329UL) + ((uint64_t)op[6] * 463160036307106100UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 463160036307106100UL) + ((uint64_t)op[1] * 6153551921839161738UL) + ((uint64_t)op[2] * 15017881964218166852UL) + ((((uint64_t)op[3] * 17015217530661285082UL) + ((uint64_t)op[4] * 13852594961215715751UL) + ((uint64_t)op[5] * 5254129526572140651UL) + ((uint64_t)op[6] * 9839873785181111329UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 9839873785181111329UL) + ((uint64_t)op[1] * 463160036307106100UL) + ((uint64_t)op[2] * 6153551921839161738UL) + ((uint64_t)op[3] * 15017881964218166852UL) + ((((uint64_t)op[4] * 17015217530661285082UL) + ((uint64_t)op[5] * 13852594961215715751UL) + ((uint64_t)op[6] * 5254129526572140651UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 5254129526572140651UL) + ((uint64_t)op[1] * 9839873785181111329UL) + ((uint64_t)op[2] * 463160036307106100UL) + ((uint64_t)op[3] * 6153551921839161738UL) + ((uint64_t)op[4] * 15017881964218166852UL) + ((((uint64_t)op[5] * 17015217530661285082UL) + ((uint64_t)op[6] * 13852594961215715751UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 13852594961215715751UL) + ((uint64_t)op[1] * 5254129526572140651UL) + ((uint64_t)op[2] * 9839873785181111329UL) + ((uint64_t)op[3] * 463160036307106100UL) + ((uint64_t)op[4] * 6153551921839161738UL) + ((uint64_t)op[5] * 15017881964218166852UL) + ((uint64_t)op[6] * 11289111358468218946UL);
	tmp_q[6] = ((uint64_t)op[0] * 17015217530661285082UL) + ((uint64_t)op[1] * 13852594961215715751UL) + ((uint64_t)op[2] * 5254129526572140651UL) + ((uint64_t)op[3] * 9839873785181111329UL) + ((uint64_t)op[4] * 463160036307106100UL) + ((uint64_t)op[5] * 6153551921839161738UL) + ((uint64_t)op[6] * 15017881964218166852UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 28416015933L) + ((((int128)tmp_q[1] * 20669594299L) + ((int128)tmp_q[2] * 12002085738L) - ((int128)tmp_q[3] * 24704833645L) + ((int128)tmp_q[4] * 9766685967L) - ((int128)tmp_q[5] * 11763023421L) - ((int128)tmp_q[6] * 7135500354L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 7135500354L) + ((int128)tmp_q[1] * 28416015933L) + ((((int128)tmp_q[2] * 20669594299L) + ((int128)tmp_q[3] * 12002085738L) - ((int128)tmp_q[4] * 24704833645L) + ((int128)tmp_q[5] * 9766685967L) - ((int128)tmp_q[6] * 11763023421L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 11763023421L) - ((int128)tmp_q[1] * 7135500354L) + ((int128)tmp_q[2] * 28416015933L) + ((((int128)tmp_q[3] * 20669594299L) + ((int128)tmp_q[4] * 12002085738L) - ((int128)tmp_q[5] * 24704833645L) + ((int128)tmp_q[6] * 9766685967L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 9766685967L) - ((int128)tmp_q[1] * 11763023421L) - ((int128)tmp_q[2] * 7135500354L) + ((int128)tmp_q[3] * 28416015933L) + ((((int128)tmp_q[4] * 20669594299L) + ((int128)tmp_q[5] * 12002085738L) - ((int128)tmp_q[6] * 24704833645L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 24704833645L) + ((int128)tmp_q[1] * 9766685967L) - ((int128)tmp_q[2] * 11763023421L) - ((int128)tmp_q[3] * 7135500354L) + ((int128)tmp_q[4] * 28416015933L) + ((((int128)tmp_q[5] * 20669594299L) + ((int128)tmp_q[6] * 12002085738L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 12002085738L) - ((int128)tmp_q[1] * 24704833645L) + ((int128)tmp_q[2] * 9766685967L) - ((int128)tmp_q[3] * 11763023421L) - ((int128)tmp_q[4] * 7135500354L) + ((int128)tmp_q[5] * 28416015933L) + ((int128)tmp_q[6] * 103347971495L);
	tmp_zero[6] = ((int128)tmp_q[0] * 20669594299L) + ((int128)tmp_q[1] * 12002085738L) - ((int128)tmp_q[2] * 24704833645L) + ((int128)tmp_q[3] * 9766685967L) - ((int128)tmp_q[4] * 11763023421L) - ((int128)tmp_q[5] * 7135500354L) + ((int128)tmp_q[6] * 28416015933L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

