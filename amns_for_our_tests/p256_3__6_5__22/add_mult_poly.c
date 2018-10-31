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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3491697155929328970UL) + ((((uint64_t)op[1] * 7648740816194539991UL) + ((uint64_t)op[2] * 16847658504861139384UL) + ((uint64_t)op[3] * 3409507938192107727UL) + ((uint64_t)op[4] * 9306499439710795869UL) + ((uint64_t)op[5] * 10392192235440260072UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 10392192235440260072UL) + ((uint64_t)op[1] * 3491697155929328970UL) + ((((uint64_t)op[2] * 7648740816194539991UL) + ((uint64_t)op[3] * 16847658504861139384UL) + ((uint64_t)op[4] * 3409507938192107727UL) + ((uint64_t)op[5] * 9306499439710795869UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 9306499439710795869UL) + ((uint64_t)op[1] * 10392192235440260072UL) + ((uint64_t)op[2] * 3491697155929328970UL) + ((((uint64_t)op[3] * 7648740816194539991UL) + ((uint64_t)op[4] * 16847658504861139384UL) + ((uint64_t)op[5] * 3409507938192107727UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 3409507938192107727UL) + ((uint64_t)op[1] * 9306499439710795869UL) + ((uint64_t)op[2] * 10392192235440260072UL) + ((uint64_t)op[3] * 3491697155929328970UL) + ((((uint64_t)op[4] * 7648740816194539991UL) + ((uint64_t)op[5] * 16847658504861139384UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 16847658504861139384UL) + ((uint64_t)op[1] * 3409507938192107727UL) + ((uint64_t)op[2] * 9306499439710795869UL) + ((uint64_t)op[3] * 10392192235440260072UL) + ((uint64_t)op[4] * 3491697155929328970UL) + ((uint64_t)op[5] * 1350215933553596723UL);
	tmp_q[5] = ((uint64_t)op[0] * 7648740816194539991UL) + ((uint64_t)op[1] * 16847658504861139384UL) + ((uint64_t)op[2] * 3409507938192107727UL) + ((uint64_t)op[3] * 9306499439710795869UL) + ((uint64_t)op[4] * 10392192235440260072UL) + ((uint64_t)op[5] * 3491697155929328970UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1337656018126L) + ((-((int128)tmp_q[1] * 4568784125893L) - ((int128)tmp_q[2] * 6108213202694L) + ((int128)tmp_q[3] * 1019392247137L) + ((int128)tmp_q[4] * 6184751287577L) - ((int128)tmp_q[5] * 6764616529262L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 6764616529262L) - ((int128)tmp_q[1] * 1337656018126L) + ((-((int128)tmp_q[2] * 4568784125893L) - ((int128)tmp_q[3] * 6108213202694L) + ((int128)tmp_q[4] * 1019392247137L) + ((int128)tmp_q[5] * 6184751287577L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 6184751287577L) - ((int128)tmp_q[1] * 6764616529262L) - ((int128)tmp_q[2] * 1337656018126L) + ((-((int128)tmp_q[3] * 4568784125893L) - ((int128)tmp_q[4] * 6108213202694L) + ((int128)tmp_q[5] * 1019392247137L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 1019392247137L) + ((int128)tmp_q[1] * 6184751287577L) - ((int128)tmp_q[2] * 6764616529262L) - ((int128)tmp_q[3] * 1337656018126L) + ((-((int128)tmp_q[4] * 4568784125893L) - ((int128)tmp_q[5] * 6108213202694L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 6108213202694L) + ((int128)tmp_q[1] * 1019392247137L) + ((int128)tmp_q[2] * 6184751287577L) - ((int128)tmp_q[3] * 6764616529262L) - ((int128)tmp_q[4] * 1337656018126L) - ((int128)tmp_q[5] * 22843920629465L);
	tmp_zero[5] = -((int128)tmp_q[0] * 4568784125893L) - ((int128)tmp_q[1] * 6108213202694L) + ((int128)tmp_q[2] * 1019392247137L) + ((int128)tmp_q[3] * 6184751287577L) - ((int128)tmp_q[4] * 6764616529262L) - ((int128)tmp_q[5] * 1337656018126L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

