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
	tmp_q[0] = ((uint64_t)op[0] * 9670776295704020863UL) + ((((uint64_t)op[1] * 13233341183492875307UL) + ((uint64_t)op[2] * 856414420260027459UL) + ((uint64_t)op[3] * 15741876947218290348UL) + ((uint64_t)op[4] * 18102560560773262278UL) + ((uint64_t)op[5] * 9807943168550545800UL) + ((uint64_t)op[6] * 13236555398348592016UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 13236555398348592016UL) + ((uint64_t)op[1] * 9670776295704020863UL) + ((((uint64_t)op[2] * 13233341183492875307UL) + ((uint64_t)op[3] * 856414420260027459UL) + ((uint64_t)op[4] * 15741876947218290348UL) + ((uint64_t)op[5] * 18102560560773262278UL) + ((uint64_t)op[6] * 9807943168550545800UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 9807943168550545800UL) + ((uint64_t)op[1] * 13236555398348592016UL) + ((uint64_t)op[2] * 9670776295704020863UL) + ((((uint64_t)op[3] * 13233341183492875307UL) + ((uint64_t)op[4] * 856414420260027459UL) + ((uint64_t)op[5] * 15741876947218290348UL) + ((uint64_t)op[6] * 18102560560773262278UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 18102560560773262278UL) + ((uint64_t)op[1] * 9807943168550545800UL) + ((uint64_t)op[2] * 13236555398348592016UL) + ((uint64_t)op[3] * 9670776295704020863UL) + ((((uint64_t)op[4] * 13233341183492875307UL) + ((uint64_t)op[5] * 856414420260027459UL) + ((uint64_t)op[6] * 15741876947218290348UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 15741876947218290348UL) + ((uint64_t)op[1] * 18102560560773262278UL) + ((uint64_t)op[2] * 9807943168550545800UL) + ((uint64_t)op[3] * 13236555398348592016UL) + ((uint64_t)op[4] * 9670776295704020863UL) + ((((uint64_t)op[5] * 13233341183492875307UL) + ((uint64_t)op[6] * 856414420260027459UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 856414420260027459UL) + ((uint64_t)op[1] * 15741876947218290348UL) + ((uint64_t)op[2] * 18102560560773262278UL) + ((uint64_t)op[3] * 9807943168550545800UL) + ((uint64_t)op[4] * 13236555398348592016UL) + ((uint64_t)op[5] * 9670776295704020863UL) + ((uint64_t)op[6] * 10826473696335721687UL);
	tmp_q[6] = ((uint64_t)op[0] * 13233341183492875307UL) + ((uint64_t)op[1] * 856414420260027459UL) + ((uint64_t)op[2] * 15741876947218290348UL) + ((uint64_t)op[3] * 18102560560773262278UL) + ((uint64_t)op[4] * 9807943168550545800UL) + ((uint64_t)op[5] * 13236555398348592016UL) + ((uint64_t)op[6] * 9670776295704020863UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 17555072741L) + ((-((int128)tmp_q[1] * 9782433752L) + ((int128)tmp_q[2] * 47715548733L) + ((int128)tmp_q[3] * 46823577201L) + ((int128)tmp_q[4] * 11480870304L) + ((int128)tmp_q[5] * 60193414899L) + ((int128)tmp_q[6] * 43047858567L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 43047858567L) + ((int128)tmp_q[1] * 17555072741L) + ((-((int128)tmp_q[2] * 9782433752L) + ((int128)tmp_q[3] * 47715548733L) + ((int128)tmp_q[4] * 46823577201L) + ((int128)tmp_q[5] * 11480870304L) + ((int128)tmp_q[6] * 60193414899L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 60193414899L) + ((int128)tmp_q[1] * 43047858567L) + ((int128)tmp_q[2] * 17555072741L) + ((-((int128)tmp_q[3] * 9782433752L) + ((int128)tmp_q[4] * 47715548733L) + ((int128)tmp_q[5] * 46823577201L) + ((int128)tmp_q[6] * 11480870304L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 11480870304L) + ((int128)tmp_q[1] * 60193414899L) + ((int128)tmp_q[2] * 43047858567L) + ((int128)tmp_q[3] * 17555072741L) + ((-((int128)tmp_q[4] * 9782433752L) + ((int128)tmp_q[5] * 47715548733L) + ((int128)tmp_q[6] * 46823577201L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 46823577201L) + ((int128)tmp_q[1] * 11480870304L) + ((int128)tmp_q[2] * 60193414899L) + ((int128)tmp_q[3] * 43047858567L) + ((int128)tmp_q[4] * 17555072741L) + ((-((int128)tmp_q[5] * 9782433752L) + ((int128)tmp_q[6] * 47715548733L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 47715548733L) + ((int128)tmp_q[1] * 46823577201L) + ((int128)tmp_q[2] * 11480870304L) + ((int128)tmp_q[3] * 60193414899L) + ((int128)tmp_q[4] * 43047858567L) + ((int128)tmp_q[5] * 17555072741L) - ((int128)tmp_q[6] * 48912168760L);
	tmp_zero[6] = -((int128)tmp_q[0] * 9782433752L) + ((int128)tmp_q[1] * 47715548733L) + ((int128)tmp_q[2] * 46823577201L) + ((int128)tmp_q[3] * 11480870304L) + ((int128)tmp_q[4] * 60193414899L) + ((int128)tmp_q[5] * 43047858567L) + ((int128)tmp_q[6] * 17555072741L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

