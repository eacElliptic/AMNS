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
	tmp_q[0] = ((uint64_t)op[0] * 4064695502864551061UL) + ((((uint64_t)op[1] * 16009587823487782172UL) + ((uint64_t)op[2] * 3411174868994677060UL) + ((uint64_t)op[3] * 17423200193717141735UL) + ((uint64_t)op[4] * 4341124398432788288UL) + ((uint64_t)op[5] * 11926766480026338668UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 11926766480026338668UL) + ((uint64_t)op[1] * 4064695502864551061UL) + ((((uint64_t)op[2] * 16009587823487782172UL) + ((uint64_t)op[3] * 3411174868994677060UL) + ((uint64_t)op[4] * 17423200193717141735UL) + ((uint64_t)op[5] * 4341124398432788288UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 4341124398432788288UL) + ((uint64_t)op[1] * 11926766480026338668UL) + ((uint64_t)op[2] * 4064695502864551061UL) + ((((uint64_t)op[3] * 16009587823487782172UL) + ((uint64_t)op[4] * 3411174868994677060UL) + ((uint64_t)op[5] * 17423200193717141735UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 17423200193717141735UL) + ((uint64_t)op[1] * 4341124398432788288UL) + ((uint64_t)op[2] * 11926766480026338668UL) + ((uint64_t)op[3] * 4064695502864551061UL) + ((((uint64_t)op[4] * 16009587823487782172UL) + ((uint64_t)op[5] * 3411174868994677060UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 3411174868994677060UL) + ((uint64_t)op[1] * 17423200193717141735UL) + ((uint64_t)op[2] * 4341124398432788288UL) + ((uint64_t)op[3] * 11926766480026338668UL) + ((uint64_t)op[4] * 4064695502864551061UL) + ((uint64_t)op[5] * 3823806572378934952UL);
	tmp_q[5] = ((uint64_t)op[0] * 16009587823487782172UL) + ((uint64_t)op[1] * 3411174868994677060UL) + ((uint64_t)op[2] * 17423200193717141735UL) + ((uint64_t)op[3] * 4341124398432788288UL) + ((uint64_t)op[4] * 11926766480026338668UL) + ((uint64_t)op[5] * 4064695502864551061UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 82149537143L) + ((-((int128)tmp_q[1] * 74721796460L) + ((int128)tmp_q[2] * 58861709108L) - ((int128)tmp_q[3] * 1065197539L) + ((int128)tmp_q[4] * 42064133856L) + ((int128)tmp_q[5] * 16977414564L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 16977414564L) - ((int128)tmp_q[1] * 82149537143L) + ((-((int128)tmp_q[2] * 74721796460L) + ((int128)tmp_q[3] * 58861709108L) - ((int128)tmp_q[4] * 1065197539L) + ((int128)tmp_q[5] * 42064133856L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 42064133856L) + ((int128)tmp_q[1] * 16977414564L) - ((int128)tmp_q[2] * 82149537143L) + ((-((int128)tmp_q[3] * 74721796460L) + ((int128)tmp_q[4] * 58861709108L) - ((int128)tmp_q[5] * 1065197539L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 1065197539L) + ((int128)tmp_q[1] * 42064133856L) + ((int128)tmp_q[2] * 16977414564L) - ((int128)tmp_q[3] * 82149537143L) + ((-((int128)tmp_q[4] * 74721796460L) + ((int128)tmp_q[5] * 58861709108L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 58861709108L) - ((int128)tmp_q[1] * 1065197539L) + ((int128)tmp_q[2] * 42064133856L) + ((int128)tmp_q[3] * 16977414564L) - ((int128)tmp_q[4] * 82149537143L) - ((int128)tmp_q[5] * 448330778760L);
	tmp_zero[5] = -((int128)tmp_q[0] * 74721796460L) + ((int128)tmp_q[1] * 58861709108L) - ((int128)tmp_q[2] * 1065197539L) + ((int128)tmp_q[3] * 42064133856L) + ((int128)tmp_q[4] * 16977414564L) - ((int128)tmp_q[5] * 82149537143L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

