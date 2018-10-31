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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) << 2);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[8]) << 2);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[8] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[8] * pa[8]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8292143105809843961UL) + ((((uint64_t)op[1] * 14143821510296452471UL) + ((uint64_t)op[2] * 16616022375094639376UL) + ((uint64_t)op[3] * 1925278081890378847UL) + ((uint64_t)op[4] * 13187538776768387256UL) + ((uint64_t)op[5] * 15396728060797277524UL) + ((uint64_t)op[6] * 4642456296233477485UL) + ((uint64_t)op[7] * 14600804472141172470UL) + ((uint64_t)op[8] * 10450597607955789708UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 10450597607955789708UL) + ((uint64_t)op[1] * 8292143105809843961UL) + ((((uint64_t)op[2] * 14143821510296452471UL) + ((uint64_t)op[3] * 16616022375094639376UL) + ((uint64_t)op[4] * 1925278081890378847UL) + ((uint64_t)op[5] * 13187538776768387256UL) + ((uint64_t)op[6] * 15396728060797277524UL) + ((uint64_t)op[7] * 4642456296233477485UL) + ((uint64_t)op[8] * 14600804472141172470UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 14600804472141172470UL) + ((uint64_t)op[1] * 10450597607955789708UL) + ((uint64_t)op[2] * 8292143105809843961UL) + ((((uint64_t)op[3] * 14143821510296452471UL) + ((uint64_t)op[4] * 16616022375094639376UL) + ((uint64_t)op[5] * 1925278081890378847UL) + ((uint64_t)op[6] * 13187538776768387256UL) + ((uint64_t)op[7] * 15396728060797277524UL) + ((uint64_t)op[8] * 4642456296233477485UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 4642456296233477485UL) + ((uint64_t)op[1] * 14600804472141172470UL) + ((uint64_t)op[2] * 10450597607955789708UL) + ((uint64_t)op[3] * 8292143105809843961UL) + ((((uint64_t)op[4] * 14143821510296452471UL) + ((uint64_t)op[5] * 16616022375094639376UL) + ((uint64_t)op[6] * 1925278081890378847UL) + ((uint64_t)op[7] * 13187538776768387256UL) + ((uint64_t)op[8] * 15396728060797277524UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 15396728060797277524UL) + ((uint64_t)op[1] * 4642456296233477485UL) + ((uint64_t)op[2] * 14600804472141172470UL) + ((uint64_t)op[3] * 10450597607955789708UL) + ((uint64_t)op[4] * 8292143105809843961UL) + ((((uint64_t)op[5] * 14143821510296452471UL) + ((uint64_t)op[6] * 16616022375094639376UL) + ((uint64_t)op[7] * 1925278081890378847UL) + ((uint64_t)op[8] * 13187538776768387256UL)) * 18446744073709551612);
	tmp_q[5] = ((uint64_t)op[0] * 13187538776768387256UL) + ((uint64_t)op[1] * 15396728060797277524UL) + ((uint64_t)op[2] * 4642456296233477485UL) + ((uint64_t)op[3] * 14600804472141172470UL) + ((uint64_t)op[4] * 10450597607955789708UL) + ((uint64_t)op[5] * 8292143105809843961UL) + ((((uint64_t)op[6] * 14143821510296452471UL) + ((uint64_t)op[7] * 16616022375094639376UL) + ((uint64_t)op[8] * 1925278081890378847UL)) * 18446744073709551612);
	tmp_q[6] = ((uint64_t)op[0] * 1925278081890378847UL) + ((uint64_t)op[1] * 13187538776768387256UL) + ((uint64_t)op[2] * 15396728060797277524UL) + ((uint64_t)op[3] * 4642456296233477485UL) + ((uint64_t)op[4] * 14600804472141172470UL) + ((uint64_t)op[5] * 10450597607955789708UL) + ((uint64_t)op[6] * 8292143105809843961UL) + ((((uint64_t)op[7] * 14143821510296452471UL) + ((uint64_t)op[8] * 16616022375094639376UL)) * 18446744073709551612);
	tmp_q[7] = ((uint64_t)op[0] * 16616022375094639376UL) + ((uint64_t)op[1] * 1925278081890378847UL) + ((uint64_t)op[2] * 13187538776768387256UL) + ((uint64_t)op[3] * 15396728060797277524UL) + ((uint64_t)op[4] * 4642456296233477485UL) + ((uint64_t)op[5] * 14600804472141172470UL) + ((uint64_t)op[6] * 10450597607955789708UL) + ((uint64_t)op[7] * 8292143105809843961UL) + ((uint64_t)op[8] * 17211690253652396580UL);
	tmp_q[8] = ((uint64_t)op[0] * 14143821510296452471UL) + ((uint64_t)op[1] * 16616022375094639376UL) + ((uint64_t)op[2] * 1925278081890378847UL) + ((uint64_t)op[3] * 13187538776768387256UL) + ((uint64_t)op[4] * 15396728060797277524UL) + ((uint64_t)op[5] * 4642456296233477485UL) + ((uint64_t)op[6] * 14600804472141172470UL) + ((uint64_t)op[7] * 10450597607955789708UL) + ((uint64_t)op[8] * 8292143105809843961UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 380872577557L) - ((((int128)tmp_q[1] * 584509749053L) - ((int128)tmp_q[2] * 1046833896412L) - ((int128)tmp_q[3] * 4218306354666L) - ((int128)tmp_q[4] * 1889351483184L) - ((int128)tmp_q[5] * 1124453604936L) + ((int128)tmp_q[6] * 1509488305121L) - ((int128)tmp_q[7] * 2057371408626L) + ((int128)tmp_q[8] * 1285671898092L)) * 4);
	tmp_zero[1] = ((int128)tmp_q[0] * 1285671898092L) - ((int128)tmp_q[1] * 380872577557L) - ((((int128)tmp_q[2] * 584509749053L) - ((int128)tmp_q[3] * 1046833896412L) - ((int128)tmp_q[4] * 4218306354666L) - ((int128)tmp_q[5] * 1889351483184L) - ((int128)tmp_q[6] * 1124453604936L) + ((int128)tmp_q[7] * 1509488305121L) - ((int128)tmp_q[8] * 2057371408626L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 2057371408626L) + ((int128)tmp_q[1] * 1285671898092L) - ((int128)tmp_q[2] * 380872577557L) - ((((int128)tmp_q[3] * 584509749053L) - ((int128)tmp_q[4] * 1046833896412L) - ((int128)tmp_q[5] * 4218306354666L) - ((int128)tmp_q[6] * 1889351483184L) - ((int128)tmp_q[7] * 1124453604936L) + ((int128)tmp_q[8] * 1509488305121L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 1509488305121L) - ((int128)tmp_q[1] * 2057371408626L) + ((int128)tmp_q[2] * 1285671898092L) - ((int128)tmp_q[3] * 380872577557L) - ((((int128)tmp_q[4] * 584509749053L) - ((int128)tmp_q[5] * 1046833896412L) - ((int128)tmp_q[6] * 4218306354666L) - ((int128)tmp_q[7] * 1889351483184L) - ((int128)tmp_q[8] * 1124453604936L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 1124453604936L) + ((int128)tmp_q[1] * 1509488305121L) - ((int128)tmp_q[2] * 2057371408626L) + ((int128)tmp_q[3] * 1285671898092L) - ((int128)tmp_q[4] * 380872577557L) - ((((int128)tmp_q[5] * 584509749053L) - ((int128)tmp_q[6] * 1046833896412L) - ((int128)tmp_q[7] * 4218306354666L) - ((int128)tmp_q[8] * 1889351483184L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 1889351483184L) - ((int128)tmp_q[1] * 1124453604936L) + ((int128)tmp_q[2] * 1509488305121L) - ((int128)tmp_q[3] * 2057371408626L) + ((int128)tmp_q[4] * 1285671898092L) - ((int128)tmp_q[5] * 380872577557L) - ((((int128)tmp_q[6] * 584509749053L) - ((int128)tmp_q[7] * 1046833896412L) - ((int128)tmp_q[8] * 4218306354666L)) * 4);
	tmp_zero[6] = -((int128)tmp_q[0] * 4218306354666L) - ((int128)tmp_q[1] * 1889351483184L) - ((int128)tmp_q[2] * 1124453604936L) + ((int128)tmp_q[3] * 1509488305121L) - ((int128)tmp_q[4] * 2057371408626L) + ((int128)tmp_q[5] * 1285671898092L) - ((int128)tmp_q[6] * 380872577557L) - ((((int128)tmp_q[7] * 584509749053L) - ((int128)tmp_q[8] * 1046833896412L)) * 4);
	tmp_zero[7] = -((int128)tmp_q[0] * 1046833896412L) - ((int128)tmp_q[1] * 4218306354666L) - ((int128)tmp_q[2] * 1889351483184L) - ((int128)tmp_q[3] * 1124453604936L) + ((int128)tmp_q[4] * 1509488305121L) - ((int128)tmp_q[5] * 2057371408626L) + ((int128)tmp_q[6] * 1285671898092L) - ((int128)tmp_q[7] * 380872577557L) - ((int128)tmp_q[8] * 2338038996212L);
	tmp_zero[8] = ((int128)tmp_q[0] * 584509749053L) - ((int128)tmp_q[1] * 1046833896412L) - ((int128)tmp_q[2] * 4218306354666L) - ((int128)tmp_q[3] * 1889351483184L) - ((int128)tmp_q[4] * 1124453604936L) + ((int128)tmp_q[5] * 1509488305121L) - ((int128)tmp_q[6] * 2057371408626L) + ((int128)tmp_q[7] * 1285671898092L) - ((int128)tmp_q[8] * 380872577557L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
}

