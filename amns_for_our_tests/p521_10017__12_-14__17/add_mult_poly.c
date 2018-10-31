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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 14);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 14);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 14);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 14);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 14);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 14);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 14);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 14);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 14);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 14);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 14);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 28);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 28);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 28);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 28);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 28);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 14);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17301829976254679583UL) + ((((uint64_t)op[1] * 3501461844349044799UL) + ((uint64_t)op[2] * 2181466012582691264UL) + ((uint64_t)op[3] * 18172868929929883528UL) + ((uint64_t)op[4] * 15360855308412574016UL) + ((uint64_t)op[5] * 12211347664417793980UL) + ((uint64_t)op[6] * 4155218134373462662UL) + ((uint64_t)op[7] * 10777683127923592730UL) + ((uint64_t)op[8] * 3375595182309003684UL) + ((uint64_t)op[9] * 14548701981487962793UL) + ((uint64_t)op[10] * 16644288787730752746UL) + ((uint64_t)op[11] * 15686842423776737731UL)) * 18446744073709551602);
	tmp_q[1] = ((uint64_t)op[0] * 15686842423776737731UL) + ((uint64_t)op[1] * 17301829976254679583UL) + ((((uint64_t)op[2] * 3501461844349044799UL) + ((uint64_t)op[3] * 2181466012582691264UL) + ((uint64_t)op[4] * 18172868929929883528UL) + ((uint64_t)op[5] * 15360855308412574016UL) + ((uint64_t)op[6] * 12211347664417793980UL) + ((uint64_t)op[7] * 4155218134373462662UL) + ((uint64_t)op[8] * 10777683127923592730UL) + ((uint64_t)op[9] * 3375595182309003684UL) + ((uint64_t)op[10] * 14548701981487962793UL) + ((uint64_t)op[11] * 16644288787730752746UL)) * 18446744073709551602);
	tmp_q[2] = ((uint64_t)op[0] * 16644288787730752746UL) + ((uint64_t)op[1] * 15686842423776737731UL) + ((uint64_t)op[2] * 17301829976254679583UL) + ((((uint64_t)op[3] * 3501461844349044799UL) + ((uint64_t)op[4] * 2181466012582691264UL) + ((uint64_t)op[5] * 18172868929929883528UL) + ((uint64_t)op[6] * 15360855308412574016UL) + ((uint64_t)op[7] * 12211347664417793980UL) + ((uint64_t)op[8] * 4155218134373462662UL) + ((uint64_t)op[9] * 10777683127923592730UL) + ((uint64_t)op[10] * 3375595182309003684UL) + ((uint64_t)op[11] * 14548701981487962793UL)) * 18446744073709551602);
	tmp_q[3] = ((uint64_t)op[0] * 14548701981487962793UL) + ((uint64_t)op[1] * 16644288787730752746UL) + ((uint64_t)op[2] * 15686842423776737731UL) + ((uint64_t)op[3] * 17301829976254679583UL) + ((((uint64_t)op[4] * 3501461844349044799UL) + ((uint64_t)op[5] * 2181466012582691264UL) + ((uint64_t)op[6] * 18172868929929883528UL) + ((uint64_t)op[7] * 15360855308412574016UL) + ((uint64_t)op[8] * 12211347664417793980UL) + ((uint64_t)op[9] * 4155218134373462662UL) + ((uint64_t)op[10] * 10777683127923592730UL) + ((uint64_t)op[11] * 3375595182309003684UL)) * 18446744073709551602);
	tmp_q[4] = ((uint64_t)op[0] * 3375595182309003684UL) + ((uint64_t)op[1] * 14548701981487962793UL) + ((uint64_t)op[2] * 16644288787730752746UL) + ((uint64_t)op[3] * 15686842423776737731UL) + ((uint64_t)op[4] * 17301829976254679583UL) + ((((uint64_t)op[5] * 3501461844349044799UL) + ((uint64_t)op[6] * 2181466012582691264UL) + ((uint64_t)op[7] * 18172868929929883528UL) + ((uint64_t)op[8] * 15360855308412574016UL) + ((uint64_t)op[9] * 12211347664417793980UL) + ((uint64_t)op[10] * 4155218134373462662UL) + ((uint64_t)op[11] * 10777683127923592730UL)) * 18446744073709551602);
	tmp_q[5] = ((uint64_t)op[0] * 10777683127923592730UL) + ((uint64_t)op[1] * 3375595182309003684UL) + ((uint64_t)op[2] * 14548701981487962793UL) + ((uint64_t)op[3] * 16644288787730752746UL) + ((uint64_t)op[4] * 15686842423776737731UL) + ((uint64_t)op[5] * 17301829976254679583UL) + ((((uint64_t)op[6] * 3501461844349044799UL) + ((uint64_t)op[7] * 2181466012582691264UL) + ((uint64_t)op[8] * 18172868929929883528UL) + ((uint64_t)op[9] * 15360855308412574016UL) + ((uint64_t)op[10] * 12211347664417793980UL) + ((uint64_t)op[11] * 4155218134373462662UL)) * 18446744073709551602);
	tmp_q[6] = ((uint64_t)op[0] * 4155218134373462662UL) + ((uint64_t)op[1] * 10777683127923592730UL) + ((uint64_t)op[2] * 3375595182309003684UL) + ((uint64_t)op[3] * 14548701981487962793UL) + ((uint64_t)op[4] * 16644288787730752746UL) + ((uint64_t)op[5] * 15686842423776737731UL) + ((uint64_t)op[6] * 17301829976254679583UL) + ((((uint64_t)op[7] * 3501461844349044799UL) + ((uint64_t)op[8] * 2181466012582691264UL) + ((uint64_t)op[9] * 18172868929929883528UL) + ((uint64_t)op[10] * 15360855308412574016UL) + ((uint64_t)op[11] * 12211347664417793980UL)) * 18446744073709551602);
	tmp_q[7] = ((uint64_t)op[0] * 12211347664417793980UL) + ((uint64_t)op[1] * 4155218134373462662UL) + ((uint64_t)op[2] * 10777683127923592730UL) + ((uint64_t)op[3] * 3375595182309003684UL) + ((uint64_t)op[4] * 14548701981487962793UL) + ((uint64_t)op[5] * 16644288787730752746UL) + ((uint64_t)op[6] * 15686842423776737731UL) + ((uint64_t)op[7] * 17301829976254679583UL) + ((((uint64_t)op[8] * 3501461844349044799UL) + ((uint64_t)op[9] * 2181466012582691264UL) + ((uint64_t)op[10] * 18172868929929883528UL) + ((uint64_t)op[11] * 15360855308412574016UL)) * 18446744073709551602);
	tmp_q[8] = ((uint64_t)op[0] * 15360855308412574016UL) + ((uint64_t)op[1] * 12211347664417793980UL) + ((uint64_t)op[2] * 4155218134373462662UL) + ((uint64_t)op[3] * 10777683127923592730UL) + ((uint64_t)op[4] * 3375595182309003684UL) + ((uint64_t)op[5] * 14548701981487962793UL) + ((uint64_t)op[6] * 16644288787730752746UL) + ((uint64_t)op[7] * 15686842423776737731UL) + ((uint64_t)op[8] * 17301829976254679583UL) + ((((uint64_t)op[9] * 3501461844349044799UL) + ((uint64_t)op[10] * 2181466012582691264UL) + ((uint64_t)op[11] * 18172868929929883528UL)) * 18446744073709551602);
	tmp_q[9] = ((uint64_t)op[0] * 18172868929929883528UL) + ((uint64_t)op[1] * 15360855308412574016UL) + ((uint64_t)op[2] * 12211347664417793980UL) + ((uint64_t)op[3] * 4155218134373462662UL) + ((uint64_t)op[4] * 10777683127923592730UL) + ((uint64_t)op[5] * 3375595182309003684UL) + ((uint64_t)op[6] * 14548701981487962793UL) + ((uint64_t)op[7] * 16644288787730752746UL) + ((uint64_t)op[8] * 15686842423776737731UL) + ((uint64_t)op[9] * 17301829976254679583UL) + ((((uint64_t)op[10] * 3501461844349044799UL) + ((uint64_t)op[11] * 2181466012582691264UL)) * 18446744073709551602);
	tmp_q[10] = ((uint64_t)op[0] * 2181466012582691264UL) + ((uint64_t)op[1] * 18172868929929883528UL) + ((uint64_t)op[2] * 15360855308412574016UL) + ((uint64_t)op[3] * 12211347664417793980UL) + ((uint64_t)op[4] * 4155218134373462662UL) + ((uint64_t)op[5] * 10777683127923592730UL) + ((uint64_t)op[6] * 3375595182309003684UL) + ((uint64_t)op[7] * 14548701981487962793UL) + ((uint64_t)op[8] * 16644288787730752746UL) + ((uint64_t)op[9] * 15686842423776737731UL) + ((uint64_t)op[10] * 17301829976254679583UL) + ((uint64_t)op[11] * 6319766400242027662UL);
	tmp_q[11] = ((uint64_t)op[0] * 3501461844349044799UL) + ((uint64_t)op[1] * 2181466012582691264UL) + ((uint64_t)op[2] * 18172868929929883528UL) + ((uint64_t)op[3] * 15360855308412574016UL) + ((uint64_t)op[4] * 12211347664417793980UL) + ((uint64_t)op[5] * 4155218134373462662UL) + ((uint64_t)op[6] * 10777683127923592730UL) + ((uint64_t)op[7] * 3375595182309003684UL) + ((uint64_t)op[8] * 14548701981487962793UL) + ((uint64_t)op[9] * 16644288787730752746UL) + ((uint64_t)op[10] * 15686842423776737731UL) + ((uint64_t)op[11] * 17301829976254679583UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4846767595125L) - ((((int128)tmp_q[1] * 2565426037200L) - ((int128)tmp_q[2] * 4436565457764L) - ((int128)tmp_q[3] * 5196053901879L) + ((int128)tmp_q[4] * 506930601087L) - ((int128)tmp_q[5] * 2652814245671L) + ((int128)tmp_q[6] * 3881935874488L) + ((int128)tmp_q[7] * 5249161693896L) + ((int128)tmp_q[8] * 2538992441587L) + ((int128)tmp_q[9] * 2840122872592L) + ((int128)tmp_q[10] * 3357863970865L) - ((int128)tmp_q[11] * 4470228205611L)) * 14);
	tmp_zero[1] = -((int128)tmp_q[0] * 4470228205611L) + ((int128)tmp_q[1] * 4846767595125L) - ((((int128)tmp_q[2] * 2565426037200L) - ((int128)tmp_q[3] * 4436565457764L) - ((int128)tmp_q[4] * 5196053901879L) + ((int128)tmp_q[5] * 506930601087L) - ((int128)tmp_q[6] * 2652814245671L) + ((int128)tmp_q[7] * 3881935874488L) + ((int128)tmp_q[8] * 5249161693896L) + ((int128)tmp_q[9] * 2538992441587L) + ((int128)tmp_q[10] * 2840122872592L) + ((int128)tmp_q[11] * 3357863970865L)) * 14);
	tmp_zero[2] = ((int128)tmp_q[0] * 3357863970865L) - ((int128)tmp_q[1] * 4470228205611L) + ((int128)tmp_q[2] * 4846767595125L) - ((((int128)tmp_q[3] * 2565426037200L) - ((int128)tmp_q[4] * 4436565457764L) - ((int128)tmp_q[5] * 5196053901879L) + ((int128)tmp_q[6] * 506930601087L) - ((int128)tmp_q[7] * 2652814245671L) + ((int128)tmp_q[8] * 3881935874488L) + ((int128)tmp_q[9] * 5249161693896L) + ((int128)tmp_q[10] * 2538992441587L) + ((int128)tmp_q[11] * 2840122872592L)) * 14);
	tmp_zero[3] = ((int128)tmp_q[0] * 2840122872592L) + ((int128)tmp_q[1] * 3357863970865L) - ((int128)tmp_q[2] * 4470228205611L) + ((int128)tmp_q[3] * 4846767595125L) - ((((int128)tmp_q[4] * 2565426037200L) - ((int128)tmp_q[5] * 4436565457764L) - ((int128)tmp_q[6] * 5196053901879L) + ((int128)tmp_q[7] * 506930601087L) - ((int128)tmp_q[8] * 2652814245671L) + ((int128)tmp_q[9] * 3881935874488L) + ((int128)tmp_q[10] * 5249161693896L) + ((int128)tmp_q[11] * 2538992441587L)) * 14);
	tmp_zero[4] = ((int128)tmp_q[0] * 2538992441587L) + ((int128)tmp_q[1] * 2840122872592L) + ((int128)tmp_q[2] * 3357863970865L) - ((int128)tmp_q[3] * 4470228205611L) + ((int128)tmp_q[4] * 4846767595125L) - ((((int128)tmp_q[5] * 2565426037200L) - ((int128)tmp_q[6] * 4436565457764L) - ((int128)tmp_q[7] * 5196053901879L) + ((int128)tmp_q[8] * 506930601087L) - ((int128)tmp_q[9] * 2652814245671L) + ((int128)tmp_q[10] * 3881935874488L) + ((int128)tmp_q[11] * 5249161693896L)) * 14);
	tmp_zero[5] = ((int128)tmp_q[0] * 5249161693896L) + ((int128)tmp_q[1] * 2538992441587L) + ((int128)tmp_q[2] * 2840122872592L) + ((int128)tmp_q[3] * 3357863970865L) - ((int128)tmp_q[4] * 4470228205611L) + ((int128)tmp_q[5] * 4846767595125L) - ((((int128)tmp_q[6] * 2565426037200L) - ((int128)tmp_q[7] * 4436565457764L) - ((int128)tmp_q[8] * 5196053901879L) + ((int128)tmp_q[9] * 506930601087L) - ((int128)tmp_q[10] * 2652814245671L) + ((int128)tmp_q[11] * 3881935874488L)) * 14);
	tmp_zero[6] = ((int128)tmp_q[0] * 3881935874488L) + ((int128)tmp_q[1] * 5249161693896L) + ((int128)tmp_q[2] * 2538992441587L) + ((int128)tmp_q[3] * 2840122872592L) + ((int128)tmp_q[4] * 3357863970865L) - ((int128)tmp_q[5] * 4470228205611L) + ((int128)tmp_q[6] * 4846767595125L) - ((((int128)tmp_q[7] * 2565426037200L) - ((int128)tmp_q[8] * 4436565457764L) - ((int128)tmp_q[9] * 5196053901879L) + ((int128)tmp_q[10] * 506930601087L) - ((int128)tmp_q[11] * 2652814245671L)) * 14);
	tmp_zero[7] = -((int128)tmp_q[0] * 2652814245671L) + ((int128)tmp_q[1] * 3881935874488L) + ((int128)tmp_q[2] * 5249161693896L) + ((int128)tmp_q[3] * 2538992441587L) + ((int128)tmp_q[4] * 2840122872592L) + ((int128)tmp_q[5] * 3357863970865L) - ((int128)tmp_q[6] * 4470228205611L) + ((int128)tmp_q[7] * 4846767595125L) - ((((int128)tmp_q[8] * 2565426037200L) - ((int128)tmp_q[9] * 4436565457764L) - ((int128)tmp_q[10] * 5196053901879L) + ((int128)tmp_q[11] * 506930601087L)) * 14);
	tmp_zero[8] = ((int128)tmp_q[0] * 506930601087L) - ((int128)tmp_q[1] * 2652814245671L) + ((int128)tmp_q[2] * 3881935874488L) + ((int128)tmp_q[3] * 5249161693896L) + ((int128)tmp_q[4] * 2538992441587L) + ((int128)tmp_q[5] * 2840122872592L) + ((int128)tmp_q[6] * 3357863970865L) - ((int128)tmp_q[7] * 4470228205611L) + ((int128)tmp_q[8] * 4846767595125L) - ((((int128)tmp_q[9] * 2565426037200L) - ((int128)tmp_q[10] * 4436565457764L) - ((int128)tmp_q[11] * 5196053901879L)) * 14);
	tmp_zero[9] = -((int128)tmp_q[0] * 5196053901879L) + ((int128)tmp_q[1] * 506930601087L) - ((int128)tmp_q[2] * 2652814245671L) + ((int128)tmp_q[3] * 3881935874488L) + ((int128)tmp_q[4] * 5249161693896L) + ((int128)tmp_q[5] * 2538992441587L) + ((int128)tmp_q[6] * 2840122872592L) + ((int128)tmp_q[7] * 3357863970865L) - ((int128)tmp_q[8] * 4470228205611L) + ((int128)tmp_q[9] * 4846767595125L) - ((((int128)tmp_q[10] * 2565426037200L) - ((int128)tmp_q[11] * 4436565457764L)) * 14);
	tmp_zero[10] = -((int128)tmp_q[0] * 4436565457764L) - ((int128)tmp_q[1] * 5196053901879L) + ((int128)tmp_q[2] * 506930601087L) - ((int128)tmp_q[3] * 2652814245671L) + ((int128)tmp_q[4] * 3881935874488L) + ((int128)tmp_q[5] * 5249161693896L) + ((int128)tmp_q[6] * 2538992441587L) + ((int128)tmp_q[7] * 2840122872592L) + ((int128)tmp_q[8] * 3357863970865L) - ((int128)tmp_q[9] * 4470228205611L) + ((int128)tmp_q[10] * 4846767595125L) - ((int128)tmp_q[11] * 35915964520800L);
	tmp_zero[11] = ((int128)tmp_q[0] * 2565426037200L) - ((int128)tmp_q[1] * 4436565457764L) - ((int128)tmp_q[2] * 5196053901879L) + ((int128)tmp_q[3] * 506930601087L) - ((int128)tmp_q[4] * 2652814245671L) + ((int128)tmp_q[5] * 3881935874488L) + ((int128)tmp_q[6] * 5249161693896L) + ((int128)tmp_q[7] * 2538992441587L) + ((int128)tmp_q[8] * 2840122872592L) + ((int128)tmp_q[9] * 3357863970865L) - ((int128)tmp_q[10] * 4470228205611L) + ((int128)tmp_q[11] * 4846767595125L);

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
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
	rop[11] = (op[11] + tmp_zero[11]) >> WORD_SIZE;
}

