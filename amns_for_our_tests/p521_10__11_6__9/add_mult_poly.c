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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 6);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 6);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 12);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16064643716726582399UL) + ((((uint64_t)op[1] * 13857023424550562548UL) + ((uint64_t)op[2] * 817959567704805977UL) + ((uint64_t)op[3] * 14248326942373012932UL) + ((uint64_t)op[4] * 9757534607716322644UL) + ((uint64_t)op[5] * 16519967568496872305UL) + ((uint64_t)op[6] * 3186450555624326980UL) + ((uint64_t)op[7] * 551340085675485405UL) + ((uint64_t)op[8] * 1419155637304691518UL) + ((uint64_t)op[9] * 13171702751325253220UL) + ((uint64_t)op[10] * 5647381080711183468UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 5647381080711183468UL) + ((uint64_t)op[1] * 16064643716726582399UL) + ((((uint64_t)op[2] * 13857023424550562548UL) + ((uint64_t)op[3] * 817959567704805977UL) + ((uint64_t)op[4] * 14248326942373012932UL) + ((uint64_t)op[5] * 9757534607716322644UL) + ((uint64_t)op[6] * 16519967568496872305UL) + ((uint64_t)op[7] * 3186450555624326980UL) + ((uint64_t)op[8] * 551340085675485405UL) + ((uint64_t)op[9] * 1419155637304691518UL) + ((uint64_t)op[10] * 13171702751325253220UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 13171702751325253220UL) + ((uint64_t)op[1] * 5647381080711183468UL) + ((uint64_t)op[2] * 16064643716726582399UL) + ((((uint64_t)op[3] * 13857023424550562548UL) + ((uint64_t)op[4] * 817959567704805977UL) + ((uint64_t)op[5] * 14248326942373012932UL) + ((uint64_t)op[6] * 9757534607716322644UL) + ((uint64_t)op[7] * 16519967568496872305UL) + ((uint64_t)op[8] * 3186450555624326980UL) + ((uint64_t)op[9] * 551340085675485405UL) + ((uint64_t)op[10] * 1419155637304691518UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 1419155637304691518UL) + ((uint64_t)op[1] * 13171702751325253220UL) + ((uint64_t)op[2] * 5647381080711183468UL) + ((uint64_t)op[3] * 16064643716726582399UL) + ((((uint64_t)op[4] * 13857023424550562548UL) + ((uint64_t)op[5] * 817959567704805977UL) + ((uint64_t)op[6] * 14248326942373012932UL) + ((uint64_t)op[7] * 9757534607716322644UL) + ((uint64_t)op[8] * 16519967568496872305UL) + ((uint64_t)op[9] * 3186450555624326980UL) + ((uint64_t)op[10] * 551340085675485405UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 551340085675485405UL) + ((uint64_t)op[1] * 1419155637304691518UL) + ((uint64_t)op[2] * 13171702751325253220UL) + ((uint64_t)op[3] * 5647381080711183468UL) + ((uint64_t)op[4] * 16064643716726582399UL) + ((((uint64_t)op[5] * 13857023424550562548UL) + ((uint64_t)op[6] * 817959567704805977UL) + ((uint64_t)op[7] * 14248326942373012932UL) + ((uint64_t)op[8] * 9757534607716322644UL) + ((uint64_t)op[9] * 16519967568496872305UL) + ((uint64_t)op[10] * 3186450555624326980UL)) * 6);
	tmp_q[5] = ((uint64_t)op[0] * 3186450555624326980UL) + ((uint64_t)op[1] * 551340085675485405UL) + ((uint64_t)op[2] * 1419155637304691518UL) + ((uint64_t)op[3] * 13171702751325253220UL) + ((uint64_t)op[4] * 5647381080711183468UL) + ((uint64_t)op[5] * 16064643716726582399UL) + ((((uint64_t)op[6] * 13857023424550562548UL) + ((uint64_t)op[7] * 817959567704805977UL) + ((uint64_t)op[8] * 14248326942373012932UL) + ((uint64_t)op[9] * 9757534607716322644UL) + ((uint64_t)op[10] * 16519967568496872305UL)) * 6);
	tmp_q[6] = ((uint64_t)op[0] * 16519967568496872305UL) + ((uint64_t)op[1] * 3186450555624326980UL) + ((uint64_t)op[2] * 551340085675485405UL) + ((uint64_t)op[3] * 1419155637304691518UL) + ((uint64_t)op[4] * 13171702751325253220UL) + ((uint64_t)op[5] * 5647381080711183468UL) + ((uint64_t)op[6] * 16064643716726582399UL) + ((((uint64_t)op[7] * 13857023424550562548UL) + ((uint64_t)op[8] * 817959567704805977UL) + ((uint64_t)op[9] * 14248326942373012932UL) + ((uint64_t)op[10] * 9757534607716322644UL)) * 6);
	tmp_q[7] = ((uint64_t)op[0] * 9757534607716322644UL) + ((uint64_t)op[1] * 16519967568496872305UL) + ((uint64_t)op[2] * 3186450555624326980UL) + ((uint64_t)op[3] * 551340085675485405UL) + ((uint64_t)op[4] * 1419155637304691518UL) + ((uint64_t)op[5] * 13171702751325253220UL) + ((uint64_t)op[6] * 5647381080711183468UL) + ((uint64_t)op[7] * 16064643716726582399UL) + ((((uint64_t)op[8] * 13857023424550562548UL) + ((uint64_t)op[9] * 817959567704805977UL) + ((uint64_t)op[10] * 14248326942373012932UL)) * 6);
	tmp_q[8] = ((uint64_t)op[0] * 14248326942373012932UL) + ((uint64_t)op[1] * 9757534607716322644UL) + ((uint64_t)op[2] * 16519967568496872305UL) + ((uint64_t)op[3] * 3186450555624326980UL) + ((uint64_t)op[4] * 551340085675485405UL) + ((uint64_t)op[5] * 1419155637304691518UL) + ((uint64_t)op[6] * 13171702751325253220UL) + ((uint64_t)op[7] * 5647381080711183468UL) + ((uint64_t)op[8] * 16064643716726582399UL) + ((((uint64_t)op[9] * 13857023424550562548UL) + ((uint64_t)op[10] * 817959567704805977UL)) * 6);
	tmp_q[9] = ((uint64_t)op[0] * 817959567704805977UL) + ((uint64_t)op[1] * 14248326942373012932UL) + ((uint64_t)op[2] * 9757534607716322644UL) + ((uint64_t)op[3] * 16519967568496872305UL) + ((uint64_t)op[4] * 3186450555624326980UL) + ((uint64_t)op[5] * 551340085675485405UL) + ((uint64_t)op[6] * 1419155637304691518UL) + ((uint64_t)op[7] * 13171702751325253220UL) + ((uint64_t)op[8] * 5647381080711183468UL) + ((uint64_t)op[9] * 16064643716726582399UL) + ((uint64_t)op[10] * 9355164252465168824UL);
	tmp_q[10] = ((uint64_t)op[0] * 13857023424550562548UL) + ((uint64_t)op[1] * 817959567704805977UL) + ((uint64_t)op[2] * 14248326942373012932UL) + ((uint64_t)op[3] * 9757534607716322644UL) + ((uint64_t)op[4] * 16519967568496872305UL) + ((uint64_t)op[5] * 3186450555624326980UL) + ((uint64_t)op[6] * 551340085675485405UL) + ((uint64_t)op[7] * 1419155637304691518UL) + ((uint64_t)op[8] * 13171702751325253220UL) + ((uint64_t)op[9] * 5647381080711183468UL) + ((uint64_t)op[10] * 16064643716726582399UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 17934360050131L) + ((-((int128)tmp_q[1] * 41958540927472L) - ((int128)tmp_q[2] * 111500408692861L) - ((int128)tmp_q[3] * 50704419087123L) + ((int128)tmp_q[4] * 14190054661940L) + ((int128)tmp_q[5] * 76026167680627L) + ((int128)tmp_q[6] * 104154394065080L) - ((int128)tmp_q[7] * 9534168169979L) + ((int128)tmp_q[8] * 35766413050492L) - ((int128)tmp_q[9] * 110327695102588L) - ((int128)tmp_q[10] * 35516174344544L)) * 6);
	tmp_zero[1] = -((int128)tmp_q[0] * 35516174344544L) - ((int128)tmp_q[1] * 17934360050131L) + ((-((int128)tmp_q[2] * 41958540927472L) - ((int128)tmp_q[3] * 111500408692861L) - ((int128)tmp_q[4] * 50704419087123L) + ((int128)tmp_q[5] * 14190054661940L) + ((int128)tmp_q[6] * 76026167680627L) + ((int128)tmp_q[7] * 104154394065080L) - ((int128)tmp_q[8] * 9534168169979L) + ((int128)tmp_q[9] * 35766413050492L) - ((int128)tmp_q[10] * 110327695102588L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 110327695102588L) - ((int128)tmp_q[1] * 35516174344544L) - ((int128)tmp_q[2] * 17934360050131L) + ((-((int128)tmp_q[3] * 41958540927472L) - ((int128)tmp_q[4] * 111500408692861L) - ((int128)tmp_q[5] * 50704419087123L) + ((int128)tmp_q[6] * 14190054661940L) + ((int128)tmp_q[7] * 76026167680627L) + ((int128)tmp_q[8] * 104154394065080L) - ((int128)tmp_q[9] * 9534168169979L) + ((int128)tmp_q[10] * 35766413050492L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 35766413050492L) - ((int128)tmp_q[1] * 110327695102588L) - ((int128)tmp_q[2] * 35516174344544L) - ((int128)tmp_q[3] * 17934360050131L) + ((-((int128)tmp_q[4] * 41958540927472L) - ((int128)tmp_q[5] * 111500408692861L) - ((int128)tmp_q[6] * 50704419087123L) + ((int128)tmp_q[7] * 14190054661940L) + ((int128)tmp_q[8] * 76026167680627L) + ((int128)tmp_q[9] * 104154394065080L) - ((int128)tmp_q[10] * 9534168169979L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 9534168169979L) + ((int128)tmp_q[1] * 35766413050492L) - ((int128)tmp_q[2] * 110327695102588L) - ((int128)tmp_q[3] * 35516174344544L) - ((int128)tmp_q[4] * 17934360050131L) + ((-((int128)tmp_q[5] * 41958540927472L) - ((int128)tmp_q[6] * 111500408692861L) - ((int128)tmp_q[7] * 50704419087123L) + ((int128)tmp_q[8] * 14190054661940L) + ((int128)tmp_q[9] * 76026167680627L) + ((int128)tmp_q[10] * 104154394065080L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 104154394065080L) - ((int128)tmp_q[1] * 9534168169979L) + ((int128)tmp_q[2] * 35766413050492L) - ((int128)tmp_q[3] * 110327695102588L) - ((int128)tmp_q[4] * 35516174344544L) - ((int128)tmp_q[5] * 17934360050131L) + ((-((int128)tmp_q[6] * 41958540927472L) - ((int128)tmp_q[7] * 111500408692861L) - ((int128)tmp_q[8] * 50704419087123L) + ((int128)tmp_q[9] * 14190054661940L) + ((int128)tmp_q[10] * 76026167680627L)) * 6);
	tmp_zero[6] = ((int128)tmp_q[0] * 76026167680627L) + ((int128)tmp_q[1] * 104154394065080L) - ((int128)tmp_q[2] * 9534168169979L) + ((int128)tmp_q[3] * 35766413050492L) - ((int128)tmp_q[4] * 110327695102588L) - ((int128)tmp_q[5] * 35516174344544L) - ((int128)tmp_q[6] * 17934360050131L) + ((-((int128)tmp_q[7] * 41958540927472L) - ((int128)tmp_q[8] * 111500408692861L) - ((int128)tmp_q[9] * 50704419087123L) + ((int128)tmp_q[10] * 14190054661940L)) * 6);
	tmp_zero[7] = ((int128)tmp_q[0] * 14190054661940L) + ((int128)tmp_q[1] * 76026167680627L) + ((int128)tmp_q[2] * 104154394065080L) - ((int128)tmp_q[3] * 9534168169979L) + ((int128)tmp_q[4] * 35766413050492L) - ((int128)tmp_q[5] * 110327695102588L) - ((int128)tmp_q[6] * 35516174344544L) - ((int128)tmp_q[7] * 17934360050131L) + ((-((int128)tmp_q[8] * 41958540927472L) - ((int128)tmp_q[9] * 111500408692861L) - ((int128)tmp_q[10] * 50704419087123L)) * 6);
	tmp_zero[8] = -((int128)tmp_q[0] * 50704419087123L) + ((int128)tmp_q[1] * 14190054661940L) + ((int128)tmp_q[2] * 76026167680627L) + ((int128)tmp_q[3] * 104154394065080L) - ((int128)tmp_q[4] * 9534168169979L) + ((int128)tmp_q[5] * 35766413050492L) - ((int128)tmp_q[6] * 110327695102588L) - ((int128)tmp_q[7] * 35516174344544L) - ((int128)tmp_q[8] * 17934360050131L) + ((-((int128)tmp_q[9] * 41958540927472L) - ((int128)tmp_q[10] * 111500408692861L)) * 6);
	tmp_zero[9] = -((int128)tmp_q[0] * 111500408692861L) - ((int128)tmp_q[1] * 50704419087123L) + ((int128)tmp_q[2] * 14190054661940L) + ((int128)tmp_q[3] * 76026167680627L) + ((int128)tmp_q[4] * 104154394065080L) - ((int128)tmp_q[5] * 9534168169979L) + ((int128)tmp_q[6] * 35766413050492L) - ((int128)tmp_q[7] * 110327695102588L) - ((int128)tmp_q[8] * 35516174344544L) - ((int128)tmp_q[9] * 17934360050131L) - ((int128)tmp_q[10] * 251751245564832L);
	tmp_zero[10] = -((int128)tmp_q[0] * 41958540927472L) - ((int128)tmp_q[1] * 111500408692861L) - ((int128)tmp_q[2] * 50704419087123L) + ((int128)tmp_q[3] * 14190054661940L) + ((int128)tmp_q[4] * 76026167680627L) + ((int128)tmp_q[5] * 104154394065080L) - ((int128)tmp_q[6] * 9534168169979L) + ((int128)tmp_q[7] * 35766413050492L) - ((int128)tmp_q[8] * 110327695102588L) - ((int128)tmp_q[9] * 35516174344544L) - ((int128)tmp_q[10] * 17934360050131L);

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
}

