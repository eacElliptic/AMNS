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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 10);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3907906721583516969UL) + ((((uint64_t)op[1] * 749902408641677508UL) + ((uint64_t)op[2] * 17818432649174562329UL) + ((uint64_t)op[3] * 12329333998969498718UL) + ((uint64_t)op[4] * 11659165117088931433UL) + ((uint64_t)op[5] * 9140299057674061060UL) + ((uint64_t)op[6] * 8281476772544062835UL) + ((uint64_t)op[7] * 14565614424064929551UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 14565614424064929551UL) + ((uint64_t)op[1] * 3907906721583516969UL) + ((((uint64_t)op[2] * 749902408641677508UL) + ((uint64_t)op[3] * 17818432649174562329UL) + ((uint64_t)op[4] * 12329333998969498718UL) + ((uint64_t)op[5] * 11659165117088931433UL) + ((uint64_t)op[6] * 9140299057674061060UL) + ((uint64_t)op[7] * 8281476772544062835UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 8281476772544062835UL) + ((uint64_t)op[1] * 14565614424064929551UL) + ((uint64_t)op[2] * 3907906721583516969UL) + ((((uint64_t)op[3] * 749902408641677508UL) + ((uint64_t)op[4] * 17818432649174562329UL) + ((uint64_t)op[5] * 12329333998969498718UL) + ((uint64_t)op[6] * 11659165117088931433UL) + ((uint64_t)op[7] * 9140299057674061060UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 9140299057674061060UL) + ((uint64_t)op[1] * 8281476772544062835UL) + ((uint64_t)op[2] * 14565614424064929551UL) + ((uint64_t)op[3] * 3907906721583516969UL) + ((((uint64_t)op[4] * 749902408641677508UL) + ((uint64_t)op[5] * 17818432649174562329UL) + ((uint64_t)op[6] * 12329333998969498718UL) + ((uint64_t)op[7] * 11659165117088931433UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 11659165117088931433UL) + ((uint64_t)op[1] * 9140299057674061060UL) + ((uint64_t)op[2] * 8281476772544062835UL) + ((uint64_t)op[3] * 14565614424064929551UL) + ((uint64_t)op[4] * 3907906721583516969UL) + ((((uint64_t)op[5] * 749902408641677508UL) + ((uint64_t)op[6] * 17818432649174562329UL) + ((uint64_t)op[7] * 12329333998969498718UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 12329333998969498718UL) + ((uint64_t)op[1] * 11659165117088931433UL) + ((uint64_t)op[2] * 9140299057674061060UL) + ((uint64_t)op[3] * 8281476772544062835UL) + ((uint64_t)op[4] * 14565614424064929551UL) + ((uint64_t)op[5] * 3907906721583516969UL) + ((((uint64_t)op[6] * 749902408641677508UL) + ((uint64_t)op[7] * 17818432649174562329UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 17818432649174562329UL) + ((uint64_t)op[1] * 12329333998969498718UL) + ((uint64_t)op[2] * 11659165117088931433UL) + ((uint64_t)op[3] * 9140299057674061060UL) + ((uint64_t)op[4] * 8281476772544062835UL) + ((uint64_t)op[5] * 14565614424064929551UL) + ((uint64_t)op[6] * 3907906721583516969UL) + ((uint64_t)op[7] * 3749512043208387540UL);
	tmp_q[7] = ((uint64_t)op[0] * 749902408641677508UL) + ((uint64_t)op[1] * 17818432649174562329UL) + ((uint64_t)op[2] * 12329333998969498718UL) + ((uint64_t)op[3] * 11659165117088931433UL) + ((uint64_t)op[4] * 9140299057674061060UL) + ((uint64_t)op[5] * 8281476772544062835UL) + ((uint64_t)op[6] * 14565614424064929551UL) + ((uint64_t)op[7] * 3907906721583516969UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 70712465573807L) + ((((int128)tmp_q[1] * 91804043609201L) - ((int128)tmp_q[2] * 37280810914111L) - ((int128)tmp_q[3] * 150413875443768L) + ((int128)tmp_q[4] * 153678404809349L) + ((int128)tmp_q[5] * 64591116018126L) + ((int128)tmp_q[6] * 97437420048269L) - ((int128)tmp_q[7] * 121907147142580L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 121907147142580L) + ((int128)tmp_q[1] * 70712465573807L) + ((((int128)tmp_q[2] * 91804043609201L) - ((int128)tmp_q[3] * 37280810914111L) - ((int128)tmp_q[4] * 150413875443768L) + ((int128)tmp_q[5] * 153678404809349L) + ((int128)tmp_q[6] * 64591116018126L) + ((int128)tmp_q[7] * 97437420048269L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 97437420048269L) - ((int128)tmp_q[1] * 121907147142580L) + ((int128)tmp_q[2] * 70712465573807L) + ((((int128)tmp_q[3] * 91804043609201L) - ((int128)tmp_q[4] * 37280810914111L) - ((int128)tmp_q[5] * 150413875443768L) + ((int128)tmp_q[6] * 153678404809349L) + ((int128)tmp_q[7] * 64591116018126L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 64591116018126L) + ((int128)tmp_q[1] * 97437420048269L) - ((int128)tmp_q[2] * 121907147142580L) + ((int128)tmp_q[3] * 70712465573807L) + ((((int128)tmp_q[4] * 91804043609201L) - ((int128)tmp_q[5] * 37280810914111L) - ((int128)tmp_q[6] * 150413875443768L) + ((int128)tmp_q[7] * 153678404809349L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 153678404809349L) + ((int128)tmp_q[1] * 64591116018126L) + ((int128)tmp_q[2] * 97437420048269L) - ((int128)tmp_q[3] * 121907147142580L) + ((int128)tmp_q[4] * 70712465573807L) + ((((int128)tmp_q[5] * 91804043609201L) - ((int128)tmp_q[6] * 37280810914111L) - ((int128)tmp_q[7] * 150413875443768L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 150413875443768L) + ((int128)tmp_q[1] * 153678404809349L) + ((int128)tmp_q[2] * 64591116018126L) + ((int128)tmp_q[3] * 97437420048269L) - ((int128)tmp_q[4] * 121907147142580L) + ((int128)tmp_q[5] * 70712465573807L) + ((((int128)tmp_q[6] * 91804043609201L) - ((int128)tmp_q[7] * 37280810914111L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 37280810914111L) - ((int128)tmp_q[1] * 150413875443768L) + ((int128)tmp_q[2] * 153678404809349L) + ((int128)tmp_q[3] * 64591116018126L) + ((int128)tmp_q[4] * 97437420048269L) - ((int128)tmp_q[5] * 121907147142580L) + ((int128)tmp_q[6] * 70712465573807L) + ((int128)tmp_q[7] * 459020218046005L);
	tmp_zero[7] = ((int128)tmp_q[0] * 91804043609201L) - ((int128)tmp_q[1] * 37280810914111L) - ((int128)tmp_q[2] * 150413875443768L) + ((int128)tmp_q[3] * 153678404809349L) + ((int128)tmp_q[4] * 64591116018126L) + ((int128)tmp_q[5] * 97437420048269L) - ((int128)tmp_q[6] * 121907147142580L) + ((int128)tmp_q[7] * 70712465573807L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

