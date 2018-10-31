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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 5);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 5);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 10);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 5);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18374822472752314011UL) + ((((uint64_t)op[1] * 14081760712615705285UL) + ((uint64_t)op[2] * 8633588027085420619UL) + ((uint64_t)op[3] * 15059449663304469479UL) + ((uint64_t)op[4] * 2164598721435742336UL) + ((uint64_t)op[5] * 10530904766646719793UL) + ((uint64_t)op[6] * 17522661386641372934UL) + ((uint64_t)op[7] * 11787433664057123403UL) + ((uint64_t)op[8] * 1766888283642219124UL) + ((uint64_t)op[9] * 12789885499916064347UL) + ((uint64_t)op[10] * 4247907704266680004UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 4247907704266680004UL) + ((uint64_t)op[1] * 18374822472752314011UL) + ((((uint64_t)op[2] * 14081760712615705285UL) + ((uint64_t)op[3] * 8633588027085420619UL) + ((uint64_t)op[4] * 15059449663304469479UL) + ((uint64_t)op[5] * 2164598721435742336UL) + ((uint64_t)op[6] * 10530904766646719793UL) + ((uint64_t)op[7] * 17522661386641372934UL) + ((uint64_t)op[8] * 11787433664057123403UL) + ((uint64_t)op[9] * 1766888283642219124UL) + ((uint64_t)op[10] * 12789885499916064347UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 12789885499916064347UL) + ((uint64_t)op[1] * 4247907704266680004UL) + ((uint64_t)op[2] * 18374822472752314011UL) + ((((uint64_t)op[3] * 14081760712615705285UL) + ((uint64_t)op[4] * 8633588027085420619UL) + ((uint64_t)op[5] * 15059449663304469479UL) + ((uint64_t)op[6] * 2164598721435742336UL) + ((uint64_t)op[7] * 10530904766646719793UL) + ((uint64_t)op[8] * 17522661386641372934UL) + ((uint64_t)op[9] * 11787433664057123403UL) + ((uint64_t)op[10] * 1766888283642219124UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 1766888283642219124UL) + ((uint64_t)op[1] * 12789885499916064347UL) + ((uint64_t)op[2] * 4247907704266680004UL) + ((uint64_t)op[3] * 18374822472752314011UL) + ((((uint64_t)op[4] * 14081760712615705285UL) + ((uint64_t)op[5] * 8633588027085420619UL) + ((uint64_t)op[6] * 15059449663304469479UL) + ((uint64_t)op[7] * 2164598721435742336UL) + ((uint64_t)op[8] * 10530904766646719793UL) + ((uint64_t)op[9] * 17522661386641372934UL) + ((uint64_t)op[10] * 11787433664057123403UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 11787433664057123403UL) + ((uint64_t)op[1] * 1766888283642219124UL) + ((uint64_t)op[2] * 12789885499916064347UL) + ((uint64_t)op[3] * 4247907704266680004UL) + ((uint64_t)op[4] * 18374822472752314011UL) + ((((uint64_t)op[5] * 14081760712615705285UL) + ((uint64_t)op[6] * 8633588027085420619UL) + ((uint64_t)op[7] * 15059449663304469479UL) + ((uint64_t)op[8] * 2164598721435742336UL) + ((uint64_t)op[9] * 10530904766646719793UL) + ((uint64_t)op[10] * 17522661386641372934UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 17522661386641372934UL) + ((uint64_t)op[1] * 11787433664057123403UL) + ((uint64_t)op[2] * 1766888283642219124UL) + ((uint64_t)op[3] * 12789885499916064347UL) + ((uint64_t)op[4] * 4247907704266680004UL) + ((uint64_t)op[5] * 18374822472752314011UL) + ((((uint64_t)op[6] * 14081760712615705285UL) + ((uint64_t)op[7] * 8633588027085420619UL) + ((uint64_t)op[8] * 15059449663304469479UL) + ((uint64_t)op[9] * 2164598721435742336UL) + ((uint64_t)op[10] * 10530904766646719793UL)) * 18446744073709551611);
	tmp_q[6] = ((uint64_t)op[0] * 10530904766646719793UL) + ((uint64_t)op[1] * 17522661386641372934UL) + ((uint64_t)op[2] * 11787433664057123403UL) + ((uint64_t)op[3] * 1766888283642219124UL) + ((uint64_t)op[4] * 12789885499916064347UL) + ((uint64_t)op[5] * 4247907704266680004UL) + ((uint64_t)op[6] * 18374822472752314011UL) + ((((uint64_t)op[7] * 14081760712615705285UL) + ((uint64_t)op[8] * 8633588027085420619UL) + ((uint64_t)op[9] * 15059449663304469479UL) + ((uint64_t)op[10] * 2164598721435742336UL)) * 18446744073709551611);
	tmp_q[7] = ((uint64_t)op[0] * 2164598721435742336UL) + ((uint64_t)op[1] * 10530904766646719793UL) + ((uint64_t)op[2] * 17522661386641372934UL) + ((uint64_t)op[3] * 11787433664057123403UL) + ((uint64_t)op[4] * 1766888283642219124UL) + ((uint64_t)op[5] * 12789885499916064347UL) + ((uint64_t)op[6] * 4247907704266680004UL) + ((uint64_t)op[7] * 18374822472752314011UL) + ((((uint64_t)op[8] * 14081760712615705285UL) + ((uint64_t)op[9] * 8633588027085420619UL) + ((uint64_t)op[10] * 15059449663304469479UL)) * 18446744073709551611);
	tmp_q[8] = ((uint64_t)op[0] * 15059449663304469479UL) + ((uint64_t)op[1] * 2164598721435742336UL) + ((uint64_t)op[2] * 10530904766646719793UL) + ((uint64_t)op[3] * 17522661386641372934UL) + ((uint64_t)op[4] * 11787433664057123403UL) + ((uint64_t)op[5] * 1766888283642219124UL) + ((uint64_t)op[6] * 12789885499916064347UL) + ((uint64_t)op[7] * 4247907704266680004UL) + ((uint64_t)op[8] * 18374822472752314011UL) + ((((uint64_t)op[9] * 14081760712615705285UL) + ((uint64_t)op[10] * 8633588027085420619UL)) * 18446744073709551611);
	tmp_q[9] = ((uint64_t)op[0] * 8633588027085420619UL) + ((uint64_t)op[1] * 15059449663304469479UL) + ((uint64_t)op[2] * 2164598721435742336UL) + ((uint64_t)op[3] * 10530904766646719793UL) + ((uint64_t)op[4] * 17522661386641372934UL) + ((uint64_t)op[5] * 11787433664057123403UL) + ((uint64_t)op[6] * 1766888283642219124UL) + ((uint64_t)op[7] * 12789885499916064347UL) + ((uint64_t)op[8] * 4247907704266680004UL) + ((uint64_t)op[9] * 18374822472752314011UL) + ((uint64_t)op[10] * 3378172731759680039UL);
	tmp_q[10] = ((uint64_t)op[0] * 14081760712615705285UL) + ((uint64_t)op[1] * 8633588027085420619UL) + ((uint64_t)op[2] * 15059449663304469479UL) + ((uint64_t)op[3] * 2164598721435742336UL) + ((uint64_t)op[4] * 10530904766646719793UL) + ((uint64_t)op[5] * 17522661386641372934UL) + ((uint64_t)op[6] * 11787433664057123403UL) + ((uint64_t)op[7] * 1766888283642219124UL) + ((uint64_t)op[8] * 12789885499916064347UL) + ((uint64_t)op[9] * 4247907704266680004UL) + ((uint64_t)op[10] * 18374822472752314011UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 71243768560058L) - ((((int128)tmp_q[1] * 54399266451725L) - ((int128)tmp_q[2] * 40485352498776L) + ((int128)tmp_q[3] * 48902816151558L) - ((int128)tmp_q[4] * 55595276287113L) + ((int128)tmp_q[5] * 79614643509380L) + ((int128)tmp_q[6] * 83734412135570L) - ((int128)tmp_q[7] * 111991004117825L) + ((int128)tmp_q[8] * 101319890617906L) + ((int128)tmp_q[9] * 79948358425910L) - ((int128)tmp_q[10] * 12169016728036L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 12169016728036L) - ((int128)tmp_q[1] * 71243768560058L) - ((((int128)tmp_q[2] * 54399266451725L) - ((int128)tmp_q[3] * 40485352498776L) + ((int128)tmp_q[4] * 48902816151558L) - ((int128)tmp_q[5] * 55595276287113L) + ((int128)tmp_q[6] * 79614643509380L) + ((int128)tmp_q[7] * 83734412135570L) - ((int128)tmp_q[8] * 111991004117825L) + ((int128)tmp_q[9] * 101319890617906L) + ((int128)tmp_q[10] * 79948358425910L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 79948358425910L) - ((int128)tmp_q[1] * 12169016728036L) - ((int128)tmp_q[2] * 71243768560058L) - ((((int128)tmp_q[3] * 54399266451725L) - ((int128)tmp_q[4] * 40485352498776L) + ((int128)tmp_q[5] * 48902816151558L) - ((int128)tmp_q[6] * 55595276287113L) + ((int128)tmp_q[7] * 79614643509380L) + ((int128)tmp_q[8] * 83734412135570L) - ((int128)tmp_q[9] * 111991004117825L) + ((int128)tmp_q[10] * 101319890617906L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 101319890617906L) + ((int128)tmp_q[1] * 79948358425910L) - ((int128)tmp_q[2] * 12169016728036L) - ((int128)tmp_q[3] * 71243768560058L) - ((((int128)tmp_q[4] * 54399266451725L) - ((int128)tmp_q[5] * 40485352498776L) + ((int128)tmp_q[6] * 48902816151558L) - ((int128)tmp_q[7] * 55595276287113L) + ((int128)tmp_q[8] * 79614643509380L) + ((int128)tmp_q[9] * 83734412135570L) - ((int128)tmp_q[10] * 111991004117825L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 111991004117825L) + ((int128)tmp_q[1] * 101319890617906L) + ((int128)tmp_q[2] * 79948358425910L) - ((int128)tmp_q[3] * 12169016728036L) - ((int128)tmp_q[4] * 71243768560058L) - ((((int128)tmp_q[5] * 54399266451725L) - ((int128)tmp_q[6] * 40485352498776L) + ((int128)tmp_q[7] * 48902816151558L) - ((int128)tmp_q[8] * 55595276287113L) + ((int128)tmp_q[9] * 79614643509380L) + ((int128)tmp_q[10] * 83734412135570L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 83734412135570L) - ((int128)tmp_q[1] * 111991004117825L) + ((int128)tmp_q[2] * 101319890617906L) + ((int128)tmp_q[3] * 79948358425910L) - ((int128)tmp_q[4] * 12169016728036L) - ((int128)tmp_q[5] * 71243768560058L) - ((((int128)tmp_q[6] * 54399266451725L) - ((int128)tmp_q[7] * 40485352498776L) + ((int128)tmp_q[8] * 48902816151558L) - ((int128)tmp_q[9] * 55595276287113L) + ((int128)tmp_q[10] * 79614643509380L)) * 5);
	tmp_zero[6] = ((int128)tmp_q[0] * 79614643509380L) + ((int128)tmp_q[1] * 83734412135570L) - ((int128)tmp_q[2] * 111991004117825L) + ((int128)tmp_q[3] * 101319890617906L) + ((int128)tmp_q[4] * 79948358425910L) - ((int128)tmp_q[5] * 12169016728036L) - ((int128)tmp_q[6] * 71243768560058L) - ((((int128)tmp_q[7] * 54399266451725L) - ((int128)tmp_q[8] * 40485352498776L) + ((int128)tmp_q[9] * 48902816151558L) - ((int128)tmp_q[10] * 55595276287113L)) * 5);
	tmp_zero[7] = -((int128)tmp_q[0] * 55595276287113L) + ((int128)tmp_q[1] * 79614643509380L) + ((int128)tmp_q[2] * 83734412135570L) - ((int128)tmp_q[3] * 111991004117825L) + ((int128)tmp_q[4] * 101319890617906L) + ((int128)tmp_q[5] * 79948358425910L) - ((int128)tmp_q[6] * 12169016728036L) - ((int128)tmp_q[7] * 71243768560058L) - ((((int128)tmp_q[8] * 54399266451725L) - ((int128)tmp_q[9] * 40485352498776L) + ((int128)tmp_q[10] * 48902816151558L)) * 5);
	tmp_zero[8] = ((int128)tmp_q[0] * 48902816151558L) - ((int128)tmp_q[1] * 55595276287113L) + ((int128)tmp_q[2] * 79614643509380L) + ((int128)tmp_q[3] * 83734412135570L) - ((int128)tmp_q[4] * 111991004117825L) + ((int128)tmp_q[5] * 101319890617906L) + ((int128)tmp_q[6] * 79948358425910L) - ((int128)tmp_q[7] * 12169016728036L) - ((int128)tmp_q[8] * 71243768560058L) - ((((int128)tmp_q[9] * 54399266451725L) - ((int128)tmp_q[10] * 40485352498776L)) * 5);
	tmp_zero[9] = -((int128)tmp_q[0] * 40485352498776L) + ((int128)tmp_q[1] * 48902816151558L) - ((int128)tmp_q[2] * 55595276287113L) + ((int128)tmp_q[3] * 79614643509380L) + ((int128)tmp_q[4] * 83734412135570L) - ((int128)tmp_q[5] * 111991004117825L) + ((int128)tmp_q[6] * 101319890617906L) + ((int128)tmp_q[7] * 79948358425910L) - ((int128)tmp_q[8] * 12169016728036L) - ((int128)tmp_q[9] * 71243768560058L) - ((int128)tmp_q[10] * 271996332258625L);
	tmp_zero[10] = ((int128)tmp_q[0] * 54399266451725L) - ((int128)tmp_q[1] * 40485352498776L) + ((int128)tmp_q[2] * 48902816151558L) - ((int128)tmp_q[3] * 55595276287113L) + ((int128)tmp_q[4] * 79614643509380L) + ((int128)tmp_q[5] * 83734412135570L) - ((int128)tmp_q[6] * 111991004117825L) + ((int128)tmp_q[7] * 101319890617906L) + ((int128)tmp_q[8] * 79948358425910L) - ((int128)tmp_q[9] * 12169016728036L) - ((int128)tmp_q[10] * 71243768560058L);

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

