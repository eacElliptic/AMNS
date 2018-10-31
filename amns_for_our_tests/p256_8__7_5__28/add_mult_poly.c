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
	tmp_q[0] = ((uint64_t)op[0] * 10045988444339870179UL) + ((((uint64_t)op[1] * 4733409949146951344UL) + ((uint64_t)op[2] * 12427639793230182082UL) + ((uint64_t)op[3] * 3801435950956238112UL) + ((uint64_t)op[4] * 623449226343587624UL) + ((uint64_t)op[5] * 15643987463610631353UL) + ((uint64_t)op[6] * 11862026474240606541UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 11862026474240606541UL) + ((uint64_t)op[1] * 10045988444339870179UL) + ((((uint64_t)op[2] * 4733409949146951344UL) + ((uint64_t)op[3] * 12427639793230182082UL) + ((uint64_t)op[4] * 3801435950956238112UL) + ((uint64_t)op[5] * 623449226343587624UL) + ((uint64_t)op[6] * 15643987463610631353UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 15643987463610631353UL) + ((uint64_t)op[1] * 11862026474240606541UL) + ((uint64_t)op[2] * 10045988444339870179UL) + ((((uint64_t)op[3] * 4733409949146951344UL) + ((uint64_t)op[4] * 12427639793230182082UL) + ((uint64_t)op[5] * 3801435950956238112UL) + ((uint64_t)op[6] * 623449226343587624UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 623449226343587624UL) + ((uint64_t)op[1] * 15643987463610631353UL) + ((uint64_t)op[2] * 11862026474240606541UL) + ((uint64_t)op[3] * 10045988444339870179UL) + ((((uint64_t)op[4] * 4733409949146951344UL) + ((uint64_t)op[5] * 12427639793230182082UL) + ((uint64_t)op[6] * 3801435950956238112UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 3801435950956238112UL) + ((uint64_t)op[1] * 623449226343587624UL) + ((uint64_t)op[2] * 15643987463610631353UL) + ((uint64_t)op[3] * 11862026474240606541UL) + ((uint64_t)op[4] * 10045988444339870179UL) + ((((uint64_t)op[5] * 4733409949146951344UL) + ((uint64_t)op[6] * 12427639793230182082UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 12427639793230182082UL) + ((uint64_t)op[1] * 3801435950956238112UL) + ((uint64_t)op[2] * 623449226343587624UL) + ((uint64_t)op[3] * 15643987463610631353UL) + ((uint64_t)op[4] * 11862026474240606541UL) + ((uint64_t)op[5] * 10045988444339870179UL) + ((uint64_t)op[6] * 5220305672025205104UL);
	tmp_q[6] = ((uint64_t)op[0] * 4733409949146951344UL) + ((uint64_t)op[1] * 12427639793230182082UL) + ((uint64_t)op[2] * 3801435950956238112UL) + ((uint64_t)op[3] * 623449226343587624UL) + ((uint64_t)op[4] * 15643987463610631353UL) + ((uint64_t)op[5] * 11862026474240606541UL) + ((uint64_t)op[6] * 10045988444339870179UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 5498856261L) + ((((int128)tmp_q[1] * 9012287211L) + ((int128)tmp_q[2] * 41873573401L) + ((int128)tmp_q[3] * 28860813282L) + ((int128)tmp_q[4] * 17304211653L) - ((int128)tmp_q[5] * 54502955085L) - ((int128)tmp_q[6] * 15380165368L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 15380165368L) - ((int128)tmp_q[1] * 5498856261L) + ((((int128)tmp_q[2] * 9012287211L) + ((int128)tmp_q[3] * 41873573401L) + ((int128)tmp_q[4] * 28860813282L) + ((int128)tmp_q[5] * 17304211653L) - ((int128)tmp_q[6] * 54502955085L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 54502955085L) - ((int128)tmp_q[1] * 15380165368L) - ((int128)tmp_q[2] * 5498856261L) + ((((int128)tmp_q[3] * 9012287211L) + ((int128)tmp_q[4] * 41873573401L) + ((int128)tmp_q[5] * 28860813282L) + ((int128)tmp_q[6] * 17304211653L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 17304211653L) - ((int128)tmp_q[1] * 54502955085L) - ((int128)tmp_q[2] * 15380165368L) - ((int128)tmp_q[3] * 5498856261L) + ((((int128)tmp_q[4] * 9012287211L) + ((int128)tmp_q[5] * 41873573401L) + ((int128)tmp_q[6] * 28860813282L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 28860813282L) + ((int128)tmp_q[1] * 17304211653L) - ((int128)tmp_q[2] * 54502955085L) - ((int128)tmp_q[3] * 15380165368L) - ((int128)tmp_q[4] * 5498856261L) + ((((int128)tmp_q[5] * 9012287211L) + ((int128)tmp_q[6] * 41873573401L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 41873573401L) + ((int128)tmp_q[1] * 28860813282L) + ((int128)tmp_q[2] * 17304211653L) - ((int128)tmp_q[3] * 54502955085L) - ((int128)tmp_q[4] * 15380165368L) - ((int128)tmp_q[5] * 5498856261L) + ((int128)tmp_q[6] * 45061436055L);
	tmp_zero[6] = ((int128)tmp_q[0] * 9012287211L) + ((int128)tmp_q[1] * 41873573401L) + ((int128)tmp_q[2] * 28860813282L) + ((int128)tmp_q[3] * 17304211653L) - ((int128)tmp_q[4] * 54502955085L) - ((int128)tmp_q[5] * 15380165368L) - ((int128)tmp_q[6] * 5498856261L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

