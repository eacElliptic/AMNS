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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16270104243975068019UL) + ((((uint64_t)op[1] * 18204478944279425880UL) + ((uint64_t)op[2] * 12508030722684870882UL) + ((uint64_t)op[3] * 3228999485287161395UL) + ((uint64_t)op[4] * 3133102657820286843UL) + ((uint64_t)op[5] * 243183058850240770UL) + ((uint64_t)op[6] * 14849891550509593550UL) + ((uint64_t)op[7] * 5779617769572003234UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 5779617769572003234UL) + ((uint64_t)op[1] * 16270104243975068019UL) + ((((uint64_t)op[2] * 18204478944279425880UL) + ((uint64_t)op[3] * 12508030722684870882UL) + ((uint64_t)op[4] * 3228999485287161395UL) + ((uint64_t)op[5] * 3133102657820286843UL) + ((uint64_t)op[6] * 243183058850240770UL) + ((uint64_t)op[7] * 14849891550509593550UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 14849891550509593550UL) + ((uint64_t)op[1] * 5779617769572003234UL) + ((uint64_t)op[2] * 16270104243975068019UL) + ((((uint64_t)op[3] * 18204478944279425880UL) + ((uint64_t)op[4] * 12508030722684870882UL) + ((uint64_t)op[5] * 3228999485287161395UL) + ((uint64_t)op[6] * 3133102657820286843UL) + ((uint64_t)op[7] * 243183058850240770UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 243183058850240770UL) + ((uint64_t)op[1] * 14849891550509593550UL) + ((uint64_t)op[2] * 5779617769572003234UL) + ((uint64_t)op[3] * 16270104243975068019UL) + ((((uint64_t)op[4] * 18204478944279425880UL) + ((uint64_t)op[5] * 12508030722684870882UL) + ((uint64_t)op[6] * 3228999485287161395UL) + ((uint64_t)op[7] * 3133102657820286843UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 3133102657820286843UL) + ((uint64_t)op[1] * 243183058850240770UL) + ((uint64_t)op[2] * 14849891550509593550UL) + ((uint64_t)op[3] * 5779617769572003234UL) + ((uint64_t)op[4] * 16270104243975068019UL) + ((((uint64_t)op[5] * 18204478944279425880UL) + ((uint64_t)op[6] * 12508030722684870882UL) + ((uint64_t)op[7] * 3228999485287161395UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 3228999485287161395UL) + ((uint64_t)op[1] * 3133102657820286843UL) + ((uint64_t)op[2] * 243183058850240770UL) + ((uint64_t)op[3] * 14849891550509593550UL) + ((uint64_t)op[4] * 5779617769572003234UL) + ((uint64_t)op[5] * 16270104243975068019UL) + ((((uint64_t)op[6] * 18204478944279425880UL) + ((uint64_t)op[7] * 12508030722684870882UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 12508030722684870882UL) + ((uint64_t)op[1] * 3228999485287161395UL) + ((uint64_t)op[2] * 3133102657820286843UL) + ((uint64_t)op[3] * 243183058850240770UL) + ((uint64_t)op[4] * 14849891550509593550UL) + ((uint64_t)op[5] * 5779617769572003234UL) + ((uint64_t)op[6] * 16270104243975068019UL) + ((uint64_t)op[7] * 726795388290377208UL);
	tmp_q[7] = ((uint64_t)op[0] * 18204478944279425880UL) + ((uint64_t)op[1] * 12508030722684870882UL) + ((uint64_t)op[2] * 3228999485287161395UL) + ((uint64_t)op[3] * 3133102657820286843UL) + ((uint64_t)op[4] * 243183058850240770UL) + ((uint64_t)op[5] * 14849891550509593550UL) + ((uint64_t)op[6] * 5779617769572003234UL) + ((uint64_t)op[7] * 16270104243975068019UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 35907038194082L) - ((((int128)tmp_q[1] * 60084633874666L) - ((int128)tmp_q[2] * 84479712586569L) - ((int128)tmp_q[3] * 109601073587410L) + ((int128)tmp_q[4] * 56169793413230L) - ((int128)tmp_q[5] * 32696317788839L) - ((int128)tmp_q[6] * 72827919233037L) + ((int128)tmp_q[7] * 68973967042796L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 68973967042796L) + ((int128)tmp_q[1] * 35907038194082L) - ((((int128)tmp_q[2] * 60084633874666L) - ((int128)tmp_q[3] * 84479712586569L) - ((int128)tmp_q[4] * 109601073587410L) + ((int128)tmp_q[5] * 56169793413230L) - ((int128)tmp_q[6] * 32696317788839L) - ((int128)tmp_q[7] * 72827919233037L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 72827919233037L) + ((int128)tmp_q[1] * 68973967042796L) + ((int128)tmp_q[2] * 35907038194082L) - ((((int128)tmp_q[3] * 60084633874666L) - ((int128)tmp_q[4] * 84479712586569L) - ((int128)tmp_q[5] * 109601073587410L) + ((int128)tmp_q[6] * 56169793413230L) - ((int128)tmp_q[7] * 32696317788839L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 32696317788839L) - ((int128)tmp_q[1] * 72827919233037L) + ((int128)tmp_q[2] * 68973967042796L) + ((int128)tmp_q[3] * 35907038194082L) - ((((int128)tmp_q[4] * 60084633874666L) - ((int128)tmp_q[5] * 84479712586569L) - ((int128)tmp_q[6] * 109601073587410L) + ((int128)tmp_q[7] * 56169793413230L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 56169793413230L) - ((int128)tmp_q[1] * 32696317788839L) - ((int128)tmp_q[2] * 72827919233037L) + ((int128)tmp_q[3] * 68973967042796L) + ((int128)tmp_q[4] * 35907038194082L) - ((((int128)tmp_q[5] * 60084633874666L) - ((int128)tmp_q[6] * 84479712586569L) - ((int128)tmp_q[7] * 109601073587410L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 109601073587410L) + ((int128)tmp_q[1] * 56169793413230L) - ((int128)tmp_q[2] * 32696317788839L) - ((int128)tmp_q[3] * 72827919233037L) + ((int128)tmp_q[4] * 68973967042796L) + ((int128)tmp_q[5] * 35907038194082L) - ((((int128)tmp_q[6] * 60084633874666L) - ((int128)tmp_q[7] * 84479712586569L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 84479712586569L) - ((int128)tmp_q[1] * 109601073587410L) + ((int128)tmp_q[2] * 56169793413230L) - ((int128)tmp_q[3] * 32696317788839L) - ((int128)tmp_q[4] * 72827919233037L) + ((int128)tmp_q[5] * 68973967042796L) + ((int128)tmp_q[6] * 35907038194082L) - ((int128)tmp_q[7] * 180253901623998L);
	tmp_zero[7] = ((int128)tmp_q[0] * 60084633874666L) - ((int128)tmp_q[1] * 84479712586569L) - ((int128)tmp_q[2] * 109601073587410L) + ((int128)tmp_q[3] * 56169793413230L) - ((int128)tmp_q[4] * 32696317788839L) - ((int128)tmp_q[5] * 72827919233037L) + ((int128)tmp_q[6] * 68973967042796L) + ((int128)tmp_q[7] * 35907038194082L);

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

