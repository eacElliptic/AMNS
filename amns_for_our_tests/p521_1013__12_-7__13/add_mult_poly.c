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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 7);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 7);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 7);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 7);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 14);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 7);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 14);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 7);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 14);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 7);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10156986006751753143UL) + ((((uint64_t)op[1] * 5462231390985352355UL) + ((uint64_t)op[2] * 10163010691078063721UL) + ((uint64_t)op[3] * 10296057139431234495UL) + ((uint64_t)op[4] * 13778009854058663889UL) + ((uint64_t)op[5] * 7409850377553130839UL) + ((uint64_t)op[6] * 16481807680229777755UL) + ((uint64_t)op[7] * 11768954815069946667UL) + ((uint64_t)op[8] * 8116730732989826917UL) + ((uint64_t)op[9] * 11989276146155319986UL) + ((uint64_t)op[10] * 962821716142879593UL) + ((uint64_t)op[11] * 2469182022666045473UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 2469182022666045473UL) + ((uint64_t)op[1] * 10156986006751753143UL) + ((((uint64_t)op[2] * 5462231390985352355UL) + ((uint64_t)op[3] * 10163010691078063721UL) + ((uint64_t)op[4] * 10296057139431234495UL) + ((uint64_t)op[5] * 13778009854058663889UL) + ((uint64_t)op[6] * 7409850377553130839UL) + ((uint64_t)op[7] * 16481807680229777755UL) + ((uint64_t)op[8] * 11768954815069946667UL) + ((uint64_t)op[9] * 8116730732989826917UL) + ((uint64_t)op[10] * 11989276146155319986UL) + ((uint64_t)op[11] * 962821716142879593UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 962821716142879593UL) + ((uint64_t)op[1] * 2469182022666045473UL) + ((uint64_t)op[2] * 10156986006751753143UL) + ((((uint64_t)op[3] * 5462231390985352355UL) + ((uint64_t)op[4] * 10163010691078063721UL) + ((uint64_t)op[5] * 10296057139431234495UL) + ((uint64_t)op[6] * 13778009854058663889UL) + ((uint64_t)op[7] * 7409850377553130839UL) + ((uint64_t)op[8] * 16481807680229777755UL) + ((uint64_t)op[9] * 11768954815069946667UL) + ((uint64_t)op[10] * 8116730732989826917UL) + ((uint64_t)op[11] * 11989276146155319986UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 11989276146155319986UL) + ((uint64_t)op[1] * 962821716142879593UL) + ((uint64_t)op[2] * 2469182022666045473UL) + ((uint64_t)op[3] * 10156986006751753143UL) + ((((uint64_t)op[4] * 5462231390985352355UL) + ((uint64_t)op[5] * 10163010691078063721UL) + ((uint64_t)op[6] * 10296057139431234495UL) + ((uint64_t)op[7] * 13778009854058663889UL) + ((uint64_t)op[8] * 7409850377553130839UL) + ((uint64_t)op[9] * 16481807680229777755UL) + ((uint64_t)op[10] * 11768954815069946667UL) + ((uint64_t)op[11] * 8116730732989826917UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 8116730732989826917UL) + ((uint64_t)op[1] * 11989276146155319986UL) + ((uint64_t)op[2] * 962821716142879593UL) + ((uint64_t)op[3] * 2469182022666045473UL) + ((uint64_t)op[4] * 10156986006751753143UL) + ((((uint64_t)op[5] * 5462231390985352355UL) + ((uint64_t)op[6] * 10163010691078063721UL) + ((uint64_t)op[7] * 10296057139431234495UL) + ((uint64_t)op[8] * 13778009854058663889UL) + ((uint64_t)op[9] * 7409850377553130839UL) + ((uint64_t)op[10] * 16481807680229777755UL) + ((uint64_t)op[11] * 11768954815069946667UL)) * 18446744073709551609);
	tmp_q[5] = ((uint64_t)op[0] * 11768954815069946667UL) + ((uint64_t)op[1] * 8116730732989826917UL) + ((uint64_t)op[2] * 11989276146155319986UL) + ((uint64_t)op[3] * 962821716142879593UL) + ((uint64_t)op[4] * 2469182022666045473UL) + ((uint64_t)op[5] * 10156986006751753143UL) + ((((uint64_t)op[6] * 5462231390985352355UL) + ((uint64_t)op[7] * 10163010691078063721UL) + ((uint64_t)op[8] * 10296057139431234495UL) + ((uint64_t)op[9] * 13778009854058663889UL) + ((uint64_t)op[10] * 7409850377553130839UL) + ((uint64_t)op[11] * 16481807680229777755UL)) * 18446744073709551609);
	tmp_q[6] = ((uint64_t)op[0] * 16481807680229777755UL) + ((uint64_t)op[1] * 11768954815069946667UL) + ((uint64_t)op[2] * 8116730732989826917UL) + ((uint64_t)op[3] * 11989276146155319986UL) + ((uint64_t)op[4] * 962821716142879593UL) + ((uint64_t)op[5] * 2469182022666045473UL) + ((uint64_t)op[6] * 10156986006751753143UL) + ((((uint64_t)op[7] * 5462231390985352355UL) + ((uint64_t)op[8] * 10163010691078063721UL) + ((uint64_t)op[9] * 10296057139431234495UL) + ((uint64_t)op[10] * 13778009854058663889UL) + ((uint64_t)op[11] * 7409850377553130839UL)) * 18446744073709551609);
	tmp_q[7] = ((uint64_t)op[0] * 7409850377553130839UL) + ((uint64_t)op[1] * 16481807680229777755UL) + ((uint64_t)op[2] * 11768954815069946667UL) + ((uint64_t)op[3] * 8116730732989826917UL) + ((uint64_t)op[4] * 11989276146155319986UL) + ((uint64_t)op[5] * 962821716142879593UL) + ((uint64_t)op[6] * 2469182022666045473UL) + ((uint64_t)op[7] * 10156986006751753143UL) + ((((uint64_t)op[8] * 5462231390985352355UL) + ((uint64_t)op[9] * 10163010691078063721UL) + ((uint64_t)op[10] * 10296057139431234495UL) + ((uint64_t)op[11] * 13778009854058663889UL)) * 18446744073709551609);
	tmp_q[8] = ((uint64_t)op[0] * 13778009854058663889UL) + ((uint64_t)op[1] * 7409850377553130839UL) + ((uint64_t)op[2] * 16481807680229777755UL) + ((uint64_t)op[3] * 11768954815069946667UL) + ((uint64_t)op[4] * 8116730732989826917UL) + ((uint64_t)op[5] * 11989276146155319986UL) + ((uint64_t)op[6] * 962821716142879593UL) + ((uint64_t)op[7] * 2469182022666045473UL) + ((uint64_t)op[8] * 10156986006751753143UL) + ((((uint64_t)op[9] * 5462231390985352355UL) + ((uint64_t)op[10] * 10163010691078063721UL) + ((uint64_t)op[11] * 10296057139431234495UL)) * 18446744073709551609);
	tmp_q[9] = ((uint64_t)op[0] * 10296057139431234495UL) + ((uint64_t)op[1] * 13778009854058663889UL) + ((uint64_t)op[2] * 7409850377553130839UL) + ((uint64_t)op[3] * 16481807680229777755UL) + ((uint64_t)op[4] * 11768954815069946667UL) + ((uint64_t)op[5] * 8116730732989826917UL) + ((uint64_t)op[6] * 11989276146155319986UL) + ((uint64_t)op[7] * 962821716142879593UL) + ((uint64_t)op[8] * 2469182022666045473UL) + ((uint64_t)op[9] * 10156986006751753143UL) + ((((uint64_t)op[10] * 5462231390985352355UL) + ((uint64_t)op[11] * 10163010691078063721UL)) * 18446744073709551609);
	tmp_q[10] = ((uint64_t)op[0] * 10163010691078063721UL) + ((uint64_t)op[1] * 10296057139431234495UL) + ((uint64_t)op[2] * 13778009854058663889UL) + ((uint64_t)op[3] * 7409850377553130839UL) + ((uint64_t)op[4] * 16481807680229777755UL) + ((uint64_t)op[5] * 11768954815069946667UL) + ((uint64_t)op[6] * 8116730732989826917UL) + ((uint64_t)op[7] * 11989276146155319986UL) + ((uint64_t)op[8] * 962821716142879593UL) + ((uint64_t)op[9] * 2469182022666045473UL) + ((uint64_t)op[10] * 10156986006751753143UL) + ((uint64_t)op[11] * 17104612484231188363UL);
	tmp_q[11] = ((uint64_t)op[0] * 5462231390985352355UL) + ((uint64_t)op[1] * 10163010691078063721UL) + ((uint64_t)op[2] * 10296057139431234495UL) + ((uint64_t)op[3] * 13778009854058663889UL) + ((uint64_t)op[4] * 7409850377553130839UL) + ((uint64_t)op[5] * 16481807680229777755UL) + ((uint64_t)op[6] * 11768954815069946667UL) + ((uint64_t)op[7] * 8116730732989826917UL) + ((uint64_t)op[8] * 11989276146155319986UL) + ((uint64_t)op[9] * 962821716142879593UL) + ((uint64_t)op[10] * 2469182022666045473UL) + ((uint64_t)op[11] * 10156986006751753143UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2470715929961L) - ((((int128)tmp_q[1] * 2819756375183L) + ((int128)tmp_q[2] * 4412459660305L) + ((int128)tmp_q[3] * 1699379088536L) + ((int128)tmp_q[4] * 7928958528837L) - ((int128)tmp_q[5] * 4726134229563L) - ((int128)tmp_q[6] * 143675493173L) + ((int128)tmp_q[7] * 4708206304131L) + ((int128)tmp_q[8] * 1331403026253L) + ((int128)tmp_q[9] * 2759590246267L) - ((int128)tmp_q[10] * 1546922990891L) - ((int128)tmp_q[11] * 6968458009853L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 6968458009853L) - ((int128)tmp_q[1] * 2470715929961L) - ((((int128)tmp_q[2] * 2819756375183L) + ((int128)tmp_q[3] * 4412459660305L) + ((int128)tmp_q[4] * 1699379088536L) + ((int128)tmp_q[5] * 7928958528837L) - ((int128)tmp_q[6] * 4726134229563L) - ((int128)tmp_q[7] * 143675493173L) + ((int128)tmp_q[8] * 4708206304131L) + ((int128)tmp_q[9] * 1331403026253L) + ((int128)tmp_q[10] * 2759590246267L) - ((int128)tmp_q[11] * 1546922990891L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 1546922990891L) - ((int128)tmp_q[1] * 6968458009853L) - ((int128)tmp_q[2] * 2470715929961L) - ((((int128)tmp_q[3] * 2819756375183L) + ((int128)tmp_q[4] * 4412459660305L) + ((int128)tmp_q[5] * 1699379088536L) + ((int128)tmp_q[6] * 7928958528837L) - ((int128)tmp_q[7] * 4726134229563L) - ((int128)tmp_q[8] * 143675493173L) + ((int128)tmp_q[9] * 4708206304131L) + ((int128)tmp_q[10] * 1331403026253L) + ((int128)tmp_q[11] * 2759590246267L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 2759590246267L) - ((int128)tmp_q[1] * 1546922990891L) - ((int128)tmp_q[2] * 6968458009853L) - ((int128)tmp_q[3] * 2470715929961L) - ((((int128)tmp_q[4] * 2819756375183L) + ((int128)tmp_q[5] * 4412459660305L) + ((int128)tmp_q[6] * 1699379088536L) + ((int128)tmp_q[7] * 7928958528837L) - ((int128)tmp_q[8] * 4726134229563L) - ((int128)tmp_q[9] * 143675493173L) + ((int128)tmp_q[10] * 4708206304131L) + ((int128)tmp_q[11] * 1331403026253L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 1331403026253L) + ((int128)tmp_q[1] * 2759590246267L) - ((int128)tmp_q[2] * 1546922990891L) - ((int128)tmp_q[3] * 6968458009853L) - ((int128)tmp_q[4] * 2470715929961L) - ((((int128)tmp_q[5] * 2819756375183L) + ((int128)tmp_q[6] * 4412459660305L) + ((int128)tmp_q[7] * 1699379088536L) + ((int128)tmp_q[8] * 7928958528837L) - ((int128)tmp_q[9] * 4726134229563L) - ((int128)tmp_q[10] * 143675493173L) + ((int128)tmp_q[11] * 4708206304131L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 4708206304131L) + ((int128)tmp_q[1] * 1331403026253L) + ((int128)tmp_q[2] * 2759590246267L) - ((int128)tmp_q[3] * 1546922990891L) - ((int128)tmp_q[4] * 6968458009853L) - ((int128)tmp_q[5] * 2470715929961L) - ((((int128)tmp_q[6] * 2819756375183L) + ((int128)tmp_q[7] * 4412459660305L) + ((int128)tmp_q[8] * 1699379088536L) + ((int128)tmp_q[9] * 7928958528837L) - ((int128)tmp_q[10] * 4726134229563L) - ((int128)tmp_q[11] * 143675493173L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 143675493173L) + ((int128)tmp_q[1] * 4708206304131L) + ((int128)tmp_q[2] * 1331403026253L) + ((int128)tmp_q[3] * 2759590246267L) - ((int128)tmp_q[4] * 1546922990891L) - ((int128)tmp_q[5] * 6968458009853L) - ((int128)tmp_q[6] * 2470715929961L) - ((((int128)tmp_q[7] * 2819756375183L) + ((int128)tmp_q[8] * 4412459660305L) + ((int128)tmp_q[9] * 1699379088536L) + ((int128)tmp_q[10] * 7928958528837L) - ((int128)tmp_q[11] * 4726134229563L)) * 7);
	tmp_zero[7] = -((int128)tmp_q[0] * 4726134229563L) - ((int128)tmp_q[1] * 143675493173L) + ((int128)tmp_q[2] * 4708206304131L) + ((int128)tmp_q[3] * 1331403026253L) + ((int128)tmp_q[4] * 2759590246267L) - ((int128)tmp_q[5] * 1546922990891L) - ((int128)tmp_q[6] * 6968458009853L) - ((int128)tmp_q[7] * 2470715929961L) - ((((int128)tmp_q[8] * 2819756375183L) + ((int128)tmp_q[9] * 4412459660305L) + ((int128)tmp_q[10] * 1699379088536L) + ((int128)tmp_q[11] * 7928958528837L)) * 7);
	tmp_zero[8] = ((int128)tmp_q[0] * 7928958528837L) - ((int128)tmp_q[1] * 4726134229563L) - ((int128)tmp_q[2] * 143675493173L) + ((int128)tmp_q[3] * 4708206304131L) + ((int128)tmp_q[4] * 1331403026253L) + ((int128)tmp_q[5] * 2759590246267L) - ((int128)tmp_q[6] * 1546922990891L) - ((int128)tmp_q[7] * 6968458009853L) - ((int128)tmp_q[8] * 2470715929961L) - ((((int128)tmp_q[9] * 2819756375183L) + ((int128)tmp_q[10] * 4412459660305L) + ((int128)tmp_q[11] * 1699379088536L)) * 7);
	tmp_zero[9] = ((int128)tmp_q[0] * 1699379088536L) + ((int128)tmp_q[1] * 7928958528837L) - ((int128)tmp_q[2] * 4726134229563L) - ((int128)tmp_q[3] * 143675493173L) + ((int128)tmp_q[4] * 4708206304131L) + ((int128)tmp_q[5] * 1331403026253L) + ((int128)tmp_q[6] * 2759590246267L) - ((int128)tmp_q[7] * 1546922990891L) - ((int128)tmp_q[8] * 6968458009853L) - ((int128)tmp_q[9] * 2470715929961L) - ((((int128)tmp_q[10] * 2819756375183L) + ((int128)tmp_q[11] * 4412459660305L)) * 7);
	tmp_zero[10] = ((int128)tmp_q[0] * 4412459660305L) + ((int128)tmp_q[1] * 1699379088536L) + ((int128)tmp_q[2] * 7928958528837L) - ((int128)tmp_q[3] * 4726134229563L) - ((int128)tmp_q[4] * 143675493173L) + ((int128)tmp_q[5] * 4708206304131L) + ((int128)tmp_q[6] * 1331403026253L) + ((int128)tmp_q[7] * 2759590246267L) - ((int128)tmp_q[8] * 1546922990891L) - ((int128)tmp_q[9] * 6968458009853L) - ((int128)tmp_q[10] * 2470715929961L) - ((int128)tmp_q[11] * 19738294626281L);
	tmp_zero[11] = ((int128)tmp_q[0] * 2819756375183L) + ((int128)tmp_q[1] * 4412459660305L) + ((int128)tmp_q[2] * 1699379088536L) + ((int128)tmp_q[3] * 7928958528837L) - ((int128)tmp_q[4] * 4726134229563L) - ((int128)tmp_q[5] * 143675493173L) + ((int128)tmp_q[6] * 4708206304131L) + ((int128)tmp_q[7] * 1331403026253L) + ((int128)tmp_q[8] * 2759590246267L) - ((int128)tmp_q[9] * 1546922990891L) - ((int128)tmp_q[10] * 6968458009853L) - ((int128)tmp_q[11] * 2470715929961L);

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

