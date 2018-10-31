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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 15736118240690741490UL) + ((((uint64_t)op[1] * 3557656452297574330UL) + ((uint64_t)op[2] * 12931517456405395594UL) + ((uint64_t)op[3] * 15841945835515786058UL) + ((uint64_t)op[4] * 16379684890418016058UL) + ((uint64_t)op[5] * 11846369034950018598UL) + ((uint64_t)op[6] * 1613212953370117717UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 1613212953370117717UL) + ((uint64_t)op[1] * 15736118240690741490UL) + ((((uint64_t)op[2] * 3557656452297574330UL) + ((uint64_t)op[3] * 12931517456405395594UL) + ((uint64_t)op[4] * 15841945835515786058UL) + ((uint64_t)op[5] * 16379684890418016058UL) + ((uint64_t)op[6] * 11846369034950018598UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 11846369034950018598UL) + ((uint64_t)op[1] * 1613212953370117717UL) + ((uint64_t)op[2] * 15736118240690741490UL) + ((((uint64_t)op[3] * 3557656452297574330UL) + ((uint64_t)op[4] * 12931517456405395594UL) + ((uint64_t)op[5] * 15841945835515786058UL) + ((uint64_t)op[6] * 16379684890418016058UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 16379684890418016058UL) + ((uint64_t)op[1] * 11846369034950018598UL) + ((uint64_t)op[2] * 1613212953370117717UL) + ((uint64_t)op[3] * 15736118240690741490UL) + ((((uint64_t)op[4] * 3557656452297574330UL) + ((uint64_t)op[5] * 12931517456405395594UL) + ((uint64_t)op[6] * 15841945835515786058UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 15841945835515786058UL) + ((uint64_t)op[1] * 16379684890418016058UL) + ((uint64_t)op[2] * 11846369034950018598UL) + ((uint64_t)op[3] * 1613212953370117717UL) + ((uint64_t)op[4] * 15736118240690741490UL) + ((((uint64_t)op[5] * 3557656452297574330UL) + ((uint64_t)op[6] * 12931517456405395594UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 12931517456405395594UL) + ((uint64_t)op[1] * 15841945835515786058UL) + ((uint64_t)op[2] * 16379684890418016058UL) + ((uint64_t)op[3] * 11846369034950018598UL) + ((uint64_t)op[4] * 1613212953370117717UL) + ((uint64_t)op[5] * 15736118240690741490UL) + ((uint64_t)op[6] * 7773774716816828626UL);
	tmp_q[6] = ((uint64_t)op[0] * 3557656452297574330UL) + ((uint64_t)op[1] * 12931517456405395594UL) + ((uint64_t)op[2] * 15841945835515786058UL) + ((uint64_t)op[3] * 16379684890418016058UL) + ((uint64_t)op[4] * 11846369034950018598UL) + ((uint64_t)op[5] * 1613212953370117717UL) + ((uint64_t)op[6] * 15736118240690741490UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 14500005994L) - ((((int128)tmp_q[1] * 96891110135L) + ((int128)tmp_q[2] * 7794027566L) - ((int128)tmp_q[3] * 14373323042L) - ((int128)tmp_q[4] * 31390896258L) + ((int128)tmp_q[5] * 56556193358L) - ((int128)tmp_q[6] * 88983529698L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 88983529698L) + ((int128)tmp_q[1] * 14500005994L) - ((((int128)tmp_q[2] * 96891110135L) + ((int128)tmp_q[3] * 7794027566L) - ((int128)tmp_q[4] * 14373323042L) - ((int128)tmp_q[5] * 31390896258L) + ((int128)tmp_q[6] * 56556193358L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 56556193358L) - ((int128)tmp_q[1] * 88983529698L) + ((int128)tmp_q[2] * 14500005994L) - ((((int128)tmp_q[3] * 96891110135L) + ((int128)tmp_q[4] * 7794027566L) - ((int128)tmp_q[5] * 14373323042L) - ((int128)tmp_q[6] * 31390896258L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 31390896258L) + ((int128)tmp_q[1] * 56556193358L) - ((int128)tmp_q[2] * 88983529698L) + ((int128)tmp_q[3] * 14500005994L) - ((((int128)tmp_q[4] * 96891110135L) + ((int128)tmp_q[5] * 7794027566L) - ((int128)tmp_q[6] * 14373323042L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 14373323042L) - ((int128)tmp_q[1] * 31390896258L) + ((int128)tmp_q[2] * 56556193358L) - ((int128)tmp_q[3] * 88983529698L) + ((int128)tmp_q[4] * 14500005994L) - ((((int128)tmp_q[5] * 96891110135L) + ((int128)tmp_q[6] * 7794027566L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 7794027566L) - ((int128)tmp_q[1] * 14373323042L) - ((int128)tmp_q[2] * 31390896258L) + ((int128)tmp_q[3] * 56556193358L) - ((int128)tmp_q[4] * 88983529698L) + ((int128)tmp_q[5] * 14500005994L) - ((int128)tmp_q[6] * 290673330405L);
	tmp_zero[6] = ((int128)tmp_q[0] * 96891110135L) + ((int128)tmp_q[1] * 7794027566L) - ((int128)tmp_q[2] * 14373323042L) - ((int128)tmp_q[3] * 31390896258L) + ((int128)tmp_q[4] * 56556193358L) - ((int128)tmp_q[5] * 88983529698L) + ((int128)tmp_q[6] * 14500005994L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

