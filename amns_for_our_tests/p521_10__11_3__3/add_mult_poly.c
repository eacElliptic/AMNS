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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 6);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12005199344540808893UL) + ((((uint64_t)op[1] * 4140514939976284599UL) + ((uint64_t)op[2] * 13711871616492941384UL) + ((uint64_t)op[3] * 6938575291076385984UL) + ((uint64_t)op[4] * 6631129054742801165UL) + ((uint64_t)op[5] * 8830362728195467902UL) + ((uint64_t)op[6] * 11939405330549907300UL) + ((uint64_t)op[7] * 2349739078527436202UL) + ((uint64_t)op[8] * 635755381247175898UL) + ((uint64_t)op[9] * 14908163069278495020UL) + ((uint64_t)op[10] * 5600022363019489930UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 5600022363019489930UL) + ((uint64_t)op[1] * 12005199344540808893UL) + ((((uint64_t)op[2] * 4140514939976284599UL) + ((uint64_t)op[3] * 13711871616492941384UL) + ((uint64_t)op[4] * 6938575291076385984UL) + ((uint64_t)op[5] * 6631129054742801165UL) + ((uint64_t)op[6] * 8830362728195467902UL) + ((uint64_t)op[7] * 11939405330549907300UL) + ((uint64_t)op[8] * 2349739078527436202UL) + ((uint64_t)op[9] * 635755381247175898UL) + ((uint64_t)op[10] * 14908163069278495020UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 14908163069278495020UL) + ((uint64_t)op[1] * 5600022363019489930UL) + ((uint64_t)op[2] * 12005199344540808893UL) + ((((uint64_t)op[3] * 4140514939976284599UL) + ((uint64_t)op[4] * 13711871616492941384UL) + ((uint64_t)op[5] * 6938575291076385984UL) + ((uint64_t)op[6] * 6631129054742801165UL) + ((uint64_t)op[7] * 8830362728195467902UL) + ((uint64_t)op[8] * 11939405330549907300UL) + ((uint64_t)op[9] * 2349739078527436202UL) + ((uint64_t)op[10] * 635755381247175898UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 635755381247175898UL) + ((uint64_t)op[1] * 14908163069278495020UL) + ((uint64_t)op[2] * 5600022363019489930UL) + ((uint64_t)op[3] * 12005199344540808893UL) + ((((uint64_t)op[4] * 4140514939976284599UL) + ((uint64_t)op[5] * 13711871616492941384UL) + ((uint64_t)op[6] * 6938575291076385984UL) + ((uint64_t)op[7] * 6631129054742801165UL) + ((uint64_t)op[8] * 8830362728195467902UL) + ((uint64_t)op[9] * 11939405330549907300UL) + ((uint64_t)op[10] * 2349739078527436202UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 2349739078527436202UL) + ((uint64_t)op[1] * 635755381247175898UL) + ((uint64_t)op[2] * 14908163069278495020UL) + ((uint64_t)op[3] * 5600022363019489930UL) + ((uint64_t)op[4] * 12005199344540808893UL) + ((((uint64_t)op[5] * 4140514939976284599UL) + ((uint64_t)op[6] * 13711871616492941384UL) + ((uint64_t)op[7] * 6938575291076385984UL) + ((uint64_t)op[8] * 6631129054742801165UL) + ((uint64_t)op[9] * 8830362728195467902UL) + ((uint64_t)op[10] * 11939405330549907300UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 11939405330549907300UL) + ((uint64_t)op[1] * 2349739078527436202UL) + ((uint64_t)op[2] * 635755381247175898UL) + ((uint64_t)op[3] * 14908163069278495020UL) + ((uint64_t)op[4] * 5600022363019489930UL) + ((uint64_t)op[5] * 12005199344540808893UL) + ((((uint64_t)op[6] * 4140514939976284599UL) + ((uint64_t)op[7] * 13711871616492941384UL) + ((uint64_t)op[8] * 6938575291076385984UL) + ((uint64_t)op[9] * 6631129054742801165UL) + ((uint64_t)op[10] * 8830362728195467902UL)) * 3);
	tmp_q[6] = ((uint64_t)op[0] * 8830362728195467902UL) + ((uint64_t)op[1] * 11939405330549907300UL) + ((uint64_t)op[2] * 2349739078527436202UL) + ((uint64_t)op[3] * 635755381247175898UL) + ((uint64_t)op[4] * 14908163069278495020UL) + ((uint64_t)op[5] * 5600022363019489930UL) + ((uint64_t)op[6] * 12005199344540808893UL) + ((((uint64_t)op[7] * 4140514939976284599UL) + ((uint64_t)op[8] * 13711871616492941384UL) + ((uint64_t)op[9] * 6938575291076385984UL) + ((uint64_t)op[10] * 6631129054742801165UL)) * 3);
	tmp_q[7] = ((uint64_t)op[0] * 6631129054742801165UL) + ((uint64_t)op[1] * 8830362728195467902UL) + ((uint64_t)op[2] * 11939405330549907300UL) + ((uint64_t)op[3] * 2349739078527436202UL) + ((uint64_t)op[4] * 635755381247175898UL) + ((uint64_t)op[5] * 14908163069278495020UL) + ((uint64_t)op[6] * 5600022363019489930UL) + ((uint64_t)op[7] * 12005199344540808893UL) + ((((uint64_t)op[8] * 4140514939976284599UL) + ((uint64_t)op[9] * 13711871616492941384UL) + ((uint64_t)op[10] * 6938575291076385984UL)) * 3);
	tmp_q[8] = ((uint64_t)op[0] * 6938575291076385984UL) + ((uint64_t)op[1] * 6631129054742801165UL) + ((uint64_t)op[2] * 8830362728195467902UL) + ((uint64_t)op[3] * 11939405330549907300UL) + ((uint64_t)op[4] * 2349739078527436202UL) + ((uint64_t)op[5] * 635755381247175898UL) + ((uint64_t)op[6] * 14908163069278495020UL) + ((uint64_t)op[7] * 5600022363019489930UL) + ((uint64_t)op[8] * 12005199344540808893UL) + ((((uint64_t)op[9] * 4140514939976284599UL) + ((uint64_t)op[10] * 13711871616492941384UL)) * 3);
	tmp_q[9] = ((uint64_t)op[0] * 13711871616492941384UL) + ((uint64_t)op[1] * 6938575291076385984UL) + ((uint64_t)op[2] * 6631129054742801165UL) + ((uint64_t)op[3] * 8830362728195467902UL) + ((uint64_t)op[4] * 11939405330549907300UL) + ((uint64_t)op[5] * 2349739078527436202UL) + ((uint64_t)op[6] * 635755381247175898UL) + ((uint64_t)op[7] * 14908163069278495020UL) + ((uint64_t)op[8] * 5600022363019489930UL) + ((uint64_t)op[9] * 12005199344540808893UL) + ((uint64_t)op[10] * 12421544819928853797UL);
	tmp_q[10] = ((uint64_t)op[0] * 4140514939976284599UL) + ((uint64_t)op[1] * 13711871616492941384UL) + ((uint64_t)op[2] * 6938575291076385984UL) + ((uint64_t)op[3] * 6631129054742801165UL) + ((uint64_t)op[4] * 8830362728195467902UL) + ((uint64_t)op[5] * 11939405330549907300UL) + ((uint64_t)op[6] * 2349739078527436202UL) + ((uint64_t)op[7] * 635755381247175898UL) + ((uint64_t)op[8] * 14908163069278495020UL) + ((uint64_t)op[9] * 5600022363019489930UL) + ((uint64_t)op[10] * 12005199344540808893UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 19022908421233L) + ((((int128)tmp_q[1] * 9704327132529L) - ((int128)tmp_q[2] * 10894263903395L) - ((int128)tmp_q[3] * 4533365462952L) - ((int128)tmp_q[4] * 54688774325005L) + ((int128)tmp_q[5] * 68393169537198L) - ((int128)tmp_q[6] * 28495547053807L) - ((int128)tmp_q[7] * 53133996878071L) - ((int128)tmp_q[8] * 78934266944114L) + ((int128)tmp_q[9] * 34153582712328L) + ((int128)tmp_q[10] * 80254157022577L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 80254157022577L) + ((int128)tmp_q[1] * 19022908421233L) + ((((int128)tmp_q[2] * 9704327132529L) - ((int128)tmp_q[3] * 10894263903395L) - ((int128)tmp_q[4] * 4533365462952L) - ((int128)tmp_q[5] * 54688774325005L) + ((int128)tmp_q[6] * 68393169537198L) - ((int128)tmp_q[7] * 28495547053807L) - ((int128)tmp_q[8] * 53133996878071L) - ((int128)tmp_q[9] * 78934266944114L) + ((int128)tmp_q[10] * 34153582712328L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 34153582712328L) + ((int128)tmp_q[1] * 80254157022577L) + ((int128)tmp_q[2] * 19022908421233L) + ((((int128)tmp_q[3] * 9704327132529L) - ((int128)tmp_q[4] * 10894263903395L) - ((int128)tmp_q[5] * 4533365462952L) - ((int128)tmp_q[6] * 54688774325005L) + ((int128)tmp_q[7] * 68393169537198L) - ((int128)tmp_q[8] * 28495547053807L) - ((int128)tmp_q[9] * 53133996878071L) - ((int128)tmp_q[10] * 78934266944114L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 78934266944114L) + ((int128)tmp_q[1] * 34153582712328L) + ((int128)tmp_q[2] * 80254157022577L) + ((int128)tmp_q[3] * 19022908421233L) + ((((int128)tmp_q[4] * 9704327132529L) - ((int128)tmp_q[5] * 10894263903395L) - ((int128)tmp_q[6] * 4533365462952L) - ((int128)tmp_q[7] * 54688774325005L) + ((int128)tmp_q[8] * 68393169537198L) - ((int128)tmp_q[9] * 28495547053807L) - ((int128)tmp_q[10] * 53133996878071L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 53133996878071L) - ((int128)tmp_q[1] * 78934266944114L) + ((int128)tmp_q[2] * 34153582712328L) + ((int128)tmp_q[3] * 80254157022577L) + ((int128)tmp_q[4] * 19022908421233L) + ((((int128)tmp_q[5] * 9704327132529L) - ((int128)tmp_q[6] * 10894263903395L) - ((int128)tmp_q[7] * 4533365462952L) - ((int128)tmp_q[8] * 54688774325005L) + ((int128)tmp_q[9] * 68393169537198L) - ((int128)tmp_q[10] * 28495547053807L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 28495547053807L) - ((int128)tmp_q[1] * 53133996878071L) - ((int128)tmp_q[2] * 78934266944114L) + ((int128)tmp_q[3] * 34153582712328L) + ((int128)tmp_q[4] * 80254157022577L) + ((int128)tmp_q[5] * 19022908421233L) + ((((int128)tmp_q[6] * 9704327132529L) - ((int128)tmp_q[7] * 10894263903395L) - ((int128)tmp_q[8] * 4533365462952L) - ((int128)tmp_q[9] * 54688774325005L) + ((int128)tmp_q[10] * 68393169537198L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 68393169537198L) - ((int128)tmp_q[1] * 28495547053807L) - ((int128)tmp_q[2] * 53133996878071L) - ((int128)tmp_q[3] * 78934266944114L) + ((int128)tmp_q[4] * 34153582712328L) + ((int128)tmp_q[5] * 80254157022577L) + ((int128)tmp_q[6] * 19022908421233L) + ((((int128)tmp_q[7] * 9704327132529L) - ((int128)tmp_q[8] * 10894263903395L) - ((int128)tmp_q[9] * 4533365462952L) - ((int128)tmp_q[10] * 54688774325005L)) * 3);
	tmp_zero[7] = -((int128)tmp_q[0] * 54688774325005L) + ((int128)tmp_q[1] * 68393169537198L) - ((int128)tmp_q[2] * 28495547053807L) - ((int128)tmp_q[3] * 53133996878071L) - ((int128)tmp_q[4] * 78934266944114L) + ((int128)tmp_q[5] * 34153582712328L) + ((int128)tmp_q[6] * 80254157022577L) + ((int128)tmp_q[7] * 19022908421233L) + ((((int128)tmp_q[8] * 9704327132529L) - ((int128)tmp_q[9] * 10894263903395L) - ((int128)tmp_q[10] * 4533365462952L)) * 3);
	tmp_zero[8] = -((int128)tmp_q[0] * 4533365462952L) - ((int128)tmp_q[1] * 54688774325005L) + ((int128)tmp_q[2] * 68393169537198L) - ((int128)tmp_q[3] * 28495547053807L) - ((int128)tmp_q[4] * 53133996878071L) - ((int128)tmp_q[5] * 78934266944114L) + ((int128)tmp_q[6] * 34153582712328L) + ((int128)tmp_q[7] * 80254157022577L) + ((int128)tmp_q[8] * 19022908421233L) + ((((int128)tmp_q[9] * 9704327132529L) - ((int128)tmp_q[10] * 10894263903395L)) * 3);
	tmp_zero[9] = -((int128)tmp_q[0] * 10894263903395L) - ((int128)tmp_q[1] * 4533365462952L) - ((int128)tmp_q[2] * 54688774325005L) + ((int128)tmp_q[3] * 68393169537198L) - ((int128)tmp_q[4] * 28495547053807L) - ((int128)tmp_q[5] * 53133996878071L) - ((int128)tmp_q[6] * 78934266944114L) + ((int128)tmp_q[7] * 34153582712328L) + ((int128)tmp_q[8] * 80254157022577L) + ((int128)tmp_q[9] * 19022908421233L) + ((int128)tmp_q[10] * 29112981397587L);
	tmp_zero[10] = ((int128)tmp_q[0] * 9704327132529L) - ((int128)tmp_q[1] * 10894263903395L) - ((int128)tmp_q[2] * 4533365462952L) - ((int128)tmp_q[3] * 54688774325005L) + ((int128)tmp_q[4] * 68393169537198L) - ((int128)tmp_q[5] * 28495547053807L) - ((int128)tmp_q[6] * 53133996878071L) - ((int128)tmp_q[7] * 78934266944114L) + ((int128)tmp_q[8] * 34153582712328L) + ((int128)tmp_q[9] * 80254157022577L) + ((int128)tmp_q[10] * 19022908421233L);

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

