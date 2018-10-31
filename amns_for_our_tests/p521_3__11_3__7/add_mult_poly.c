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
	tmp_q[0] = ((uint64_t)op[0] * 5510310906050008464UL) + ((((uint64_t)op[1] * 10431918495173763022UL) + ((uint64_t)op[2] * 17931579735822029944UL) + ((uint64_t)op[3] * 1961756715986181945UL) + ((uint64_t)op[4] * 954676183017573513UL) + ((uint64_t)op[5] * 2888905511320666775UL) + ((uint64_t)op[6] * 2015417584680868154UL) + ((uint64_t)op[7] * 5540982362796196841UL) + ((uint64_t)op[8] * 6789352237611447427UL) + ((uint64_t)op[9] * 5971117534308177417UL) + ((uint64_t)op[10] * 3076127628672531163UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 3076127628672531163UL) + ((uint64_t)op[1] * 5510310906050008464UL) + ((((uint64_t)op[2] * 10431918495173763022UL) + ((uint64_t)op[3] * 17931579735822029944UL) + ((uint64_t)op[4] * 1961756715986181945UL) + ((uint64_t)op[5] * 954676183017573513UL) + ((uint64_t)op[6] * 2888905511320666775UL) + ((uint64_t)op[7] * 2015417584680868154UL) + ((uint64_t)op[8] * 5540982362796196841UL) + ((uint64_t)op[9] * 6789352237611447427UL) + ((uint64_t)op[10] * 5971117534308177417UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 5971117534308177417UL) + ((uint64_t)op[1] * 3076127628672531163UL) + ((uint64_t)op[2] * 5510310906050008464UL) + ((((uint64_t)op[3] * 10431918495173763022UL) + ((uint64_t)op[4] * 17931579735822029944UL) + ((uint64_t)op[5] * 1961756715986181945UL) + ((uint64_t)op[6] * 954676183017573513UL) + ((uint64_t)op[7] * 2888905511320666775UL) + ((uint64_t)op[8] * 2015417584680868154UL) + ((uint64_t)op[9] * 5540982362796196841UL) + ((uint64_t)op[10] * 6789352237611447427UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 6789352237611447427UL) + ((uint64_t)op[1] * 5971117534308177417UL) + ((uint64_t)op[2] * 3076127628672531163UL) + ((uint64_t)op[3] * 5510310906050008464UL) + ((((uint64_t)op[4] * 10431918495173763022UL) + ((uint64_t)op[5] * 17931579735822029944UL) + ((uint64_t)op[6] * 1961756715986181945UL) + ((uint64_t)op[7] * 954676183017573513UL) + ((uint64_t)op[8] * 2888905511320666775UL) + ((uint64_t)op[9] * 2015417584680868154UL) + ((uint64_t)op[10] * 5540982362796196841UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 5540982362796196841UL) + ((uint64_t)op[1] * 6789352237611447427UL) + ((uint64_t)op[2] * 5971117534308177417UL) + ((uint64_t)op[3] * 3076127628672531163UL) + ((uint64_t)op[4] * 5510310906050008464UL) + ((((uint64_t)op[5] * 10431918495173763022UL) + ((uint64_t)op[6] * 17931579735822029944UL) + ((uint64_t)op[7] * 1961756715986181945UL) + ((uint64_t)op[8] * 954676183017573513UL) + ((uint64_t)op[9] * 2888905511320666775UL) + ((uint64_t)op[10] * 2015417584680868154UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 2015417584680868154UL) + ((uint64_t)op[1] * 5540982362796196841UL) + ((uint64_t)op[2] * 6789352237611447427UL) + ((uint64_t)op[3] * 5971117534308177417UL) + ((uint64_t)op[4] * 3076127628672531163UL) + ((uint64_t)op[5] * 5510310906050008464UL) + ((((uint64_t)op[6] * 10431918495173763022UL) + ((uint64_t)op[7] * 17931579735822029944UL) + ((uint64_t)op[8] * 1961756715986181945UL) + ((uint64_t)op[9] * 954676183017573513UL) + ((uint64_t)op[10] * 2888905511320666775UL)) * 3);
	tmp_q[6] = ((uint64_t)op[0] * 2888905511320666775UL) + ((uint64_t)op[1] * 2015417584680868154UL) + ((uint64_t)op[2] * 5540982362796196841UL) + ((uint64_t)op[3] * 6789352237611447427UL) + ((uint64_t)op[4] * 5971117534308177417UL) + ((uint64_t)op[5] * 3076127628672531163UL) + ((uint64_t)op[6] * 5510310906050008464UL) + ((((uint64_t)op[7] * 10431918495173763022UL) + ((uint64_t)op[8] * 17931579735822029944UL) + ((uint64_t)op[9] * 1961756715986181945UL) + ((uint64_t)op[10] * 954676183017573513UL)) * 3);
	tmp_q[7] = ((uint64_t)op[0] * 954676183017573513UL) + ((uint64_t)op[1] * 2888905511320666775UL) + ((uint64_t)op[2] * 2015417584680868154UL) + ((uint64_t)op[3] * 5540982362796196841UL) + ((uint64_t)op[4] * 6789352237611447427UL) + ((uint64_t)op[5] * 5971117534308177417UL) + ((uint64_t)op[6] * 3076127628672531163UL) + ((uint64_t)op[7] * 5510310906050008464UL) + ((((uint64_t)op[8] * 10431918495173763022UL) + ((uint64_t)op[9] * 17931579735822029944UL) + ((uint64_t)op[10] * 1961756715986181945UL)) * 3);
	tmp_q[8] = ((uint64_t)op[0] * 1961756715986181945UL) + ((uint64_t)op[1] * 954676183017573513UL) + ((uint64_t)op[2] * 2888905511320666775UL) + ((uint64_t)op[3] * 2015417584680868154UL) + ((uint64_t)op[4] * 5540982362796196841UL) + ((uint64_t)op[5] * 6789352237611447427UL) + ((uint64_t)op[6] * 5971117534308177417UL) + ((uint64_t)op[7] * 3076127628672531163UL) + ((uint64_t)op[8] * 5510310906050008464UL) + ((((uint64_t)op[9] * 10431918495173763022UL) + ((uint64_t)op[10] * 17931579735822029944UL)) * 3);
	tmp_q[9] = ((uint64_t)op[0] * 17931579735822029944UL) + ((uint64_t)op[1] * 1961756715986181945UL) + ((uint64_t)op[2] * 954676183017573513UL) + ((uint64_t)op[3] * 2888905511320666775UL) + ((uint64_t)op[4] * 2015417584680868154UL) + ((uint64_t)op[5] * 5540982362796196841UL) + ((uint64_t)op[6] * 6789352237611447427UL) + ((uint64_t)op[7] * 5971117534308177417UL) + ((uint64_t)op[8] * 3076127628672531163UL) + ((uint64_t)op[9] * 5510310906050008464UL) + ((uint64_t)op[10] * 12849011411811737450UL);
	tmp_q[10] = ((uint64_t)op[0] * 10431918495173763022UL) + ((uint64_t)op[1] * 17931579735822029944UL) + ((uint64_t)op[2] * 1961756715986181945UL) + ((uint64_t)op[3] * 954676183017573513UL) + ((uint64_t)op[4] * 2888905511320666775UL) + ((uint64_t)op[5] * 2015417584680868154UL) + ((uint64_t)op[6] * 5540982362796196841UL) + ((uint64_t)op[7] * 6789352237611447427UL) + ((uint64_t)op[8] * 5971117534308177417UL) + ((uint64_t)op[9] * 3076127628672531163UL) + ((uint64_t)op[10] * 5510310906050008464UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 121168971387998L) + ((((int128)tmp_q[1] * 64555460360566L) - ((int128)tmp_q[2] * 25939610498852L) + ((int128)tmp_q[3] * 58881609217007L) - ((int128)tmp_q[4] * 56508713578851L) - ((int128)tmp_q[5] * 11026755104694L) - ((int128)tmp_q[6] * 14319793846230L) - ((int128)tmp_q[7] * 95968030528387L) - ((int128)tmp_q[8] * 38068969464516L) - ((int128)tmp_q[9] * 88349226804803L) + ((int128)tmp_q[10] * 18050480051135L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 18050480051135L) - ((int128)tmp_q[1] * 121168971387998L) + ((((int128)tmp_q[2] * 64555460360566L) - ((int128)tmp_q[3] * 25939610498852L) + ((int128)tmp_q[4] * 58881609217007L) - ((int128)tmp_q[5] * 56508713578851L) - ((int128)tmp_q[6] * 11026755104694L) - ((int128)tmp_q[7] * 14319793846230L) - ((int128)tmp_q[8] * 95968030528387L) - ((int128)tmp_q[9] * 38068969464516L) - ((int128)tmp_q[10] * 88349226804803L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 88349226804803L) + ((int128)tmp_q[1] * 18050480051135L) - ((int128)tmp_q[2] * 121168971387998L) + ((((int128)tmp_q[3] * 64555460360566L) - ((int128)tmp_q[4] * 25939610498852L) + ((int128)tmp_q[5] * 58881609217007L) - ((int128)tmp_q[6] * 56508713578851L) - ((int128)tmp_q[7] * 11026755104694L) - ((int128)tmp_q[8] * 14319793846230L) - ((int128)tmp_q[9] * 95968030528387L) - ((int128)tmp_q[10] * 38068969464516L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 38068969464516L) - ((int128)tmp_q[1] * 88349226804803L) + ((int128)tmp_q[2] * 18050480051135L) - ((int128)tmp_q[3] * 121168971387998L) + ((((int128)tmp_q[4] * 64555460360566L) - ((int128)tmp_q[5] * 25939610498852L) + ((int128)tmp_q[6] * 58881609217007L) - ((int128)tmp_q[7] * 56508713578851L) - ((int128)tmp_q[8] * 11026755104694L) - ((int128)tmp_q[9] * 14319793846230L) - ((int128)tmp_q[10] * 95968030528387L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 95968030528387L) - ((int128)tmp_q[1] * 38068969464516L) - ((int128)tmp_q[2] * 88349226804803L) + ((int128)tmp_q[3] * 18050480051135L) - ((int128)tmp_q[4] * 121168971387998L) + ((((int128)tmp_q[5] * 64555460360566L) - ((int128)tmp_q[6] * 25939610498852L) + ((int128)tmp_q[7] * 58881609217007L) - ((int128)tmp_q[8] * 56508713578851L) - ((int128)tmp_q[9] * 11026755104694L) - ((int128)tmp_q[10] * 14319793846230L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 14319793846230L) - ((int128)tmp_q[1] * 95968030528387L) - ((int128)tmp_q[2] * 38068969464516L) - ((int128)tmp_q[3] * 88349226804803L) + ((int128)tmp_q[4] * 18050480051135L) - ((int128)tmp_q[5] * 121168971387998L) + ((((int128)tmp_q[6] * 64555460360566L) - ((int128)tmp_q[7] * 25939610498852L) + ((int128)tmp_q[8] * 58881609217007L) - ((int128)tmp_q[9] * 56508713578851L) - ((int128)tmp_q[10] * 11026755104694L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 11026755104694L) - ((int128)tmp_q[1] * 14319793846230L) - ((int128)tmp_q[2] * 95968030528387L) - ((int128)tmp_q[3] * 38068969464516L) - ((int128)tmp_q[4] * 88349226804803L) + ((int128)tmp_q[5] * 18050480051135L) - ((int128)tmp_q[6] * 121168971387998L) + ((((int128)tmp_q[7] * 64555460360566L) - ((int128)tmp_q[8] * 25939610498852L) + ((int128)tmp_q[9] * 58881609217007L) - ((int128)tmp_q[10] * 56508713578851L)) * 3);
	tmp_zero[7] = -((int128)tmp_q[0] * 56508713578851L) - ((int128)tmp_q[1] * 11026755104694L) - ((int128)tmp_q[2] * 14319793846230L) - ((int128)tmp_q[3] * 95968030528387L) - ((int128)tmp_q[4] * 38068969464516L) - ((int128)tmp_q[5] * 88349226804803L) + ((int128)tmp_q[6] * 18050480051135L) - ((int128)tmp_q[7] * 121168971387998L) + ((((int128)tmp_q[8] * 64555460360566L) - ((int128)tmp_q[9] * 25939610498852L) + ((int128)tmp_q[10] * 58881609217007L)) * 3);
	tmp_zero[8] = ((int128)tmp_q[0] * 58881609217007L) - ((int128)tmp_q[1] * 56508713578851L) - ((int128)tmp_q[2] * 11026755104694L) - ((int128)tmp_q[3] * 14319793846230L) - ((int128)tmp_q[4] * 95968030528387L) - ((int128)tmp_q[5] * 38068969464516L) - ((int128)tmp_q[6] * 88349226804803L) + ((int128)tmp_q[7] * 18050480051135L) - ((int128)tmp_q[8] * 121168971387998L) + ((((int128)tmp_q[9] * 64555460360566L) - ((int128)tmp_q[10] * 25939610498852L)) * 3);
	tmp_zero[9] = -((int128)tmp_q[0] * 25939610498852L) + ((int128)tmp_q[1] * 58881609217007L) - ((int128)tmp_q[2] * 56508713578851L) - ((int128)tmp_q[3] * 11026755104694L) - ((int128)tmp_q[4] * 14319793846230L) - ((int128)tmp_q[5] * 95968030528387L) - ((int128)tmp_q[6] * 38068969464516L) - ((int128)tmp_q[7] * 88349226804803L) + ((int128)tmp_q[8] * 18050480051135L) - ((int128)tmp_q[9] * 121168971387998L) + ((int128)tmp_q[10] * 193666381081698L);
	tmp_zero[10] = ((int128)tmp_q[0] * 64555460360566L) - ((int128)tmp_q[1] * 25939610498852L) + ((int128)tmp_q[2] * 58881609217007L) - ((int128)tmp_q[3] * 56508713578851L) - ((int128)tmp_q[4] * 11026755104694L) - ((int128)tmp_q[5] * 14319793846230L) - ((int128)tmp_q[6] * 95968030528387L) - ((int128)tmp_q[7] * 38068969464516L) - ((int128)tmp_q[8] * 88349226804803L) + ((int128)tmp_q[9] * 18050480051135L) - ((int128)tmp_q[10] * 121168971387998L);

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

