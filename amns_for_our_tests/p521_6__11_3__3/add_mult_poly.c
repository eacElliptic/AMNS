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
	tmp_q[0] = ((uint64_t)op[0] * 12694316800453417886UL) + ((((uint64_t)op[1] * 5246118601987035746UL) + ((uint64_t)op[2] * 14964024854082554203UL) + ((uint64_t)op[3] * 11489664896712655265UL) + ((uint64_t)op[4] * 12363459581221491832UL) + ((uint64_t)op[5] * 16370095129477163223UL) + ((uint64_t)op[6] * 12625624663942593045UL) + ((uint64_t)op[7] * 5731018458855753572UL) + ((uint64_t)op[8] * 15731122408470736355UL) + ((uint64_t)op[9] * 13232851297713076765UL) + ((uint64_t)op[10] * 6669782978401248703UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 6669782978401248703UL) + ((uint64_t)op[1] * 12694316800453417886UL) + ((((uint64_t)op[2] * 5246118601987035746UL) + ((uint64_t)op[3] * 14964024854082554203UL) + ((uint64_t)op[4] * 11489664896712655265UL) + ((uint64_t)op[5] * 12363459581221491832UL) + ((uint64_t)op[6] * 16370095129477163223UL) + ((uint64_t)op[7] * 12625624663942593045UL) + ((uint64_t)op[8] * 5731018458855753572UL) + ((uint64_t)op[9] * 15731122408470736355UL) + ((uint64_t)op[10] * 13232851297713076765UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 13232851297713076765UL) + ((uint64_t)op[1] * 6669782978401248703UL) + ((uint64_t)op[2] * 12694316800453417886UL) + ((((uint64_t)op[3] * 5246118601987035746UL) + ((uint64_t)op[4] * 14964024854082554203UL) + ((uint64_t)op[5] * 11489664896712655265UL) + ((uint64_t)op[6] * 12363459581221491832UL) + ((uint64_t)op[7] * 16370095129477163223UL) + ((uint64_t)op[8] * 12625624663942593045UL) + ((uint64_t)op[9] * 5731018458855753572UL) + ((uint64_t)op[10] * 15731122408470736355UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 15731122408470736355UL) + ((uint64_t)op[1] * 13232851297713076765UL) + ((uint64_t)op[2] * 6669782978401248703UL) + ((uint64_t)op[3] * 12694316800453417886UL) + ((((uint64_t)op[4] * 5246118601987035746UL) + ((uint64_t)op[5] * 14964024854082554203UL) + ((uint64_t)op[6] * 11489664896712655265UL) + ((uint64_t)op[7] * 12363459581221491832UL) + ((uint64_t)op[8] * 16370095129477163223UL) + ((uint64_t)op[9] * 12625624663942593045UL) + ((uint64_t)op[10] * 5731018458855753572UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 5731018458855753572UL) + ((uint64_t)op[1] * 15731122408470736355UL) + ((uint64_t)op[2] * 13232851297713076765UL) + ((uint64_t)op[3] * 6669782978401248703UL) + ((uint64_t)op[4] * 12694316800453417886UL) + ((((uint64_t)op[5] * 5246118601987035746UL) + ((uint64_t)op[6] * 14964024854082554203UL) + ((uint64_t)op[7] * 11489664896712655265UL) + ((uint64_t)op[8] * 12363459581221491832UL) + ((uint64_t)op[9] * 16370095129477163223UL) + ((uint64_t)op[10] * 12625624663942593045UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 12625624663942593045UL) + ((uint64_t)op[1] * 5731018458855753572UL) + ((uint64_t)op[2] * 15731122408470736355UL) + ((uint64_t)op[3] * 13232851297713076765UL) + ((uint64_t)op[4] * 6669782978401248703UL) + ((uint64_t)op[5] * 12694316800453417886UL) + ((((uint64_t)op[6] * 5246118601987035746UL) + ((uint64_t)op[7] * 14964024854082554203UL) + ((uint64_t)op[8] * 11489664896712655265UL) + ((uint64_t)op[9] * 12363459581221491832UL) + ((uint64_t)op[10] * 16370095129477163223UL)) * 3);
	tmp_q[6] = ((uint64_t)op[0] * 16370095129477163223UL) + ((uint64_t)op[1] * 12625624663942593045UL) + ((uint64_t)op[2] * 5731018458855753572UL) + ((uint64_t)op[3] * 15731122408470736355UL) + ((uint64_t)op[4] * 13232851297713076765UL) + ((uint64_t)op[5] * 6669782978401248703UL) + ((uint64_t)op[6] * 12694316800453417886UL) + ((((uint64_t)op[7] * 5246118601987035746UL) + ((uint64_t)op[8] * 14964024854082554203UL) + ((uint64_t)op[9] * 11489664896712655265UL) + ((uint64_t)op[10] * 12363459581221491832UL)) * 3);
	tmp_q[7] = ((uint64_t)op[0] * 12363459581221491832UL) + ((uint64_t)op[1] * 16370095129477163223UL) + ((uint64_t)op[2] * 12625624663942593045UL) + ((uint64_t)op[3] * 5731018458855753572UL) + ((uint64_t)op[4] * 15731122408470736355UL) + ((uint64_t)op[5] * 13232851297713076765UL) + ((uint64_t)op[6] * 6669782978401248703UL) + ((uint64_t)op[7] * 12694316800453417886UL) + ((((uint64_t)op[8] * 5246118601987035746UL) + ((uint64_t)op[9] * 14964024854082554203UL) + ((uint64_t)op[10] * 11489664896712655265UL)) * 3);
	tmp_q[8] = ((uint64_t)op[0] * 11489664896712655265UL) + ((uint64_t)op[1] * 12363459581221491832UL) + ((uint64_t)op[2] * 16370095129477163223UL) + ((uint64_t)op[3] * 12625624663942593045UL) + ((uint64_t)op[4] * 5731018458855753572UL) + ((uint64_t)op[5] * 15731122408470736355UL) + ((uint64_t)op[6] * 13232851297713076765UL) + ((uint64_t)op[7] * 6669782978401248703UL) + ((uint64_t)op[8] * 12694316800453417886UL) + ((((uint64_t)op[9] * 5246118601987035746UL) + ((uint64_t)op[10] * 14964024854082554203UL)) * 3);
	tmp_q[9] = ((uint64_t)op[0] * 14964024854082554203UL) + ((uint64_t)op[1] * 11489664896712655265UL) + ((uint64_t)op[2] * 12363459581221491832UL) + ((uint64_t)op[3] * 16370095129477163223UL) + ((uint64_t)op[4] * 12625624663942593045UL) + ((uint64_t)op[5] * 5731018458855753572UL) + ((uint64_t)op[6] * 15731122408470736355UL) + ((uint64_t)op[7] * 13232851297713076765UL) + ((uint64_t)op[8] * 6669782978401248703UL) + ((uint64_t)op[9] * 12694316800453417886UL) + ((uint64_t)op[10] * 15738355805961107238UL);
	tmp_q[10] = ((uint64_t)op[0] * 5246118601987035746UL) + ((uint64_t)op[1] * 14964024854082554203UL) + ((uint64_t)op[2] * 11489664896712655265UL) + ((uint64_t)op[3] * 12363459581221491832UL) + ((uint64_t)op[4] * 16370095129477163223UL) + ((uint64_t)op[5] * 12625624663942593045UL) + ((uint64_t)op[6] * 5731018458855753572UL) + ((uint64_t)op[7] * 15731122408470736355UL) + ((uint64_t)op[8] * 13232851297713076765UL) + ((uint64_t)op[9] * 6669782978401248703UL) + ((uint64_t)op[10] * 12694316800453417886UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 120549136349366L) + ((-((int128)tmp_q[1] * 59506204673116L) + ((int128)tmp_q[2] * 30966918800194L) - ((int128)tmp_q[3] * 2127715356577L) + ((int128)tmp_q[4] * 97115222088408L) - ((int128)tmp_q[5] * 1522270830161L) - ((int128)tmp_q[6] * 44517189046579L) - ((int128)tmp_q[7] * 3667720066427L) - ((int128)tmp_q[8] * 71985634957564L) - ((int128)tmp_q[9] * 118327328156310L) + ((int128)tmp_q[10] * 60523056091851L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 60523056091851L) + ((int128)tmp_q[1] * 120549136349366L) + ((-((int128)tmp_q[2] * 59506204673116L) + ((int128)tmp_q[3] * 30966918800194L) - ((int128)tmp_q[4] * 2127715356577L) + ((int128)tmp_q[5] * 97115222088408L) - ((int128)tmp_q[6] * 1522270830161L) - ((int128)tmp_q[7] * 44517189046579L) - ((int128)tmp_q[8] * 3667720066427L) - ((int128)tmp_q[9] * 71985634957564L) - ((int128)tmp_q[10] * 118327328156310L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 118327328156310L) + ((int128)tmp_q[1] * 60523056091851L) + ((int128)tmp_q[2] * 120549136349366L) + ((-((int128)tmp_q[3] * 59506204673116L) + ((int128)tmp_q[4] * 30966918800194L) - ((int128)tmp_q[5] * 2127715356577L) + ((int128)tmp_q[6] * 97115222088408L) - ((int128)tmp_q[7] * 1522270830161L) - ((int128)tmp_q[8] * 44517189046579L) - ((int128)tmp_q[9] * 3667720066427L) - ((int128)tmp_q[10] * 71985634957564L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 71985634957564L) - ((int128)tmp_q[1] * 118327328156310L) + ((int128)tmp_q[2] * 60523056091851L) + ((int128)tmp_q[3] * 120549136349366L) + ((-((int128)tmp_q[4] * 59506204673116L) + ((int128)tmp_q[5] * 30966918800194L) - ((int128)tmp_q[6] * 2127715356577L) + ((int128)tmp_q[7] * 97115222088408L) - ((int128)tmp_q[8] * 1522270830161L) - ((int128)tmp_q[9] * 44517189046579L) - ((int128)tmp_q[10] * 3667720066427L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 3667720066427L) - ((int128)tmp_q[1] * 71985634957564L) - ((int128)tmp_q[2] * 118327328156310L) + ((int128)tmp_q[3] * 60523056091851L) + ((int128)tmp_q[4] * 120549136349366L) + ((-((int128)tmp_q[5] * 59506204673116L) + ((int128)tmp_q[6] * 30966918800194L) - ((int128)tmp_q[7] * 2127715356577L) + ((int128)tmp_q[8] * 97115222088408L) - ((int128)tmp_q[9] * 1522270830161L) - ((int128)tmp_q[10] * 44517189046579L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 44517189046579L) - ((int128)tmp_q[1] * 3667720066427L) - ((int128)tmp_q[2] * 71985634957564L) - ((int128)tmp_q[3] * 118327328156310L) + ((int128)tmp_q[4] * 60523056091851L) + ((int128)tmp_q[5] * 120549136349366L) + ((-((int128)tmp_q[6] * 59506204673116L) + ((int128)tmp_q[7] * 30966918800194L) - ((int128)tmp_q[8] * 2127715356577L) + ((int128)tmp_q[9] * 97115222088408L) - ((int128)tmp_q[10] * 1522270830161L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 1522270830161L) - ((int128)tmp_q[1] * 44517189046579L) - ((int128)tmp_q[2] * 3667720066427L) - ((int128)tmp_q[3] * 71985634957564L) - ((int128)tmp_q[4] * 118327328156310L) + ((int128)tmp_q[5] * 60523056091851L) + ((int128)tmp_q[6] * 120549136349366L) + ((-((int128)tmp_q[7] * 59506204673116L) + ((int128)tmp_q[8] * 30966918800194L) - ((int128)tmp_q[9] * 2127715356577L) + ((int128)tmp_q[10] * 97115222088408L)) * 3);
	tmp_zero[7] = ((int128)tmp_q[0] * 97115222088408L) - ((int128)tmp_q[1] * 1522270830161L) - ((int128)tmp_q[2] * 44517189046579L) - ((int128)tmp_q[3] * 3667720066427L) - ((int128)tmp_q[4] * 71985634957564L) - ((int128)tmp_q[5] * 118327328156310L) + ((int128)tmp_q[6] * 60523056091851L) + ((int128)tmp_q[7] * 120549136349366L) + ((-((int128)tmp_q[8] * 59506204673116L) + ((int128)tmp_q[9] * 30966918800194L) - ((int128)tmp_q[10] * 2127715356577L)) * 3);
	tmp_zero[8] = -((int128)tmp_q[0] * 2127715356577L) + ((int128)tmp_q[1] * 97115222088408L) - ((int128)tmp_q[2] * 1522270830161L) - ((int128)tmp_q[3] * 44517189046579L) - ((int128)tmp_q[4] * 3667720066427L) - ((int128)tmp_q[5] * 71985634957564L) - ((int128)tmp_q[6] * 118327328156310L) + ((int128)tmp_q[7] * 60523056091851L) + ((int128)tmp_q[8] * 120549136349366L) + ((-((int128)tmp_q[9] * 59506204673116L) + ((int128)tmp_q[10] * 30966918800194L)) * 3);
	tmp_zero[9] = ((int128)tmp_q[0] * 30966918800194L) - ((int128)tmp_q[1] * 2127715356577L) + ((int128)tmp_q[2] * 97115222088408L) - ((int128)tmp_q[3] * 1522270830161L) - ((int128)tmp_q[4] * 44517189046579L) - ((int128)tmp_q[5] * 3667720066427L) - ((int128)tmp_q[6] * 71985634957564L) - ((int128)tmp_q[7] * 118327328156310L) + ((int128)tmp_q[8] * 60523056091851L) + ((int128)tmp_q[9] * 120549136349366L) - ((int128)tmp_q[10] * 178518614019348L);
	tmp_zero[10] = -((int128)tmp_q[0] * 59506204673116L) + ((int128)tmp_q[1] * 30966918800194L) - ((int128)tmp_q[2] * 2127715356577L) + ((int128)tmp_q[3] * 97115222088408L) - ((int128)tmp_q[4] * 1522270830161L) - ((int128)tmp_q[5] * 44517189046579L) - ((int128)tmp_q[6] * 3667720066427L) - ((int128)tmp_q[7] * 71985634957564L) - ((int128)tmp_q[8] * 118327328156310L) + ((int128)tmp_q[9] * 60523056091851L) + ((int128)tmp_q[10] * 120549136349366L);

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

