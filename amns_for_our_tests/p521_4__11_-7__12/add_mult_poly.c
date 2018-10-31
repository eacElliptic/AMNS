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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 7);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 7);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 7);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 7);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 7);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11280122342308212616UL) + ((((uint64_t)op[1] * 15019669605021546041UL) + ((uint64_t)op[2] * 3377335617350371296UL) + ((uint64_t)op[3] * 12806228966935637524UL) + ((uint64_t)op[4] * 15673472094271095791UL) + ((uint64_t)op[5] * 18347310934792866405UL) + ((uint64_t)op[6] * 50516686209208426UL) + ((uint64_t)op[7] * 2386740410197409589UL) + ((uint64_t)op[8] * 13901190488604541461UL) + ((uint64_t)op[9] * 4226845870281841013UL) + ((uint64_t)op[10] * 12548713084972334701UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 12548713084972334701UL) + ((uint64_t)op[1] * 11280122342308212616UL) + ((((uint64_t)op[2] * 15019669605021546041UL) + ((uint64_t)op[3] * 3377335617350371296UL) + ((uint64_t)op[4] * 12806228966935637524UL) + ((uint64_t)op[5] * 15673472094271095791UL) + ((uint64_t)op[6] * 18347310934792866405UL) + ((uint64_t)op[7] * 50516686209208426UL) + ((uint64_t)op[8] * 2386740410197409589UL) + ((uint64_t)op[9] * 13901190488604541461UL) + ((uint64_t)op[10] * 4226845870281841013UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 4226845870281841013UL) + ((uint64_t)op[1] * 12548713084972334701UL) + ((uint64_t)op[2] * 11280122342308212616UL) + ((((uint64_t)op[3] * 15019669605021546041UL) + ((uint64_t)op[4] * 3377335617350371296UL) + ((uint64_t)op[5] * 12806228966935637524UL) + ((uint64_t)op[6] * 15673472094271095791UL) + ((uint64_t)op[7] * 18347310934792866405UL) + ((uint64_t)op[8] * 50516686209208426UL) + ((uint64_t)op[9] * 2386740410197409589UL) + ((uint64_t)op[10] * 13901190488604541461UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 13901190488604541461UL) + ((uint64_t)op[1] * 4226845870281841013UL) + ((uint64_t)op[2] * 12548713084972334701UL) + ((uint64_t)op[3] * 11280122342308212616UL) + ((((uint64_t)op[4] * 15019669605021546041UL) + ((uint64_t)op[5] * 3377335617350371296UL) + ((uint64_t)op[6] * 12806228966935637524UL) + ((uint64_t)op[7] * 15673472094271095791UL) + ((uint64_t)op[8] * 18347310934792866405UL) + ((uint64_t)op[9] * 50516686209208426UL) + ((uint64_t)op[10] * 2386740410197409589UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 2386740410197409589UL) + ((uint64_t)op[1] * 13901190488604541461UL) + ((uint64_t)op[2] * 4226845870281841013UL) + ((uint64_t)op[3] * 12548713084972334701UL) + ((uint64_t)op[4] * 11280122342308212616UL) + ((((uint64_t)op[5] * 15019669605021546041UL) + ((uint64_t)op[6] * 3377335617350371296UL) + ((uint64_t)op[7] * 12806228966935637524UL) + ((uint64_t)op[8] * 15673472094271095791UL) + ((uint64_t)op[9] * 18347310934792866405UL) + ((uint64_t)op[10] * 50516686209208426UL)) * 18446744073709551609);
	tmp_q[5] = ((uint64_t)op[0] * 50516686209208426UL) + ((uint64_t)op[1] * 2386740410197409589UL) + ((uint64_t)op[2] * 13901190488604541461UL) + ((uint64_t)op[3] * 4226845870281841013UL) + ((uint64_t)op[4] * 12548713084972334701UL) + ((uint64_t)op[5] * 11280122342308212616UL) + ((((uint64_t)op[6] * 15019669605021546041UL) + ((uint64_t)op[7] * 3377335617350371296UL) + ((uint64_t)op[8] * 12806228966935637524UL) + ((uint64_t)op[9] * 15673472094271095791UL) + ((uint64_t)op[10] * 18347310934792866405UL)) * 18446744073709551609);
	tmp_q[6] = ((uint64_t)op[0] * 18347310934792866405UL) + ((uint64_t)op[1] * 50516686209208426UL) + ((uint64_t)op[2] * 2386740410197409589UL) + ((uint64_t)op[3] * 13901190488604541461UL) + ((uint64_t)op[4] * 4226845870281841013UL) + ((uint64_t)op[5] * 12548713084972334701UL) + ((uint64_t)op[6] * 11280122342308212616UL) + ((((uint64_t)op[7] * 15019669605021546041UL) + ((uint64_t)op[8] * 3377335617350371296UL) + ((uint64_t)op[9] * 12806228966935637524UL) + ((uint64_t)op[10] * 15673472094271095791UL)) * 18446744073709551609);
	tmp_q[7] = ((uint64_t)op[0] * 15673472094271095791UL) + ((uint64_t)op[1] * 18347310934792866405UL) + ((uint64_t)op[2] * 50516686209208426UL) + ((uint64_t)op[3] * 2386740410197409589UL) + ((uint64_t)op[4] * 13901190488604541461UL) + ((uint64_t)op[5] * 4226845870281841013UL) + ((uint64_t)op[6] * 12548713084972334701UL) + ((uint64_t)op[7] * 11280122342308212616UL) + ((((uint64_t)op[8] * 15019669605021546041UL) + ((uint64_t)op[9] * 3377335617350371296UL) + ((uint64_t)op[10] * 12806228966935637524UL)) * 18446744073709551609);
	tmp_q[8] = ((uint64_t)op[0] * 12806228966935637524UL) + ((uint64_t)op[1] * 15673472094271095791UL) + ((uint64_t)op[2] * 18347310934792866405UL) + ((uint64_t)op[3] * 50516686209208426UL) + ((uint64_t)op[4] * 2386740410197409589UL) + ((uint64_t)op[5] * 13901190488604541461UL) + ((uint64_t)op[6] * 4226845870281841013UL) + ((uint64_t)op[7] * 12548713084972334701UL) + ((uint64_t)op[8] * 11280122342308212616UL) + ((((uint64_t)op[9] * 15019669605021546041UL) + ((uint64_t)op[10] * 3377335617350371296UL)) * 18446744073709551609);
	tmp_q[9] = ((uint64_t)op[0] * 3377335617350371296UL) + ((uint64_t)op[1] * 12806228966935637524UL) + ((uint64_t)op[2] * 15673472094271095791UL) + ((uint64_t)op[3] * 18347310934792866405UL) + ((uint64_t)op[4] * 50516686209208426UL) + ((uint64_t)op[5] * 2386740410197409589UL) + ((uint64_t)op[6] * 13901190488604541461UL) + ((uint64_t)op[7] * 4226845870281841013UL) + ((uint64_t)op[8] * 12548713084972334701UL) + ((uint64_t)op[9] * 11280122342308212616UL) + ((uint64_t)op[10] * 5542777207106487409UL);
	tmp_q[10] = ((uint64_t)op[0] * 15019669605021546041UL) + ((uint64_t)op[1] * 3377335617350371296UL) + ((uint64_t)op[2] * 12806228966935637524UL) + ((uint64_t)op[3] * 15673472094271095791UL) + ((uint64_t)op[4] * 18347310934792866405UL) + ((uint64_t)op[5] * 50516686209208426UL) + ((uint64_t)op[6] * 2386740410197409589UL) + ((uint64_t)op[7] * 13901190488604541461UL) + ((uint64_t)op[8] * 4226845870281841013UL) + ((uint64_t)op[9] * 12548713084972334701UL) + ((uint64_t)op[10] * 11280122342308212616UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 39773243233408L) - ((-((int128)tmp_q[1] * 11931660794139L) + ((int128)tmp_q[2] * 62099515842006L) - ((int128)tmp_q[3] * 51586997215882L) + ((int128)tmp_q[4] * 96984901017004L) - ((int128)tmp_q[5] * 44455680663492L) + ((int128)tmp_q[6] * 37119613791963L) - ((int128)tmp_q[7] * 37613167190874L) + ((int128)tmp_q[8] * 30011776261424L) + ((int128)tmp_q[9] * 59345631955314L) - ((int128)tmp_q[10] * 340591115891L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 340591115891L) - ((int128)tmp_q[1] * 39773243233408L) - ((-((int128)tmp_q[2] * 11931660794139L) + ((int128)tmp_q[3] * 62099515842006L) - ((int128)tmp_q[4] * 51586997215882L) + ((int128)tmp_q[5] * 96984901017004L) - ((int128)tmp_q[6] * 44455680663492L) + ((int128)tmp_q[7] * 37119613791963L) - ((int128)tmp_q[8] * 37613167190874L) + ((int128)tmp_q[9] * 30011776261424L) + ((int128)tmp_q[10] * 59345631955314L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 59345631955314L) - ((int128)tmp_q[1] * 340591115891L) - ((int128)tmp_q[2] * 39773243233408L) - ((-((int128)tmp_q[3] * 11931660794139L) + ((int128)tmp_q[4] * 62099515842006L) - ((int128)tmp_q[5] * 51586997215882L) + ((int128)tmp_q[6] * 96984901017004L) - ((int128)tmp_q[7] * 44455680663492L) + ((int128)tmp_q[8] * 37119613791963L) - ((int128)tmp_q[9] * 37613167190874L) + ((int128)tmp_q[10] * 30011776261424L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 30011776261424L) + ((int128)tmp_q[1] * 59345631955314L) - ((int128)tmp_q[2] * 340591115891L) - ((int128)tmp_q[3] * 39773243233408L) - ((-((int128)tmp_q[4] * 11931660794139L) + ((int128)tmp_q[5] * 62099515842006L) - ((int128)tmp_q[6] * 51586997215882L) + ((int128)tmp_q[7] * 96984901017004L) - ((int128)tmp_q[8] * 44455680663492L) + ((int128)tmp_q[9] * 37119613791963L) - ((int128)tmp_q[10] * 37613167190874L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 37613167190874L) + ((int128)tmp_q[1] * 30011776261424L) + ((int128)tmp_q[2] * 59345631955314L) - ((int128)tmp_q[3] * 340591115891L) - ((int128)tmp_q[4] * 39773243233408L) - ((-((int128)tmp_q[5] * 11931660794139L) + ((int128)tmp_q[6] * 62099515842006L) - ((int128)tmp_q[7] * 51586997215882L) + ((int128)tmp_q[8] * 96984901017004L) - ((int128)tmp_q[9] * 44455680663492L) + ((int128)tmp_q[10] * 37119613791963L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 37119613791963L) - ((int128)tmp_q[1] * 37613167190874L) + ((int128)tmp_q[2] * 30011776261424L) + ((int128)tmp_q[3] * 59345631955314L) - ((int128)tmp_q[4] * 340591115891L) - ((int128)tmp_q[5] * 39773243233408L) - ((-((int128)tmp_q[6] * 11931660794139L) + ((int128)tmp_q[7] * 62099515842006L) - ((int128)tmp_q[8] * 51586997215882L) + ((int128)tmp_q[9] * 96984901017004L) - ((int128)tmp_q[10] * 44455680663492L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 44455680663492L) + ((int128)tmp_q[1] * 37119613791963L) - ((int128)tmp_q[2] * 37613167190874L) + ((int128)tmp_q[3] * 30011776261424L) + ((int128)tmp_q[4] * 59345631955314L) - ((int128)tmp_q[5] * 340591115891L) - ((int128)tmp_q[6] * 39773243233408L) - ((-((int128)tmp_q[7] * 11931660794139L) + ((int128)tmp_q[8] * 62099515842006L) - ((int128)tmp_q[9] * 51586997215882L) + ((int128)tmp_q[10] * 96984901017004L)) * 7);
	tmp_zero[7] = ((int128)tmp_q[0] * 96984901017004L) - ((int128)tmp_q[1] * 44455680663492L) + ((int128)tmp_q[2] * 37119613791963L) - ((int128)tmp_q[3] * 37613167190874L) + ((int128)tmp_q[4] * 30011776261424L) + ((int128)tmp_q[5] * 59345631955314L) - ((int128)tmp_q[6] * 340591115891L) - ((int128)tmp_q[7] * 39773243233408L) - ((-((int128)tmp_q[8] * 11931660794139L) + ((int128)tmp_q[9] * 62099515842006L) - ((int128)tmp_q[10] * 51586997215882L)) * 7);
	tmp_zero[8] = -((int128)tmp_q[0] * 51586997215882L) + ((int128)tmp_q[1] * 96984901017004L) - ((int128)tmp_q[2] * 44455680663492L) + ((int128)tmp_q[3] * 37119613791963L) - ((int128)tmp_q[4] * 37613167190874L) + ((int128)tmp_q[5] * 30011776261424L) + ((int128)tmp_q[6] * 59345631955314L) - ((int128)tmp_q[7] * 340591115891L) - ((int128)tmp_q[8] * 39773243233408L) - ((-((int128)tmp_q[9] * 11931660794139L) + ((int128)tmp_q[10] * 62099515842006L)) * 7);
	tmp_zero[9] = ((int128)tmp_q[0] * 62099515842006L) - ((int128)tmp_q[1] * 51586997215882L) + ((int128)tmp_q[2] * 96984901017004L) - ((int128)tmp_q[3] * 44455680663492L) + ((int128)tmp_q[4] * 37119613791963L) - ((int128)tmp_q[5] * 37613167190874L) + ((int128)tmp_q[6] * 30011776261424L) + ((int128)tmp_q[7] * 59345631955314L) - ((int128)tmp_q[8] * 340591115891L) - ((int128)tmp_q[9] * 39773243233408L) + ((int128)tmp_q[10] * 83521625558973L);
	tmp_zero[10] = -((int128)tmp_q[0] * 11931660794139L) + ((int128)tmp_q[1] * 62099515842006L) - ((int128)tmp_q[2] * 51586997215882L) + ((int128)tmp_q[3] * 96984901017004L) - ((int128)tmp_q[4] * 44455680663492L) + ((int128)tmp_q[5] * 37119613791963L) - ((int128)tmp_q[6] * 37613167190874L) + ((int128)tmp_q[7] * 30011776261424L) + ((int128)tmp_q[8] * 59345631955314L) - ((int128)tmp_q[9] * 340591115891L) - ((int128)tmp_q[10] * 39773243233408L);

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

