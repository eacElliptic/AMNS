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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 7);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 7);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 7);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 7);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 7);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12138843064315026658UL) + ((((uint64_t)op[1] * 3951147444060466469UL) + ((uint64_t)op[2] * 8139484458311945695UL) + ((uint64_t)op[3] * 8334986941968672927UL) + ((uint64_t)op[4] * 1289423320853126629UL) + ((uint64_t)op[5] * 15537302282930645955UL) + ((uint64_t)op[6] * 11085811501740726250UL) + ((uint64_t)op[7] * 15537574914040770804UL) + ((uint64_t)op[8] * 8604190431724035639UL) + ((uint64_t)op[9] * 16023725082466963326UL) + ((uint64_t)op[10] * 11532446015081824051UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 11532446015081824051UL) + ((uint64_t)op[1] * 12138843064315026658UL) + ((((uint64_t)op[2] * 3951147444060466469UL) + ((uint64_t)op[3] * 8139484458311945695UL) + ((uint64_t)op[4] * 8334986941968672927UL) + ((uint64_t)op[5] * 1289423320853126629UL) + ((uint64_t)op[6] * 15537302282930645955UL) + ((uint64_t)op[7] * 11085811501740726250UL) + ((uint64_t)op[8] * 15537574914040770804UL) + ((uint64_t)op[9] * 8604190431724035639UL) + ((uint64_t)op[10] * 16023725082466963326UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 16023725082466963326UL) + ((uint64_t)op[1] * 11532446015081824051UL) + ((uint64_t)op[2] * 12138843064315026658UL) + ((((uint64_t)op[3] * 3951147444060466469UL) + ((uint64_t)op[4] * 8139484458311945695UL) + ((uint64_t)op[5] * 8334986941968672927UL) + ((uint64_t)op[6] * 1289423320853126629UL) + ((uint64_t)op[7] * 15537302282930645955UL) + ((uint64_t)op[8] * 11085811501740726250UL) + ((uint64_t)op[9] * 15537574914040770804UL) + ((uint64_t)op[10] * 8604190431724035639UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 8604190431724035639UL) + ((uint64_t)op[1] * 16023725082466963326UL) + ((uint64_t)op[2] * 11532446015081824051UL) + ((uint64_t)op[3] * 12138843064315026658UL) + ((((uint64_t)op[4] * 3951147444060466469UL) + ((uint64_t)op[5] * 8139484458311945695UL) + ((uint64_t)op[6] * 8334986941968672927UL) + ((uint64_t)op[7] * 1289423320853126629UL) + ((uint64_t)op[8] * 15537302282930645955UL) + ((uint64_t)op[9] * 11085811501740726250UL) + ((uint64_t)op[10] * 15537574914040770804UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 15537574914040770804UL) + ((uint64_t)op[1] * 8604190431724035639UL) + ((uint64_t)op[2] * 16023725082466963326UL) + ((uint64_t)op[3] * 11532446015081824051UL) + ((uint64_t)op[4] * 12138843064315026658UL) + ((((uint64_t)op[5] * 3951147444060466469UL) + ((uint64_t)op[6] * 8139484458311945695UL) + ((uint64_t)op[7] * 8334986941968672927UL) + ((uint64_t)op[8] * 1289423320853126629UL) + ((uint64_t)op[9] * 15537302282930645955UL) + ((uint64_t)op[10] * 11085811501740726250UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 11085811501740726250UL) + ((uint64_t)op[1] * 15537574914040770804UL) + ((uint64_t)op[2] * 8604190431724035639UL) + ((uint64_t)op[3] * 16023725082466963326UL) + ((uint64_t)op[4] * 11532446015081824051UL) + ((uint64_t)op[5] * 12138843064315026658UL) + ((((uint64_t)op[6] * 3951147444060466469UL) + ((uint64_t)op[7] * 8139484458311945695UL) + ((uint64_t)op[8] * 8334986941968672927UL) + ((uint64_t)op[9] * 1289423320853126629UL) + ((uint64_t)op[10] * 15537302282930645955UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 15537302282930645955UL) + ((uint64_t)op[1] * 11085811501740726250UL) + ((uint64_t)op[2] * 15537574914040770804UL) + ((uint64_t)op[3] * 8604190431724035639UL) + ((uint64_t)op[4] * 16023725082466963326UL) + ((uint64_t)op[5] * 11532446015081824051UL) + ((uint64_t)op[6] * 12138843064315026658UL) + ((((uint64_t)op[7] * 3951147444060466469UL) + ((uint64_t)op[8] * 8139484458311945695UL) + ((uint64_t)op[9] * 8334986941968672927UL) + ((uint64_t)op[10] * 1289423320853126629UL)) * 7);
	tmp_q[7] = ((uint64_t)op[0] * 1289423320853126629UL) + ((uint64_t)op[1] * 15537302282930645955UL) + ((uint64_t)op[2] * 11085811501740726250UL) + ((uint64_t)op[3] * 15537574914040770804UL) + ((uint64_t)op[4] * 8604190431724035639UL) + ((uint64_t)op[5] * 16023725082466963326UL) + ((uint64_t)op[6] * 11532446015081824051UL) + ((uint64_t)op[7] * 12138843064315026658UL) + ((((uint64_t)op[8] * 3951147444060466469UL) + ((uint64_t)op[9] * 8139484458311945695UL) + ((uint64_t)op[10] * 8334986941968672927UL)) * 7);
	tmp_q[8] = ((uint64_t)op[0] * 8334986941968672927UL) + ((uint64_t)op[1] * 1289423320853126629UL) + ((uint64_t)op[2] * 15537302282930645955UL) + ((uint64_t)op[3] * 11085811501740726250UL) + ((uint64_t)op[4] * 15537574914040770804UL) + ((uint64_t)op[5] * 8604190431724035639UL) + ((uint64_t)op[6] * 16023725082466963326UL) + ((uint64_t)op[7] * 11532446015081824051UL) + ((uint64_t)op[8] * 12138843064315026658UL) + ((((uint64_t)op[9] * 3951147444060466469UL) + ((uint64_t)op[10] * 8139484458311945695UL)) * 7);
	tmp_q[9] = ((uint64_t)op[0] * 8139484458311945695UL) + ((uint64_t)op[1] * 8334986941968672927UL) + ((uint64_t)op[2] * 1289423320853126629UL) + ((uint64_t)op[3] * 15537302282930645955UL) + ((uint64_t)op[4] * 11085811501740726250UL) + ((uint64_t)op[5] * 15537574914040770804UL) + ((uint64_t)op[6] * 8604190431724035639UL) + ((uint64_t)op[7] * 16023725082466963326UL) + ((uint64_t)op[8] * 11532446015081824051UL) + ((uint64_t)op[9] * 12138843064315026658UL) + ((uint64_t)op[10] * 9211288034713713667UL);
	tmp_q[10] = ((uint64_t)op[0] * 3951147444060466469UL) + ((uint64_t)op[1] * 8139484458311945695UL) + ((uint64_t)op[2] * 8334986941968672927UL) + ((uint64_t)op[3] * 1289423320853126629UL) + ((uint64_t)op[4] * 15537302282930645955UL) + ((uint64_t)op[5] * 11085811501740726250UL) + ((uint64_t)op[6] * 15537574914040770804UL) + ((uint64_t)op[7] * 8604190431724035639UL) + ((uint64_t)op[8] * 16023725082466963326UL) + ((uint64_t)op[9] * 11532446015081824051UL) + ((uint64_t)op[10] * 12138843064315026658UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 23664034743236L) + ((-((int128)tmp_q[1] * 75298483926170L) + ((int128)tmp_q[2] * 40129352224870L) + ((int128)tmp_q[3] * 83613586901033L) - ((int128)tmp_q[4] * 71625886439504L) - ((int128)tmp_q[5] * 73095569950280L) + ((int128)tmp_q[6] * 16212499142000L) - ((int128)tmp_q[7] * 39120728066137L) - ((int128)tmp_q[8] * 122570424992969L) - ((int128)tmp_q[9] * 7265070031016L) + ((int128)tmp_q[10] * 41106866914776L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 41106866914776L) - ((int128)tmp_q[1] * 23664034743236L) + ((-((int128)tmp_q[2] * 75298483926170L) + ((int128)tmp_q[3] * 40129352224870L) + ((int128)tmp_q[4] * 83613586901033L) - ((int128)tmp_q[5] * 71625886439504L) - ((int128)tmp_q[6] * 73095569950280L) + ((int128)tmp_q[7] * 16212499142000L) - ((int128)tmp_q[8] * 39120728066137L) - ((int128)tmp_q[9] * 122570424992969L) - ((int128)tmp_q[10] * 7265070031016L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 7265070031016L) + ((int128)tmp_q[1] * 41106866914776L) - ((int128)tmp_q[2] * 23664034743236L) + ((-((int128)tmp_q[3] * 75298483926170L) + ((int128)tmp_q[4] * 40129352224870L) + ((int128)tmp_q[5] * 83613586901033L) - ((int128)tmp_q[6] * 71625886439504L) - ((int128)tmp_q[7] * 73095569950280L) + ((int128)tmp_q[8] * 16212499142000L) - ((int128)tmp_q[9] * 39120728066137L) - ((int128)tmp_q[10] * 122570424992969L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 122570424992969L) - ((int128)tmp_q[1] * 7265070031016L) + ((int128)tmp_q[2] * 41106866914776L) - ((int128)tmp_q[3] * 23664034743236L) + ((-((int128)tmp_q[4] * 75298483926170L) + ((int128)tmp_q[5] * 40129352224870L) + ((int128)tmp_q[6] * 83613586901033L) - ((int128)tmp_q[7] * 71625886439504L) - ((int128)tmp_q[8] * 73095569950280L) + ((int128)tmp_q[9] * 16212499142000L) - ((int128)tmp_q[10] * 39120728066137L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 39120728066137L) - ((int128)tmp_q[1] * 122570424992969L) - ((int128)tmp_q[2] * 7265070031016L) + ((int128)tmp_q[3] * 41106866914776L) - ((int128)tmp_q[4] * 23664034743236L) + ((-((int128)tmp_q[5] * 75298483926170L) + ((int128)tmp_q[6] * 40129352224870L) + ((int128)tmp_q[7] * 83613586901033L) - ((int128)tmp_q[8] * 71625886439504L) - ((int128)tmp_q[9] * 73095569950280L) + ((int128)tmp_q[10] * 16212499142000L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 16212499142000L) - ((int128)tmp_q[1] * 39120728066137L) - ((int128)tmp_q[2] * 122570424992969L) - ((int128)tmp_q[3] * 7265070031016L) + ((int128)tmp_q[4] * 41106866914776L) - ((int128)tmp_q[5] * 23664034743236L) + ((-((int128)tmp_q[6] * 75298483926170L) + ((int128)tmp_q[7] * 40129352224870L) + ((int128)tmp_q[8] * 83613586901033L) - ((int128)tmp_q[9] * 71625886439504L) - ((int128)tmp_q[10] * 73095569950280L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 73095569950280L) + ((int128)tmp_q[1] * 16212499142000L) - ((int128)tmp_q[2] * 39120728066137L) - ((int128)tmp_q[3] * 122570424992969L) - ((int128)tmp_q[4] * 7265070031016L) + ((int128)tmp_q[5] * 41106866914776L) - ((int128)tmp_q[6] * 23664034743236L) + ((-((int128)tmp_q[7] * 75298483926170L) + ((int128)tmp_q[8] * 40129352224870L) + ((int128)tmp_q[9] * 83613586901033L) - ((int128)tmp_q[10] * 71625886439504L)) * 7);
	tmp_zero[7] = -((int128)tmp_q[0] * 71625886439504L) - ((int128)tmp_q[1] * 73095569950280L) + ((int128)tmp_q[2] * 16212499142000L) - ((int128)tmp_q[3] * 39120728066137L) - ((int128)tmp_q[4] * 122570424992969L) - ((int128)tmp_q[5] * 7265070031016L) + ((int128)tmp_q[6] * 41106866914776L) - ((int128)tmp_q[7] * 23664034743236L) + ((-((int128)tmp_q[8] * 75298483926170L) + ((int128)tmp_q[9] * 40129352224870L) + ((int128)tmp_q[10] * 83613586901033L)) * 7);
	tmp_zero[8] = ((int128)tmp_q[0] * 83613586901033L) - ((int128)tmp_q[1] * 71625886439504L) - ((int128)tmp_q[2] * 73095569950280L) + ((int128)tmp_q[3] * 16212499142000L) - ((int128)tmp_q[4] * 39120728066137L) - ((int128)tmp_q[5] * 122570424992969L) - ((int128)tmp_q[6] * 7265070031016L) + ((int128)tmp_q[7] * 41106866914776L) - ((int128)tmp_q[8] * 23664034743236L) + ((-((int128)tmp_q[9] * 75298483926170L) + ((int128)tmp_q[10] * 40129352224870L)) * 7);
	tmp_zero[9] = ((int128)tmp_q[0] * 40129352224870L) + ((int128)tmp_q[1] * 83613586901033L) - ((int128)tmp_q[2] * 71625886439504L) - ((int128)tmp_q[3] * 73095569950280L) + ((int128)tmp_q[4] * 16212499142000L) - ((int128)tmp_q[5] * 39120728066137L) - ((int128)tmp_q[6] * 122570424992969L) - ((int128)tmp_q[7] * 7265070031016L) + ((int128)tmp_q[8] * 41106866914776L) - ((int128)tmp_q[9] * 23664034743236L) - ((int128)tmp_q[10] * 527089387483190L);
	tmp_zero[10] = -((int128)tmp_q[0] * 75298483926170L) + ((int128)tmp_q[1] * 40129352224870L) + ((int128)tmp_q[2] * 83613586901033L) - ((int128)tmp_q[3] * 71625886439504L) - ((int128)tmp_q[4] * 73095569950280L) + ((int128)tmp_q[5] * 16212499142000L) - ((int128)tmp_q[6] * 39120728066137L) - ((int128)tmp_q[7] * 122570424992969L) - ((int128)tmp_q[8] * 7265070031016L) + ((int128)tmp_q[9] * 41106866914776L) - ((int128)tmp_q[10] * 23664034743236L);

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

