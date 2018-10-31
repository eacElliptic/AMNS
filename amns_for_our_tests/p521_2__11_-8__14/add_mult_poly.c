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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 4);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 4);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13351237599665602403UL) + ((((uint64_t)op[1] * 7845323693149871297UL) + ((uint64_t)op[2] * 9915429646857274653UL) + ((uint64_t)op[3] * 5375262232209213100UL) + ((uint64_t)op[4] * 18031177138167795442UL) + ((uint64_t)op[5] * 9574629338144349931UL) + ((uint64_t)op[6] * 6373181085038296657UL) + ((uint64_t)op[7] * 6320522728187506924UL) + ((uint64_t)op[8] * 16530307123512761827UL) + ((uint64_t)op[9] * 10282907103168763805UL) + ((uint64_t)op[10] * 16579370566975886809UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 16579370566975886809UL) + ((uint64_t)op[1] * 13351237599665602403UL) + ((((uint64_t)op[2] * 7845323693149871297UL) + ((uint64_t)op[3] * 9915429646857274653UL) + ((uint64_t)op[4] * 5375262232209213100UL) + ((uint64_t)op[5] * 18031177138167795442UL) + ((uint64_t)op[6] * 9574629338144349931UL) + ((uint64_t)op[7] * 6373181085038296657UL) + ((uint64_t)op[8] * 6320522728187506924UL) + ((uint64_t)op[9] * 16530307123512761827UL) + ((uint64_t)op[10] * 10282907103168763805UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 10282907103168763805UL) + ((uint64_t)op[1] * 16579370566975886809UL) + ((uint64_t)op[2] * 13351237599665602403UL) + ((((uint64_t)op[3] * 7845323693149871297UL) + ((uint64_t)op[4] * 9915429646857274653UL) + ((uint64_t)op[5] * 5375262232209213100UL) + ((uint64_t)op[6] * 18031177138167795442UL) + ((uint64_t)op[7] * 9574629338144349931UL) + ((uint64_t)op[8] * 6373181085038296657UL) + ((uint64_t)op[9] * 6320522728187506924UL) + ((uint64_t)op[10] * 16530307123512761827UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 16530307123512761827UL) + ((uint64_t)op[1] * 10282907103168763805UL) + ((uint64_t)op[2] * 16579370566975886809UL) + ((uint64_t)op[3] * 13351237599665602403UL) + ((((uint64_t)op[4] * 7845323693149871297UL) + ((uint64_t)op[5] * 9915429646857274653UL) + ((uint64_t)op[6] * 5375262232209213100UL) + ((uint64_t)op[7] * 18031177138167795442UL) + ((uint64_t)op[8] * 9574629338144349931UL) + ((uint64_t)op[9] * 6373181085038296657UL) + ((uint64_t)op[10] * 6320522728187506924UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 6320522728187506924UL) + ((uint64_t)op[1] * 16530307123512761827UL) + ((uint64_t)op[2] * 10282907103168763805UL) + ((uint64_t)op[3] * 16579370566975886809UL) + ((uint64_t)op[4] * 13351237599665602403UL) + ((((uint64_t)op[5] * 7845323693149871297UL) + ((uint64_t)op[6] * 9915429646857274653UL) + ((uint64_t)op[7] * 5375262232209213100UL) + ((uint64_t)op[8] * 18031177138167795442UL) + ((uint64_t)op[9] * 9574629338144349931UL) + ((uint64_t)op[10] * 6373181085038296657UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 6373181085038296657UL) + ((uint64_t)op[1] * 6320522728187506924UL) + ((uint64_t)op[2] * 16530307123512761827UL) + ((uint64_t)op[3] * 10282907103168763805UL) + ((uint64_t)op[4] * 16579370566975886809UL) + ((uint64_t)op[5] * 13351237599665602403UL) + ((((uint64_t)op[6] * 7845323693149871297UL) + ((uint64_t)op[7] * 9915429646857274653UL) + ((uint64_t)op[8] * 5375262232209213100UL) + ((uint64_t)op[9] * 18031177138167795442UL) + ((uint64_t)op[10] * 9574629338144349931UL)) * 18446744073709551608);
	tmp_q[6] = ((uint64_t)op[0] * 9574629338144349931UL) + ((uint64_t)op[1] * 6373181085038296657UL) + ((uint64_t)op[2] * 6320522728187506924UL) + ((uint64_t)op[3] * 16530307123512761827UL) + ((uint64_t)op[4] * 10282907103168763805UL) + ((uint64_t)op[5] * 16579370566975886809UL) + ((uint64_t)op[6] * 13351237599665602403UL) + ((((uint64_t)op[7] * 7845323693149871297UL) + ((uint64_t)op[8] * 9915429646857274653UL) + ((uint64_t)op[9] * 5375262232209213100UL) + ((uint64_t)op[10] * 18031177138167795442UL)) * 18446744073709551608);
	tmp_q[7] = ((uint64_t)op[0] * 18031177138167795442UL) + ((uint64_t)op[1] * 9574629338144349931UL) + ((uint64_t)op[2] * 6373181085038296657UL) + ((uint64_t)op[3] * 6320522728187506924UL) + ((uint64_t)op[4] * 16530307123512761827UL) + ((uint64_t)op[5] * 10282907103168763805UL) + ((uint64_t)op[6] * 16579370566975886809UL) + ((uint64_t)op[7] * 13351237599665602403UL) + ((((uint64_t)op[8] * 7845323693149871297UL) + ((uint64_t)op[9] * 9915429646857274653UL) + ((uint64_t)op[10] * 5375262232209213100UL)) * 18446744073709551608);
	tmp_q[8] = ((uint64_t)op[0] * 5375262232209213100UL) + ((uint64_t)op[1] * 18031177138167795442UL) + ((uint64_t)op[2] * 9574629338144349931UL) + ((uint64_t)op[3] * 6373181085038296657UL) + ((uint64_t)op[4] * 6320522728187506924UL) + ((uint64_t)op[5] * 16530307123512761827UL) + ((uint64_t)op[6] * 10282907103168763805UL) + ((uint64_t)op[7] * 16579370566975886809UL) + ((uint64_t)op[8] * 13351237599665602403UL) + ((((uint64_t)op[9] * 7845323693149871297UL) + ((uint64_t)op[10] * 9915429646857274653UL)) * 18446744073709551608);
	tmp_q[9] = ((uint64_t)op[0] * 9915429646857274653UL) + ((uint64_t)op[1] * 5375262232209213100UL) + ((uint64_t)op[2] * 18031177138167795442UL) + ((uint64_t)op[3] * 9574629338144349931UL) + ((uint64_t)op[4] * 6373181085038296657UL) + ((uint64_t)op[5] * 6320522728187506924UL) + ((uint64_t)op[6] * 16530307123512761827UL) + ((uint64_t)op[7] * 10282907103168763805UL) + ((uint64_t)op[8] * 16579370566975886809UL) + ((uint64_t)op[9] * 13351237599665602403UL) + ((uint64_t)op[10] * 11024386749639236088UL);
	tmp_q[10] = ((uint64_t)op[0] * 7845323693149871297UL) + ((uint64_t)op[1] * 9915429646857274653UL) + ((uint64_t)op[2] * 5375262232209213100UL) + ((uint64_t)op[3] * 18031177138167795442UL) + ((uint64_t)op[4] * 9574629338144349931UL) + ((uint64_t)op[5] * 6373181085038296657UL) + ((uint64_t)op[6] * 6320522728187506924UL) + ((uint64_t)op[7] * 16530307123512761827UL) + ((uint64_t)op[8] * 10282907103168763805UL) + ((uint64_t)op[9] * 16579370566975886809UL) + ((uint64_t)op[10] * 13351237599665602403UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 34365976605291L) - ((-((int128)tmp_q[1] * 18261437793912L) - ((int128)tmp_q[2] * 64292896159510L) + ((int128)tmp_q[3] * 84356735431328L) + ((int128)tmp_q[4] * 53830333670049L) + ((int128)tmp_q[5] * 7644529351525L) + ((int128)tmp_q[6] * 99669386618192L) - ((int128)tmp_q[7] * 50319401881621L) + ((int128)tmp_q[8] * 89940726687870L) + ((int128)tmp_q[9] * 117163052545170L) - ((int128)tmp_q[10] * 81516603527199L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 81516603527199L) - ((int128)tmp_q[1] * 34365976605291L) - ((-((int128)tmp_q[2] * 18261437793912L) - ((int128)tmp_q[3] * 64292896159510L) + ((int128)tmp_q[4] * 84356735431328L) + ((int128)tmp_q[5] * 53830333670049L) + ((int128)tmp_q[6] * 7644529351525L) + ((int128)tmp_q[7] * 99669386618192L) - ((int128)tmp_q[8] * 50319401881621L) + ((int128)tmp_q[9] * 89940726687870L) + ((int128)tmp_q[10] * 117163052545170L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 117163052545170L) - ((int128)tmp_q[1] * 81516603527199L) - ((int128)tmp_q[2] * 34365976605291L) - ((-((int128)tmp_q[3] * 18261437793912L) - ((int128)tmp_q[4] * 64292896159510L) + ((int128)tmp_q[5] * 84356735431328L) + ((int128)tmp_q[6] * 53830333670049L) + ((int128)tmp_q[7] * 7644529351525L) + ((int128)tmp_q[8] * 99669386618192L) - ((int128)tmp_q[9] * 50319401881621L) + ((int128)tmp_q[10] * 89940726687870L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 89940726687870L) + ((int128)tmp_q[1] * 117163052545170L) - ((int128)tmp_q[2] * 81516603527199L) - ((int128)tmp_q[3] * 34365976605291L) - ((-((int128)tmp_q[4] * 18261437793912L) - ((int128)tmp_q[5] * 64292896159510L) + ((int128)tmp_q[6] * 84356735431328L) + ((int128)tmp_q[7] * 53830333670049L) + ((int128)tmp_q[8] * 7644529351525L) + ((int128)tmp_q[9] * 99669386618192L) - ((int128)tmp_q[10] * 50319401881621L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 50319401881621L) + ((int128)tmp_q[1] * 89940726687870L) + ((int128)tmp_q[2] * 117163052545170L) - ((int128)tmp_q[3] * 81516603527199L) - ((int128)tmp_q[4] * 34365976605291L) - ((-((int128)tmp_q[5] * 18261437793912L) - ((int128)tmp_q[6] * 64292896159510L) + ((int128)tmp_q[7] * 84356735431328L) + ((int128)tmp_q[8] * 53830333670049L) + ((int128)tmp_q[9] * 7644529351525L) + ((int128)tmp_q[10] * 99669386618192L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 99669386618192L) - ((int128)tmp_q[1] * 50319401881621L) + ((int128)tmp_q[2] * 89940726687870L) + ((int128)tmp_q[3] * 117163052545170L) - ((int128)tmp_q[4] * 81516603527199L) - ((int128)tmp_q[5] * 34365976605291L) - ((-((int128)tmp_q[6] * 18261437793912L) - ((int128)tmp_q[7] * 64292896159510L) + ((int128)tmp_q[8] * 84356735431328L) + ((int128)tmp_q[9] * 53830333670049L) + ((int128)tmp_q[10] * 7644529351525L)) * 8);
	tmp_zero[6] = ((int128)tmp_q[0] * 7644529351525L) + ((int128)tmp_q[1] * 99669386618192L) - ((int128)tmp_q[2] * 50319401881621L) + ((int128)tmp_q[3] * 89940726687870L) + ((int128)tmp_q[4] * 117163052545170L) - ((int128)tmp_q[5] * 81516603527199L) - ((int128)tmp_q[6] * 34365976605291L) - ((-((int128)tmp_q[7] * 18261437793912L) - ((int128)tmp_q[8] * 64292896159510L) + ((int128)tmp_q[9] * 84356735431328L) + ((int128)tmp_q[10] * 53830333670049L)) * 8);
	tmp_zero[7] = ((int128)tmp_q[0] * 53830333670049L) + ((int128)tmp_q[1] * 7644529351525L) + ((int128)tmp_q[2] * 99669386618192L) - ((int128)tmp_q[3] * 50319401881621L) + ((int128)tmp_q[4] * 89940726687870L) + ((int128)tmp_q[5] * 117163052545170L) - ((int128)tmp_q[6] * 81516603527199L) - ((int128)tmp_q[7] * 34365976605291L) - ((-((int128)tmp_q[8] * 18261437793912L) - ((int128)tmp_q[9] * 64292896159510L) + ((int128)tmp_q[10] * 84356735431328L)) * 8);
	tmp_zero[8] = ((int128)tmp_q[0] * 84356735431328L) + ((int128)tmp_q[1] * 53830333670049L) + ((int128)tmp_q[2] * 7644529351525L) + ((int128)tmp_q[3] * 99669386618192L) - ((int128)tmp_q[4] * 50319401881621L) + ((int128)tmp_q[5] * 89940726687870L) + ((int128)tmp_q[6] * 117163052545170L) - ((int128)tmp_q[7] * 81516603527199L) - ((int128)tmp_q[8] * 34365976605291L) - ((-((int128)tmp_q[9] * 18261437793912L) - ((int128)tmp_q[10] * 64292896159510L)) * 8);
	tmp_zero[9] = -((int128)tmp_q[0] * 64292896159510L) + ((int128)tmp_q[1] * 84356735431328L) + ((int128)tmp_q[2] * 53830333670049L) + ((int128)tmp_q[3] * 7644529351525L) + ((int128)tmp_q[4] * 99669386618192L) - ((int128)tmp_q[5] * 50319401881621L) + ((int128)tmp_q[6] * 89940726687870L) + ((int128)tmp_q[7] * 117163052545170L) - ((int128)tmp_q[8] * 81516603527199L) - ((int128)tmp_q[9] * 34365976605291L) + ((int128)tmp_q[10] * 146091502351296L);
	tmp_zero[10] = -((int128)tmp_q[0] * 18261437793912L) - ((int128)tmp_q[1] * 64292896159510L) + ((int128)tmp_q[2] * 84356735431328L) + ((int128)tmp_q[3] * 53830333670049L) + ((int128)tmp_q[4] * 7644529351525L) + ((int128)tmp_q[5] * 99669386618192L) - ((int128)tmp_q[6] * 50319401881621L) + ((int128)tmp_q[7] * 89940726687870L) + ((int128)tmp_q[8] * 117163052545170L) - ((int128)tmp_q[9] * 81516603527199L) - ((int128)tmp_q[10] * 34365976605291L);

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

