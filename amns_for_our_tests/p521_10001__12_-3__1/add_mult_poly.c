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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 3);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 3);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7583782571699760194UL) + ((((uint64_t)op[1] * 7268319623568262330UL) + ((uint64_t)op[2] * 15160106442157393186UL) + ((uint64_t)op[3] * 12957797612956778390UL) + ((uint64_t)op[4] * 6455467546671129108UL) + ((uint64_t)op[5] * 1560317081909984160UL) + ((uint64_t)op[6] * 843625904122372817UL) + ((uint64_t)op[7] * 17403553707027694272UL) + ((uint64_t)op[8] * 7596478329713939526UL) + ((uint64_t)op[9] * 9365511004969277521UL) + ((uint64_t)op[10] * 14358077871891965259UL) + ((uint64_t)op[11] * 12020541618962656418UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 12020541618962656418UL) + ((uint64_t)op[1] * 7583782571699760194UL) + ((((uint64_t)op[2] * 7268319623568262330UL) + ((uint64_t)op[3] * 15160106442157393186UL) + ((uint64_t)op[4] * 12957797612956778390UL) + ((uint64_t)op[5] * 6455467546671129108UL) + ((uint64_t)op[6] * 1560317081909984160UL) + ((uint64_t)op[7] * 843625904122372817UL) + ((uint64_t)op[8] * 17403553707027694272UL) + ((uint64_t)op[9] * 7596478329713939526UL) + ((uint64_t)op[10] * 9365511004969277521UL) + ((uint64_t)op[11] * 14358077871891965259UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 14358077871891965259UL) + ((uint64_t)op[1] * 12020541618962656418UL) + ((uint64_t)op[2] * 7583782571699760194UL) + ((((uint64_t)op[3] * 7268319623568262330UL) + ((uint64_t)op[4] * 15160106442157393186UL) + ((uint64_t)op[5] * 12957797612956778390UL) + ((uint64_t)op[6] * 6455467546671129108UL) + ((uint64_t)op[7] * 1560317081909984160UL) + ((uint64_t)op[8] * 843625904122372817UL) + ((uint64_t)op[9] * 17403553707027694272UL) + ((uint64_t)op[10] * 7596478329713939526UL) + ((uint64_t)op[11] * 9365511004969277521UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 9365511004969277521UL) + ((uint64_t)op[1] * 14358077871891965259UL) + ((uint64_t)op[2] * 12020541618962656418UL) + ((uint64_t)op[3] * 7583782571699760194UL) + ((((uint64_t)op[4] * 7268319623568262330UL) + ((uint64_t)op[5] * 15160106442157393186UL) + ((uint64_t)op[6] * 12957797612956778390UL) + ((uint64_t)op[7] * 6455467546671129108UL) + ((uint64_t)op[8] * 1560317081909984160UL) + ((uint64_t)op[9] * 843625904122372817UL) + ((uint64_t)op[10] * 17403553707027694272UL) + ((uint64_t)op[11] * 7596478329713939526UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 7596478329713939526UL) + ((uint64_t)op[1] * 9365511004969277521UL) + ((uint64_t)op[2] * 14358077871891965259UL) + ((uint64_t)op[3] * 12020541618962656418UL) + ((uint64_t)op[4] * 7583782571699760194UL) + ((((uint64_t)op[5] * 7268319623568262330UL) + ((uint64_t)op[6] * 15160106442157393186UL) + ((uint64_t)op[7] * 12957797612956778390UL) + ((uint64_t)op[8] * 6455467546671129108UL) + ((uint64_t)op[9] * 1560317081909984160UL) + ((uint64_t)op[10] * 843625904122372817UL) + ((uint64_t)op[11] * 17403553707027694272UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 17403553707027694272UL) + ((uint64_t)op[1] * 7596478329713939526UL) + ((uint64_t)op[2] * 9365511004969277521UL) + ((uint64_t)op[3] * 14358077871891965259UL) + ((uint64_t)op[4] * 12020541618962656418UL) + ((uint64_t)op[5] * 7583782571699760194UL) + ((((uint64_t)op[6] * 7268319623568262330UL) + ((uint64_t)op[7] * 15160106442157393186UL) + ((uint64_t)op[8] * 12957797612956778390UL) + ((uint64_t)op[9] * 6455467546671129108UL) + ((uint64_t)op[10] * 1560317081909984160UL) + ((uint64_t)op[11] * 843625904122372817UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 843625904122372817UL) + ((uint64_t)op[1] * 17403553707027694272UL) + ((uint64_t)op[2] * 7596478329713939526UL) + ((uint64_t)op[3] * 9365511004969277521UL) + ((uint64_t)op[4] * 14358077871891965259UL) + ((uint64_t)op[5] * 12020541618962656418UL) + ((uint64_t)op[6] * 7583782571699760194UL) + ((((uint64_t)op[7] * 7268319623568262330UL) + ((uint64_t)op[8] * 15160106442157393186UL) + ((uint64_t)op[9] * 12957797612956778390UL) + ((uint64_t)op[10] * 6455467546671129108UL) + ((uint64_t)op[11] * 1560317081909984160UL)) * 18446744073709551613);
	tmp_q[7] = ((uint64_t)op[0] * 1560317081909984160UL) + ((uint64_t)op[1] * 843625904122372817UL) + ((uint64_t)op[2] * 17403553707027694272UL) + ((uint64_t)op[3] * 7596478329713939526UL) + ((uint64_t)op[4] * 9365511004969277521UL) + ((uint64_t)op[5] * 14358077871891965259UL) + ((uint64_t)op[6] * 12020541618962656418UL) + ((uint64_t)op[7] * 7583782571699760194UL) + ((((uint64_t)op[8] * 7268319623568262330UL) + ((uint64_t)op[9] * 15160106442157393186UL) + ((uint64_t)op[10] * 12957797612956778390UL) + ((uint64_t)op[11] * 6455467546671129108UL)) * 18446744073709551613);
	tmp_q[8] = ((uint64_t)op[0] * 6455467546671129108UL) + ((uint64_t)op[1] * 1560317081909984160UL) + ((uint64_t)op[2] * 843625904122372817UL) + ((uint64_t)op[3] * 17403553707027694272UL) + ((uint64_t)op[4] * 7596478329713939526UL) + ((uint64_t)op[5] * 9365511004969277521UL) + ((uint64_t)op[6] * 14358077871891965259UL) + ((uint64_t)op[7] * 12020541618962656418UL) + ((uint64_t)op[8] * 7583782571699760194UL) + ((((uint64_t)op[9] * 7268319623568262330UL) + ((uint64_t)op[10] * 15160106442157393186UL) + ((uint64_t)op[11] * 12957797612956778390UL)) * 18446744073709551613);
	tmp_q[9] = ((uint64_t)op[0] * 12957797612956778390UL) + ((uint64_t)op[1] * 6455467546671129108UL) + ((uint64_t)op[2] * 1560317081909984160UL) + ((uint64_t)op[3] * 843625904122372817UL) + ((uint64_t)op[4] * 17403553707027694272UL) + ((uint64_t)op[5] * 7596478329713939526UL) + ((uint64_t)op[6] * 9365511004969277521UL) + ((uint64_t)op[7] * 14358077871891965259UL) + ((uint64_t)op[8] * 12020541618962656418UL) + ((uint64_t)op[9] * 7583782571699760194UL) + ((((uint64_t)op[10] * 7268319623568262330UL) + ((uint64_t)op[11] * 15160106442157393186UL)) * 18446744073709551613);
	tmp_q[10] = ((uint64_t)op[0] * 15160106442157393186UL) + ((uint64_t)op[1] * 12957797612956778390UL) + ((uint64_t)op[2] * 6455467546671129108UL) + ((uint64_t)op[3] * 1560317081909984160UL) + ((uint64_t)op[4] * 843625904122372817UL) + ((uint64_t)op[5] * 17403553707027694272UL) + ((uint64_t)op[6] * 7596478329713939526UL) + ((uint64_t)op[7] * 9365511004969277521UL) + ((uint64_t)op[8] * 14358077871891965259UL) + ((uint64_t)op[9] * 12020541618962656418UL) + ((uint64_t)op[10] * 7583782571699760194UL) + ((uint64_t)op[11] * 15088529276714316242UL);
	tmp_q[11] = ((uint64_t)op[0] * 7268319623568262330UL) + ((uint64_t)op[1] * 15160106442157393186UL) + ((uint64_t)op[2] * 12957797612956778390UL) + ((uint64_t)op[3] * 6455467546671129108UL) + ((uint64_t)op[4] * 1560317081909984160UL) + ((uint64_t)op[5] * 843625904122372817UL) + ((uint64_t)op[6] * 17403553707027694272UL) + ((uint64_t)op[7] * 7596478329713939526UL) + ((uint64_t)op[8] * 9365511004969277521UL) + ((uint64_t)op[9] * 14358077871891965259UL) + ((uint64_t)op[10] * 12020541618962656418UL) + ((uint64_t)op[11] * 7583782571699760194UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1968618642035L) - ((-((int128)tmp_q[1] * 4904829847845L) + ((int128)tmp_q[2] * 1315854976020L) + ((int128)tmp_q[3] * 1175337582122L) - ((int128)tmp_q[4] * 2783542768572L) - ((int128)tmp_q[5] * 4116825652081L) + ((int128)tmp_q[6] * 4441010233241L) - ((int128)tmp_q[7] * 7039244533820L) - ((int128)tmp_q[8] * 4635694603419L) - ((int128)tmp_q[9] * 4904218397858L) - ((int128)tmp_q[10] * 4365856906325L) + ((int128)tmp_q[11] * 559136315739L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 559136315739L) - ((int128)tmp_q[1] * 1968618642035L) - ((-((int128)tmp_q[2] * 4904829847845L) + ((int128)tmp_q[3] * 1315854976020L) + ((int128)tmp_q[4] * 1175337582122L) - ((int128)tmp_q[5] * 2783542768572L) - ((int128)tmp_q[6] * 4116825652081L) + ((int128)tmp_q[7] * 4441010233241L) - ((int128)tmp_q[8] * 7039244533820L) - ((int128)tmp_q[9] * 4635694603419L) - ((int128)tmp_q[10] * 4904218397858L) - ((int128)tmp_q[11] * 4365856906325L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 4365856906325L) + ((int128)tmp_q[1] * 559136315739L) - ((int128)tmp_q[2] * 1968618642035L) - ((-((int128)tmp_q[3] * 4904829847845L) + ((int128)tmp_q[4] * 1315854976020L) + ((int128)tmp_q[5] * 1175337582122L) - ((int128)tmp_q[6] * 2783542768572L) - ((int128)tmp_q[7] * 4116825652081L) + ((int128)tmp_q[8] * 4441010233241L) - ((int128)tmp_q[9] * 7039244533820L) - ((int128)tmp_q[10] * 4635694603419L) - ((int128)tmp_q[11] * 4904218397858L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 4904218397858L) - ((int128)tmp_q[1] * 4365856906325L) + ((int128)tmp_q[2] * 559136315739L) - ((int128)tmp_q[3] * 1968618642035L) - ((-((int128)tmp_q[4] * 4904829847845L) + ((int128)tmp_q[5] * 1315854976020L) + ((int128)tmp_q[6] * 1175337582122L) - ((int128)tmp_q[7] * 2783542768572L) - ((int128)tmp_q[8] * 4116825652081L) + ((int128)tmp_q[9] * 4441010233241L) - ((int128)tmp_q[10] * 7039244533820L) - ((int128)tmp_q[11] * 4635694603419L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 4635694603419L) - ((int128)tmp_q[1] * 4904218397858L) - ((int128)tmp_q[2] * 4365856906325L) + ((int128)tmp_q[3] * 559136315739L) - ((int128)tmp_q[4] * 1968618642035L) - ((-((int128)tmp_q[5] * 4904829847845L) + ((int128)tmp_q[6] * 1315854976020L) + ((int128)tmp_q[7] * 1175337582122L) - ((int128)tmp_q[8] * 2783542768572L) - ((int128)tmp_q[9] * 4116825652081L) + ((int128)tmp_q[10] * 4441010233241L) - ((int128)tmp_q[11] * 7039244533820L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 7039244533820L) - ((int128)tmp_q[1] * 4635694603419L) - ((int128)tmp_q[2] * 4904218397858L) - ((int128)tmp_q[3] * 4365856906325L) + ((int128)tmp_q[4] * 559136315739L) - ((int128)tmp_q[5] * 1968618642035L) - ((-((int128)tmp_q[6] * 4904829847845L) + ((int128)tmp_q[7] * 1315854976020L) + ((int128)tmp_q[8] * 1175337582122L) - ((int128)tmp_q[9] * 2783542768572L) - ((int128)tmp_q[10] * 4116825652081L) + ((int128)tmp_q[11] * 4441010233241L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 4441010233241L) - ((int128)tmp_q[1] * 7039244533820L) - ((int128)tmp_q[2] * 4635694603419L) - ((int128)tmp_q[3] * 4904218397858L) - ((int128)tmp_q[4] * 4365856906325L) + ((int128)tmp_q[5] * 559136315739L) - ((int128)tmp_q[6] * 1968618642035L) - ((-((int128)tmp_q[7] * 4904829847845L) + ((int128)tmp_q[8] * 1315854976020L) + ((int128)tmp_q[9] * 1175337582122L) - ((int128)tmp_q[10] * 2783542768572L) - ((int128)tmp_q[11] * 4116825652081L)) * 3);
	tmp_zero[7] = -((int128)tmp_q[0] * 4116825652081L) + ((int128)tmp_q[1] * 4441010233241L) - ((int128)tmp_q[2] * 7039244533820L) - ((int128)tmp_q[3] * 4635694603419L) - ((int128)tmp_q[4] * 4904218397858L) - ((int128)tmp_q[5] * 4365856906325L) + ((int128)tmp_q[6] * 559136315739L) - ((int128)tmp_q[7] * 1968618642035L) - ((-((int128)tmp_q[8] * 4904829847845L) + ((int128)tmp_q[9] * 1315854976020L) + ((int128)tmp_q[10] * 1175337582122L) - ((int128)tmp_q[11] * 2783542768572L)) * 3);
	tmp_zero[8] = -((int128)tmp_q[0] * 2783542768572L) - ((int128)tmp_q[1] * 4116825652081L) + ((int128)tmp_q[2] * 4441010233241L) - ((int128)tmp_q[3] * 7039244533820L) - ((int128)tmp_q[4] * 4635694603419L) - ((int128)tmp_q[5] * 4904218397858L) - ((int128)tmp_q[6] * 4365856906325L) + ((int128)tmp_q[7] * 559136315739L) - ((int128)tmp_q[8] * 1968618642035L) - ((-((int128)tmp_q[9] * 4904829847845L) + ((int128)tmp_q[10] * 1315854976020L) + ((int128)tmp_q[11] * 1175337582122L)) * 3);
	tmp_zero[9] = ((int128)tmp_q[0] * 1175337582122L) - ((int128)tmp_q[1] * 2783542768572L) - ((int128)tmp_q[2] * 4116825652081L) + ((int128)tmp_q[3] * 4441010233241L) - ((int128)tmp_q[4] * 7039244533820L) - ((int128)tmp_q[5] * 4635694603419L) - ((int128)tmp_q[6] * 4904218397858L) - ((int128)tmp_q[7] * 4365856906325L) + ((int128)tmp_q[8] * 559136315739L) - ((int128)tmp_q[9] * 1968618642035L) - ((-((int128)tmp_q[10] * 4904829847845L) + ((int128)tmp_q[11] * 1315854976020L)) * 3);
	tmp_zero[10] = ((int128)tmp_q[0] * 1315854976020L) + ((int128)tmp_q[1] * 1175337582122L) - ((int128)tmp_q[2] * 2783542768572L) - ((int128)tmp_q[3] * 4116825652081L) + ((int128)tmp_q[4] * 4441010233241L) - ((int128)tmp_q[5] * 7039244533820L) - ((int128)tmp_q[6] * 4635694603419L) - ((int128)tmp_q[7] * 4904218397858L) - ((int128)tmp_q[8] * 4365856906325L) + ((int128)tmp_q[9] * 559136315739L) - ((int128)tmp_q[10] * 1968618642035L) + ((int128)tmp_q[11] * 14714489543535L);
	tmp_zero[11] = -((int128)tmp_q[0] * 4904829847845L) + ((int128)tmp_q[1] * 1315854976020L) + ((int128)tmp_q[2] * 1175337582122L) - ((int128)tmp_q[3] * 2783542768572L) - ((int128)tmp_q[4] * 4116825652081L) + ((int128)tmp_q[5] * 4441010233241L) - ((int128)tmp_q[6] * 7039244533820L) - ((int128)tmp_q[7] * 4635694603419L) - ((int128)tmp_q[8] * 4904218397858L) - ((int128)tmp_q[9] * 4365856906325L) + ((int128)tmp_q[10] * 559136315739L) - ((int128)tmp_q[11] * 1968618642035L);

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

