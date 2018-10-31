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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 2);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8903936968432633433UL) + ((((uint64_t)op[1] * 5810291276659875253UL) + ((uint64_t)op[2] * 9048799147118507318UL) + ((uint64_t)op[3] * 14185911293813305429UL) + ((uint64_t)op[4] * 8721465246765925815UL) + ((uint64_t)op[5] * 8955200321095641414UL) + ((uint64_t)op[6] * 2785156593613506705UL) + ((uint64_t)op[7] * 17013342952758375884UL) + ((uint64_t)op[8] * 18436196810728490178UL) + ((uint64_t)op[9] * 17371775266185246550UL) + ((uint64_t)op[10] * 17277613494661190017UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 17277613494661190017UL) + ((uint64_t)op[1] * 8903936968432633433UL) + ((((uint64_t)op[2] * 5810291276659875253UL) + ((uint64_t)op[3] * 9048799147118507318UL) + ((uint64_t)op[4] * 14185911293813305429UL) + ((uint64_t)op[5] * 8721465246765925815UL) + ((uint64_t)op[6] * 8955200321095641414UL) + ((uint64_t)op[7] * 2785156593613506705UL) + ((uint64_t)op[8] * 17013342952758375884UL) + ((uint64_t)op[9] * 18436196810728490178UL) + ((uint64_t)op[10] * 17371775266185246550UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 17371775266185246550UL) + ((uint64_t)op[1] * 17277613494661190017UL) + ((uint64_t)op[2] * 8903936968432633433UL) + ((((uint64_t)op[3] * 5810291276659875253UL) + ((uint64_t)op[4] * 9048799147118507318UL) + ((uint64_t)op[5] * 14185911293813305429UL) + ((uint64_t)op[6] * 8721465246765925815UL) + ((uint64_t)op[7] * 8955200321095641414UL) + ((uint64_t)op[8] * 2785156593613506705UL) + ((uint64_t)op[9] * 17013342952758375884UL) + ((uint64_t)op[10] * 18436196810728490178UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 18436196810728490178UL) + ((uint64_t)op[1] * 17371775266185246550UL) + ((uint64_t)op[2] * 17277613494661190017UL) + ((uint64_t)op[3] * 8903936968432633433UL) + ((((uint64_t)op[4] * 5810291276659875253UL) + ((uint64_t)op[5] * 9048799147118507318UL) + ((uint64_t)op[6] * 14185911293813305429UL) + ((uint64_t)op[7] * 8721465246765925815UL) + ((uint64_t)op[8] * 8955200321095641414UL) + ((uint64_t)op[9] * 2785156593613506705UL) + ((uint64_t)op[10] * 17013342952758375884UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 17013342952758375884UL) + ((uint64_t)op[1] * 18436196810728490178UL) + ((uint64_t)op[2] * 17371775266185246550UL) + ((uint64_t)op[3] * 17277613494661190017UL) + ((uint64_t)op[4] * 8903936968432633433UL) + ((((uint64_t)op[5] * 5810291276659875253UL) + ((uint64_t)op[6] * 9048799147118507318UL) + ((uint64_t)op[7] * 14185911293813305429UL) + ((uint64_t)op[8] * 8721465246765925815UL) + ((uint64_t)op[9] * 8955200321095641414UL) + ((uint64_t)op[10] * 2785156593613506705UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 2785156593613506705UL) + ((uint64_t)op[1] * 17013342952758375884UL) + ((uint64_t)op[2] * 18436196810728490178UL) + ((uint64_t)op[3] * 17371775266185246550UL) + ((uint64_t)op[4] * 17277613494661190017UL) + ((uint64_t)op[5] * 8903936968432633433UL) + ((((uint64_t)op[6] * 5810291276659875253UL) + ((uint64_t)op[7] * 9048799147118507318UL) + ((uint64_t)op[8] * 14185911293813305429UL) + ((uint64_t)op[9] * 8721465246765925815UL) + ((uint64_t)op[10] * 8955200321095641414UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 8955200321095641414UL) + ((uint64_t)op[1] * 2785156593613506705UL) + ((uint64_t)op[2] * 17013342952758375884UL) + ((uint64_t)op[3] * 18436196810728490178UL) + ((uint64_t)op[4] * 17371775266185246550UL) + ((uint64_t)op[5] * 17277613494661190017UL) + ((uint64_t)op[6] * 8903936968432633433UL) + ((((uint64_t)op[7] * 5810291276659875253UL) + ((uint64_t)op[8] * 9048799147118507318UL) + ((uint64_t)op[9] * 14185911293813305429UL) + ((uint64_t)op[10] * 8721465246765925815UL)) * 2);
	tmp_q[7] = ((uint64_t)op[0] * 8721465246765925815UL) + ((uint64_t)op[1] * 8955200321095641414UL) + ((uint64_t)op[2] * 2785156593613506705UL) + ((uint64_t)op[3] * 17013342952758375884UL) + ((uint64_t)op[4] * 18436196810728490178UL) + ((uint64_t)op[5] * 17371775266185246550UL) + ((uint64_t)op[6] * 17277613494661190017UL) + ((uint64_t)op[7] * 8903936968432633433UL) + ((((uint64_t)op[8] * 5810291276659875253UL) + ((uint64_t)op[9] * 9048799147118507318UL) + ((uint64_t)op[10] * 14185911293813305429UL)) * 2);
	tmp_q[8] = ((uint64_t)op[0] * 14185911293813305429UL) + ((uint64_t)op[1] * 8721465246765925815UL) + ((uint64_t)op[2] * 8955200321095641414UL) + ((uint64_t)op[3] * 2785156593613506705UL) + ((uint64_t)op[4] * 17013342952758375884UL) + ((uint64_t)op[5] * 18436196810728490178UL) + ((uint64_t)op[6] * 17371775266185246550UL) + ((uint64_t)op[7] * 17277613494661190017UL) + ((uint64_t)op[8] * 8903936968432633433UL) + ((((uint64_t)op[9] * 5810291276659875253UL) + ((uint64_t)op[10] * 9048799147118507318UL)) * 2);
	tmp_q[9] = ((uint64_t)op[0] * 9048799147118507318UL) + ((uint64_t)op[1] * 14185911293813305429UL) + ((uint64_t)op[2] * 8721465246765925815UL) + ((uint64_t)op[3] * 8955200321095641414UL) + ((uint64_t)op[4] * 2785156593613506705UL) + ((uint64_t)op[5] * 17013342952758375884UL) + ((uint64_t)op[6] * 18436196810728490178UL) + ((uint64_t)op[7] * 17371775266185246550UL) + ((uint64_t)op[8] * 17277613494661190017UL) + ((uint64_t)op[9] * 8903936968432633433UL) + ((uint64_t)op[10] * 11620582553319750506UL);
	tmp_q[10] = ((uint64_t)op[0] * 5810291276659875253UL) + ((uint64_t)op[1] * 9048799147118507318UL) + ((uint64_t)op[2] * 14185911293813305429UL) + ((uint64_t)op[3] * 8721465246765925815UL) + ((uint64_t)op[4] * 8955200321095641414UL) + ((uint64_t)op[5] * 2785156593613506705UL) + ((uint64_t)op[6] * 17013342952758375884UL) + ((uint64_t)op[7] * 18436196810728490178UL) + ((uint64_t)op[8] * 17371775266185246550UL) + ((uint64_t)op[9] * 17277613494661190017UL) + ((uint64_t)op[10] * 8903936968432633433UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 31210372450601L) + ((((int128)tmp_q[1] * 97032295113618L) + ((int128)tmp_q[2] * 17337589542291L) - ((int128)tmp_q[3] * 9297891094844L) + ((int128)tmp_q[4] * 73600006440081L) - ((int128)tmp_q[5] * 19220109595585L) + ((int128)tmp_q[6] * 23177491134596L) + ((int128)tmp_q[7] * 76975980984217L) - ((int128)tmp_q[8] * 34875575895355L) - ((int128)tmp_q[9] * 84796780952641L) - ((int128)tmp_q[10] * 25250884237533L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 25250884237533L) - ((int128)tmp_q[1] * 31210372450601L) + ((((int128)tmp_q[2] * 97032295113618L) + ((int128)tmp_q[3] * 17337589542291L) - ((int128)tmp_q[4] * 9297891094844L) + ((int128)tmp_q[5] * 73600006440081L) - ((int128)tmp_q[6] * 19220109595585L) + ((int128)tmp_q[7] * 23177491134596L) + ((int128)tmp_q[8] * 76975980984217L) - ((int128)tmp_q[9] * 34875575895355L) - ((int128)tmp_q[10] * 84796780952641L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 84796780952641L) - ((int128)tmp_q[1] * 25250884237533L) - ((int128)tmp_q[2] * 31210372450601L) + ((((int128)tmp_q[3] * 97032295113618L) + ((int128)tmp_q[4] * 17337589542291L) - ((int128)tmp_q[5] * 9297891094844L) + ((int128)tmp_q[6] * 73600006440081L) - ((int128)tmp_q[7] * 19220109595585L) + ((int128)tmp_q[8] * 23177491134596L) + ((int128)tmp_q[9] * 76975980984217L) - ((int128)tmp_q[10] * 34875575895355L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 34875575895355L) - ((int128)tmp_q[1] * 84796780952641L) - ((int128)tmp_q[2] * 25250884237533L) - ((int128)tmp_q[3] * 31210372450601L) + ((((int128)tmp_q[4] * 97032295113618L) + ((int128)tmp_q[5] * 17337589542291L) - ((int128)tmp_q[6] * 9297891094844L) + ((int128)tmp_q[7] * 73600006440081L) - ((int128)tmp_q[8] * 19220109595585L) + ((int128)tmp_q[9] * 23177491134596L) + ((int128)tmp_q[10] * 76975980984217L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 76975980984217L) - ((int128)tmp_q[1] * 34875575895355L) - ((int128)tmp_q[2] * 84796780952641L) - ((int128)tmp_q[3] * 25250884237533L) - ((int128)tmp_q[4] * 31210372450601L) + ((((int128)tmp_q[5] * 97032295113618L) + ((int128)tmp_q[6] * 17337589542291L) - ((int128)tmp_q[7] * 9297891094844L) + ((int128)tmp_q[8] * 73600006440081L) - ((int128)tmp_q[9] * 19220109595585L) + ((int128)tmp_q[10] * 23177491134596L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 23177491134596L) + ((int128)tmp_q[1] * 76975980984217L) - ((int128)tmp_q[2] * 34875575895355L) - ((int128)tmp_q[3] * 84796780952641L) - ((int128)tmp_q[4] * 25250884237533L) - ((int128)tmp_q[5] * 31210372450601L) + ((((int128)tmp_q[6] * 97032295113618L) + ((int128)tmp_q[7] * 17337589542291L) - ((int128)tmp_q[8] * 9297891094844L) + ((int128)tmp_q[9] * 73600006440081L) - ((int128)tmp_q[10] * 19220109595585L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 19220109595585L) + ((int128)tmp_q[1] * 23177491134596L) + ((int128)tmp_q[2] * 76975980984217L) - ((int128)tmp_q[3] * 34875575895355L) - ((int128)tmp_q[4] * 84796780952641L) - ((int128)tmp_q[5] * 25250884237533L) - ((int128)tmp_q[6] * 31210372450601L) + ((((int128)tmp_q[7] * 97032295113618L) + ((int128)tmp_q[8] * 17337589542291L) - ((int128)tmp_q[9] * 9297891094844L) + ((int128)tmp_q[10] * 73600006440081L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 73600006440081L) - ((int128)tmp_q[1] * 19220109595585L) + ((int128)tmp_q[2] * 23177491134596L) + ((int128)tmp_q[3] * 76975980984217L) - ((int128)tmp_q[4] * 34875575895355L) - ((int128)tmp_q[5] * 84796780952641L) - ((int128)tmp_q[6] * 25250884237533L) - ((int128)tmp_q[7] * 31210372450601L) + ((((int128)tmp_q[8] * 97032295113618L) + ((int128)tmp_q[9] * 17337589542291L) - ((int128)tmp_q[10] * 9297891094844L)) * 2);
	tmp_zero[8] = -((int128)tmp_q[0] * 9297891094844L) + ((int128)tmp_q[1] * 73600006440081L) - ((int128)tmp_q[2] * 19220109595585L) + ((int128)tmp_q[3] * 23177491134596L) + ((int128)tmp_q[4] * 76975980984217L) - ((int128)tmp_q[5] * 34875575895355L) - ((int128)tmp_q[6] * 84796780952641L) - ((int128)tmp_q[7] * 25250884237533L) - ((int128)tmp_q[8] * 31210372450601L) + ((((int128)tmp_q[9] * 97032295113618L) + ((int128)tmp_q[10] * 17337589542291L)) * 2);
	tmp_zero[9] = ((int128)tmp_q[0] * 17337589542291L) - ((int128)tmp_q[1] * 9297891094844L) + ((int128)tmp_q[2] * 73600006440081L) - ((int128)tmp_q[3] * 19220109595585L) + ((int128)tmp_q[4] * 23177491134596L) + ((int128)tmp_q[5] * 76975980984217L) - ((int128)tmp_q[6] * 34875575895355L) - ((int128)tmp_q[7] * 84796780952641L) - ((int128)tmp_q[8] * 25250884237533L) - ((int128)tmp_q[9] * 31210372450601L) + ((int128)tmp_q[10] * 194064590227236L);
	tmp_zero[10] = ((int128)tmp_q[0] * 97032295113618L) + ((int128)tmp_q[1] * 17337589542291L) - ((int128)tmp_q[2] * 9297891094844L) + ((int128)tmp_q[3] * 73600006440081L) - ((int128)tmp_q[4] * 19220109595585L) + ((int128)tmp_q[5] * 23177491134596L) + ((int128)tmp_q[6] * 76975980984217L) - ((int128)tmp_q[7] * 34875575895355L) - ((int128)tmp_q[8] * 84796780952641L) - ((int128)tmp_q[9] * 25250884237533L) - ((int128)tmp_q[10] * 31210372450601L);

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

