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
	tmp_q[0] = ((uint64_t)op[0] * 7864426315193501600UL) + ((((uint64_t)op[1] * 551922171309873464UL) + ((uint64_t)op[2] * 1444030011803709839UL) + ((uint64_t)op[3] * 10697592599300877827UL) + ((uint64_t)op[4] * 2986680434646315904UL) + ((uint64_t)op[5] * 3352462246150110384UL) + ((uint64_t)op[6] * 14288431535442607209UL) + ((uint64_t)op[7] * 15761886103529735425UL) + ((uint64_t)op[8] * 7727674483658416977UL) + ((uint64_t)op[9] * 1231326365958178985UL) + ((uint64_t)op[10] * 10309490330565444625UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 10309490330565444625UL) + ((uint64_t)op[1] * 7864426315193501600UL) + ((((uint64_t)op[2] * 551922171309873464UL) + ((uint64_t)op[3] * 1444030011803709839UL) + ((uint64_t)op[4] * 10697592599300877827UL) + ((uint64_t)op[5] * 2986680434646315904UL) + ((uint64_t)op[6] * 3352462246150110384UL) + ((uint64_t)op[7] * 14288431535442607209UL) + ((uint64_t)op[8] * 15761886103529735425UL) + ((uint64_t)op[9] * 7727674483658416977UL) + ((uint64_t)op[10] * 1231326365958178985UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 1231326365958178985UL) + ((uint64_t)op[1] * 10309490330565444625UL) + ((uint64_t)op[2] * 7864426315193501600UL) + ((((uint64_t)op[3] * 551922171309873464UL) + ((uint64_t)op[4] * 1444030011803709839UL) + ((uint64_t)op[5] * 10697592599300877827UL) + ((uint64_t)op[6] * 2986680434646315904UL) + ((uint64_t)op[7] * 3352462246150110384UL) + ((uint64_t)op[8] * 14288431535442607209UL) + ((uint64_t)op[9] * 15761886103529735425UL) + ((uint64_t)op[10] * 7727674483658416977UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 7727674483658416977UL) + ((uint64_t)op[1] * 1231326365958178985UL) + ((uint64_t)op[2] * 10309490330565444625UL) + ((uint64_t)op[3] * 7864426315193501600UL) + ((((uint64_t)op[4] * 551922171309873464UL) + ((uint64_t)op[5] * 1444030011803709839UL) + ((uint64_t)op[6] * 10697592599300877827UL) + ((uint64_t)op[7] * 2986680434646315904UL) + ((uint64_t)op[8] * 3352462246150110384UL) + ((uint64_t)op[9] * 14288431535442607209UL) + ((uint64_t)op[10] * 15761886103529735425UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 15761886103529735425UL) + ((uint64_t)op[1] * 7727674483658416977UL) + ((uint64_t)op[2] * 1231326365958178985UL) + ((uint64_t)op[3] * 10309490330565444625UL) + ((uint64_t)op[4] * 7864426315193501600UL) + ((((uint64_t)op[5] * 551922171309873464UL) + ((uint64_t)op[6] * 1444030011803709839UL) + ((uint64_t)op[7] * 10697592599300877827UL) + ((uint64_t)op[8] * 2986680434646315904UL) + ((uint64_t)op[9] * 3352462246150110384UL) + ((uint64_t)op[10] * 14288431535442607209UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 14288431535442607209UL) + ((uint64_t)op[1] * 15761886103529735425UL) + ((uint64_t)op[2] * 7727674483658416977UL) + ((uint64_t)op[3] * 1231326365958178985UL) + ((uint64_t)op[4] * 10309490330565444625UL) + ((uint64_t)op[5] * 7864426315193501600UL) + ((((uint64_t)op[6] * 551922171309873464UL) + ((uint64_t)op[7] * 1444030011803709839UL) + ((uint64_t)op[8] * 10697592599300877827UL) + ((uint64_t)op[9] * 2986680434646315904UL) + ((uint64_t)op[10] * 3352462246150110384UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 3352462246150110384UL) + ((uint64_t)op[1] * 14288431535442607209UL) + ((uint64_t)op[2] * 15761886103529735425UL) + ((uint64_t)op[3] * 7727674483658416977UL) + ((uint64_t)op[4] * 1231326365958178985UL) + ((uint64_t)op[5] * 10309490330565444625UL) + ((uint64_t)op[6] * 7864426315193501600UL) + ((((uint64_t)op[7] * 551922171309873464UL) + ((uint64_t)op[8] * 1444030011803709839UL) + ((uint64_t)op[9] * 10697592599300877827UL) + ((uint64_t)op[10] * 2986680434646315904UL)) * 7);
	tmp_q[7] = ((uint64_t)op[0] * 2986680434646315904UL) + ((uint64_t)op[1] * 3352462246150110384UL) + ((uint64_t)op[2] * 14288431535442607209UL) + ((uint64_t)op[3] * 15761886103529735425UL) + ((uint64_t)op[4] * 7727674483658416977UL) + ((uint64_t)op[5] * 1231326365958178985UL) + ((uint64_t)op[6] * 10309490330565444625UL) + ((uint64_t)op[7] * 7864426315193501600UL) + ((((uint64_t)op[8] * 551922171309873464UL) + ((uint64_t)op[9] * 1444030011803709839UL) + ((uint64_t)op[10] * 10697592599300877827UL)) * 7);
	tmp_q[8] = ((uint64_t)op[0] * 10697592599300877827UL) + ((uint64_t)op[1] * 2986680434646315904UL) + ((uint64_t)op[2] * 3352462246150110384UL) + ((uint64_t)op[3] * 14288431535442607209UL) + ((uint64_t)op[4] * 15761886103529735425UL) + ((uint64_t)op[5] * 7727674483658416977UL) + ((uint64_t)op[6] * 1231326365958178985UL) + ((uint64_t)op[7] * 10309490330565444625UL) + ((uint64_t)op[8] * 7864426315193501600UL) + ((((uint64_t)op[9] * 551922171309873464UL) + ((uint64_t)op[10] * 1444030011803709839UL)) * 7);
	tmp_q[9] = ((uint64_t)op[0] * 1444030011803709839UL) + ((uint64_t)op[1] * 10697592599300877827UL) + ((uint64_t)op[2] * 2986680434646315904UL) + ((uint64_t)op[3] * 3352462246150110384UL) + ((uint64_t)op[4] * 14288431535442607209UL) + ((uint64_t)op[5] * 15761886103529735425UL) + ((uint64_t)op[6] * 7727674483658416977UL) + ((uint64_t)op[7] * 1231326365958178985UL) + ((uint64_t)op[8] * 10309490330565444625UL) + ((uint64_t)op[9] * 7864426315193501600UL) + ((uint64_t)op[10] * 3863455199169114248UL);
	tmp_q[10] = ((uint64_t)op[0] * 551922171309873464UL) + ((uint64_t)op[1] * 1444030011803709839UL) + ((uint64_t)op[2] * 10697592599300877827UL) + ((uint64_t)op[3] * 2986680434646315904UL) + ((uint64_t)op[4] * 3352462246150110384UL) + ((uint64_t)op[5] * 14288431535442607209UL) + ((uint64_t)op[6] * 15761886103529735425UL) + ((uint64_t)op[7] * 7727674483658416977UL) + ((uint64_t)op[8] * 1231326365958178985UL) + ((uint64_t)op[9] * 10309490330565444625UL) + ((uint64_t)op[10] * 7864426315193501600UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 43448070718530L) + ((((int128)tmp_q[1] * 71150873568084L) + ((int128)tmp_q[2] * 27986509043603L) + ((int128)tmp_q[3] * 77979715585343L) + ((int128)tmp_q[4] * 53332178474619L) - ((int128)tmp_q[5] * 2540940244994L) + ((int128)tmp_q[6] * 45439720877410L) + ((int128)tmp_q[7] * 52372751634947L) + ((int128)tmp_q[8] * 32015356200528L) + ((int128)tmp_q[9] * 82302910548966L) - ((int128)tmp_q[10] * 53972515285225L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 53972515285225L) + ((int128)tmp_q[1] * 43448070718530L) + ((((int128)tmp_q[2] * 71150873568084L) + ((int128)tmp_q[3] * 27986509043603L) + ((int128)tmp_q[4] * 77979715585343L) + ((int128)tmp_q[5] * 53332178474619L) - ((int128)tmp_q[6] * 2540940244994L) + ((int128)tmp_q[7] * 45439720877410L) + ((int128)tmp_q[8] * 52372751634947L) + ((int128)tmp_q[9] * 32015356200528L) + ((int128)tmp_q[10] * 82302910548966L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 82302910548966L) - ((int128)tmp_q[1] * 53972515285225L) + ((int128)tmp_q[2] * 43448070718530L) + ((((int128)tmp_q[3] * 71150873568084L) + ((int128)tmp_q[4] * 27986509043603L) + ((int128)tmp_q[5] * 77979715585343L) + ((int128)tmp_q[6] * 53332178474619L) - ((int128)tmp_q[7] * 2540940244994L) + ((int128)tmp_q[8] * 45439720877410L) + ((int128)tmp_q[9] * 52372751634947L) + ((int128)tmp_q[10] * 32015356200528L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 32015356200528L) + ((int128)tmp_q[1] * 82302910548966L) - ((int128)tmp_q[2] * 53972515285225L) + ((int128)tmp_q[3] * 43448070718530L) + ((((int128)tmp_q[4] * 71150873568084L) + ((int128)tmp_q[5] * 27986509043603L) + ((int128)tmp_q[6] * 77979715585343L) + ((int128)tmp_q[7] * 53332178474619L) - ((int128)tmp_q[8] * 2540940244994L) + ((int128)tmp_q[9] * 45439720877410L) + ((int128)tmp_q[10] * 52372751634947L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 52372751634947L) + ((int128)tmp_q[1] * 32015356200528L) + ((int128)tmp_q[2] * 82302910548966L) - ((int128)tmp_q[3] * 53972515285225L) + ((int128)tmp_q[4] * 43448070718530L) + ((((int128)tmp_q[5] * 71150873568084L) + ((int128)tmp_q[6] * 27986509043603L) + ((int128)tmp_q[7] * 77979715585343L) + ((int128)tmp_q[8] * 53332178474619L) - ((int128)tmp_q[9] * 2540940244994L) + ((int128)tmp_q[10] * 45439720877410L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 45439720877410L) + ((int128)tmp_q[1] * 52372751634947L) + ((int128)tmp_q[2] * 32015356200528L) + ((int128)tmp_q[3] * 82302910548966L) - ((int128)tmp_q[4] * 53972515285225L) + ((int128)tmp_q[5] * 43448070718530L) + ((((int128)tmp_q[6] * 71150873568084L) + ((int128)tmp_q[7] * 27986509043603L) + ((int128)tmp_q[8] * 77979715585343L) + ((int128)tmp_q[9] * 53332178474619L) - ((int128)tmp_q[10] * 2540940244994L)) * 7);
	tmp_zero[6] = -((int128)tmp_q[0] * 2540940244994L) + ((int128)tmp_q[1] * 45439720877410L) + ((int128)tmp_q[2] * 52372751634947L) + ((int128)tmp_q[3] * 32015356200528L) + ((int128)tmp_q[4] * 82302910548966L) - ((int128)tmp_q[5] * 53972515285225L) + ((int128)tmp_q[6] * 43448070718530L) + ((((int128)tmp_q[7] * 71150873568084L) + ((int128)tmp_q[8] * 27986509043603L) + ((int128)tmp_q[9] * 77979715585343L) + ((int128)tmp_q[10] * 53332178474619L)) * 7);
	tmp_zero[7] = ((int128)tmp_q[0] * 53332178474619L) - ((int128)tmp_q[1] * 2540940244994L) + ((int128)tmp_q[2] * 45439720877410L) + ((int128)tmp_q[3] * 52372751634947L) + ((int128)tmp_q[4] * 32015356200528L) + ((int128)tmp_q[5] * 82302910548966L) - ((int128)tmp_q[6] * 53972515285225L) + ((int128)tmp_q[7] * 43448070718530L) + ((((int128)tmp_q[8] * 71150873568084L) + ((int128)tmp_q[9] * 27986509043603L) + ((int128)tmp_q[10] * 77979715585343L)) * 7);
	tmp_zero[8] = ((int128)tmp_q[0] * 77979715585343L) + ((int128)tmp_q[1] * 53332178474619L) - ((int128)tmp_q[2] * 2540940244994L) + ((int128)tmp_q[3] * 45439720877410L) + ((int128)tmp_q[4] * 52372751634947L) + ((int128)tmp_q[5] * 32015356200528L) + ((int128)tmp_q[6] * 82302910548966L) - ((int128)tmp_q[7] * 53972515285225L) + ((int128)tmp_q[8] * 43448070718530L) + ((((int128)tmp_q[9] * 71150873568084L) + ((int128)tmp_q[10] * 27986509043603L)) * 7);
	tmp_zero[9] = ((int128)tmp_q[0] * 27986509043603L) + ((int128)tmp_q[1] * 77979715585343L) + ((int128)tmp_q[2] * 53332178474619L) - ((int128)tmp_q[3] * 2540940244994L) + ((int128)tmp_q[4] * 45439720877410L) + ((int128)tmp_q[5] * 52372751634947L) + ((int128)tmp_q[6] * 32015356200528L) + ((int128)tmp_q[7] * 82302910548966L) - ((int128)tmp_q[8] * 53972515285225L) + ((int128)tmp_q[9] * 43448070718530L) + ((int128)tmp_q[10] * 498056114976588L);
	tmp_zero[10] = ((int128)tmp_q[0] * 71150873568084L) + ((int128)tmp_q[1] * 27986509043603L) + ((int128)tmp_q[2] * 77979715585343L) + ((int128)tmp_q[3] * 53332178474619L) - ((int128)tmp_q[4] * 2540940244994L) + ((int128)tmp_q[5] * 45439720877410L) + ((int128)tmp_q[6] * 52372751634947L) + ((int128)tmp_q[7] * 32015356200528L) + ((int128)tmp_q[8] * 82302910548966L) - ((int128)tmp_q[9] * 53972515285225L) + ((int128)tmp_q[10] * 43448070718530L);

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

