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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 5);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 5);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] + (((int128)pa[11] * pb[11]) * 5);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 10);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 5);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 10);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 5);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[11] * pa[10]) * 10);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] + (((int128)pa[11] * pa[11]) * 5);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17892405478078705354UL) + ((((uint64_t)op[1] * 10340692951313638744UL) + ((uint64_t)op[2] * 11245754194899759681UL) + ((uint64_t)op[3] * 11548460736881947698UL) + ((uint64_t)op[4] * 10836810805779995005UL) + ((uint64_t)op[5] * 5773617361995173452UL) + ((uint64_t)op[6] * 18266796174982462578UL) + ((uint64_t)op[7] * 18442144638809822964UL) + ((uint64_t)op[8] * 12873464325429654447UL) + ((uint64_t)op[9] * 7094750280415094423UL) + ((uint64_t)op[10] * 14340053468723864591UL) + ((uint64_t)op[11] * 12810105075866897950UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 12810105075866897950UL) + ((uint64_t)op[1] * 17892405478078705354UL) + ((((uint64_t)op[2] * 10340692951313638744UL) + ((uint64_t)op[3] * 11245754194899759681UL) + ((uint64_t)op[4] * 11548460736881947698UL) + ((uint64_t)op[5] * 10836810805779995005UL) + ((uint64_t)op[6] * 5773617361995173452UL) + ((uint64_t)op[7] * 18266796174982462578UL) + ((uint64_t)op[8] * 18442144638809822964UL) + ((uint64_t)op[9] * 12873464325429654447UL) + ((uint64_t)op[10] * 7094750280415094423UL) + ((uint64_t)op[11] * 14340053468723864591UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 14340053468723864591UL) + ((uint64_t)op[1] * 12810105075866897950UL) + ((uint64_t)op[2] * 17892405478078705354UL) + ((((uint64_t)op[3] * 10340692951313638744UL) + ((uint64_t)op[4] * 11245754194899759681UL) + ((uint64_t)op[5] * 11548460736881947698UL) + ((uint64_t)op[6] * 10836810805779995005UL) + ((uint64_t)op[7] * 5773617361995173452UL) + ((uint64_t)op[8] * 18266796174982462578UL) + ((uint64_t)op[9] * 18442144638809822964UL) + ((uint64_t)op[10] * 12873464325429654447UL) + ((uint64_t)op[11] * 7094750280415094423UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 7094750280415094423UL) + ((uint64_t)op[1] * 14340053468723864591UL) + ((uint64_t)op[2] * 12810105075866897950UL) + ((uint64_t)op[3] * 17892405478078705354UL) + ((((uint64_t)op[4] * 10340692951313638744UL) + ((uint64_t)op[5] * 11245754194899759681UL) + ((uint64_t)op[6] * 11548460736881947698UL) + ((uint64_t)op[7] * 10836810805779995005UL) + ((uint64_t)op[8] * 5773617361995173452UL) + ((uint64_t)op[9] * 18266796174982462578UL) + ((uint64_t)op[10] * 18442144638809822964UL) + ((uint64_t)op[11] * 12873464325429654447UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 12873464325429654447UL) + ((uint64_t)op[1] * 7094750280415094423UL) + ((uint64_t)op[2] * 14340053468723864591UL) + ((uint64_t)op[3] * 12810105075866897950UL) + ((uint64_t)op[4] * 17892405478078705354UL) + ((((uint64_t)op[5] * 10340692951313638744UL) + ((uint64_t)op[6] * 11245754194899759681UL) + ((uint64_t)op[7] * 11548460736881947698UL) + ((uint64_t)op[8] * 10836810805779995005UL) + ((uint64_t)op[9] * 5773617361995173452UL) + ((uint64_t)op[10] * 18266796174982462578UL) + ((uint64_t)op[11] * 18442144638809822964UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 18442144638809822964UL) + ((uint64_t)op[1] * 12873464325429654447UL) + ((uint64_t)op[2] * 7094750280415094423UL) + ((uint64_t)op[3] * 14340053468723864591UL) + ((uint64_t)op[4] * 12810105075866897950UL) + ((uint64_t)op[5] * 17892405478078705354UL) + ((((uint64_t)op[6] * 10340692951313638744UL) + ((uint64_t)op[7] * 11245754194899759681UL) + ((uint64_t)op[8] * 11548460736881947698UL) + ((uint64_t)op[9] * 10836810805779995005UL) + ((uint64_t)op[10] * 5773617361995173452UL) + ((uint64_t)op[11] * 18266796174982462578UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 18266796174982462578UL) + ((uint64_t)op[1] * 18442144638809822964UL) + ((uint64_t)op[2] * 12873464325429654447UL) + ((uint64_t)op[3] * 7094750280415094423UL) + ((uint64_t)op[4] * 14340053468723864591UL) + ((uint64_t)op[5] * 12810105075866897950UL) + ((uint64_t)op[6] * 17892405478078705354UL) + ((((uint64_t)op[7] * 10340692951313638744UL) + ((uint64_t)op[8] * 11245754194899759681UL) + ((uint64_t)op[9] * 11548460736881947698UL) + ((uint64_t)op[10] * 10836810805779995005UL) + ((uint64_t)op[11] * 5773617361995173452UL)) * 5);
	tmp_q[7] = ((uint64_t)op[0] * 5773617361995173452UL) + ((uint64_t)op[1] * 18266796174982462578UL) + ((uint64_t)op[2] * 18442144638809822964UL) + ((uint64_t)op[3] * 12873464325429654447UL) + ((uint64_t)op[4] * 7094750280415094423UL) + ((uint64_t)op[5] * 14340053468723864591UL) + ((uint64_t)op[6] * 12810105075866897950UL) + ((uint64_t)op[7] * 17892405478078705354UL) + ((((uint64_t)op[8] * 10340692951313638744UL) + ((uint64_t)op[9] * 11245754194899759681UL) + ((uint64_t)op[10] * 11548460736881947698UL) + ((uint64_t)op[11] * 10836810805779995005UL)) * 5);
	tmp_q[8] = ((uint64_t)op[0] * 10836810805779995005UL) + ((uint64_t)op[1] * 5773617361995173452UL) + ((uint64_t)op[2] * 18266796174982462578UL) + ((uint64_t)op[3] * 18442144638809822964UL) + ((uint64_t)op[4] * 12873464325429654447UL) + ((uint64_t)op[5] * 7094750280415094423UL) + ((uint64_t)op[6] * 14340053468723864591UL) + ((uint64_t)op[7] * 12810105075866897950UL) + ((uint64_t)op[8] * 17892405478078705354UL) + ((((uint64_t)op[9] * 10340692951313638744UL) + ((uint64_t)op[10] * 11245754194899759681UL) + ((uint64_t)op[11] * 11548460736881947698UL)) * 5);
	tmp_q[9] = ((uint64_t)op[0] * 11548460736881947698UL) + ((uint64_t)op[1] * 10836810805779995005UL) + ((uint64_t)op[2] * 5773617361995173452UL) + ((uint64_t)op[3] * 18266796174982462578UL) + ((uint64_t)op[4] * 18442144638809822964UL) + ((uint64_t)op[5] * 12873464325429654447UL) + ((uint64_t)op[6] * 7094750280415094423UL) + ((uint64_t)op[7] * 14340053468723864591UL) + ((uint64_t)op[8] * 12810105075866897950UL) + ((uint64_t)op[9] * 17892405478078705354UL) + ((((uint64_t)op[10] * 10340692951313638744UL) + ((uint64_t)op[11] * 11245754194899759681UL)) * 5);
	tmp_q[10] = ((uint64_t)op[0] * 11245754194899759681UL) + ((uint64_t)op[1] * 11548460736881947698UL) + ((uint64_t)op[2] * 10836810805779995005UL) + ((uint64_t)op[3] * 5773617361995173452UL) + ((uint64_t)op[4] * 18266796174982462578UL) + ((uint64_t)op[5] * 18442144638809822964UL) + ((uint64_t)op[6] * 12873464325429654447UL) + ((uint64_t)op[7] * 7094750280415094423UL) + ((uint64_t)op[8] * 14340053468723864591UL) + ((uint64_t)op[9] * 12810105075866897950UL) + ((uint64_t)op[10] * 17892405478078705354UL) + ((uint64_t)op[11] * 14809976609149090488UL);
	tmp_q[11] = ((uint64_t)op[0] * 10340692951313638744UL) + ((uint64_t)op[1] * 11245754194899759681UL) + ((uint64_t)op[2] * 11548460736881947698UL) + ((uint64_t)op[3] * 10836810805779995005UL) + ((uint64_t)op[4] * 5773617361995173452UL) + ((uint64_t)op[5] * 18266796174982462578UL) + ((uint64_t)op[6] * 18442144638809822964UL) + ((uint64_t)op[7] * 12873464325429654447UL) + ((uint64_t)op[8] * 7094750280415094423UL) + ((uint64_t)op[9] * 14340053468723864591UL) + ((uint64_t)op[10] * 12810105075866897950UL) + ((uint64_t)op[11] * 17892405478078705354UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 6173072770014L) + ((-((int128)tmp_q[1] * 3351462963846L) + ((int128)tmp_q[2] * 3864706493367L) - ((int128)tmp_q[3] * 1142377339511L) + ((int128)tmp_q[4] * 6246735182063L) - ((int128)tmp_q[5] * 900093356720L) - ((int128)tmp_q[6] * 178266331166L) - ((int128)tmp_q[7] * 779550961030L) + ((int128)tmp_q[8] * 1915534336725L) + ((int128)tmp_q[9] * 334585225538L) - ((int128)tmp_q[10] * 747273006231L) - ((int128)tmp_q[11] * 3175005397834L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 3175005397834L) - ((int128)tmp_q[1] * 6173072770014L) + ((-((int128)tmp_q[2] * 3351462963846L) + ((int128)tmp_q[3] * 3864706493367L) - ((int128)tmp_q[4] * 1142377339511L) + ((int128)tmp_q[5] * 6246735182063L) - ((int128)tmp_q[6] * 900093356720L) - ((int128)tmp_q[7] * 178266331166L) - ((int128)tmp_q[8] * 779550961030L) + ((int128)tmp_q[9] * 1915534336725L) + ((int128)tmp_q[10] * 334585225538L) - ((int128)tmp_q[11] * 747273006231L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 747273006231L) - ((int128)tmp_q[1] * 3175005397834L) - ((int128)tmp_q[2] * 6173072770014L) + ((-((int128)tmp_q[3] * 3351462963846L) + ((int128)tmp_q[4] * 3864706493367L) - ((int128)tmp_q[5] * 1142377339511L) + ((int128)tmp_q[6] * 6246735182063L) - ((int128)tmp_q[7] * 900093356720L) - ((int128)tmp_q[8] * 178266331166L) - ((int128)tmp_q[9] * 779550961030L) + ((int128)tmp_q[10] * 1915534336725L) + ((int128)tmp_q[11] * 334585225538L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 334585225538L) - ((int128)tmp_q[1] * 747273006231L) - ((int128)tmp_q[2] * 3175005397834L) - ((int128)tmp_q[3] * 6173072770014L) + ((-((int128)tmp_q[4] * 3351462963846L) + ((int128)tmp_q[5] * 3864706493367L) - ((int128)tmp_q[6] * 1142377339511L) + ((int128)tmp_q[7] * 6246735182063L) - ((int128)tmp_q[8] * 900093356720L) - ((int128)tmp_q[9] * 178266331166L) - ((int128)tmp_q[10] * 779550961030L) + ((int128)tmp_q[11] * 1915534336725L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 1915534336725L) + ((int128)tmp_q[1] * 334585225538L) - ((int128)tmp_q[2] * 747273006231L) - ((int128)tmp_q[3] * 3175005397834L) - ((int128)tmp_q[4] * 6173072770014L) + ((-((int128)tmp_q[5] * 3351462963846L) + ((int128)tmp_q[6] * 3864706493367L) - ((int128)tmp_q[7] * 1142377339511L) + ((int128)tmp_q[8] * 6246735182063L) - ((int128)tmp_q[9] * 900093356720L) - ((int128)tmp_q[10] * 178266331166L) - ((int128)tmp_q[11] * 779550961030L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 779550961030L) + ((int128)tmp_q[1] * 1915534336725L) + ((int128)tmp_q[2] * 334585225538L) - ((int128)tmp_q[3] * 747273006231L) - ((int128)tmp_q[4] * 3175005397834L) - ((int128)tmp_q[5] * 6173072770014L) + ((-((int128)tmp_q[6] * 3351462963846L) + ((int128)tmp_q[7] * 3864706493367L) - ((int128)tmp_q[8] * 1142377339511L) + ((int128)tmp_q[9] * 6246735182063L) - ((int128)tmp_q[10] * 900093356720L) - ((int128)tmp_q[11] * 178266331166L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 178266331166L) - ((int128)tmp_q[1] * 779550961030L) + ((int128)tmp_q[2] * 1915534336725L) + ((int128)tmp_q[3] * 334585225538L) - ((int128)tmp_q[4] * 747273006231L) - ((int128)tmp_q[5] * 3175005397834L) - ((int128)tmp_q[6] * 6173072770014L) + ((-((int128)tmp_q[7] * 3351462963846L) + ((int128)tmp_q[8] * 3864706493367L) - ((int128)tmp_q[9] * 1142377339511L) + ((int128)tmp_q[10] * 6246735182063L) - ((int128)tmp_q[11] * 900093356720L)) * 5);
	tmp_zero[7] = -((int128)tmp_q[0] * 900093356720L) - ((int128)tmp_q[1] * 178266331166L) - ((int128)tmp_q[2] * 779550961030L) + ((int128)tmp_q[3] * 1915534336725L) + ((int128)tmp_q[4] * 334585225538L) - ((int128)tmp_q[5] * 747273006231L) - ((int128)tmp_q[6] * 3175005397834L) - ((int128)tmp_q[7] * 6173072770014L) + ((-((int128)tmp_q[8] * 3351462963846L) + ((int128)tmp_q[9] * 3864706493367L) - ((int128)tmp_q[10] * 1142377339511L) + ((int128)tmp_q[11] * 6246735182063L)) * 5);
	tmp_zero[8] = ((int128)tmp_q[0] * 6246735182063L) - ((int128)tmp_q[1] * 900093356720L) - ((int128)tmp_q[2] * 178266331166L) - ((int128)tmp_q[3] * 779550961030L) + ((int128)tmp_q[4] * 1915534336725L) + ((int128)tmp_q[5] * 334585225538L) - ((int128)tmp_q[6] * 747273006231L) - ((int128)tmp_q[7] * 3175005397834L) - ((int128)tmp_q[8] * 6173072770014L) + ((-((int128)tmp_q[9] * 3351462963846L) + ((int128)tmp_q[10] * 3864706493367L) - ((int128)tmp_q[11] * 1142377339511L)) * 5);
	tmp_zero[9] = -((int128)tmp_q[0] * 1142377339511L) + ((int128)tmp_q[1] * 6246735182063L) - ((int128)tmp_q[2] * 900093356720L) - ((int128)tmp_q[3] * 178266331166L) - ((int128)tmp_q[4] * 779550961030L) + ((int128)tmp_q[5] * 1915534336725L) + ((int128)tmp_q[6] * 334585225538L) - ((int128)tmp_q[7] * 747273006231L) - ((int128)tmp_q[8] * 3175005397834L) - ((int128)tmp_q[9] * 6173072770014L) + ((-((int128)tmp_q[10] * 3351462963846L) + ((int128)tmp_q[11] * 3864706493367L)) * 5);
	tmp_zero[10] = ((int128)tmp_q[0] * 3864706493367L) - ((int128)tmp_q[1] * 1142377339511L) + ((int128)tmp_q[2] * 6246735182063L) - ((int128)tmp_q[3] * 900093356720L) - ((int128)tmp_q[4] * 178266331166L) - ((int128)tmp_q[5] * 779550961030L) + ((int128)tmp_q[6] * 1915534336725L) + ((int128)tmp_q[7] * 334585225538L) - ((int128)tmp_q[8] * 747273006231L) - ((int128)tmp_q[9] * 3175005397834L) - ((int128)tmp_q[10] * 6173072770014L) - ((int128)tmp_q[11] * 16757314819230L);
	tmp_zero[11] = -((int128)tmp_q[0] * 3351462963846L) + ((int128)tmp_q[1] * 3864706493367L) - ((int128)tmp_q[2] * 1142377339511L) + ((int128)tmp_q[3] * 6246735182063L) - ((int128)tmp_q[4] * 900093356720L) - ((int128)tmp_q[5] * 178266331166L) - ((int128)tmp_q[6] * 779550961030L) + ((int128)tmp_q[7] * 1915534336725L) + ((int128)tmp_q[8] * 334585225538L) - ((int128)tmp_q[9] * 747273006231L) - ((int128)tmp_q[10] * 3175005397834L) - ((int128)tmp_q[11] * 6173072770014L);

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

