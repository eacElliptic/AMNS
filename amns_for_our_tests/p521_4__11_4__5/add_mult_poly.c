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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 2);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 2);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 2);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 2);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 2);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 8505313128859288941UL) + ((((uint64_t)op[1] * 1491847016890945692UL) + ((uint64_t)op[2] * 1381244999413868328UL) + ((uint64_t)op[3] * 6106244022919691923UL) + ((uint64_t)op[4] * 11760725813402448903UL) + ((uint64_t)op[5] * 5261299762833598517UL) + ((uint64_t)op[6] * 4732436334796450951UL) + ((uint64_t)op[7] * 15062502108896822889UL) + ((uint64_t)op[8] * 5048225074864028234UL) + ((uint64_t)op[9] * 7195801850613405060UL) + ((uint64_t)op[10] * 14966728221943731428UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 14966728221943731428UL) + ((uint64_t)op[1] * 8505313128859288941UL) + ((((uint64_t)op[2] * 1491847016890945692UL) + ((uint64_t)op[3] * 1381244999413868328UL) + ((uint64_t)op[4] * 6106244022919691923UL) + ((uint64_t)op[5] * 11760725813402448903UL) + ((uint64_t)op[6] * 5261299762833598517UL) + ((uint64_t)op[7] * 4732436334796450951UL) + ((uint64_t)op[8] * 15062502108896822889UL) + ((uint64_t)op[9] * 5048225074864028234UL) + ((uint64_t)op[10] * 7195801850613405060UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 7195801850613405060UL) + ((uint64_t)op[1] * 14966728221943731428UL) + ((uint64_t)op[2] * 8505313128859288941UL) + ((((uint64_t)op[3] * 1491847016890945692UL) + ((uint64_t)op[4] * 1381244999413868328UL) + ((uint64_t)op[5] * 6106244022919691923UL) + ((uint64_t)op[6] * 11760725813402448903UL) + ((uint64_t)op[7] * 5261299762833598517UL) + ((uint64_t)op[8] * 4732436334796450951UL) + ((uint64_t)op[9] * 15062502108896822889UL) + ((uint64_t)op[10] * 5048225074864028234UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 5048225074864028234UL) + ((uint64_t)op[1] * 7195801850613405060UL) + ((uint64_t)op[2] * 14966728221943731428UL) + ((uint64_t)op[3] * 8505313128859288941UL) + ((((uint64_t)op[4] * 1491847016890945692UL) + ((uint64_t)op[5] * 1381244999413868328UL) + ((uint64_t)op[6] * 6106244022919691923UL) + ((uint64_t)op[7] * 11760725813402448903UL) + ((uint64_t)op[8] * 5261299762833598517UL) + ((uint64_t)op[9] * 4732436334796450951UL) + ((uint64_t)op[10] * 15062502108896822889UL)) * 4);
	tmp_q[4] = ((uint64_t)op[0] * 15062502108896822889UL) + ((uint64_t)op[1] * 5048225074864028234UL) + ((uint64_t)op[2] * 7195801850613405060UL) + ((uint64_t)op[3] * 14966728221943731428UL) + ((uint64_t)op[4] * 8505313128859288941UL) + ((((uint64_t)op[5] * 1491847016890945692UL) + ((uint64_t)op[6] * 1381244999413868328UL) + ((uint64_t)op[7] * 6106244022919691923UL) + ((uint64_t)op[8] * 11760725813402448903UL) + ((uint64_t)op[9] * 5261299762833598517UL) + ((uint64_t)op[10] * 4732436334796450951UL)) * 4);
	tmp_q[5] = ((uint64_t)op[0] * 4732436334796450951UL) + ((uint64_t)op[1] * 15062502108896822889UL) + ((uint64_t)op[2] * 5048225074864028234UL) + ((uint64_t)op[3] * 7195801850613405060UL) + ((uint64_t)op[4] * 14966728221943731428UL) + ((uint64_t)op[5] * 8505313128859288941UL) + ((((uint64_t)op[6] * 1491847016890945692UL) + ((uint64_t)op[7] * 1381244999413868328UL) + ((uint64_t)op[8] * 6106244022919691923UL) + ((uint64_t)op[9] * 11760725813402448903UL) + ((uint64_t)op[10] * 5261299762833598517UL)) * 4);
	tmp_q[6] = ((uint64_t)op[0] * 5261299762833598517UL) + ((uint64_t)op[1] * 4732436334796450951UL) + ((uint64_t)op[2] * 15062502108896822889UL) + ((uint64_t)op[3] * 5048225074864028234UL) + ((uint64_t)op[4] * 7195801850613405060UL) + ((uint64_t)op[5] * 14966728221943731428UL) + ((uint64_t)op[6] * 8505313128859288941UL) + ((((uint64_t)op[7] * 1491847016890945692UL) + ((uint64_t)op[8] * 1381244999413868328UL) + ((uint64_t)op[9] * 6106244022919691923UL) + ((uint64_t)op[10] * 11760725813402448903UL)) * 4);
	tmp_q[7] = ((uint64_t)op[0] * 11760725813402448903UL) + ((uint64_t)op[1] * 5261299762833598517UL) + ((uint64_t)op[2] * 4732436334796450951UL) + ((uint64_t)op[3] * 15062502108896822889UL) + ((uint64_t)op[4] * 5048225074864028234UL) + ((uint64_t)op[5] * 7195801850613405060UL) + ((uint64_t)op[6] * 14966728221943731428UL) + ((uint64_t)op[7] * 8505313128859288941UL) + ((((uint64_t)op[8] * 1491847016890945692UL) + ((uint64_t)op[9] * 1381244999413868328UL) + ((uint64_t)op[10] * 6106244022919691923UL)) * 4);
	tmp_q[8] = ((uint64_t)op[0] * 6106244022919691923UL) + ((uint64_t)op[1] * 11760725813402448903UL) + ((uint64_t)op[2] * 5261299762833598517UL) + ((uint64_t)op[3] * 4732436334796450951UL) + ((uint64_t)op[4] * 15062502108896822889UL) + ((uint64_t)op[5] * 5048225074864028234UL) + ((uint64_t)op[6] * 7195801850613405060UL) + ((uint64_t)op[7] * 14966728221943731428UL) + ((uint64_t)op[8] * 8505313128859288941UL) + ((((uint64_t)op[9] * 1491847016890945692UL) + ((uint64_t)op[10] * 1381244999413868328UL)) * 4);
	tmp_q[9] = ((uint64_t)op[0] * 1381244999413868328UL) + ((uint64_t)op[1] * 6106244022919691923UL) + ((uint64_t)op[2] * 11760725813402448903UL) + ((uint64_t)op[3] * 5261299762833598517UL) + ((uint64_t)op[4] * 4732436334796450951UL) + ((uint64_t)op[5] * 15062502108896822889UL) + ((uint64_t)op[6] * 5048225074864028234UL) + ((uint64_t)op[7] * 7195801850613405060UL) + ((uint64_t)op[8] * 14966728221943731428UL) + ((uint64_t)op[9] * 8505313128859288941UL) + ((uint64_t)op[10] * 5967388067563782768UL);
	tmp_q[10] = ((uint64_t)op[0] * 1491847016890945692UL) + ((uint64_t)op[1] * 1381244999413868328UL) + ((uint64_t)op[2] * 6106244022919691923UL) + ((uint64_t)op[3] * 11760725813402448903UL) + ((uint64_t)op[4] * 5261299762833598517UL) + ((uint64_t)op[5] * 4732436334796450951UL) + ((uint64_t)op[6] * 15062502108896822889UL) + ((uint64_t)op[7] * 5048225074864028234UL) + ((uint64_t)op[8] * 7195801850613405060UL) + ((uint64_t)op[9] * 14966728221943731428UL) + ((uint64_t)op[10] * 8505313128859288941UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 30527467821571L) + ((((int128)tmp_q[1] * 22777134166141L) + ((int128)tmp_q[2] * 72612318877134L) - ((int128)tmp_q[3] * 23835179751766L) - ((int128)tmp_q[4] * 9478304829081L) + ((int128)tmp_q[5] * 28628200333609L) + ((int128)tmp_q[6] * 30962081056915L) + ((int128)tmp_q[7] * 26604763100633L) + ((int128)tmp_q[8] * 63517352159718L) - ((int128)tmp_q[9] * 98030402288992L) - ((int128)tmp_q[10] * 8463480129980L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 8463480129980L) + ((int128)tmp_q[1] * 30527467821571L) + ((((int128)tmp_q[2] * 22777134166141L) + ((int128)tmp_q[3] * 72612318877134L) - ((int128)tmp_q[4] * 23835179751766L) - ((int128)tmp_q[5] * 9478304829081L) + ((int128)tmp_q[6] * 28628200333609L) + ((int128)tmp_q[7] * 30962081056915L) + ((int128)tmp_q[8] * 26604763100633L) + ((int128)tmp_q[9] * 63517352159718L) - ((int128)tmp_q[10] * 98030402288992L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 98030402288992L) - ((int128)tmp_q[1] * 8463480129980L) + ((int128)tmp_q[2] * 30527467821571L) + ((((int128)tmp_q[3] * 22777134166141L) + ((int128)tmp_q[4] * 72612318877134L) - ((int128)tmp_q[5] * 23835179751766L) - ((int128)tmp_q[6] * 9478304829081L) + ((int128)tmp_q[7] * 28628200333609L) + ((int128)tmp_q[8] * 30962081056915L) + ((int128)tmp_q[9] * 26604763100633L) + ((int128)tmp_q[10] * 63517352159718L)) * 4);
	tmp_zero[3] = ((int128)tmp_q[0] * 63517352159718L) - ((int128)tmp_q[1] * 98030402288992L) - ((int128)tmp_q[2] * 8463480129980L) + ((int128)tmp_q[3] * 30527467821571L) + ((((int128)tmp_q[4] * 22777134166141L) + ((int128)tmp_q[5] * 72612318877134L) - ((int128)tmp_q[6] * 23835179751766L) - ((int128)tmp_q[7] * 9478304829081L) + ((int128)tmp_q[8] * 28628200333609L) + ((int128)tmp_q[9] * 30962081056915L) + ((int128)tmp_q[10] * 26604763100633L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 26604763100633L) + ((int128)tmp_q[1] * 63517352159718L) - ((int128)tmp_q[2] * 98030402288992L) - ((int128)tmp_q[3] * 8463480129980L) + ((int128)tmp_q[4] * 30527467821571L) + ((((int128)tmp_q[5] * 22777134166141L) + ((int128)tmp_q[6] * 72612318877134L) - ((int128)tmp_q[7] * 23835179751766L) - ((int128)tmp_q[8] * 9478304829081L) + ((int128)tmp_q[9] * 28628200333609L) + ((int128)tmp_q[10] * 30962081056915L)) * 4);
	tmp_zero[5] = ((int128)tmp_q[0] * 30962081056915L) + ((int128)tmp_q[1] * 26604763100633L) + ((int128)tmp_q[2] * 63517352159718L) - ((int128)tmp_q[3] * 98030402288992L) - ((int128)tmp_q[4] * 8463480129980L) + ((int128)tmp_q[5] * 30527467821571L) + ((((int128)tmp_q[6] * 22777134166141L) + ((int128)tmp_q[7] * 72612318877134L) - ((int128)tmp_q[8] * 23835179751766L) - ((int128)tmp_q[9] * 9478304829081L) + ((int128)tmp_q[10] * 28628200333609L)) * 4);
	tmp_zero[6] = ((int128)tmp_q[0] * 28628200333609L) + ((int128)tmp_q[1] * 30962081056915L) + ((int128)tmp_q[2] * 26604763100633L) + ((int128)tmp_q[3] * 63517352159718L) - ((int128)tmp_q[4] * 98030402288992L) - ((int128)tmp_q[5] * 8463480129980L) + ((int128)tmp_q[6] * 30527467821571L) + ((((int128)tmp_q[7] * 22777134166141L) + ((int128)tmp_q[8] * 72612318877134L) - ((int128)tmp_q[9] * 23835179751766L) - ((int128)tmp_q[10] * 9478304829081L)) * 4);
	tmp_zero[7] = -((int128)tmp_q[0] * 9478304829081L) + ((int128)tmp_q[1] * 28628200333609L) + ((int128)tmp_q[2] * 30962081056915L) + ((int128)tmp_q[3] * 26604763100633L) + ((int128)tmp_q[4] * 63517352159718L) - ((int128)tmp_q[5] * 98030402288992L) - ((int128)tmp_q[6] * 8463480129980L) + ((int128)tmp_q[7] * 30527467821571L) + ((((int128)tmp_q[8] * 22777134166141L) + ((int128)tmp_q[9] * 72612318877134L) - ((int128)tmp_q[10] * 23835179751766L)) * 4);
	tmp_zero[8] = -((int128)tmp_q[0] * 23835179751766L) - ((int128)tmp_q[1] * 9478304829081L) + ((int128)tmp_q[2] * 28628200333609L) + ((int128)tmp_q[3] * 30962081056915L) + ((int128)tmp_q[4] * 26604763100633L) + ((int128)tmp_q[5] * 63517352159718L) - ((int128)tmp_q[6] * 98030402288992L) - ((int128)tmp_q[7] * 8463480129980L) + ((int128)tmp_q[8] * 30527467821571L) + ((((int128)tmp_q[9] * 22777134166141L) + ((int128)tmp_q[10] * 72612318877134L)) * 4);
	tmp_zero[9] = ((int128)tmp_q[0] * 72612318877134L) - ((int128)tmp_q[1] * 23835179751766L) - ((int128)tmp_q[2] * 9478304829081L) + ((int128)tmp_q[3] * 28628200333609L) + ((int128)tmp_q[4] * 30962081056915L) + ((int128)tmp_q[5] * 26604763100633L) + ((int128)tmp_q[6] * 63517352159718L) - ((int128)tmp_q[7] * 98030402288992L) - ((int128)tmp_q[8] * 8463480129980L) + ((int128)tmp_q[9] * 30527467821571L) + ((int128)tmp_q[10] * 91108536664564L);
	tmp_zero[10] = ((int128)tmp_q[0] * 22777134166141L) + ((int128)tmp_q[1] * 72612318877134L) - ((int128)tmp_q[2] * 23835179751766L) - ((int128)tmp_q[3] * 9478304829081L) + ((int128)tmp_q[4] * 28628200333609L) + ((int128)tmp_q[5] * 30962081056915L) + ((int128)tmp_q[6] * 26604763100633L) + ((int128)tmp_q[7] * 63517352159718L) - ((int128)tmp_q[8] * 98030402288992L) - ((int128)tmp_q[9] * 8463480129980L) + ((int128)tmp_q[10] * 30527467821571L);

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
