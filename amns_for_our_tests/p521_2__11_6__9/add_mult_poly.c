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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 6);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 6);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 12);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3218963048509694375UL) + ((((uint64_t)op[1] * 1984705834937581380UL) + ((uint64_t)op[2] * 10048368719627705526UL) + ((uint64_t)op[3] * 17317612611506448436UL) + ((uint64_t)op[4] * 17869404735284978751UL) + ((uint64_t)op[5] * 15614354907554336032UL) + ((uint64_t)op[6] * 4412557358008408995UL) + ((uint64_t)op[7] * 4629173112344707480UL) + ((uint64_t)op[8] * 8493843549362610491UL) + ((uint64_t)op[9] * 6241635919442698342UL) + ((uint64_t)op[10] * 5796583550015554858UL)) * 6);
	tmp_q[1] = ((uint64_t)op[0] * 5796583550015554858UL) + ((uint64_t)op[1] * 3218963048509694375UL) + ((((uint64_t)op[2] * 1984705834937581380UL) + ((uint64_t)op[3] * 10048368719627705526UL) + ((uint64_t)op[4] * 17317612611506448436UL) + ((uint64_t)op[5] * 17869404735284978751UL) + ((uint64_t)op[6] * 15614354907554336032UL) + ((uint64_t)op[7] * 4412557358008408995UL) + ((uint64_t)op[8] * 4629173112344707480UL) + ((uint64_t)op[9] * 8493843549362610491UL) + ((uint64_t)op[10] * 6241635919442698342UL)) * 6);
	tmp_q[2] = ((uint64_t)op[0] * 6241635919442698342UL) + ((uint64_t)op[1] * 5796583550015554858UL) + ((uint64_t)op[2] * 3218963048509694375UL) + ((((uint64_t)op[3] * 1984705834937581380UL) + ((uint64_t)op[4] * 10048368719627705526UL) + ((uint64_t)op[5] * 17317612611506448436UL) + ((uint64_t)op[6] * 17869404735284978751UL) + ((uint64_t)op[7] * 15614354907554336032UL) + ((uint64_t)op[8] * 4412557358008408995UL) + ((uint64_t)op[9] * 4629173112344707480UL) + ((uint64_t)op[10] * 8493843549362610491UL)) * 6);
	tmp_q[3] = ((uint64_t)op[0] * 8493843549362610491UL) + ((uint64_t)op[1] * 6241635919442698342UL) + ((uint64_t)op[2] * 5796583550015554858UL) + ((uint64_t)op[3] * 3218963048509694375UL) + ((((uint64_t)op[4] * 1984705834937581380UL) + ((uint64_t)op[5] * 10048368719627705526UL) + ((uint64_t)op[6] * 17317612611506448436UL) + ((uint64_t)op[7] * 17869404735284978751UL) + ((uint64_t)op[8] * 15614354907554336032UL) + ((uint64_t)op[9] * 4412557358008408995UL) + ((uint64_t)op[10] * 4629173112344707480UL)) * 6);
	tmp_q[4] = ((uint64_t)op[0] * 4629173112344707480UL) + ((uint64_t)op[1] * 8493843549362610491UL) + ((uint64_t)op[2] * 6241635919442698342UL) + ((uint64_t)op[3] * 5796583550015554858UL) + ((uint64_t)op[4] * 3218963048509694375UL) + ((((uint64_t)op[5] * 1984705834937581380UL) + ((uint64_t)op[6] * 10048368719627705526UL) + ((uint64_t)op[7] * 17317612611506448436UL) + ((uint64_t)op[8] * 17869404735284978751UL) + ((uint64_t)op[9] * 15614354907554336032UL) + ((uint64_t)op[10] * 4412557358008408995UL)) * 6);
	tmp_q[5] = ((uint64_t)op[0] * 4412557358008408995UL) + ((uint64_t)op[1] * 4629173112344707480UL) + ((uint64_t)op[2] * 8493843549362610491UL) + ((uint64_t)op[3] * 6241635919442698342UL) + ((uint64_t)op[4] * 5796583550015554858UL) + ((uint64_t)op[5] * 3218963048509694375UL) + ((((uint64_t)op[6] * 1984705834937581380UL) + ((uint64_t)op[7] * 10048368719627705526UL) + ((uint64_t)op[8] * 17317612611506448436UL) + ((uint64_t)op[9] * 17869404735284978751UL) + ((uint64_t)op[10] * 15614354907554336032UL)) * 6);
	tmp_q[6] = ((uint64_t)op[0] * 15614354907554336032UL) + ((uint64_t)op[1] * 4412557358008408995UL) + ((uint64_t)op[2] * 4629173112344707480UL) + ((uint64_t)op[3] * 8493843549362610491UL) + ((uint64_t)op[4] * 6241635919442698342UL) + ((uint64_t)op[5] * 5796583550015554858UL) + ((uint64_t)op[6] * 3218963048509694375UL) + ((((uint64_t)op[7] * 1984705834937581380UL) + ((uint64_t)op[8] * 10048368719627705526UL) + ((uint64_t)op[9] * 17317612611506448436UL) + ((uint64_t)op[10] * 17869404735284978751UL)) * 6);
	tmp_q[7] = ((uint64_t)op[0] * 17869404735284978751UL) + ((uint64_t)op[1] * 15614354907554336032UL) + ((uint64_t)op[2] * 4412557358008408995UL) + ((uint64_t)op[3] * 4629173112344707480UL) + ((uint64_t)op[4] * 8493843549362610491UL) + ((uint64_t)op[5] * 6241635919442698342UL) + ((uint64_t)op[6] * 5796583550015554858UL) + ((uint64_t)op[7] * 3218963048509694375UL) + ((((uint64_t)op[8] * 1984705834937581380UL) + ((uint64_t)op[9] * 10048368719627705526UL) + ((uint64_t)op[10] * 17317612611506448436UL)) * 6);
	tmp_q[8] = ((uint64_t)op[0] * 17317612611506448436UL) + ((uint64_t)op[1] * 17869404735284978751UL) + ((uint64_t)op[2] * 15614354907554336032UL) + ((uint64_t)op[3] * 4412557358008408995UL) + ((uint64_t)op[4] * 4629173112344707480UL) + ((uint64_t)op[5] * 8493843549362610491UL) + ((uint64_t)op[6] * 6241635919442698342UL) + ((uint64_t)op[7] * 5796583550015554858UL) + ((uint64_t)op[8] * 3218963048509694375UL) + ((((uint64_t)op[9] * 1984705834937581380UL) + ((uint64_t)op[10] * 10048368719627705526UL)) * 6);
	tmp_q[9] = ((uint64_t)op[0] * 10048368719627705526UL) + ((uint64_t)op[1] * 17317612611506448436UL) + ((uint64_t)op[2] * 17869404735284978751UL) + ((uint64_t)op[3] * 15614354907554336032UL) + ((uint64_t)op[4] * 4412557358008408995UL) + ((uint64_t)op[5] * 4629173112344707480UL) + ((uint64_t)op[6] * 8493843549362610491UL) + ((uint64_t)op[7] * 6241635919442698342UL) + ((uint64_t)op[8] * 5796583550015554858UL) + ((uint64_t)op[9] * 3218963048509694375UL) + ((uint64_t)op[10] * 11908235009625488280UL);
	tmp_q[10] = ((uint64_t)op[0] * 1984705834937581380UL) + ((uint64_t)op[1] * 10048368719627705526UL) + ((uint64_t)op[2] * 17317612611506448436UL) + ((uint64_t)op[3] * 17869404735284978751UL) + ((uint64_t)op[4] * 15614354907554336032UL) + ((uint64_t)op[5] * 4412557358008408995UL) + ((uint64_t)op[6] * 4629173112344707480UL) + ((uint64_t)op[7] * 8493843549362610491UL) + ((uint64_t)op[8] * 6241635919442698342UL) + ((uint64_t)op[9] * 5796583550015554858UL) + ((uint64_t)op[10] * 3218963048509694375UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 47152739198567L) + ((((int128)tmp_q[1] * 85385572776723L) - ((int128)tmp_q[2] * 80739975872373L) + ((int128)tmp_q[3] * 79715087220688L) + ((int128)tmp_q[4] * 71861881608951L) + ((int128)tmp_q[5] * 12431342723415L) + ((int128)tmp_q[6] * 36165464101959L) - ((int128)tmp_q[7] * 97662454639928L) - ((int128)tmp_q[8] * 11658231501935L) + ((int128)tmp_q[9] * 46526520283174L) + ((int128)tmp_q[10] * 36595499839376L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 36595499839376L) + ((int128)tmp_q[1] * 47152739198567L) + ((((int128)tmp_q[2] * 85385572776723L) - ((int128)tmp_q[3] * 80739975872373L) + ((int128)tmp_q[4] * 79715087220688L) + ((int128)tmp_q[5] * 71861881608951L) + ((int128)tmp_q[6] * 12431342723415L) + ((int128)tmp_q[7] * 36165464101959L) - ((int128)tmp_q[8] * 97662454639928L) - ((int128)tmp_q[9] * 11658231501935L) + ((int128)tmp_q[10] * 46526520283174L)) * 6);
	tmp_zero[2] = ((int128)tmp_q[0] * 46526520283174L) + ((int128)tmp_q[1] * 36595499839376L) + ((int128)tmp_q[2] * 47152739198567L) + ((((int128)tmp_q[3] * 85385572776723L) - ((int128)tmp_q[4] * 80739975872373L) + ((int128)tmp_q[5] * 79715087220688L) + ((int128)tmp_q[6] * 71861881608951L) + ((int128)tmp_q[7] * 12431342723415L) + ((int128)tmp_q[8] * 36165464101959L) - ((int128)tmp_q[9] * 97662454639928L) - ((int128)tmp_q[10] * 11658231501935L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 11658231501935L) + ((int128)tmp_q[1] * 46526520283174L) + ((int128)tmp_q[2] * 36595499839376L) + ((int128)tmp_q[3] * 47152739198567L) + ((((int128)tmp_q[4] * 85385572776723L) - ((int128)tmp_q[5] * 80739975872373L) + ((int128)tmp_q[6] * 79715087220688L) + ((int128)tmp_q[7] * 71861881608951L) + ((int128)tmp_q[8] * 12431342723415L) + ((int128)tmp_q[9] * 36165464101959L) - ((int128)tmp_q[10] * 97662454639928L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 97662454639928L) - ((int128)tmp_q[1] * 11658231501935L) + ((int128)tmp_q[2] * 46526520283174L) + ((int128)tmp_q[3] * 36595499839376L) + ((int128)tmp_q[4] * 47152739198567L) + ((((int128)tmp_q[5] * 85385572776723L) - ((int128)tmp_q[6] * 80739975872373L) + ((int128)tmp_q[7] * 79715087220688L) + ((int128)tmp_q[8] * 71861881608951L) + ((int128)tmp_q[9] * 12431342723415L) + ((int128)tmp_q[10] * 36165464101959L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 36165464101959L) - ((int128)tmp_q[1] * 97662454639928L) - ((int128)tmp_q[2] * 11658231501935L) + ((int128)tmp_q[3] * 46526520283174L) + ((int128)tmp_q[4] * 36595499839376L) + ((int128)tmp_q[5] * 47152739198567L) + ((((int128)tmp_q[6] * 85385572776723L) - ((int128)tmp_q[7] * 80739975872373L) + ((int128)tmp_q[8] * 79715087220688L) + ((int128)tmp_q[9] * 71861881608951L) + ((int128)tmp_q[10] * 12431342723415L)) * 6);
	tmp_zero[6] = ((int128)tmp_q[0] * 12431342723415L) + ((int128)tmp_q[1] * 36165464101959L) - ((int128)tmp_q[2] * 97662454639928L) - ((int128)tmp_q[3] * 11658231501935L) + ((int128)tmp_q[4] * 46526520283174L) + ((int128)tmp_q[5] * 36595499839376L) + ((int128)tmp_q[6] * 47152739198567L) + ((((int128)tmp_q[7] * 85385572776723L) - ((int128)tmp_q[8] * 80739975872373L) + ((int128)tmp_q[9] * 79715087220688L) + ((int128)tmp_q[10] * 71861881608951L)) * 6);
	tmp_zero[7] = ((int128)tmp_q[0] * 71861881608951L) + ((int128)tmp_q[1] * 12431342723415L) + ((int128)tmp_q[2] * 36165464101959L) - ((int128)tmp_q[3] * 97662454639928L) - ((int128)tmp_q[4] * 11658231501935L) + ((int128)tmp_q[5] * 46526520283174L) + ((int128)tmp_q[6] * 36595499839376L) + ((int128)tmp_q[7] * 47152739198567L) + ((((int128)tmp_q[8] * 85385572776723L) - ((int128)tmp_q[9] * 80739975872373L) + ((int128)tmp_q[10] * 79715087220688L)) * 6);
	tmp_zero[8] = ((int128)tmp_q[0] * 79715087220688L) + ((int128)tmp_q[1] * 71861881608951L) + ((int128)tmp_q[2] * 12431342723415L) + ((int128)tmp_q[3] * 36165464101959L) - ((int128)tmp_q[4] * 97662454639928L) - ((int128)tmp_q[5] * 11658231501935L) + ((int128)tmp_q[6] * 46526520283174L) + ((int128)tmp_q[7] * 36595499839376L) + ((int128)tmp_q[8] * 47152739198567L) + ((((int128)tmp_q[9] * 85385572776723L) - ((int128)tmp_q[10] * 80739975872373L)) * 6);
	tmp_zero[9] = -((int128)tmp_q[0] * 80739975872373L) + ((int128)tmp_q[1] * 79715087220688L) + ((int128)tmp_q[2] * 71861881608951L) + ((int128)tmp_q[3] * 12431342723415L) + ((int128)tmp_q[4] * 36165464101959L) - ((int128)tmp_q[5] * 97662454639928L) - ((int128)tmp_q[6] * 11658231501935L) + ((int128)tmp_q[7] * 46526520283174L) + ((int128)tmp_q[8] * 36595499839376L) + ((int128)tmp_q[9] * 47152739198567L) + ((int128)tmp_q[10] * 512313436660338L);
	tmp_zero[10] = ((int128)tmp_q[0] * 85385572776723L) - ((int128)tmp_q[1] * 80739975872373L) + ((int128)tmp_q[2] * 79715087220688L) + ((int128)tmp_q[3] * 71861881608951L) + ((int128)tmp_q[4] * 12431342723415L) + ((int128)tmp_q[5] * 36165464101959L) - ((int128)tmp_q[6] * 97662454639928L) - ((int128)tmp_q[7] * 11658231501935L) + ((int128)tmp_q[8] * 46526520283174L) + ((int128)tmp_q[9] * 36595499839376L) + ((int128)tmp_q[10] * 47152739198567L);

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

