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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[9] + (int128)pa[9] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[9] * pa[7]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[9] * pa[8]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18320529210780447357UL) + ((((uint64_t)op[1] * 16034886971845528088UL) + ((uint64_t)op[2] * 15590941804916051280UL) + ((uint64_t)op[3] * 2584551659114682868UL) + ((uint64_t)op[4] * 14494164424885512764UL) + ((uint64_t)op[5] * 15740146837941994911UL) + ((uint64_t)op[6] * 15045817954949259965UL) + ((uint64_t)op[7] * 3774698560784541302UL) + ((uint64_t)op[8] * 2450953248114263220UL) + ((uint64_t)op[9] * 10333695591267346857UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 10333695591267346857UL) + ((uint64_t)op[1] * 18320529210780447357UL) + ((((uint64_t)op[2] * 16034886971845528088UL) + ((uint64_t)op[3] * 15590941804916051280UL) + ((uint64_t)op[4] * 2584551659114682868UL) + ((uint64_t)op[5] * 14494164424885512764UL) + ((uint64_t)op[6] * 15740146837941994911UL) + ((uint64_t)op[7] * 15045817954949259965UL) + ((uint64_t)op[8] * 3774698560784541302UL) + ((uint64_t)op[9] * 2450953248114263220UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 2450953248114263220UL) + ((uint64_t)op[1] * 10333695591267346857UL) + ((uint64_t)op[2] * 18320529210780447357UL) + ((((uint64_t)op[3] * 16034886971845528088UL) + ((uint64_t)op[4] * 15590941804916051280UL) + ((uint64_t)op[5] * 2584551659114682868UL) + ((uint64_t)op[6] * 14494164424885512764UL) + ((uint64_t)op[7] * 15740146837941994911UL) + ((uint64_t)op[8] * 15045817954949259965UL) + ((uint64_t)op[9] * 3774698560784541302UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 3774698560784541302UL) + ((uint64_t)op[1] * 2450953248114263220UL) + ((uint64_t)op[2] * 10333695591267346857UL) + ((uint64_t)op[3] * 18320529210780447357UL) + ((((uint64_t)op[4] * 16034886971845528088UL) + ((uint64_t)op[5] * 15590941804916051280UL) + ((uint64_t)op[6] * 2584551659114682868UL) + ((uint64_t)op[7] * 14494164424885512764UL) + ((uint64_t)op[8] * 15740146837941994911UL) + ((uint64_t)op[9] * 15045817954949259965UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 15045817954949259965UL) + ((uint64_t)op[1] * 3774698560784541302UL) + ((uint64_t)op[2] * 2450953248114263220UL) + ((uint64_t)op[3] * 10333695591267346857UL) + ((uint64_t)op[4] * 18320529210780447357UL) + ((((uint64_t)op[5] * 16034886971845528088UL) + ((uint64_t)op[6] * 15590941804916051280UL) + ((uint64_t)op[7] * 2584551659114682868UL) + ((uint64_t)op[8] * 14494164424885512764UL) + ((uint64_t)op[9] * 15740146837941994911UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 15740146837941994911UL) + ((uint64_t)op[1] * 15045817954949259965UL) + ((uint64_t)op[2] * 3774698560784541302UL) + ((uint64_t)op[3] * 2450953248114263220UL) + ((uint64_t)op[4] * 10333695591267346857UL) + ((uint64_t)op[5] * 18320529210780447357UL) + ((((uint64_t)op[6] * 16034886971845528088UL) + ((uint64_t)op[7] * 15590941804916051280UL) + ((uint64_t)op[8] * 2584551659114682868UL) + ((uint64_t)op[9] * 14494164424885512764UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 14494164424885512764UL) + ((uint64_t)op[1] * 15740146837941994911UL) + ((uint64_t)op[2] * 15045817954949259965UL) + ((uint64_t)op[3] * 3774698560784541302UL) + ((uint64_t)op[4] * 2450953248114263220UL) + ((uint64_t)op[5] * 10333695591267346857UL) + ((uint64_t)op[6] * 18320529210780447357UL) + ((((uint64_t)op[7] * 16034886971845528088UL) + ((uint64_t)op[8] * 15590941804916051280UL) + ((uint64_t)op[9] * 2584551659114682868UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 2584551659114682868UL) + ((uint64_t)op[1] * 14494164424885512764UL) + ((uint64_t)op[2] * 15740146837941994911UL) + ((uint64_t)op[3] * 15045817954949259965UL) + ((uint64_t)op[4] * 3774698560784541302UL) + ((uint64_t)op[5] * 2450953248114263220UL) + ((uint64_t)op[6] * 10333695591267346857UL) + ((uint64_t)op[7] * 18320529210780447357UL) + ((((uint64_t)op[8] * 16034886971845528088UL) + ((uint64_t)op[9] * 15590941804916051280UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 15590941804916051280UL) + ((uint64_t)op[1] * 2584551659114682868UL) + ((uint64_t)op[2] * 14494164424885512764UL) + ((uint64_t)op[3] * 15740146837941994911UL) + ((uint64_t)op[4] * 15045817954949259965UL) + ((uint64_t)op[5] * 3774698560784541302UL) + ((uint64_t)op[6] * 2450953248114263220UL) + ((uint64_t)op[7] * 10333695591267346857UL) + ((uint64_t)op[8] * 18320529210780447357UL) + ((uint64_t)op[9] * 4823714203728047056UL);
	tmp_q[9] = ((uint64_t)op[0] * 16034886971845528088UL) + ((uint64_t)op[1] * 15590941804916051280UL) + ((uint64_t)op[2] * 2584551659114682868UL) + ((uint64_t)op[3] * 14494164424885512764UL) + ((uint64_t)op[4] * 15740146837941994911UL) + ((uint64_t)op[5] * 15045817954949259965UL) + ((uint64_t)op[6] * 3774698560784541302UL) + ((uint64_t)op[7] * 2450953248114263220UL) + ((uint64_t)op[8] * 10333695591267346857UL) + ((uint64_t)op[9] * 18320529210780447357UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 2441139206479357L) - ((-((int128)tmp_q[1] * 759297382984217L) + ((int128)tmp_q[2] * 4013611926011657L) - ((int128)tmp_q[3] * 2188532926095474L) + ((int128)tmp_q[4] * 656125764664006L) + ((int128)tmp_q[5] * 1598805414498436L) - ((int128)tmp_q[6] * 1347537404766780L) + ((int128)tmp_q[7] * 12319161511271L) + ((int128)tmp_q[8] * 2102589105442311L) - ((int128)tmp_q[9] * 115308423024421L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 115308423024421L) + ((int128)tmp_q[1] * 2441139206479357L) - ((-((int128)tmp_q[2] * 759297382984217L) + ((int128)tmp_q[3] * 4013611926011657L) - ((int128)tmp_q[4] * 2188532926095474L) + ((int128)tmp_q[5] * 656125764664006L) + ((int128)tmp_q[6] * 1598805414498436L) - ((int128)tmp_q[7] * 1347537404766780L) + ((int128)tmp_q[8] * 12319161511271L) + ((int128)tmp_q[9] * 2102589105442311L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 2102589105442311L) - ((int128)tmp_q[1] * 115308423024421L) + ((int128)tmp_q[2] * 2441139206479357L) - ((-((int128)tmp_q[3] * 759297382984217L) + ((int128)tmp_q[4] * 4013611926011657L) - ((int128)tmp_q[5] * 2188532926095474L) + ((int128)tmp_q[6] * 656125764664006L) + ((int128)tmp_q[7] * 1598805414498436L) - ((int128)tmp_q[8] * 1347537404766780L) + ((int128)tmp_q[9] * 12319161511271L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 12319161511271L) + ((int128)tmp_q[1] * 2102589105442311L) - ((int128)tmp_q[2] * 115308423024421L) + ((int128)tmp_q[3] * 2441139206479357L) - ((-((int128)tmp_q[4] * 759297382984217L) + ((int128)tmp_q[5] * 4013611926011657L) - ((int128)tmp_q[6] * 2188532926095474L) + ((int128)tmp_q[7] * 656125764664006L) + ((int128)tmp_q[8] * 1598805414498436L) - ((int128)tmp_q[9] * 1347537404766780L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 1347537404766780L) + ((int128)tmp_q[1] * 12319161511271L) + ((int128)tmp_q[2] * 2102589105442311L) - ((int128)tmp_q[3] * 115308423024421L) + ((int128)tmp_q[4] * 2441139206479357L) - ((-((int128)tmp_q[5] * 759297382984217L) + ((int128)tmp_q[6] * 4013611926011657L) - ((int128)tmp_q[7] * 2188532926095474L) + ((int128)tmp_q[8] * 656125764664006L) + ((int128)tmp_q[9] * 1598805414498436L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 1598805414498436L) - ((int128)tmp_q[1] * 1347537404766780L) + ((int128)tmp_q[2] * 12319161511271L) + ((int128)tmp_q[3] * 2102589105442311L) - ((int128)tmp_q[4] * 115308423024421L) + ((int128)tmp_q[5] * 2441139206479357L) - ((-((int128)tmp_q[6] * 759297382984217L) + ((int128)tmp_q[7] * 4013611926011657L) - ((int128)tmp_q[8] * 2188532926095474L) + ((int128)tmp_q[9] * 656125764664006L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 656125764664006L) + ((int128)tmp_q[1] * 1598805414498436L) - ((int128)tmp_q[2] * 1347537404766780L) + ((int128)tmp_q[3] * 12319161511271L) + ((int128)tmp_q[4] * 2102589105442311L) - ((int128)tmp_q[5] * 115308423024421L) + ((int128)tmp_q[6] * 2441139206479357L) - ((-((int128)tmp_q[7] * 759297382984217L) + ((int128)tmp_q[8] * 4013611926011657L) - ((int128)tmp_q[9] * 2188532926095474L)) * 2);
	tmp_zero[7] = -((int128)tmp_q[0] * 2188532926095474L) + ((int128)tmp_q[1] * 656125764664006L) + ((int128)tmp_q[2] * 1598805414498436L) - ((int128)tmp_q[3] * 1347537404766780L) + ((int128)tmp_q[4] * 12319161511271L) + ((int128)tmp_q[5] * 2102589105442311L) - ((int128)tmp_q[6] * 115308423024421L) + ((int128)tmp_q[7] * 2441139206479357L) - ((-((int128)tmp_q[8] * 759297382984217L) + ((int128)tmp_q[9] * 4013611926011657L)) * 2);
	tmp_zero[8] = ((int128)tmp_q[0] * 4013611926011657L) - ((int128)tmp_q[1] * 2188532926095474L) + ((int128)tmp_q[2] * 656125764664006L) + ((int128)tmp_q[3] * 1598805414498436L) - ((int128)tmp_q[4] * 1347537404766780L) + ((int128)tmp_q[5] * 12319161511271L) + ((int128)tmp_q[6] * 2102589105442311L) - ((int128)tmp_q[7] * 115308423024421L) + ((int128)tmp_q[8] * 2441139206479357L) + ((int128)tmp_q[9] * 1518594765968434L);
	tmp_zero[9] = -((int128)tmp_q[0] * 759297382984217L) + ((int128)tmp_q[1] * 4013611926011657L) - ((int128)tmp_q[2] * 2188532926095474L) + ((int128)tmp_q[3] * 656125764664006L) + ((int128)tmp_q[4] * 1598805414498436L) - ((int128)tmp_q[5] * 1347537404766780L) + ((int128)tmp_q[6] * 12319161511271L) + ((int128)tmp_q[7] * 2102589105442311L) - ((int128)tmp_q[8] * 115308423024421L) + ((int128)tmp_q[9] * 2441139206479357L);

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
}

