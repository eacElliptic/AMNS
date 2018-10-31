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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18184940075136867810UL) + ((((uint64_t)op[1] * 7134202182633307359UL) + ((uint64_t)op[2] * 10820901318892532685UL) + ((uint64_t)op[3] * 6040626948313468055UL) + ((uint64_t)op[4] * 16386756582142359467UL) + ((uint64_t)op[5] * 9610133863819972183UL) + ((uint64_t)op[6] * 72502929421284068UL) + ((uint64_t)op[7] * 10899878568120147014UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 10899878568120147014UL) + ((uint64_t)op[1] * 18184940075136867810UL) + ((((uint64_t)op[2] * 7134202182633307359UL) + ((uint64_t)op[3] * 10820901318892532685UL) + ((uint64_t)op[4] * 6040626948313468055UL) + ((uint64_t)op[5] * 16386756582142359467UL) + ((uint64_t)op[6] * 9610133863819972183UL) + ((uint64_t)op[7] * 72502929421284068UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 72502929421284068UL) + ((uint64_t)op[1] * 10899878568120147014UL) + ((uint64_t)op[2] * 18184940075136867810UL) + ((((uint64_t)op[3] * 7134202182633307359UL) + ((uint64_t)op[4] * 10820901318892532685UL) + ((uint64_t)op[5] * 6040626948313468055UL) + ((uint64_t)op[6] * 16386756582142359467UL) + ((uint64_t)op[7] * 9610133863819972183UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 9610133863819972183UL) + ((uint64_t)op[1] * 72502929421284068UL) + ((uint64_t)op[2] * 10899878568120147014UL) + ((uint64_t)op[3] * 18184940075136867810UL) + ((((uint64_t)op[4] * 7134202182633307359UL) + ((uint64_t)op[5] * 10820901318892532685UL) + ((uint64_t)op[6] * 6040626948313468055UL) + ((uint64_t)op[7] * 16386756582142359467UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 16386756582142359467UL) + ((uint64_t)op[1] * 9610133863819972183UL) + ((uint64_t)op[2] * 72502929421284068UL) + ((uint64_t)op[3] * 10899878568120147014UL) + ((uint64_t)op[4] * 18184940075136867810UL) + ((((uint64_t)op[5] * 7134202182633307359UL) + ((uint64_t)op[6] * 10820901318892532685UL) + ((uint64_t)op[7] * 6040626948313468055UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 6040626948313468055UL) + ((uint64_t)op[1] * 16386756582142359467UL) + ((uint64_t)op[2] * 9610133863819972183UL) + ((uint64_t)op[3] * 72502929421284068UL) + ((uint64_t)op[4] * 10899878568120147014UL) + ((uint64_t)op[5] * 18184940075136867810UL) + ((((uint64_t)op[6] * 7134202182633307359UL) + ((uint64_t)op[7] * 10820901318892532685UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 10820901318892532685UL) + ((uint64_t)op[1] * 6040626948313468055UL) + ((uint64_t)op[2] * 16386756582142359467UL) + ((uint64_t)op[3] * 9610133863819972183UL) + ((uint64_t)op[4] * 72502929421284068UL) + ((uint64_t)op[5] * 10899878568120147014UL) + ((uint64_t)op[6] * 18184940075136867810UL) + ((uint64_t)op[7] * 15490881599519181155UL);
	tmp_q[7] = ((uint64_t)op[0] * 7134202182633307359UL) + ((uint64_t)op[1] * 10820901318892532685UL) + ((uint64_t)op[2] * 6040626948313468055UL) + ((uint64_t)op[3] * 16386756582142359467UL) + ((uint64_t)op[4] * 9610133863819972183UL) + ((uint64_t)op[5] * 72502929421284068UL) + ((uint64_t)op[6] * 10899878568120147014UL) + ((uint64_t)op[7] * 18184940075136867810UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 15581726171473L) - ((((int128)tmp_q[1] * 122403907074404L) - ((int128)tmp_q[2] * 39148902052091L) - ((int128)tmp_q[3] * 65287096788708L) + ((int128)tmp_q[4] * 48431393442890L) + ((int128)tmp_q[5] * 28555699510717L) - ((int128)tmp_q[6] * 70805152811712L) + ((int128)tmp_q[7] * 108127136314404L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 108127136314404L) - ((int128)tmp_q[1] * 15581726171473L) - ((((int128)tmp_q[2] * 122403907074404L) - ((int128)tmp_q[3] * 39148902052091L) - ((int128)tmp_q[4] * 65287096788708L) + ((int128)tmp_q[5] * 48431393442890L) + ((int128)tmp_q[6] * 28555699510717L) - ((int128)tmp_q[7] * 70805152811712L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 70805152811712L) + ((int128)tmp_q[1] * 108127136314404L) - ((int128)tmp_q[2] * 15581726171473L) - ((((int128)tmp_q[3] * 122403907074404L) - ((int128)tmp_q[4] * 39148902052091L) - ((int128)tmp_q[5] * 65287096788708L) + ((int128)tmp_q[6] * 48431393442890L) + ((int128)tmp_q[7] * 28555699510717L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 28555699510717L) - ((int128)tmp_q[1] * 70805152811712L) + ((int128)tmp_q[2] * 108127136314404L) - ((int128)tmp_q[3] * 15581726171473L) - ((((int128)tmp_q[4] * 122403907074404L) - ((int128)tmp_q[5] * 39148902052091L) - ((int128)tmp_q[6] * 65287096788708L) + ((int128)tmp_q[7] * 48431393442890L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 48431393442890L) + ((int128)tmp_q[1] * 28555699510717L) - ((int128)tmp_q[2] * 70805152811712L) + ((int128)tmp_q[3] * 108127136314404L) - ((int128)tmp_q[4] * 15581726171473L) - ((((int128)tmp_q[5] * 122403907074404L) - ((int128)tmp_q[6] * 39148902052091L) - ((int128)tmp_q[7] * 65287096788708L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 65287096788708L) + ((int128)tmp_q[1] * 48431393442890L) + ((int128)tmp_q[2] * 28555699510717L) - ((int128)tmp_q[3] * 70805152811712L) + ((int128)tmp_q[4] * 108127136314404L) - ((int128)tmp_q[5] * 15581726171473L) - ((((int128)tmp_q[6] * 122403907074404L) - ((int128)tmp_q[7] * 39148902052091L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 39148902052091L) - ((int128)tmp_q[1] * 65287096788708L) + ((int128)tmp_q[2] * 48431393442890L) + ((int128)tmp_q[3] * 28555699510717L) - ((int128)tmp_q[4] * 70805152811712L) + ((int128)tmp_q[5] * 108127136314404L) - ((int128)tmp_q[6] * 15581726171473L) - ((int128)tmp_q[7] * 367211721223212L);
	tmp_zero[7] = ((int128)tmp_q[0] * 122403907074404L) - ((int128)tmp_q[1] * 39148902052091L) - ((int128)tmp_q[2] * 65287096788708L) + ((int128)tmp_q[3] * 48431393442890L) + ((int128)tmp_q[4] * 28555699510717L) - ((int128)tmp_q[5] * 70805152811712L) + ((int128)tmp_q[6] * 108127136314404L) - ((int128)tmp_q[7] * 15581726171473L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

