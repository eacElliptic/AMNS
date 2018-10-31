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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) << 4);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10442919707002766693UL) + ((((uint64_t)op[1] * 14947164525477042353UL) + ((uint64_t)op[2] * 18062204397794232885UL) + ((uint64_t)op[3] * 15104184364375046874UL) + ((uint64_t)op[4] * 15158839697286561885UL) + ((uint64_t)op[5] * 2954669803412582023UL) + ((uint64_t)op[6] * 10105219247070817780UL) + ((uint64_t)op[7] * 3984435849483622649UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 3984435849483622649UL) + ((uint64_t)op[1] * 10442919707002766693UL) + ((((uint64_t)op[2] * 14947164525477042353UL) + ((uint64_t)op[3] * 18062204397794232885UL) + ((uint64_t)op[4] * 15104184364375046874UL) + ((uint64_t)op[5] * 15158839697286561885UL) + ((uint64_t)op[6] * 2954669803412582023UL) + ((uint64_t)op[7] * 10105219247070817780UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 10105219247070817780UL) + ((uint64_t)op[1] * 3984435849483622649UL) + ((uint64_t)op[2] * 10442919707002766693UL) + ((((uint64_t)op[3] * 14947164525477042353UL) + ((uint64_t)op[4] * 18062204397794232885UL) + ((uint64_t)op[5] * 15104184364375046874UL) + ((uint64_t)op[6] * 15158839697286561885UL) + ((uint64_t)op[7] * 2954669803412582023UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 2954669803412582023UL) + ((uint64_t)op[1] * 10105219247070817780UL) + ((uint64_t)op[2] * 3984435849483622649UL) + ((uint64_t)op[3] * 10442919707002766693UL) + ((((uint64_t)op[4] * 14947164525477042353UL) + ((uint64_t)op[5] * 18062204397794232885UL) + ((uint64_t)op[6] * 15104184364375046874UL) + ((uint64_t)op[7] * 15158839697286561885UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 15158839697286561885UL) + ((uint64_t)op[1] * 2954669803412582023UL) + ((uint64_t)op[2] * 10105219247070817780UL) + ((uint64_t)op[3] * 3984435849483622649UL) + ((uint64_t)op[4] * 10442919707002766693UL) + ((((uint64_t)op[5] * 14947164525477042353UL) + ((uint64_t)op[6] * 18062204397794232885UL) + ((uint64_t)op[7] * 15104184364375046874UL)) * 8);
	tmp_q[5] = ((uint64_t)op[0] * 15104184364375046874UL) + ((uint64_t)op[1] * 15158839697286561885UL) + ((uint64_t)op[2] * 2954669803412582023UL) + ((uint64_t)op[3] * 10105219247070817780UL) + ((uint64_t)op[4] * 3984435849483622649UL) + ((uint64_t)op[5] * 10442919707002766693UL) + ((((uint64_t)op[6] * 14947164525477042353UL) + ((uint64_t)op[7] * 18062204397794232885UL)) * 8);
	tmp_q[6] = ((uint64_t)op[0] * 18062204397794232885UL) + ((uint64_t)op[1] * 15104184364375046874UL) + ((uint64_t)op[2] * 15158839697286561885UL) + ((uint64_t)op[3] * 2954669803412582023UL) + ((uint64_t)op[4] * 10105219247070817780UL) + ((uint64_t)op[5] * 3984435849483622649UL) + ((uint64_t)op[6] * 10442919707002766693UL) + ((uint64_t)op[7] * 8896851761559029128UL);
	tmp_q[7] = ((uint64_t)op[0] * 14947164525477042353UL) + ((uint64_t)op[1] * 18062204397794232885UL) + ((uint64_t)op[2] * 15104184364375046874UL) + ((uint64_t)op[3] * 15158839697286561885UL) + ((uint64_t)op[4] * 2954669803412582023UL) + ((uint64_t)op[5] * 10105219247070817780UL) + ((uint64_t)op[6] * 3984435849483622649UL) + ((uint64_t)op[7] * 10442919707002766693UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 111675757435171L) + ((((int128)tmp_q[1] * 70218156552194L) - ((int128)tmp_q[2] * 80576502773874L) + ((int128)tmp_q[3] * 85251778344038L) - ((int128)tmp_q[4] * 6942477871586L) + ((int128)tmp_q[5] * 103206056842328L) + ((int128)tmp_q[6] * 73642965066023L) - ((int128)tmp_q[7] * 29275290993767L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 29275290993767L) + ((int128)tmp_q[1] * 111675757435171L) + ((((int128)tmp_q[2] * 70218156552194L) - ((int128)tmp_q[3] * 80576502773874L) + ((int128)tmp_q[4] * 85251778344038L) - ((int128)tmp_q[5] * 6942477871586L) + ((int128)tmp_q[6] * 103206056842328L) + ((int128)tmp_q[7] * 73642965066023L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 73642965066023L) - ((int128)tmp_q[1] * 29275290993767L) + ((int128)tmp_q[2] * 111675757435171L) + ((((int128)tmp_q[3] * 70218156552194L) - ((int128)tmp_q[4] * 80576502773874L) + ((int128)tmp_q[5] * 85251778344038L) - ((int128)tmp_q[6] * 6942477871586L) + ((int128)tmp_q[7] * 103206056842328L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 103206056842328L) + ((int128)tmp_q[1] * 73642965066023L) - ((int128)tmp_q[2] * 29275290993767L) + ((int128)tmp_q[3] * 111675757435171L) + ((((int128)tmp_q[4] * 70218156552194L) - ((int128)tmp_q[5] * 80576502773874L) + ((int128)tmp_q[6] * 85251778344038L) - ((int128)tmp_q[7] * 6942477871586L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 6942477871586L) + ((int128)tmp_q[1] * 103206056842328L) + ((int128)tmp_q[2] * 73642965066023L) - ((int128)tmp_q[3] * 29275290993767L) + ((int128)tmp_q[4] * 111675757435171L) + ((((int128)tmp_q[5] * 70218156552194L) - ((int128)tmp_q[6] * 80576502773874L) + ((int128)tmp_q[7] * 85251778344038L)) * 8);
	tmp_zero[5] = ((int128)tmp_q[0] * 85251778344038L) - ((int128)tmp_q[1] * 6942477871586L) + ((int128)tmp_q[2] * 103206056842328L) + ((int128)tmp_q[3] * 73642965066023L) - ((int128)tmp_q[4] * 29275290993767L) + ((int128)tmp_q[5] * 111675757435171L) + ((((int128)tmp_q[6] * 70218156552194L) - ((int128)tmp_q[7] * 80576502773874L)) * 8);
	tmp_zero[6] = -((int128)tmp_q[0] * 80576502773874L) + ((int128)tmp_q[1] * 85251778344038L) - ((int128)tmp_q[2] * 6942477871586L) + ((int128)tmp_q[3] * 103206056842328L) + ((int128)tmp_q[4] * 73642965066023L) - ((int128)tmp_q[5] * 29275290993767L) + ((int128)tmp_q[6] * 111675757435171L) + ((int128)tmp_q[7] * 561745252417552L);
	tmp_zero[7] = ((int128)tmp_q[0] * 70218156552194L) - ((int128)tmp_q[1] * 80576502773874L) + ((int128)tmp_q[2] * 85251778344038L) - ((int128)tmp_q[3] * 6942477871586L) + ((int128)tmp_q[4] * 103206056842328L) + ((int128)tmp_q[5] * 73642965066023L) - ((int128)tmp_q[6] * 29275290993767L) + ((int128)tmp_q[7] * 111675757435171L);

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

