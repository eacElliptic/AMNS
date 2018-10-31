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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3960051668190589333UL) + ((((uint64_t)op[1] * 13607290401169728001UL) + ((uint64_t)op[2] * 14728504804155588835UL) + ((uint64_t)op[3] * 5352988330204528503UL) + ((uint64_t)op[4] * 3825831201505954966UL) + ((uint64_t)op[5] * 15674275639148873508UL) + ((uint64_t)op[6] * 526124909351925477UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 526124909351925477UL) + ((uint64_t)op[1] * 3960051668190589333UL) + ((((uint64_t)op[2] * 13607290401169728001UL) + ((uint64_t)op[3] * 14728504804155588835UL) + ((uint64_t)op[4] * 5352988330204528503UL) + ((uint64_t)op[5] * 3825831201505954966UL) + ((uint64_t)op[6] * 15674275639148873508UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 15674275639148873508UL) + ((uint64_t)op[1] * 526124909351925477UL) + ((uint64_t)op[2] * 3960051668190589333UL) + ((((uint64_t)op[3] * 13607290401169728001UL) + ((uint64_t)op[4] * 14728504804155588835UL) + ((uint64_t)op[5] * 5352988330204528503UL) + ((uint64_t)op[6] * 3825831201505954966UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 3825831201505954966UL) + ((uint64_t)op[1] * 15674275639148873508UL) + ((uint64_t)op[2] * 526124909351925477UL) + ((uint64_t)op[3] * 3960051668190589333UL) + ((((uint64_t)op[4] * 13607290401169728001UL) + ((uint64_t)op[5] * 14728504804155588835UL) + ((uint64_t)op[6] * 5352988330204528503UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 5352988330204528503UL) + ((uint64_t)op[1] * 3825831201505954966UL) + ((uint64_t)op[2] * 15674275639148873508UL) + ((uint64_t)op[3] * 526124909351925477UL) + ((uint64_t)op[4] * 3960051668190589333UL) + ((((uint64_t)op[5] * 13607290401169728001UL) + ((uint64_t)op[6] * 14728504804155588835UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 14728504804155588835UL) + ((uint64_t)op[1] * 5352988330204528503UL) + ((uint64_t)op[2] * 3825831201505954966UL) + ((uint64_t)op[3] * 15674275639148873508UL) + ((uint64_t)op[4] * 526124909351925477UL) + ((uint64_t)op[5] * 3960051668190589333UL) + ((uint64_t)op[6] * 8767836728629904386UL);
	tmp_q[6] = ((uint64_t)op[0] * 13607290401169728001UL) + ((uint64_t)op[1] * 14728504804155588835UL) + ((uint64_t)op[2] * 5352988330204528503UL) + ((uint64_t)op[3] * 3825831201505954966UL) + ((uint64_t)op[4] * 15674275639148873508UL) + ((uint64_t)op[5] * 526124909351925477UL) + ((uint64_t)op[6] * 3960051668190589333UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 11583220573030033L) + ((((int128)tmp_q[1] * 8380177963329737L) + ((int128)tmp_q[2] * 2722631030593256L) + ((int128)tmp_q[3] * 12072269700435860L) - ((int128)tmp_q[4] * 14218384034243185L) + ((int128)tmp_q[5] * 19952194539596893L) + ((int128)tmp_q[6] * 25213882820509L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 25213882820509L) - ((int128)tmp_q[1] * 11583220573030033L) + ((((int128)tmp_q[2] * 8380177963329737L) + ((int128)tmp_q[3] * 2722631030593256L) + ((int128)tmp_q[4] * 12072269700435860L) - ((int128)tmp_q[5] * 14218384034243185L) + ((int128)tmp_q[6] * 19952194539596893L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 19952194539596893L) + ((int128)tmp_q[1] * 25213882820509L) - ((int128)tmp_q[2] * 11583220573030033L) + ((((int128)tmp_q[3] * 8380177963329737L) + ((int128)tmp_q[4] * 2722631030593256L) + ((int128)tmp_q[5] * 12072269700435860L) - ((int128)tmp_q[6] * 14218384034243185L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 14218384034243185L) + ((int128)tmp_q[1] * 19952194539596893L) + ((int128)tmp_q[2] * 25213882820509L) - ((int128)tmp_q[3] * 11583220573030033L) + ((((int128)tmp_q[4] * 8380177963329737L) + ((int128)tmp_q[5] * 2722631030593256L) + ((int128)tmp_q[6] * 12072269700435860L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 12072269700435860L) - ((int128)tmp_q[1] * 14218384034243185L) + ((int128)tmp_q[2] * 19952194539596893L) + ((int128)tmp_q[3] * 25213882820509L) - ((int128)tmp_q[4] * 11583220573030033L) + ((((int128)tmp_q[5] * 8380177963329737L) + ((int128)tmp_q[6] * 2722631030593256L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 2722631030593256L) + ((int128)tmp_q[1] * 12072269700435860L) - ((int128)tmp_q[2] * 14218384034243185L) + ((int128)tmp_q[3] * 19952194539596893L) + ((int128)tmp_q[4] * 25213882820509L) - ((int128)tmp_q[5] * 11583220573030033L) + ((int128)tmp_q[6] * 16760355926659474L);
	tmp_zero[6] = ((int128)tmp_q[0] * 8380177963329737L) + ((int128)tmp_q[1] * 2722631030593256L) + ((int128)tmp_q[2] * 12072269700435860L) - ((int128)tmp_q[3] * 14218384034243185L) + ((int128)tmp_q[4] * 19952194539596893L) + ((int128)tmp_q[5] * 25213882820509L) - ((int128)tmp_q[6] * 11583220573030033L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

