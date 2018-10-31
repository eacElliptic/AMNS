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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 14789674838022521054UL) + ((((uint64_t)op[1] * 16104516398327938680UL) + ((uint64_t)op[2] * 9744049796979293571UL) + ((uint64_t)op[3] * 7152085864099246027UL) + ((uint64_t)op[4] * 16724600013994344178UL) + ((uint64_t)op[5] * 5595590847197166008UL) + ((uint64_t)op[6] * 17489786964445625165UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 17489786964445625165UL) + ((uint64_t)op[1] * 14789674838022521054UL) + ((((uint64_t)op[2] * 16104516398327938680UL) + ((uint64_t)op[3] * 9744049796979293571UL) + ((uint64_t)op[4] * 7152085864099246027UL) + ((uint64_t)op[5] * 16724600013994344178UL) + ((uint64_t)op[6] * 5595590847197166008UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 5595590847197166008UL) + ((uint64_t)op[1] * 17489786964445625165UL) + ((uint64_t)op[2] * 14789674838022521054UL) + ((((uint64_t)op[3] * 16104516398327938680UL) + ((uint64_t)op[4] * 9744049796979293571UL) + ((uint64_t)op[5] * 7152085864099246027UL) + ((uint64_t)op[6] * 16724600013994344178UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 16724600013994344178UL) + ((uint64_t)op[1] * 5595590847197166008UL) + ((uint64_t)op[2] * 17489786964445625165UL) + ((uint64_t)op[3] * 14789674838022521054UL) + ((((uint64_t)op[4] * 16104516398327938680UL) + ((uint64_t)op[5] * 9744049796979293571UL) + ((uint64_t)op[6] * 7152085864099246027UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 7152085864099246027UL) + ((uint64_t)op[1] * 16724600013994344178UL) + ((uint64_t)op[2] * 5595590847197166008UL) + ((uint64_t)op[3] * 17489786964445625165UL) + ((uint64_t)op[4] * 14789674838022521054UL) + ((((uint64_t)op[5] * 16104516398327938680UL) + ((uint64_t)op[6] * 9744049796979293571UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 9744049796979293571UL) + ((uint64_t)op[1] * 7152085864099246027UL) + ((uint64_t)op[2] * 16724600013994344178UL) + ((uint64_t)op[3] * 5595590847197166008UL) + ((uint64_t)op[4] * 17489786964445625165UL) + ((uint64_t)op[5] * 14789674838022521054UL) + ((uint64_t)op[6] * 11711138376908064680UL);
	tmp_q[6] = ((uint64_t)op[0] * 16104516398327938680UL) + ((uint64_t)op[1] * 9744049796979293571UL) + ((uint64_t)op[2] * 7152085864099246027UL) + ((uint64_t)op[3] * 16724600013994344178UL) + ((uint64_t)op[4] * 5595590847197166008UL) + ((uint64_t)op[5] * 17489786964445625165UL) + ((uint64_t)op[6] * 14789674838022521054UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 30209025138L) - ((-((int128)tmp_q[1] * 74053871389L) + ((int128)tmp_q[2] * 42659805242L) - ((int128)tmp_q[3] * 19651489807L) + ((int128)tmp_q[4] * 8402924103L) + ((int128)tmp_q[5] * 74112181681L) - ((int128)tmp_q[6] * 12522457263L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 12522457263L) - ((int128)tmp_q[1] * 30209025138L) - ((-((int128)tmp_q[2] * 74053871389L) + ((int128)tmp_q[3] * 42659805242L) - ((int128)tmp_q[4] * 19651489807L) + ((int128)tmp_q[5] * 8402924103L) + ((int128)tmp_q[6] * 74112181681L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 74112181681L) - ((int128)tmp_q[1] * 12522457263L) - ((int128)tmp_q[2] * 30209025138L) - ((-((int128)tmp_q[3] * 74053871389L) + ((int128)tmp_q[4] * 42659805242L) - ((int128)tmp_q[5] * 19651489807L) + ((int128)tmp_q[6] * 8402924103L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 8402924103L) + ((int128)tmp_q[1] * 74112181681L) - ((int128)tmp_q[2] * 12522457263L) - ((int128)tmp_q[3] * 30209025138L) - ((-((int128)tmp_q[4] * 74053871389L) + ((int128)tmp_q[5] * 42659805242L) - ((int128)tmp_q[6] * 19651489807L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 19651489807L) + ((int128)tmp_q[1] * 8402924103L) + ((int128)tmp_q[2] * 74112181681L) - ((int128)tmp_q[3] * 12522457263L) - ((int128)tmp_q[4] * 30209025138L) - ((-((int128)tmp_q[5] * 74053871389L) + ((int128)tmp_q[6] * 42659805242L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 42659805242L) - ((int128)tmp_q[1] * 19651489807L) + ((int128)tmp_q[2] * 8402924103L) + ((int128)tmp_q[3] * 74112181681L) - ((int128)tmp_q[4] * 12522457263L) - ((int128)tmp_q[5] * 30209025138L) + ((int128)tmp_q[6] * 370269356945L);
	tmp_zero[6] = -((int128)tmp_q[0] * 74053871389L) + ((int128)tmp_q[1] * 42659805242L) - ((int128)tmp_q[2] * 19651489807L) + ((int128)tmp_q[3] * 8402924103L) + ((int128)tmp_q[4] * 74112181681L) - ((int128)tmp_q[5] * 12522457263L) - ((int128)tmp_q[6] * 30209025138L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

