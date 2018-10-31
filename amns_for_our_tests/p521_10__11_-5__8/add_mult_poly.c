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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 5);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 5);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 5);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 5);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 10);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 5);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4959817874137510342UL) + ((((uint64_t)op[1] * 18193573600994438424UL) + ((uint64_t)op[2] * 10699650705715634880UL) + ((uint64_t)op[3] * 18129161245843848764UL) + ((uint64_t)op[4] * 14619657837089910189UL) + ((uint64_t)op[5] * 10250757443112781298UL) + ((uint64_t)op[6] * 11679079052174087203UL) + ((uint64_t)op[7] * 14691207031793242146UL) + ((uint64_t)op[8] * 12183690091653718862UL) + ((uint64_t)op[9] * 1555475769190180746UL) + ((uint64_t)op[10] * 5252935754210701341UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 5252935754210701341UL) + ((uint64_t)op[1] * 4959817874137510342UL) + ((((uint64_t)op[2] * 18193573600994438424UL) + ((uint64_t)op[3] * 10699650705715634880UL) + ((uint64_t)op[4] * 18129161245843848764UL) + ((uint64_t)op[5] * 14619657837089910189UL) + ((uint64_t)op[6] * 10250757443112781298UL) + ((uint64_t)op[7] * 11679079052174087203UL) + ((uint64_t)op[8] * 14691207031793242146UL) + ((uint64_t)op[9] * 12183690091653718862UL) + ((uint64_t)op[10] * 1555475769190180746UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 1555475769190180746UL) + ((uint64_t)op[1] * 5252935754210701341UL) + ((uint64_t)op[2] * 4959817874137510342UL) + ((((uint64_t)op[3] * 18193573600994438424UL) + ((uint64_t)op[4] * 10699650705715634880UL) + ((uint64_t)op[5] * 18129161245843848764UL) + ((uint64_t)op[6] * 14619657837089910189UL) + ((uint64_t)op[7] * 10250757443112781298UL) + ((uint64_t)op[8] * 11679079052174087203UL) + ((uint64_t)op[9] * 14691207031793242146UL) + ((uint64_t)op[10] * 12183690091653718862UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 12183690091653718862UL) + ((uint64_t)op[1] * 1555475769190180746UL) + ((uint64_t)op[2] * 5252935754210701341UL) + ((uint64_t)op[3] * 4959817874137510342UL) + ((((uint64_t)op[4] * 18193573600994438424UL) + ((uint64_t)op[5] * 10699650705715634880UL) + ((uint64_t)op[6] * 18129161245843848764UL) + ((uint64_t)op[7] * 14619657837089910189UL) + ((uint64_t)op[8] * 10250757443112781298UL) + ((uint64_t)op[9] * 11679079052174087203UL) + ((uint64_t)op[10] * 14691207031793242146UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 14691207031793242146UL) + ((uint64_t)op[1] * 12183690091653718862UL) + ((uint64_t)op[2] * 1555475769190180746UL) + ((uint64_t)op[3] * 5252935754210701341UL) + ((uint64_t)op[4] * 4959817874137510342UL) + ((((uint64_t)op[5] * 18193573600994438424UL) + ((uint64_t)op[6] * 10699650705715634880UL) + ((uint64_t)op[7] * 18129161245843848764UL) + ((uint64_t)op[8] * 14619657837089910189UL) + ((uint64_t)op[9] * 10250757443112781298UL) + ((uint64_t)op[10] * 11679079052174087203UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 11679079052174087203UL) + ((uint64_t)op[1] * 14691207031793242146UL) + ((uint64_t)op[2] * 12183690091653718862UL) + ((uint64_t)op[3] * 1555475769190180746UL) + ((uint64_t)op[4] * 5252935754210701341UL) + ((uint64_t)op[5] * 4959817874137510342UL) + ((((uint64_t)op[6] * 18193573600994438424UL) + ((uint64_t)op[7] * 10699650705715634880UL) + ((uint64_t)op[8] * 18129161245843848764UL) + ((uint64_t)op[9] * 14619657837089910189UL) + ((uint64_t)op[10] * 10250757443112781298UL)) * 18446744073709551611);
	tmp_q[6] = ((uint64_t)op[0] * 10250757443112781298UL) + ((uint64_t)op[1] * 11679079052174087203UL) + ((uint64_t)op[2] * 14691207031793242146UL) + ((uint64_t)op[3] * 12183690091653718862UL) + ((uint64_t)op[4] * 1555475769190180746UL) + ((uint64_t)op[5] * 5252935754210701341UL) + ((uint64_t)op[6] * 4959817874137510342UL) + ((((uint64_t)op[7] * 18193573600994438424UL) + ((uint64_t)op[8] * 10699650705715634880UL) + ((uint64_t)op[9] * 18129161245843848764UL) + ((uint64_t)op[10] * 14619657837089910189UL)) * 18446744073709551611);
	tmp_q[7] = ((uint64_t)op[0] * 14619657837089910189UL) + ((uint64_t)op[1] * 10250757443112781298UL) + ((uint64_t)op[2] * 11679079052174087203UL) + ((uint64_t)op[3] * 14691207031793242146UL) + ((uint64_t)op[4] * 12183690091653718862UL) + ((uint64_t)op[5] * 1555475769190180746UL) + ((uint64_t)op[6] * 5252935754210701341UL) + ((uint64_t)op[7] * 4959817874137510342UL) + ((((uint64_t)op[8] * 18193573600994438424UL) + ((uint64_t)op[9] * 10699650705715634880UL) + ((uint64_t)op[10] * 18129161245843848764UL)) * 18446744073709551611);
	tmp_q[8] = ((uint64_t)op[0] * 18129161245843848764UL) + ((uint64_t)op[1] * 14619657837089910189UL) + ((uint64_t)op[2] * 10250757443112781298UL) + ((uint64_t)op[3] * 11679079052174087203UL) + ((uint64_t)op[4] * 14691207031793242146UL) + ((uint64_t)op[5] * 12183690091653718862UL) + ((uint64_t)op[6] * 1555475769190180746UL) + ((uint64_t)op[7] * 5252935754210701341UL) + ((uint64_t)op[8] * 4959817874137510342UL) + ((((uint64_t)op[9] * 18193573600994438424UL) + ((uint64_t)op[10] * 10699650705715634880UL)) * 18446744073709551611);
	tmp_q[9] = ((uint64_t)op[0] * 10699650705715634880UL) + ((uint64_t)op[1] * 18129161245843848764UL) + ((uint64_t)op[2] * 14619657837089910189UL) + ((uint64_t)op[3] * 10250757443112781298UL) + ((uint64_t)op[4] * 11679079052174087203UL) + ((uint64_t)op[5] * 14691207031793242146UL) + ((uint64_t)op[6] * 12183690091653718862UL) + ((uint64_t)op[7] * 1555475769190180746UL) + ((uint64_t)op[8] * 5252935754210701341UL) + ((uint64_t)op[9] * 4959817874137510342UL) + ((uint64_t)op[10] * 1265852363575565960UL);
	tmp_q[10] = ((uint64_t)op[0] * 18193573600994438424UL) + ((uint64_t)op[1] * 10699650705715634880UL) + ((uint64_t)op[2] * 18129161245843848764UL) + ((uint64_t)op[3] * 14619657837089910189UL) + ((uint64_t)op[4] * 10250757443112781298UL) + ((uint64_t)op[5] * 11679079052174087203UL) + ((uint64_t)op[6] * 14691207031793242146UL) + ((uint64_t)op[7] * 12183690091653718862UL) + ((uint64_t)op[8] * 1555475769190180746UL) + ((uint64_t)op[9] * 5252935754210701341UL) + ((uint64_t)op[10] * 4959817874137510342UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 408079240564L) - ((((int128)tmp_q[1] * 9676229518693L) - ((int128)tmp_q[2] * 67336914501467L) - ((int128)tmp_q[3] * 74588753833730L) + ((int128)tmp_q[4] * 51035157264224L) + ((int128)tmp_q[5] * 7085185465783L) + ((int128)tmp_q[6] * 51481702309648L) - ((int128)tmp_q[7] * 119270724689543L) + ((int128)tmp_q[8] * 49298227078589L) - ((int128)tmp_q[9] * 80112063509891L) - ((int128)tmp_q[10] * 31060975047097L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 31060975047097L) + ((int128)tmp_q[1] * 408079240564L) - ((((int128)tmp_q[2] * 9676229518693L) - ((int128)tmp_q[3] * 67336914501467L) - ((int128)tmp_q[4] * 74588753833730L) + ((int128)tmp_q[5] * 51035157264224L) + ((int128)tmp_q[6] * 7085185465783L) + ((int128)tmp_q[7] * 51481702309648L) - ((int128)tmp_q[8] * 119270724689543L) + ((int128)tmp_q[9] * 49298227078589L) - ((int128)tmp_q[10] * 80112063509891L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 80112063509891L) - ((int128)tmp_q[1] * 31060975047097L) + ((int128)tmp_q[2] * 408079240564L) - ((((int128)tmp_q[3] * 9676229518693L) - ((int128)tmp_q[4] * 67336914501467L) - ((int128)tmp_q[5] * 74588753833730L) + ((int128)tmp_q[6] * 51035157264224L) + ((int128)tmp_q[7] * 7085185465783L) + ((int128)tmp_q[8] * 51481702309648L) - ((int128)tmp_q[9] * 119270724689543L) + ((int128)tmp_q[10] * 49298227078589L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 49298227078589L) - ((int128)tmp_q[1] * 80112063509891L) - ((int128)tmp_q[2] * 31060975047097L) + ((int128)tmp_q[3] * 408079240564L) - ((((int128)tmp_q[4] * 9676229518693L) - ((int128)tmp_q[5] * 67336914501467L) - ((int128)tmp_q[6] * 74588753833730L) + ((int128)tmp_q[7] * 51035157264224L) + ((int128)tmp_q[8] * 7085185465783L) + ((int128)tmp_q[9] * 51481702309648L) - ((int128)tmp_q[10] * 119270724689543L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 119270724689543L) + ((int128)tmp_q[1] * 49298227078589L) - ((int128)tmp_q[2] * 80112063509891L) - ((int128)tmp_q[3] * 31060975047097L) + ((int128)tmp_q[4] * 408079240564L) - ((((int128)tmp_q[5] * 9676229518693L) - ((int128)tmp_q[6] * 67336914501467L) - ((int128)tmp_q[7] * 74588753833730L) + ((int128)tmp_q[8] * 51035157264224L) + ((int128)tmp_q[9] * 7085185465783L) + ((int128)tmp_q[10] * 51481702309648L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 51481702309648L) - ((int128)tmp_q[1] * 119270724689543L) + ((int128)tmp_q[2] * 49298227078589L) - ((int128)tmp_q[3] * 80112063509891L) - ((int128)tmp_q[4] * 31060975047097L) + ((int128)tmp_q[5] * 408079240564L) - ((((int128)tmp_q[6] * 9676229518693L) - ((int128)tmp_q[7] * 67336914501467L) - ((int128)tmp_q[8] * 74588753833730L) + ((int128)tmp_q[9] * 51035157264224L) + ((int128)tmp_q[10] * 7085185465783L)) * 5);
	tmp_zero[6] = ((int128)tmp_q[0] * 7085185465783L) + ((int128)tmp_q[1] * 51481702309648L) - ((int128)tmp_q[2] * 119270724689543L) + ((int128)tmp_q[3] * 49298227078589L) - ((int128)tmp_q[4] * 80112063509891L) - ((int128)tmp_q[5] * 31060975047097L) + ((int128)tmp_q[6] * 408079240564L) - ((((int128)tmp_q[7] * 9676229518693L) - ((int128)tmp_q[8] * 67336914501467L) - ((int128)tmp_q[9] * 74588753833730L) + ((int128)tmp_q[10] * 51035157264224L)) * 5);
	tmp_zero[7] = ((int128)tmp_q[0] * 51035157264224L) + ((int128)tmp_q[1] * 7085185465783L) + ((int128)tmp_q[2] * 51481702309648L) - ((int128)tmp_q[3] * 119270724689543L) + ((int128)tmp_q[4] * 49298227078589L) - ((int128)tmp_q[5] * 80112063509891L) - ((int128)tmp_q[6] * 31060975047097L) + ((int128)tmp_q[7] * 408079240564L) - ((((int128)tmp_q[8] * 9676229518693L) - ((int128)tmp_q[9] * 67336914501467L) - ((int128)tmp_q[10] * 74588753833730L)) * 5);
	tmp_zero[8] = -((int128)tmp_q[0] * 74588753833730L) + ((int128)tmp_q[1] * 51035157264224L) + ((int128)tmp_q[2] * 7085185465783L) + ((int128)tmp_q[3] * 51481702309648L) - ((int128)tmp_q[4] * 119270724689543L) + ((int128)tmp_q[5] * 49298227078589L) - ((int128)tmp_q[6] * 80112063509891L) - ((int128)tmp_q[7] * 31060975047097L) + ((int128)tmp_q[8] * 408079240564L) - ((((int128)tmp_q[9] * 9676229518693L) - ((int128)tmp_q[10] * 67336914501467L)) * 5);
	tmp_zero[9] = -((int128)tmp_q[0] * 67336914501467L) - ((int128)tmp_q[1] * 74588753833730L) + ((int128)tmp_q[2] * 51035157264224L) + ((int128)tmp_q[3] * 7085185465783L) + ((int128)tmp_q[4] * 51481702309648L) - ((int128)tmp_q[5] * 119270724689543L) + ((int128)tmp_q[6] * 49298227078589L) - ((int128)tmp_q[7] * 80112063509891L) - ((int128)tmp_q[8] * 31060975047097L) + ((int128)tmp_q[9] * 408079240564L) - ((int128)tmp_q[10] * 48381147593465L);
	tmp_zero[10] = ((int128)tmp_q[0] * 9676229518693L) - ((int128)tmp_q[1] * 67336914501467L) - ((int128)tmp_q[2] * 74588753833730L) + ((int128)tmp_q[3] * 51035157264224L) + ((int128)tmp_q[4] * 7085185465783L) + ((int128)tmp_q[5] * 51481702309648L) - ((int128)tmp_q[6] * 119270724689543L) + ((int128)tmp_q[7] * 49298227078589L) - ((int128)tmp_q[8] * 80112063509891L) - ((int128)tmp_q[9] * 31060975047097L) + ((int128)tmp_q[10] * 408079240564L);

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

