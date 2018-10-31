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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) << 3);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 4);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 4);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 4);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 4);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) << 4);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) << 3);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17947023287055257507UL) + ((((uint64_t)op[1] * 16311789314014081027UL) + ((uint64_t)op[2] * 5783112667331440672UL) + ((uint64_t)op[3] * 3399182276153008025UL) + ((uint64_t)op[4] * 3988111934782238085UL) + ((uint64_t)op[5] * 8245687144552033108UL) + ((uint64_t)op[6] * 889251288853899086UL) + ((uint64_t)op[7] * 12110918440994315591UL) + ((uint64_t)op[8] * 2125018091213515981UL) + ((uint64_t)op[9] * 12443158838592292133UL) + ((uint64_t)op[10] * 2873667334930326986UL)) * 18446744073709551608);
	tmp_q[1] = ((uint64_t)op[0] * 2873667334930326986UL) + ((uint64_t)op[1] * 17947023287055257507UL) + ((((uint64_t)op[2] * 16311789314014081027UL) + ((uint64_t)op[3] * 5783112667331440672UL) + ((uint64_t)op[4] * 3399182276153008025UL) + ((uint64_t)op[5] * 3988111934782238085UL) + ((uint64_t)op[6] * 8245687144552033108UL) + ((uint64_t)op[7] * 889251288853899086UL) + ((uint64_t)op[8] * 12110918440994315591UL) + ((uint64_t)op[9] * 2125018091213515981UL) + ((uint64_t)op[10] * 12443158838592292133UL)) * 18446744073709551608);
	tmp_q[2] = ((uint64_t)op[0] * 12443158838592292133UL) + ((uint64_t)op[1] * 2873667334930326986UL) + ((uint64_t)op[2] * 17947023287055257507UL) + ((((uint64_t)op[3] * 16311789314014081027UL) + ((uint64_t)op[4] * 5783112667331440672UL) + ((uint64_t)op[5] * 3399182276153008025UL) + ((uint64_t)op[6] * 3988111934782238085UL) + ((uint64_t)op[7] * 8245687144552033108UL) + ((uint64_t)op[8] * 889251288853899086UL) + ((uint64_t)op[9] * 12110918440994315591UL) + ((uint64_t)op[10] * 2125018091213515981UL)) * 18446744073709551608);
	tmp_q[3] = ((uint64_t)op[0] * 2125018091213515981UL) + ((uint64_t)op[1] * 12443158838592292133UL) + ((uint64_t)op[2] * 2873667334930326986UL) + ((uint64_t)op[3] * 17947023287055257507UL) + ((((uint64_t)op[4] * 16311789314014081027UL) + ((uint64_t)op[5] * 5783112667331440672UL) + ((uint64_t)op[6] * 3399182276153008025UL) + ((uint64_t)op[7] * 3988111934782238085UL) + ((uint64_t)op[8] * 8245687144552033108UL) + ((uint64_t)op[9] * 889251288853899086UL) + ((uint64_t)op[10] * 12110918440994315591UL)) * 18446744073709551608);
	tmp_q[4] = ((uint64_t)op[0] * 12110918440994315591UL) + ((uint64_t)op[1] * 2125018091213515981UL) + ((uint64_t)op[2] * 12443158838592292133UL) + ((uint64_t)op[3] * 2873667334930326986UL) + ((uint64_t)op[4] * 17947023287055257507UL) + ((((uint64_t)op[5] * 16311789314014081027UL) + ((uint64_t)op[6] * 5783112667331440672UL) + ((uint64_t)op[7] * 3399182276153008025UL) + ((uint64_t)op[8] * 3988111934782238085UL) + ((uint64_t)op[9] * 8245687144552033108UL) + ((uint64_t)op[10] * 889251288853899086UL)) * 18446744073709551608);
	tmp_q[5] = ((uint64_t)op[0] * 889251288853899086UL) + ((uint64_t)op[1] * 12110918440994315591UL) + ((uint64_t)op[2] * 2125018091213515981UL) + ((uint64_t)op[3] * 12443158838592292133UL) + ((uint64_t)op[4] * 2873667334930326986UL) + ((uint64_t)op[5] * 17947023287055257507UL) + ((((uint64_t)op[6] * 16311789314014081027UL) + ((uint64_t)op[7] * 5783112667331440672UL) + ((uint64_t)op[8] * 3399182276153008025UL) + ((uint64_t)op[9] * 3988111934782238085UL) + ((uint64_t)op[10] * 8245687144552033108UL)) * 18446744073709551608);
	tmp_q[6] = ((uint64_t)op[0] * 8245687144552033108UL) + ((uint64_t)op[1] * 889251288853899086UL) + ((uint64_t)op[2] * 12110918440994315591UL) + ((uint64_t)op[3] * 2125018091213515981UL) + ((uint64_t)op[4] * 12443158838592292133UL) + ((uint64_t)op[5] * 2873667334930326986UL) + ((uint64_t)op[6] * 17947023287055257507UL) + ((((uint64_t)op[7] * 16311789314014081027UL) + ((uint64_t)op[8] * 5783112667331440672UL) + ((uint64_t)op[9] * 3399182276153008025UL) + ((uint64_t)op[10] * 3988111934782238085UL)) * 18446744073709551608);
	tmp_q[7] = ((uint64_t)op[0] * 3988111934782238085UL) + ((uint64_t)op[1] * 8245687144552033108UL) + ((uint64_t)op[2] * 889251288853899086UL) + ((uint64_t)op[3] * 12110918440994315591UL) + ((uint64_t)op[4] * 2125018091213515981UL) + ((uint64_t)op[5] * 12443158838592292133UL) + ((uint64_t)op[6] * 2873667334930326986UL) + ((uint64_t)op[7] * 17947023287055257507UL) + ((((uint64_t)op[8] * 16311789314014081027UL) + ((uint64_t)op[9] * 5783112667331440672UL) + ((uint64_t)op[10] * 3399182276153008025UL)) * 18446744073709551608);
	tmp_q[8] = ((uint64_t)op[0] * 3399182276153008025UL) + ((uint64_t)op[1] * 3988111934782238085UL) + ((uint64_t)op[2] * 8245687144552033108UL) + ((uint64_t)op[3] * 889251288853899086UL) + ((uint64_t)op[4] * 12110918440994315591UL) + ((uint64_t)op[5] * 2125018091213515981UL) + ((uint64_t)op[6] * 12443158838592292133UL) + ((uint64_t)op[7] * 2873667334930326986UL) + ((uint64_t)op[8] * 17947023287055257507UL) + ((((uint64_t)op[9] * 16311789314014081027UL) + ((uint64_t)op[10] * 5783112667331440672UL)) * 18446744073709551608);
	tmp_q[9] = ((uint64_t)op[0] * 5783112667331440672UL) + ((uint64_t)op[1] * 3399182276153008025UL) + ((uint64_t)op[2] * 3988111934782238085UL) + ((uint64_t)op[3] * 8245687144552033108UL) + ((uint64_t)op[4] * 889251288853899086UL) + ((uint64_t)op[5] * 12110918440994315591UL) + ((uint64_t)op[6] * 2125018091213515981UL) + ((uint64_t)op[7] * 12443158838592292133UL) + ((uint64_t)op[8] * 2873667334930326986UL) + ((uint64_t)op[9] * 17947023287055257507UL) + ((uint64_t)op[10] * 17079638077563764712UL);
	tmp_q[10] = ((uint64_t)op[0] * 16311789314014081027UL) + ((uint64_t)op[1] * 5783112667331440672UL) + ((uint64_t)op[2] * 3399182276153008025UL) + ((uint64_t)op[3] * 3988111934782238085UL) + ((uint64_t)op[4] * 8245687144552033108UL) + ((uint64_t)op[5] * 889251288853899086UL) + ((uint64_t)op[6] * 12110918440994315591UL) + ((uint64_t)op[7] * 2125018091213515981UL) + ((uint64_t)op[8] * 12443158838592292133UL) + ((uint64_t)op[9] * 2873667334930326986UL) + ((uint64_t)op[10] * 17947023287055257507UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4639737410413L) - ((-((int128)tmp_q[1] * 19737558950238L) - ((int128)tmp_q[2] * 90775189087549L) - ((int128)tmp_q[3] * 51935287322981L) - ((int128)tmp_q[4] * 44085576686072L) - ((int128)tmp_q[5] * 90888659537356L) - ((int128)tmp_q[6] * 79663036808130L) - ((int128)tmp_q[7] * 21567196713756L) + ((int128)tmp_q[8] * 18244532947617L) - ((int128)tmp_q[9] * 30205871604599L) + ((int128)tmp_q[10] * 550537383922L)) * 8);
	tmp_zero[1] = ((int128)tmp_q[0] * 550537383922L) + ((int128)tmp_q[1] * 4639737410413L) - ((-((int128)tmp_q[2] * 19737558950238L) - ((int128)tmp_q[3] * 90775189087549L) - ((int128)tmp_q[4] * 51935287322981L) - ((int128)tmp_q[5] * 44085576686072L) - ((int128)tmp_q[6] * 90888659537356L) - ((int128)tmp_q[7] * 79663036808130L) - ((int128)tmp_q[8] * 21567196713756L) + ((int128)tmp_q[9] * 18244532947617L) - ((int128)tmp_q[10] * 30205871604599L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 30205871604599L) + ((int128)tmp_q[1] * 550537383922L) + ((int128)tmp_q[2] * 4639737410413L) - ((-((int128)tmp_q[3] * 19737558950238L) - ((int128)tmp_q[4] * 90775189087549L) - ((int128)tmp_q[5] * 51935287322981L) - ((int128)tmp_q[6] * 44085576686072L) - ((int128)tmp_q[7] * 90888659537356L) - ((int128)tmp_q[8] * 79663036808130L) - ((int128)tmp_q[9] * 21567196713756L) + ((int128)tmp_q[10] * 18244532947617L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 18244532947617L) - ((int128)tmp_q[1] * 30205871604599L) + ((int128)tmp_q[2] * 550537383922L) + ((int128)tmp_q[3] * 4639737410413L) - ((-((int128)tmp_q[4] * 19737558950238L) - ((int128)tmp_q[5] * 90775189087549L) - ((int128)tmp_q[6] * 51935287322981L) - ((int128)tmp_q[7] * 44085576686072L) - ((int128)tmp_q[8] * 90888659537356L) - ((int128)tmp_q[9] * 79663036808130L) - ((int128)tmp_q[10] * 21567196713756L)) * 8);
	tmp_zero[4] = -((int128)tmp_q[0] * 21567196713756L) + ((int128)tmp_q[1] * 18244532947617L) - ((int128)tmp_q[2] * 30205871604599L) + ((int128)tmp_q[3] * 550537383922L) + ((int128)tmp_q[4] * 4639737410413L) - ((-((int128)tmp_q[5] * 19737558950238L) - ((int128)tmp_q[6] * 90775189087549L) - ((int128)tmp_q[7] * 51935287322981L) - ((int128)tmp_q[8] * 44085576686072L) - ((int128)tmp_q[9] * 90888659537356L) - ((int128)tmp_q[10] * 79663036808130L)) * 8);
	tmp_zero[5] = -((int128)tmp_q[0] * 79663036808130L) - ((int128)tmp_q[1] * 21567196713756L) + ((int128)tmp_q[2] * 18244532947617L) - ((int128)tmp_q[3] * 30205871604599L) + ((int128)tmp_q[4] * 550537383922L) + ((int128)tmp_q[5] * 4639737410413L) - ((-((int128)tmp_q[6] * 19737558950238L) - ((int128)tmp_q[7] * 90775189087549L) - ((int128)tmp_q[8] * 51935287322981L) - ((int128)tmp_q[9] * 44085576686072L) - ((int128)tmp_q[10] * 90888659537356L)) * 8);
	tmp_zero[6] = -((int128)tmp_q[0] * 90888659537356L) - ((int128)tmp_q[1] * 79663036808130L) - ((int128)tmp_q[2] * 21567196713756L) + ((int128)tmp_q[3] * 18244532947617L) - ((int128)tmp_q[4] * 30205871604599L) + ((int128)tmp_q[5] * 550537383922L) + ((int128)tmp_q[6] * 4639737410413L) - ((-((int128)tmp_q[7] * 19737558950238L) - ((int128)tmp_q[8] * 90775189087549L) - ((int128)tmp_q[9] * 51935287322981L) - ((int128)tmp_q[10] * 44085576686072L)) * 8);
	tmp_zero[7] = -((int128)tmp_q[0] * 44085576686072L) - ((int128)tmp_q[1] * 90888659537356L) - ((int128)tmp_q[2] * 79663036808130L) - ((int128)tmp_q[3] * 21567196713756L) + ((int128)tmp_q[4] * 18244532947617L) - ((int128)tmp_q[5] * 30205871604599L) + ((int128)tmp_q[6] * 550537383922L) + ((int128)tmp_q[7] * 4639737410413L) - ((-((int128)tmp_q[8] * 19737558950238L) - ((int128)tmp_q[9] * 90775189087549L) - ((int128)tmp_q[10] * 51935287322981L)) * 8);
	tmp_zero[8] = -((int128)tmp_q[0] * 51935287322981L) - ((int128)tmp_q[1] * 44085576686072L) - ((int128)tmp_q[2] * 90888659537356L) - ((int128)tmp_q[3] * 79663036808130L) - ((int128)tmp_q[4] * 21567196713756L) + ((int128)tmp_q[5] * 18244532947617L) - ((int128)tmp_q[6] * 30205871604599L) + ((int128)tmp_q[7] * 550537383922L) + ((int128)tmp_q[8] * 4639737410413L) - ((-((int128)tmp_q[9] * 19737558950238L) - ((int128)tmp_q[10] * 90775189087549L)) * 8);
	tmp_zero[9] = -((int128)tmp_q[0] * 90775189087549L) - ((int128)tmp_q[1] * 51935287322981L) - ((int128)tmp_q[2] * 44085576686072L) - ((int128)tmp_q[3] * 90888659537356L) - ((int128)tmp_q[4] * 79663036808130L) - ((int128)tmp_q[5] * 21567196713756L) + ((int128)tmp_q[6] * 18244532947617L) - ((int128)tmp_q[7] * 30205871604599L) + ((int128)tmp_q[8] * 550537383922L) + ((int128)tmp_q[9] * 4639737410413L) + ((int128)tmp_q[10] * 157900471601904L);
	tmp_zero[10] = -((int128)tmp_q[0] * 19737558950238L) - ((int128)tmp_q[1] * 90775189087549L) - ((int128)tmp_q[2] * 51935287322981L) - ((int128)tmp_q[3] * 44085576686072L) - ((int128)tmp_q[4] * 90888659537356L) - ((int128)tmp_q[5] * 79663036808130L) - ((int128)tmp_q[6] * 21567196713756L) + ((int128)tmp_q[7] * 18244532947617L) - ((int128)tmp_q[8] * 30205871604599L) + ((int128)tmp_q[9] * 550537383922L) + ((int128)tmp_q[10] * 4639737410413L);

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

