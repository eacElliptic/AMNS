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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 6);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 6);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[10]) * 6);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 12);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 12);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 12);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 12);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[10] * pa[9]) * 12);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[10] * pa[10]) * 6);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 5658328008036602655UL) + ((((uint64_t)op[1] * 4867231046790035161UL) + ((uint64_t)op[2] * 9306810177458241629UL) + ((uint64_t)op[3] * 4675441231982960697UL) + ((uint64_t)op[4] * 4763906738431177706UL) + ((uint64_t)op[5] * 683745108849402065UL) + ((uint64_t)op[6] * 11502197098927065429UL) + ((uint64_t)op[7] * 1502077965055206585UL) + ((uint64_t)op[8] * 12293104920482715914UL) + ((uint64_t)op[9] * 12824144572736731398UL) + ((uint64_t)op[10] * 17068822549925918684UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 17068822549925918684UL) + ((uint64_t)op[1] * 5658328008036602655UL) + ((((uint64_t)op[2] * 4867231046790035161UL) + ((uint64_t)op[3] * 9306810177458241629UL) + ((uint64_t)op[4] * 4675441231982960697UL) + ((uint64_t)op[5] * 4763906738431177706UL) + ((uint64_t)op[6] * 683745108849402065UL) + ((uint64_t)op[7] * 11502197098927065429UL) + ((uint64_t)op[8] * 1502077965055206585UL) + ((uint64_t)op[9] * 12293104920482715914UL) + ((uint64_t)op[10] * 12824144572736731398UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 12824144572736731398UL) + ((uint64_t)op[1] * 17068822549925918684UL) + ((uint64_t)op[2] * 5658328008036602655UL) + ((((uint64_t)op[3] * 4867231046790035161UL) + ((uint64_t)op[4] * 9306810177458241629UL) + ((uint64_t)op[5] * 4675441231982960697UL) + ((uint64_t)op[6] * 4763906738431177706UL) + ((uint64_t)op[7] * 683745108849402065UL) + ((uint64_t)op[8] * 11502197098927065429UL) + ((uint64_t)op[9] * 1502077965055206585UL) + ((uint64_t)op[10] * 12293104920482715914UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 12293104920482715914UL) + ((uint64_t)op[1] * 12824144572736731398UL) + ((uint64_t)op[2] * 17068822549925918684UL) + ((uint64_t)op[3] * 5658328008036602655UL) + ((((uint64_t)op[4] * 4867231046790035161UL) + ((uint64_t)op[5] * 9306810177458241629UL) + ((uint64_t)op[6] * 4675441231982960697UL) + ((uint64_t)op[7] * 4763906738431177706UL) + ((uint64_t)op[8] * 683745108849402065UL) + ((uint64_t)op[9] * 11502197098927065429UL) + ((uint64_t)op[10] * 1502077965055206585UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 1502077965055206585UL) + ((uint64_t)op[1] * 12293104920482715914UL) + ((uint64_t)op[2] * 12824144572736731398UL) + ((uint64_t)op[3] * 17068822549925918684UL) + ((uint64_t)op[4] * 5658328008036602655UL) + ((((uint64_t)op[5] * 4867231046790035161UL) + ((uint64_t)op[6] * 9306810177458241629UL) + ((uint64_t)op[7] * 4675441231982960697UL) + ((uint64_t)op[8] * 4763906738431177706UL) + ((uint64_t)op[9] * 683745108849402065UL) + ((uint64_t)op[10] * 11502197098927065429UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 11502197098927065429UL) + ((uint64_t)op[1] * 1502077965055206585UL) + ((uint64_t)op[2] * 12293104920482715914UL) + ((uint64_t)op[3] * 12824144572736731398UL) + ((uint64_t)op[4] * 17068822549925918684UL) + ((uint64_t)op[5] * 5658328008036602655UL) + ((((uint64_t)op[6] * 4867231046790035161UL) + ((uint64_t)op[7] * 9306810177458241629UL) + ((uint64_t)op[8] * 4675441231982960697UL) + ((uint64_t)op[9] * 4763906738431177706UL) + ((uint64_t)op[10] * 683745108849402065UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 683745108849402065UL) + ((uint64_t)op[1] * 11502197098927065429UL) + ((uint64_t)op[2] * 1502077965055206585UL) + ((uint64_t)op[3] * 12293104920482715914UL) + ((uint64_t)op[4] * 12824144572736731398UL) + ((uint64_t)op[5] * 17068822549925918684UL) + ((uint64_t)op[6] * 5658328008036602655UL) + ((((uint64_t)op[7] * 4867231046790035161UL) + ((uint64_t)op[8] * 9306810177458241629UL) + ((uint64_t)op[9] * 4675441231982960697UL) + ((uint64_t)op[10] * 4763906738431177706UL)) * 18446744073709551610);
	tmp_q[7] = ((uint64_t)op[0] * 4763906738431177706UL) + ((uint64_t)op[1] * 683745108849402065UL) + ((uint64_t)op[2] * 11502197098927065429UL) + ((uint64_t)op[3] * 1502077965055206585UL) + ((uint64_t)op[4] * 12293104920482715914UL) + ((uint64_t)op[5] * 12824144572736731398UL) + ((uint64_t)op[6] * 17068822549925918684UL) + ((uint64_t)op[7] * 5658328008036602655UL) + ((((uint64_t)op[8] * 4867231046790035161UL) + ((uint64_t)op[9] * 9306810177458241629UL) + ((uint64_t)op[10] * 4675441231982960697UL)) * 18446744073709551610);
	tmp_q[8] = ((uint64_t)op[0] * 4675441231982960697UL) + ((uint64_t)op[1] * 4763906738431177706UL) + ((uint64_t)op[2] * 683745108849402065UL) + ((uint64_t)op[3] * 11502197098927065429UL) + ((uint64_t)op[4] * 1502077965055206585UL) + ((uint64_t)op[5] * 12293104920482715914UL) + ((uint64_t)op[6] * 12824144572736731398UL) + ((uint64_t)op[7] * 17068822549925918684UL) + ((uint64_t)op[8] * 5658328008036602655UL) + ((((uint64_t)op[9] * 4867231046790035161UL) + ((uint64_t)op[10] * 9306810177458241629UL)) * 18446744073709551610);
	tmp_q[9] = ((uint64_t)op[0] * 9306810177458241629UL) + ((uint64_t)op[1] * 4675441231982960697UL) + ((uint64_t)op[2] * 4763906738431177706UL) + ((uint64_t)op[3] * 683745108849402065UL) + ((uint64_t)op[4] * 11502197098927065429UL) + ((uint64_t)op[5] * 1502077965055206585UL) + ((uint64_t)op[6] * 12293104920482715914UL) + ((uint64_t)op[7] * 12824144572736731398UL) + ((uint64_t)op[8] * 17068822549925918684UL) + ((uint64_t)op[9] * 5658328008036602655UL) + ((uint64_t)op[10] * 7690101866678892266UL);
	tmp_q[10] = ((uint64_t)op[0] * 4867231046790035161UL) + ((uint64_t)op[1] * 9306810177458241629UL) + ((uint64_t)op[2] * 4675441231982960697UL) + ((uint64_t)op[3] * 4763906738431177706UL) + ((uint64_t)op[4] * 683745108849402065UL) + ((uint64_t)op[5] * 11502197098927065429UL) + ((uint64_t)op[6] * 1502077965055206585UL) + ((uint64_t)op[7] * 12293104920482715914UL) + ((uint64_t)op[8] * 12824144572736731398UL) + ((uint64_t)op[9] * 17068822549925918684UL) + ((uint64_t)op[10] * 5658328008036602655UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 109026022523815L) - ((-((int128)tmp_q[1] * 23903006611900L) - ((int128)tmp_q[2] * 39564097074449L) + ((int128)tmp_q[3] * 69972061077800L) - ((int128)tmp_q[4] * 61071605434878L) + ((int128)tmp_q[5] * 29422606755893L) - ((int128)tmp_q[6] * 120090050419917L) - ((int128)tmp_q[7] * 1273795524733L) + ((int128)tmp_q[8] * 45126736733366L) - ((int128)tmp_q[9] * 49780031885708L) + ((int128)tmp_q[10] * 38951125524220L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 38951125524220L) - ((int128)tmp_q[1] * 109026022523815L) - ((-((int128)tmp_q[2] * 23903006611900L) - ((int128)tmp_q[3] * 39564097074449L) + ((int128)tmp_q[4] * 69972061077800L) - ((int128)tmp_q[5] * 61071605434878L) + ((int128)tmp_q[6] * 29422606755893L) - ((int128)tmp_q[7] * 120090050419917L) - ((int128)tmp_q[8] * 1273795524733L) + ((int128)tmp_q[9] * 45126736733366L) - ((int128)tmp_q[10] * 49780031885708L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 49780031885708L) + ((int128)tmp_q[1] * 38951125524220L) - ((int128)tmp_q[2] * 109026022523815L) - ((-((int128)tmp_q[3] * 23903006611900L) - ((int128)tmp_q[4] * 39564097074449L) + ((int128)tmp_q[5] * 69972061077800L) - ((int128)tmp_q[6] * 61071605434878L) + ((int128)tmp_q[7] * 29422606755893L) - ((int128)tmp_q[8] * 120090050419917L) - ((int128)tmp_q[9] * 1273795524733L) + ((int128)tmp_q[10] * 45126736733366L)) * 6);
	tmp_zero[3] = ((int128)tmp_q[0] * 45126736733366L) - ((int128)tmp_q[1] * 49780031885708L) + ((int128)tmp_q[2] * 38951125524220L) - ((int128)tmp_q[3] * 109026022523815L) - ((-((int128)tmp_q[4] * 23903006611900L) - ((int128)tmp_q[5] * 39564097074449L) + ((int128)tmp_q[6] * 69972061077800L) - ((int128)tmp_q[7] * 61071605434878L) + ((int128)tmp_q[8] * 29422606755893L) - ((int128)tmp_q[9] * 120090050419917L) - ((int128)tmp_q[10] * 1273795524733L)) * 6);
	tmp_zero[4] = -((int128)tmp_q[0] * 1273795524733L) + ((int128)tmp_q[1] * 45126736733366L) - ((int128)tmp_q[2] * 49780031885708L) + ((int128)tmp_q[3] * 38951125524220L) - ((int128)tmp_q[4] * 109026022523815L) - ((-((int128)tmp_q[5] * 23903006611900L) - ((int128)tmp_q[6] * 39564097074449L) + ((int128)tmp_q[7] * 69972061077800L) - ((int128)tmp_q[8] * 61071605434878L) + ((int128)tmp_q[9] * 29422606755893L) - ((int128)tmp_q[10] * 120090050419917L)) * 6);
	tmp_zero[5] = -((int128)tmp_q[0] * 120090050419917L) - ((int128)tmp_q[1] * 1273795524733L) + ((int128)tmp_q[2] * 45126736733366L) - ((int128)tmp_q[3] * 49780031885708L) + ((int128)tmp_q[4] * 38951125524220L) - ((int128)tmp_q[5] * 109026022523815L) - ((-((int128)tmp_q[6] * 23903006611900L) - ((int128)tmp_q[7] * 39564097074449L) + ((int128)tmp_q[8] * 69972061077800L) - ((int128)tmp_q[9] * 61071605434878L) + ((int128)tmp_q[10] * 29422606755893L)) * 6);
	tmp_zero[6] = ((int128)tmp_q[0] * 29422606755893L) - ((int128)tmp_q[1] * 120090050419917L) - ((int128)tmp_q[2] * 1273795524733L) + ((int128)tmp_q[3] * 45126736733366L) - ((int128)tmp_q[4] * 49780031885708L) + ((int128)tmp_q[5] * 38951125524220L) - ((int128)tmp_q[6] * 109026022523815L) - ((-((int128)tmp_q[7] * 23903006611900L) - ((int128)tmp_q[8] * 39564097074449L) + ((int128)tmp_q[9] * 69972061077800L) - ((int128)tmp_q[10] * 61071605434878L)) * 6);
	tmp_zero[7] = -((int128)tmp_q[0] * 61071605434878L) + ((int128)tmp_q[1] * 29422606755893L) - ((int128)tmp_q[2] * 120090050419917L) - ((int128)tmp_q[3] * 1273795524733L) + ((int128)tmp_q[4] * 45126736733366L) - ((int128)tmp_q[5] * 49780031885708L) + ((int128)tmp_q[6] * 38951125524220L) - ((int128)tmp_q[7] * 109026022523815L) - ((-((int128)tmp_q[8] * 23903006611900L) - ((int128)tmp_q[9] * 39564097074449L) + ((int128)tmp_q[10] * 69972061077800L)) * 6);
	tmp_zero[8] = ((int128)tmp_q[0] * 69972061077800L) - ((int128)tmp_q[1] * 61071605434878L) + ((int128)tmp_q[2] * 29422606755893L) - ((int128)tmp_q[3] * 120090050419917L) - ((int128)tmp_q[4] * 1273795524733L) + ((int128)tmp_q[5] * 45126736733366L) - ((int128)tmp_q[6] * 49780031885708L) + ((int128)tmp_q[7] * 38951125524220L) - ((int128)tmp_q[8] * 109026022523815L) - ((-((int128)tmp_q[9] * 23903006611900L) - ((int128)tmp_q[10] * 39564097074449L)) * 6);
	tmp_zero[9] = -((int128)tmp_q[0] * 39564097074449L) + ((int128)tmp_q[1] * 69972061077800L) - ((int128)tmp_q[2] * 61071605434878L) + ((int128)tmp_q[3] * 29422606755893L) - ((int128)tmp_q[4] * 120090050419917L) - ((int128)tmp_q[5] * 1273795524733L) + ((int128)tmp_q[6] * 45126736733366L) - ((int128)tmp_q[7] * 49780031885708L) + ((int128)tmp_q[8] * 38951125524220L) - ((int128)tmp_q[9] * 109026022523815L) + ((int128)tmp_q[10] * 143418039671400L);
	tmp_zero[10] = -((int128)tmp_q[0] * 23903006611900L) - ((int128)tmp_q[1] * 39564097074449L) + ((int128)tmp_q[2] * 69972061077800L) - ((int128)tmp_q[3] * 61071605434878L) + ((int128)tmp_q[4] * 29422606755893L) - ((int128)tmp_q[5] * 120090050419917L) - ((int128)tmp_q[6] * 1273795524733L) + ((int128)tmp_q[7] * 45126736733366L) - ((int128)tmp_q[8] * 49780031885708L) + ((int128)tmp_q[9] * 38951125524220L) - ((int128)tmp_q[10] * 109026022523815L);

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

