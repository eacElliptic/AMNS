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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17550477667521893759UL) + ((((uint64_t)op[1] * 7967654103451169970UL) + ((uint64_t)op[2] * 10914039799526154063UL) + ((uint64_t)op[3] * 18105761530571767533UL) + ((uint64_t)op[4] * 8640085630230484516UL) + ((uint64_t)op[5] * 12039854930332121270UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 12039854930332121270UL) + ((uint64_t)op[1] * 17550477667521893759UL) + ((((uint64_t)op[2] * 7967654103451169970UL) + ((uint64_t)op[3] * 10914039799526154063UL) + ((uint64_t)op[4] * 18105761530571767533UL) + ((uint64_t)op[5] * 8640085630230484516UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 8640085630230484516UL) + ((uint64_t)op[1] * 12039854930332121270UL) + ((uint64_t)op[2] * 17550477667521893759UL) + ((((uint64_t)op[3] * 7967654103451169970UL) + ((uint64_t)op[4] * 10914039799526154063UL) + ((uint64_t)op[5] * 18105761530571767533UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 18105761530571767533UL) + ((uint64_t)op[1] * 8640085630230484516UL) + ((uint64_t)op[2] * 12039854930332121270UL) + ((uint64_t)op[3] * 17550477667521893759UL) + ((((uint64_t)op[4] * 7967654103451169970UL) + ((uint64_t)op[5] * 10914039799526154063UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 10914039799526154063UL) + ((uint64_t)op[1] * 18105761530571767533UL) + ((uint64_t)op[2] * 8640085630230484516UL) + ((uint64_t)op[3] * 12039854930332121270UL) + ((uint64_t)op[4] * 17550477667521893759UL) + ((uint64_t)op[5] * 8401000606480704912UL);
	tmp_q[5] = ((uint64_t)op[0] * 7967654103451169970UL) + ((uint64_t)op[1] * 10914039799526154063UL) + ((uint64_t)op[2] * 18105761530571767533UL) + ((uint64_t)op[3] * 8640085630230484516UL) + ((uint64_t)op[4] * 12039854930332121270UL) + ((uint64_t)op[5] * 17550477667521893759UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 12602208012905L) + ((-((int128)tmp_q[1] * 7345757879278L) + ((int128)tmp_q[2] * 37612602706435L) + ((int128)tmp_q[3] * 49892670819213L) - ((int128)tmp_q[4] * 44010717379024L) - ((int128)tmp_q[5] * 23109895939786L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 23109895939786L) + ((int128)tmp_q[1] * 12602208012905L) + ((-((int128)tmp_q[2] * 7345757879278L) + ((int128)tmp_q[3] * 37612602706435L) + ((int128)tmp_q[4] * 49892670819213L) - ((int128)tmp_q[5] * 44010717379024L)) * 8);
	tmp_zero[2] = -((int128)tmp_q[0] * 44010717379024L) - ((int128)tmp_q[1] * 23109895939786L) + ((int128)tmp_q[2] * 12602208012905L) + ((-((int128)tmp_q[3] * 7345757879278L) + ((int128)tmp_q[4] * 37612602706435L) + ((int128)tmp_q[5] * 49892670819213L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 49892670819213L) - ((int128)tmp_q[1] * 44010717379024L) - ((int128)tmp_q[2] * 23109895939786L) + ((int128)tmp_q[3] * 12602208012905L) + ((-((int128)tmp_q[4] * 7345757879278L) + ((int128)tmp_q[5] * 37612602706435L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 37612602706435L) + ((int128)tmp_q[1] * 49892670819213L) - ((int128)tmp_q[2] * 44010717379024L) - ((int128)tmp_q[3] * 23109895939786L) + ((int128)tmp_q[4] * 12602208012905L) - ((int128)tmp_q[5] * 58766063034224L);
	tmp_zero[5] = -((int128)tmp_q[0] * 7345757879278L) + ((int128)tmp_q[1] * 37612602706435L) + ((int128)tmp_q[2] * 49892670819213L) - ((int128)tmp_q[3] * 44010717379024L) - ((int128)tmp_q[4] * 23109895939786L) + ((int128)tmp_q[5] * 12602208012905L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

