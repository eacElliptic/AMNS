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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11414601438657336255UL) + ((((uint64_t)op[1] * 6056880519251800030UL) + ((uint64_t)op[2] * 6709136824034728427UL) + ((uint64_t)op[3] * 13245475494167697889UL) + ((uint64_t)op[4] * 4598283202487634850UL) + ((uint64_t)op[5] * 10896161970100450052UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 10896161970100450052UL) + ((uint64_t)op[1] * 11414601438657336255UL) + ((((uint64_t)op[2] * 6056880519251800030UL) + ((uint64_t)op[3] * 6709136824034728427UL) + ((uint64_t)op[4] * 13245475494167697889UL) + ((uint64_t)op[5] * 4598283202487634850UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 4598283202487634850UL) + ((uint64_t)op[1] * 10896161970100450052UL) + ((uint64_t)op[2] * 11414601438657336255UL) + ((((uint64_t)op[3] * 6056880519251800030UL) + ((uint64_t)op[4] * 6709136824034728427UL) + ((uint64_t)op[5] * 13245475494167697889UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 13245475494167697889UL) + ((uint64_t)op[1] * 4598283202487634850UL) + ((uint64_t)op[2] * 10896161970100450052UL) + ((uint64_t)op[3] * 11414601438657336255UL) + ((((uint64_t)op[4] * 6056880519251800030UL) + ((uint64_t)op[5] * 6709136824034728427UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 6709136824034728427UL) + ((uint64_t)op[1] * 13245475494167697889UL) + ((uint64_t)op[2] * 4598283202487634850UL) + ((uint64_t)op[3] * 10896161970100450052UL) + ((uint64_t)op[4] * 11414601438657336255UL) + ((uint64_t)op[5] * 6609085551160103082UL);
	tmp_q[5] = ((uint64_t)op[0] * 6056880519251800030UL) + ((uint64_t)op[1] * 6709136824034728427UL) + ((uint64_t)op[2] * 13245475494167697889UL) + ((uint64_t)op[3] * 4598283202487634850UL) + ((uint64_t)op[4] * 10896161970100450052UL) + ((uint64_t)op[5] * 11414601438657336255UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 350717552394L) - ((((int128)tmp_q[1] * 3013191814774L) - ((int128)tmp_q[2] * 2441969940785L) - ((int128)tmp_q[3] * 4032010917100L) - ((int128)tmp_q[4] * 2750339934455L) + ((int128)tmp_q[5] * 1417764849471L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 1417764849471L) - ((int128)tmp_q[1] * 350717552394L) - ((((int128)tmp_q[2] * 3013191814774L) - ((int128)tmp_q[3] * 2441969940785L) - ((int128)tmp_q[4] * 4032010917100L) - ((int128)tmp_q[5] * 2750339934455L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 2750339934455L) + ((int128)tmp_q[1] * 1417764849471L) - ((int128)tmp_q[2] * 350717552394L) - ((((int128)tmp_q[3] * 3013191814774L) - ((int128)tmp_q[4] * 2441969940785L) - ((int128)tmp_q[5] * 4032010917100L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 4032010917100L) - ((int128)tmp_q[1] * 2750339934455L) + ((int128)tmp_q[2] * 1417764849471L) - ((int128)tmp_q[3] * 350717552394L) - ((((int128)tmp_q[4] * 3013191814774L) - ((int128)tmp_q[5] * 2441969940785L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 2441969940785L) - ((int128)tmp_q[1] * 4032010917100L) - ((int128)tmp_q[2] * 2750339934455L) + ((int128)tmp_q[3] * 1417764849471L) - ((int128)tmp_q[4] * 350717552394L) - ((int128)tmp_q[5] * 15065959073870L);
	tmp_zero[5] = ((int128)tmp_q[0] * 3013191814774L) - ((int128)tmp_q[1] * 2441969940785L) - ((int128)tmp_q[2] * 4032010917100L) - ((int128)tmp_q[3] * 2750339934455L) + ((int128)tmp_q[4] * 1417764849471L) - ((int128)tmp_q[5] * 350717552394L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

