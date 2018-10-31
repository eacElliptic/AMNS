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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 7773707919896081497UL) + ((((uint64_t)op[1] * 2742786669726423538UL) + ((uint64_t)op[2] * 6745813010471695254UL) + ((uint64_t)op[3] * 15156954002309502868UL) + ((uint64_t)op[4] * 4634679969328621604UL) + ((uint64_t)op[5] * 13439234607408819884UL) + ((uint64_t)op[6] * 2468789025586649500UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 2468789025586649500UL) + ((uint64_t)op[1] * 7773707919896081497UL) + ((((uint64_t)op[2] * 2742786669726423538UL) + ((uint64_t)op[3] * 6745813010471695254UL) + ((uint64_t)op[4] * 15156954002309502868UL) + ((uint64_t)op[5] * 4634679969328621604UL) + ((uint64_t)op[6] * 13439234607408819884UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 13439234607408819884UL) + ((uint64_t)op[1] * 2468789025586649500UL) + ((uint64_t)op[2] * 7773707919896081497UL) + ((((uint64_t)op[3] * 2742786669726423538UL) + ((uint64_t)op[4] * 6745813010471695254UL) + ((uint64_t)op[5] * 15156954002309502868UL) + ((uint64_t)op[6] * 4634679969328621604UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 4634679969328621604UL) + ((uint64_t)op[1] * 13439234607408819884UL) + ((uint64_t)op[2] * 2468789025586649500UL) + ((uint64_t)op[3] * 7773707919896081497UL) + ((((uint64_t)op[4] * 2742786669726423538UL) + ((uint64_t)op[5] * 6745813010471695254UL) + ((uint64_t)op[6] * 15156954002309502868UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 15156954002309502868UL) + ((uint64_t)op[1] * 4634679969328621604UL) + ((uint64_t)op[2] * 13439234607408819884UL) + ((uint64_t)op[3] * 2468789025586649500UL) + ((uint64_t)op[4] * 7773707919896081497UL) + ((((uint64_t)op[5] * 2742786669726423538UL) + ((uint64_t)op[6] * 6745813010471695254UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 6745813010471695254UL) + ((uint64_t)op[1] * 15156954002309502868UL) + ((uint64_t)op[2] * 4634679969328621604UL) + ((uint64_t)op[3] * 13439234607408819884UL) + ((uint64_t)op[4] * 2468789025586649500UL) + ((uint64_t)op[5] * 7773707919896081497UL) + ((uint64_t)op[6] * 13713933348632117690UL);
	tmp_q[6] = ((uint64_t)op[0] * 2742786669726423538UL) + ((uint64_t)op[1] * 6745813010471695254UL) + ((uint64_t)op[2] * 15156954002309502868UL) + ((uint64_t)op[3] * 4634679969328621604UL) + ((uint64_t)op[4] * 13439234607408819884UL) + ((uint64_t)op[5] * 2468789025586649500UL) + ((uint64_t)op[6] * 7773707919896081497UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 37440974953L) + ((((int128)tmp_q[1] * 30604007874L) + ((int128)tmp_q[2] * 33664586546L) + ((int128)tmp_q[3] * 11727875076L) + ((int128)tmp_q[4] * 17942070648L) - ((int128)tmp_q[5] * 22781143132L) - ((int128)tmp_q[6] * 35794170892L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 35794170892L) - ((int128)tmp_q[1] * 37440974953L) + ((((int128)tmp_q[2] * 30604007874L) + ((int128)tmp_q[3] * 33664586546L) + ((int128)tmp_q[4] * 11727875076L) + ((int128)tmp_q[5] * 17942070648L) - ((int128)tmp_q[6] * 22781143132L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 22781143132L) - ((int128)tmp_q[1] * 35794170892L) - ((int128)tmp_q[2] * 37440974953L) + ((((int128)tmp_q[3] * 30604007874L) + ((int128)tmp_q[4] * 33664586546L) + ((int128)tmp_q[5] * 11727875076L) + ((int128)tmp_q[6] * 17942070648L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 17942070648L) - ((int128)tmp_q[1] * 22781143132L) - ((int128)tmp_q[2] * 35794170892L) - ((int128)tmp_q[3] * 37440974953L) + ((((int128)tmp_q[4] * 30604007874L) + ((int128)tmp_q[5] * 33664586546L) + ((int128)tmp_q[6] * 11727875076L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 11727875076L) + ((int128)tmp_q[1] * 17942070648L) - ((int128)tmp_q[2] * 22781143132L) - ((int128)tmp_q[3] * 35794170892L) - ((int128)tmp_q[4] * 37440974953L) + ((((int128)tmp_q[5] * 30604007874L) + ((int128)tmp_q[6] * 33664586546L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 33664586546L) + ((int128)tmp_q[1] * 11727875076L) + ((int128)tmp_q[2] * 17942070648L) - ((int128)tmp_q[3] * 22781143132L) - ((int128)tmp_q[4] * 35794170892L) - ((int128)tmp_q[5] * 37440974953L) + ((int128)tmp_q[6] * 153020039370L);
	tmp_zero[6] = ((int128)tmp_q[0] * 30604007874L) + ((int128)tmp_q[1] * 33664586546L) + ((int128)tmp_q[2] * 11727875076L) + ((int128)tmp_q[3] * 17942070648L) - ((int128)tmp_q[4] * 22781143132L) - ((int128)tmp_q[5] * 35794170892L) - ((int128)tmp_q[6] * 37440974953L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

