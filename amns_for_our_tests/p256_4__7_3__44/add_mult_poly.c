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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2404199250268664676UL) + ((((uint64_t)op[1] * 11700397478150561449UL) + ((uint64_t)op[2] * 16074104109414192404UL) + ((uint64_t)op[3] * 1140464794597479168UL) + ((uint64_t)op[4] * 8214489931626463245UL) + ((uint64_t)op[5] * 10168178384642489423UL) + ((uint64_t)op[6] * 9547228024202927656UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 9547228024202927656UL) + ((uint64_t)op[1] * 2404199250268664676UL) + ((((uint64_t)op[2] * 11700397478150561449UL) + ((uint64_t)op[3] * 16074104109414192404UL) + ((uint64_t)op[4] * 1140464794597479168UL) + ((uint64_t)op[5] * 8214489931626463245UL) + ((uint64_t)op[6] * 10168178384642489423UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 10168178384642489423UL) + ((uint64_t)op[1] * 9547228024202927656UL) + ((uint64_t)op[2] * 2404199250268664676UL) + ((((uint64_t)op[3] * 11700397478150561449UL) + ((uint64_t)op[4] * 16074104109414192404UL) + ((uint64_t)op[5] * 1140464794597479168UL) + ((uint64_t)op[6] * 8214489931626463245UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 8214489931626463245UL) + ((uint64_t)op[1] * 10168178384642489423UL) + ((uint64_t)op[2] * 9547228024202927656UL) + ((uint64_t)op[3] * 2404199250268664676UL) + ((((uint64_t)op[4] * 11700397478150561449UL) + ((uint64_t)op[5] * 16074104109414192404UL) + ((uint64_t)op[6] * 1140464794597479168UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 1140464794597479168UL) + ((uint64_t)op[1] * 8214489931626463245UL) + ((uint64_t)op[2] * 10168178384642489423UL) + ((uint64_t)op[3] * 9547228024202927656UL) + ((uint64_t)op[4] * 2404199250268664676UL) + ((((uint64_t)op[5] * 11700397478150561449UL) + ((uint64_t)op[6] * 16074104109414192404UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 16074104109414192404UL) + ((uint64_t)op[1] * 1140464794597479168UL) + ((uint64_t)op[2] * 8214489931626463245UL) + ((uint64_t)op[3] * 10168178384642489423UL) + ((uint64_t)op[4] * 9547228024202927656UL) + ((uint64_t)op[5] * 2404199250268664676UL) + ((uint64_t)op[6] * 16654448360742132731UL);
	tmp_q[6] = ((uint64_t)op[0] * 11700397478150561449UL) + ((uint64_t)op[1] * 16074104109414192404UL) + ((uint64_t)op[2] * 1140464794597479168UL) + ((uint64_t)op[3] * 8214489931626463245UL) + ((uint64_t)op[4] * 10168178384642489423UL) + ((uint64_t)op[5] * 9547228024202927656UL) + ((uint64_t)op[6] * 2404199250268664676UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 14048195004L) + ((((int128)tmp_q[1] * 32262515511L) + ((int128)tmp_q[2] * 18363607335L) - ((int128)tmp_q[3] * 19374039587L) + ((int128)tmp_q[4] * 47690867365L) + ((int128)tmp_q[5] * 4394786948L) + ((int128)tmp_q[6] * 67640504675L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 67640504675L) + ((int128)tmp_q[1] * 14048195004L) + ((((int128)tmp_q[2] * 32262515511L) + ((int128)tmp_q[3] * 18363607335L) - ((int128)tmp_q[4] * 19374039587L) + ((int128)tmp_q[5] * 47690867365L) + ((int128)tmp_q[6] * 4394786948L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 4394786948L) + ((int128)tmp_q[1] * 67640504675L) + ((int128)tmp_q[2] * 14048195004L) + ((((int128)tmp_q[3] * 32262515511L) + ((int128)tmp_q[4] * 18363607335L) - ((int128)tmp_q[5] * 19374039587L) + ((int128)tmp_q[6] * 47690867365L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 47690867365L) + ((int128)tmp_q[1] * 4394786948L) + ((int128)tmp_q[2] * 67640504675L) + ((int128)tmp_q[3] * 14048195004L) + ((((int128)tmp_q[4] * 32262515511L) + ((int128)tmp_q[5] * 18363607335L) - ((int128)tmp_q[6] * 19374039587L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 19374039587L) + ((int128)tmp_q[1] * 47690867365L) + ((int128)tmp_q[2] * 4394786948L) + ((int128)tmp_q[3] * 67640504675L) + ((int128)tmp_q[4] * 14048195004L) + ((((int128)tmp_q[5] * 32262515511L) + ((int128)tmp_q[6] * 18363607335L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 18363607335L) - ((int128)tmp_q[1] * 19374039587L) + ((int128)tmp_q[2] * 47690867365L) + ((int128)tmp_q[3] * 4394786948L) + ((int128)tmp_q[4] * 67640504675L) + ((int128)tmp_q[5] * 14048195004L) + ((int128)tmp_q[6] * 96787546533L);
	tmp_zero[6] = ((int128)tmp_q[0] * 32262515511L) + ((int128)tmp_q[1] * 18363607335L) - ((int128)tmp_q[2] * 19374039587L) + ((int128)tmp_q[3] * 47690867365L) + ((int128)tmp_q[4] * 4394786948L) + ((int128)tmp_q[5] * 67640504675L) + ((int128)tmp_q[6] * 14048195004L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

