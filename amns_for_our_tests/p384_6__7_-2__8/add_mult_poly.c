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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4021458740018842705UL) + ((((uint64_t)op[1] * 4913103495224804980UL) + ((uint64_t)op[2] * 1390977055019084803UL) + ((uint64_t)op[3] * 381246198788864593UL) + ((uint64_t)op[4] * 16427877941362806753UL) + ((uint64_t)op[5] * 4932841301236407558UL) + ((uint64_t)op[6] * 8261383235025501465UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 8261383235025501465UL) + ((uint64_t)op[1] * 4021458740018842705UL) + ((((uint64_t)op[2] * 4913103495224804980UL) + ((uint64_t)op[3] * 1390977055019084803UL) + ((uint64_t)op[4] * 381246198788864593UL) + ((uint64_t)op[5] * 16427877941362806753UL) + ((uint64_t)op[6] * 4932841301236407558UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 4932841301236407558UL) + ((uint64_t)op[1] * 8261383235025501465UL) + ((uint64_t)op[2] * 4021458740018842705UL) + ((((uint64_t)op[3] * 4913103495224804980UL) + ((uint64_t)op[4] * 1390977055019084803UL) + ((uint64_t)op[5] * 381246198788864593UL) + ((uint64_t)op[6] * 16427877941362806753UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 16427877941362806753UL) + ((uint64_t)op[1] * 4932841301236407558UL) + ((uint64_t)op[2] * 8261383235025501465UL) + ((uint64_t)op[3] * 4021458740018842705UL) + ((((uint64_t)op[4] * 4913103495224804980UL) + ((uint64_t)op[5] * 1390977055019084803UL) + ((uint64_t)op[6] * 381246198788864593UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 381246198788864593UL) + ((uint64_t)op[1] * 16427877941362806753UL) + ((uint64_t)op[2] * 4932841301236407558UL) + ((uint64_t)op[3] * 8261383235025501465UL) + ((uint64_t)op[4] * 4021458740018842705UL) + ((((uint64_t)op[5] * 4913103495224804980UL) + ((uint64_t)op[6] * 1390977055019084803UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 1390977055019084803UL) + ((uint64_t)op[1] * 381246198788864593UL) + ((uint64_t)op[2] * 16427877941362806753UL) + ((uint64_t)op[3] * 4932841301236407558UL) + ((uint64_t)op[4] * 8261383235025501465UL) + ((uint64_t)op[5] * 4021458740018842705UL) + ((uint64_t)op[6] * 8620537083259941656UL);
	tmp_q[6] = ((uint64_t)op[0] * 4913103495224804980UL) + ((uint64_t)op[1] * 1390977055019084803UL) + ((uint64_t)op[2] * 381246198788864593UL) + ((uint64_t)op[3] * 16427877941362806753UL) + ((uint64_t)op[4] * 4932841301236407558UL) + ((uint64_t)op[5] * 8261383235025501465UL) + ((uint64_t)op[6] * 4021458740018842705UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 11098019207229363L) - ((((int128)tmp_q[1] * 14319232059732561L) + ((int128)tmp_q[2] * 18253715618052879L) + ((int128)tmp_q[3] * 17496820430955910L) - ((int128)tmp_q[4] * 1526490923146376L) + ((int128)tmp_q[5] * 3291135611941607L) - ((int128)tmp_q[6] * 13786229621256781L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 13786229621256781L) + ((int128)tmp_q[1] * 11098019207229363L) - ((((int128)tmp_q[2] * 14319232059732561L) + ((int128)tmp_q[3] * 18253715618052879L) + ((int128)tmp_q[4] * 17496820430955910L) - ((int128)tmp_q[5] * 1526490923146376L) + ((int128)tmp_q[6] * 3291135611941607L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 3291135611941607L) - ((int128)tmp_q[1] * 13786229621256781L) + ((int128)tmp_q[2] * 11098019207229363L) - ((((int128)tmp_q[3] * 14319232059732561L) + ((int128)tmp_q[4] * 18253715618052879L) + ((int128)tmp_q[5] * 17496820430955910L) - ((int128)tmp_q[6] * 1526490923146376L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 1526490923146376L) + ((int128)tmp_q[1] * 3291135611941607L) - ((int128)tmp_q[2] * 13786229621256781L) + ((int128)tmp_q[3] * 11098019207229363L) - ((((int128)tmp_q[4] * 14319232059732561L) + ((int128)tmp_q[5] * 18253715618052879L) + ((int128)tmp_q[6] * 17496820430955910L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 17496820430955910L) - ((int128)tmp_q[1] * 1526490923146376L) + ((int128)tmp_q[2] * 3291135611941607L) - ((int128)tmp_q[3] * 13786229621256781L) + ((int128)tmp_q[4] * 11098019207229363L) - ((((int128)tmp_q[5] * 14319232059732561L) + ((int128)tmp_q[6] * 18253715618052879L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 18253715618052879L) + ((int128)tmp_q[1] * 17496820430955910L) - ((int128)tmp_q[2] * 1526490923146376L) + ((int128)tmp_q[3] * 3291135611941607L) - ((int128)tmp_q[4] * 13786229621256781L) + ((int128)tmp_q[5] * 11098019207229363L) - ((int128)tmp_q[6] * 28638464119465122L);
	tmp_zero[6] = ((int128)tmp_q[0] * 14319232059732561L) + ((int128)tmp_q[1] * 18253715618052879L) + ((int128)tmp_q[2] * 17496820430955910L) - ((int128)tmp_q[3] * 1526490923146376L) + ((int128)tmp_q[4] * 3291135611941607L) - ((int128)tmp_q[5] * 13786229621256781L) + ((int128)tmp_q[6] * 11098019207229363L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

