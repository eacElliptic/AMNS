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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 2);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 2);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 2);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 2);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 2);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 2);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17593516296464130625UL) + ((((uint64_t)op[1] * 18146359713994388565UL) + ((uint64_t)op[2] * 15421109477468608762UL) + ((uint64_t)op[3] * 13968831632164098843UL) + ((uint64_t)op[4] * 14743637921695531250UL) + ((uint64_t)op[5] * 8460680700205068700UL) + ((uint64_t)op[6] * 8773400919491583440UL) + ((uint64_t)op[7] * 12944079383504853579UL) + ((uint64_t)op[8] * 9123448522558398123UL) + ((uint64_t)op[9] * 15947939803200273476UL) + ((uint64_t)op[10] * 9902392609903306416UL)) * 4);
	tmp_q[1] = ((uint64_t)op[0] * 9902392609903306416UL) + ((uint64_t)op[1] * 17593516296464130625UL) + ((((uint64_t)op[2] * 18146359713994388565UL) + ((uint64_t)op[3] * 15421109477468608762UL) + ((uint64_t)op[4] * 13968831632164098843UL) + ((uint64_t)op[5] * 14743637921695531250UL) + ((uint64_t)op[6] * 8460680700205068700UL) + ((uint64_t)op[7] * 8773400919491583440UL) + ((uint64_t)op[8] * 12944079383504853579UL) + ((uint64_t)op[9] * 9123448522558398123UL) + ((uint64_t)op[10] * 15947939803200273476UL)) * 4);
	tmp_q[2] = ((uint64_t)op[0] * 15947939803200273476UL) + ((uint64_t)op[1] * 9902392609903306416UL) + ((uint64_t)op[2] * 17593516296464130625UL) + ((((uint64_t)op[3] * 18146359713994388565UL) + ((uint64_t)op[4] * 15421109477468608762UL) + ((uint64_t)op[5] * 13968831632164098843UL) + ((uint64_t)op[6] * 14743637921695531250UL) + ((uint64_t)op[7] * 8460680700205068700UL) + ((uint64_t)op[8] * 8773400919491583440UL) + ((uint64_t)op[9] * 12944079383504853579UL) + ((uint64_t)op[10] * 9123448522558398123UL)) * 4);
	tmp_q[3] = ((uint64_t)op[0] * 9123448522558398123UL) + ((uint64_t)op[1] * 15947939803200273476UL) + ((uint64_t)op[2] * 9902392609903306416UL) + ((uint64_t)op[3] * 17593516296464130625UL) + ((((uint64_t)op[4] * 18146359713994388565UL) + ((uint64_t)op[5] * 15421109477468608762UL) + ((uint64_t)op[6] * 13968831632164098843UL) + ((uint64_t)op[7] * 14743637921695531250UL) + ((uint64_t)op[8] * 8460680700205068700UL) + ((uint64_t)op[9] * 8773400919491583440UL) + ((uint64_t)op[10] * 12944079383504853579UL)) * 4);
	tmp_q[4] = ((uint64_t)op[0] * 12944079383504853579UL) + ((uint64_t)op[1] * 9123448522558398123UL) + ((uint64_t)op[2] * 15947939803200273476UL) + ((uint64_t)op[3] * 9902392609903306416UL) + ((uint64_t)op[4] * 17593516296464130625UL) + ((((uint64_t)op[5] * 18146359713994388565UL) + ((uint64_t)op[6] * 15421109477468608762UL) + ((uint64_t)op[7] * 13968831632164098843UL) + ((uint64_t)op[8] * 14743637921695531250UL) + ((uint64_t)op[9] * 8460680700205068700UL) + ((uint64_t)op[10] * 8773400919491583440UL)) * 4);
	tmp_q[5] = ((uint64_t)op[0] * 8773400919491583440UL) + ((uint64_t)op[1] * 12944079383504853579UL) + ((uint64_t)op[2] * 9123448522558398123UL) + ((uint64_t)op[3] * 15947939803200273476UL) + ((uint64_t)op[4] * 9902392609903306416UL) + ((uint64_t)op[5] * 17593516296464130625UL) + ((((uint64_t)op[6] * 18146359713994388565UL) + ((uint64_t)op[7] * 15421109477468608762UL) + ((uint64_t)op[8] * 13968831632164098843UL) + ((uint64_t)op[9] * 14743637921695531250UL) + ((uint64_t)op[10] * 8460680700205068700UL)) * 4);
	tmp_q[6] = ((uint64_t)op[0] * 8460680700205068700UL) + ((uint64_t)op[1] * 8773400919491583440UL) + ((uint64_t)op[2] * 12944079383504853579UL) + ((uint64_t)op[3] * 9123448522558398123UL) + ((uint64_t)op[4] * 15947939803200273476UL) + ((uint64_t)op[5] * 9902392609903306416UL) + ((uint64_t)op[6] * 17593516296464130625UL) + ((((uint64_t)op[7] * 18146359713994388565UL) + ((uint64_t)op[8] * 15421109477468608762UL) + ((uint64_t)op[9] * 13968831632164098843UL) + ((uint64_t)op[10] * 14743637921695531250UL)) * 4);
	tmp_q[7] = ((uint64_t)op[0] * 14743637921695531250UL) + ((uint64_t)op[1] * 8460680700205068700UL) + ((uint64_t)op[2] * 8773400919491583440UL) + ((uint64_t)op[3] * 12944079383504853579UL) + ((uint64_t)op[4] * 9123448522558398123UL) + ((uint64_t)op[5] * 15947939803200273476UL) + ((uint64_t)op[6] * 9902392609903306416UL) + ((uint64_t)op[7] * 17593516296464130625UL) + ((((uint64_t)op[8] * 18146359713994388565UL) + ((uint64_t)op[9] * 15421109477468608762UL) + ((uint64_t)op[10] * 13968831632164098843UL)) * 4);
	tmp_q[8] = ((uint64_t)op[0] * 13968831632164098843UL) + ((uint64_t)op[1] * 14743637921695531250UL) + ((uint64_t)op[2] * 8460680700205068700UL) + ((uint64_t)op[3] * 8773400919491583440UL) + ((uint64_t)op[4] * 12944079383504853579UL) + ((uint64_t)op[5] * 9123448522558398123UL) + ((uint64_t)op[6] * 15947939803200273476UL) + ((uint64_t)op[7] * 9902392609903306416UL) + ((uint64_t)op[8] * 17593516296464130625UL) + ((((uint64_t)op[9] * 18146359713994388565UL) + ((uint64_t)op[10] * 15421109477468608762UL)) * 4);
	tmp_q[9] = ((uint64_t)op[0] * 15421109477468608762UL) + ((uint64_t)op[1] * 13968831632164098843UL) + ((uint64_t)op[2] * 14743637921695531250UL) + ((uint64_t)op[3] * 8460680700205068700UL) + ((uint64_t)op[4] * 8773400919491583440UL) + ((uint64_t)op[5] * 12944079383504853579UL) + ((uint64_t)op[6] * 9123448522558398123UL) + ((uint64_t)op[7] * 15947939803200273476UL) + ((uint64_t)op[8] * 9902392609903306416UL) + ((uint64_t)op[9] * 17593516296464130625UL) + ((uint64_t)op[10] * 17245206634848899412UL);
	tmp_q[10] = ((uint64_t)op[0] * 18146359713994388565UL) + ((uint64_t)op[1] * 15421109477468608762UL) + ((uint64_t)op[2] * 13968831632164098843UL) + ((uint64_t)op[3] * 14743637921695531250UL) + ((uint64_t)op[4] * 8460680700205068700UL) + ((uint64_t)op[5] * 8773400919491583440UL) + ((uint64_t)op[6] * 12944079383504853579UL) + ((uint64_t)op[7] * 9123448522558398123UL) + ((uint64_t)op[8] * 15947939803200273476UL) + ((uint64_t)op[9] * 9902392609903306416UL) + ((uint64_t)op[10] * 17593516296464130625UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 89700212801157L) + ((((int128)tmp_q[1] * 77730827842554L) - ((int128)tmp_q[2] * 7203936091799L) - ((int128)tmp_q[3] * 49138734452626L) - ((int128)tmp_q[4] * 18508991835704L) - ((int128)tmp_q[5] * 63662580311885L) - ((int128)tmp_q[6] * 69107638563876L) - ((int128)tmp_q[7] * 51789182099809L) - ((int128)tmp_q[8] * 94578140134673L) - ((int128)tmp_q[9] * 31968933187292L) - ((int128)tmp_q[10] * 63079512773712L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 63079512773712L) - ((int128)tmp_q[1] * 89700212801157L) + ((((int128)tmp_q[2] * 77730827842554L) - ((int128)tmp_q[3] * 7203936091799L) - ((int128)tmp_q[4] * 49138734452626L) - ((int128)tmp_q[5] * 18508991835704L) - ((int128)tmp_q[6] * 63662580311885L) - ((int128)tmp_q[7] * 69107638563876L) - ((int128)tmp_q[8] * 51789182099809L) - ((int128)tmp_q[9] * 94578140134673L) - ((int128)tmp_q[10] * 31968933187292L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 31968933187292L) - ((int128)tmp_q[1] * 63079512773712L) - ((int128)tmp_q[2] * 89700212801157L) + ((((int128)tmp_q[3] * 77730827842554L) - ((int128)tmp_q[4] * 7203936091799L) - ((int128)tmp_q[5] * 49138734452626L) - ((int128)tmp_q[6] * 18508991835704L) - ((int128)tmp_q[7] * 63662580311885L) - ((int128)tmp_q[8] * 69107638563876L) - ((int128)tmp_q[9] * 51789182099809L) - ((int128)tmp_q[10] * 94578140134673L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 94578140134673L) - ((int128)tmp_q[1] * 31968933187292L) - ((int128)tmp_q[2] * 63079512773712L) - ((int128)tmp_q[3] * 89700212801157L) + ((((int128)tmp_q[4] * 77730827842554L) - ((int128)tmp_q[5] * 7203936091799L) - ((int128)tmp_q[6] * 49138734452626L) - ((int128)tmp_q[7] * 18508991835704L) - ((int128)tmp_q[8] * 63662580311885L) - ((int128)tmp_q[9] * 69107638563876L) - ((int128)tmp_q[10] * 51789182099809L)) * 4);
	tmp_zero[4] = -((int128)tmp_q[0] * 51789182099809L) - ((int128)tmp_q[1] * 94578140134673L) - ((int128)tmp_q[2] * 31968933187292L) - ((int128)tmp_q[3] * 63079512773712L) - ((int128)tmp_q[4] * 89700212801157L) + ((((int128)tmp_q[5] * 77730827842554L) - ((int128)tmp_q[6] * 7203936091799L) - ((int128)tmp_q[7] * 49138734452626L) - ((int128)tmp_q[8] * 18508991835704L) - ((int128)tmp_q[9] * 63662580311885L) - ((int128)tmp_q[10] * 69107638563876L)) * 4);
	tmp_zero[5] = -((int128)tmp_q[0] * 69107638563876L) - ((int128)tmp_q[1] * 51789182099809L) - ((int128)tmp_q[2] * 94578140134673L) - ((int128)tmp_q[3] * 31968933187292L) - ((int128)tmp_q[4] * 63079512773712L) - ((int128)tmp_q[5] * 89700212801157L) + ((((int128)tmp_q[6] * 77730827842554L) - ((int128)tmp_q[7] * 7203936091799L) - ((int128)tmp_q[8] * 49138734452626L) - ((int128)tmp_q[9] * 18508991835704L) - ((int128)tmp_q[10] * 63662580311885L)) * 4);
	tmp_zero[6] = -((int128)tmp_q[0] * 63662580311885L) - ((int128)tmp_q[1] * 69107638563876L) - ((int128)tmp_q[2] * 51789182099809L) - ((int128)tmp_q[3] * 94578140134673L) - ((int128)tmp_q[4] * 31968933187292L) - ((int128)tmp_q[5] * 63079512773712L) - ((int128)tmp_q[6] * 89700212801157L) + ((((int128)tmp_q[7] * 77730827842554L) - ((int128)tmp_q[8] * 7203936091799L) - ((int128)tmp_q[9] * 49138734452626L) - ((int128)tmp_q[10] * 18508991835704L)) * 4);
	tmp_zero[7] = -((int128)tmp_q[0] * 18508991835704L) - ((int128)tmp_q[1] * 63662580311885L) - ((int128)tmp_q[2] * 69107638563876L) - ((int128)tmp_q[3] * 51789182099809L) - ((int128)tmp_q[4] * 94578140134673L) - ((int128)tmp_q[5] * 31968933187292L) - ((int128)tmp_q[6] * 63079512773712L) - ((int128)tmp_q[7] * 89700212801157L) + ((((int128)tmp_q[8] * 77730827842554L) - ((int128)tmp_q[9] * 7203936091799L) - ((int128)tmp_q[10] * 49138734452626L)) * 4);
	tmp_zero[8] = -((int128)tmp_q[0] * 49138734452626L) - ((int128)tmp_q[1] * 18508991835704L) - ((int128)tmp_q[2] * 63662580311885L) - ((int128)tmp_q[3] * 69107638563876L) - ((int128)tmp_q[4] * 51789182099809L) - ((int128)tmp_q[5] * 94578140134673L) - ((int128)tmp_q[6] * 31968933187292L) - ((int128)tmp_q[7] * 63079512773712L) - ((int128)tmp_q[8] * 89700212801157L) + ((((int128)tmp_q[9] * 77730827842554L) - ((int128)tmp_q[10] * 7203936091799L)) * 4);
	tmp_zero[9] = -((int128)tmp_q[0] * 7203936091799L) - ((int128)tmp_q[1] * 49138734452626L) - ((int128)tmp_q[2] * 18508991835704L) - ((int128)tmp_q[3] * 63662580311885L) - ((int128)tmp_q[4] * 69107638563876L) - ((int128)tmp_q[5] * 51789182099809L) - ((int128)tmp_q[6] * 94578140134673L) - ((int128)tmp_q[7] * 31968933187292L) - ((int128)tmp_q[8] * 63079512773712L) - ((int128)tmp_q[9] * 89700212801157L) + ((int128)tmp_q[10] * 310923311370216L);
	tmp_zero[10] = ((int128)tmp_q[0] * 77730827842554L) - ((int128)tmp_q[1] * 7203936091799L) - ((int128)tmp_q[2] * 49138734452626L) - ((int128)tmp_q[3] * 18508991835704L) - ((int128)tmp_q[4] * 63662580311885L) - ((int128)tmp_q[5] * 69107638563876L) - ((int128)tmp_q[6] * 51789182099809L) - ((int128)tmp_q[7] * 94578140134673L) - ((int128)tmp_q[8] * 31968933187292L) - ((int128)tmp_q[9] * 63079512773712L) - ((int128)tmp_q[10] * 89700212801157L);

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

