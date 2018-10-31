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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12047730704390324920UL) + ((((uint64_t)op[1] * 13451043304859541295UL) + ((uint64_t)op[2] * 15339297955215758435UL) + ((uint64_t)op[3] * 5902455129883066159UL) + ((uint64_t)op[4] * 15516448467868261596UL) + ((uint64_t)op[5] * 11085501154934285090UL) + ((uint64_t)op[6] * 10737793803226303324UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 10737793803226303324UL) + ((uint64_t)op[1] * 12047730704390324920UL) + ((((uint64_t)op[2] * 13451043304859541295UL) + ((uint64_t)op[3] * 15339297955215758435UL) + ((uint64_t)op[4] * 5902455129883066159UL) + ((uint64_t)op[5] * 15516448467868261596UL) + ((uint64_t)op[6] * 11085501154934285090UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 11085501154934285090UL) + ((uint64_t)op[1] * 10737793803226303324UL) + ((uint64_t)op[2] * 12047730704390324920UL) + ((((uint64_t)op[3] * 13451043304859541295UL) + ((uint64_t)op[4] * 15339297955215758435UL) + ((uint64_t)op[5] * 5902455129883066159UL) + ((uint64_t)op[6] * 15516448467868261596UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 15516448467868261596UL) + ((uint64_t)op[1] * 11085501154934285090UL) + ((uint64_t)op[2] * 10737793803226303324UL) + ((uint64_t)op[3] * 12047730704390324920UL) + ((((uint64_t)op[4] * 13451043304859541295UL) + ((uint64_t)op[5] * 15339297955215758435UL) + ((uint64_t)op[6] * 5902455129883066159UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 5902455129883066159UL) + ((uint64_t)op[1] * 15516448467868261596UL) + ((uint64_t)op[2] * 11085501154934285090UL) + ((uint64_t)op[3] * 10737793803226303324UL) + ((uint64_t)op[4] * 12047730704390324920UL) + ((((uint64_t)op[5] * 13451043304859541295UL) + ((uint64_t)op[6] * 15339297955215758435UL)) * 18446744073709551609);
	tmp_q[5] = ((uint64_t)op[0] * 15339297955215758435UL) + ((uint64_t)op[1] * 5902455129883066159UL) + ((uint64_t)op[2] * 15516448467868261596UL) + ((uint64_t)op[3] * 11085501154934285090UL) + ((uint64_t)op[4] * 10737793803226303324UL) + ((uint64_t)op[5] * 12047730704390324920UL) + ((uint64_t)op[6] * 16523161308240520631UL);
	tmp_q[6] = ((uint64_t)op[0] * 13451043304859541295UL) + ((uint64_t)op[1] * 15339297955215758435UL) + ((uint64_t)op[2] * 5902455129883066159UL) + ((uint64_t)op[3] * 15516448467868261596UL) + ((uint64_t)op[4] * 11085501154934285090UL) + ((uint64_t)op[5] * 10737793803226303324UL) + ((uint64_t)op[6] * 12047730704390324920UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 42331789584L) - ((-((int128)tmp_q[1] * 20867233703L) - ((int128)tmp_q[2] * 44130107199L) - ((int128)tmp_q[3] * 22348983976L) + ((int128)tmp_q[4] * 49202499879L) + ((int128)tmp_q[5] * 86873588723L) - ((int128)tmp_q[6] * 923460903L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 923460903L) + ((int128)tmp_q[1] * 42331789584L) - ((-((int128)tmp_q[2] * 20867233703L) - ((int128)tmp_q[3] * 44130107199L) - ((int128)tmp_q[4] * 22348983976L) + ((int128)tmp_q[5] * 49202499879L) + ((int128)tmp_q[6] * 86873588723L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 86873588723L) - ((int128)tmp_q[1] * 923460903L) + ((int128)tmp_q[2] * 42331789584L) - ((-((int128)tmp_q[3] * 20867233703L) - ((int128)tmp_q[4] * 44130107199L) - ((int128)tmp_q[5] * 22348983976L) + ((int128)tmp_q[6] * 49202499879L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 49202499879L) + ((int128)tmp_q[1] * 86873588723L) - ((int128)tmp_q[2] * 923460903L) + ((int128)tmp_q[3] * 42331789584L) - ((-((int128)tmp_q[4] * 20867233703L) - ((int128)tmp_q[5] * 44130107199L) - ((int128)tmp_q[6] * 22348983976L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 22348983976L) + ((int128)tmp_q[1] * 49202499879L) + ((int128)tmp_q[2] * 86873588723L) - ((int128)tmp_q[3] * 923460903L) + ((int128)tmp_q[4] * 42331789584L) - ((-((int128)tmp_q[5] * 20867233703L) - ((int128)tmp_q[6] * 44130107199L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 44130107199L) - ((int128)tmp_q[1] * 22348983976L) + ((int128)tmp_q[2] * 49202499879L) + ((int128)tmp_q[3] * 86873588723L) - ((int128)tmp_q[4] * 923460903L) + ((int128)tmp_q[5] * 42331789584L) + ((int128)tmp_q[6] * 146070635921L);
	tmp_zero[6] = -((int128)tmp_q[0] * 20867233703L) - ((int128)tmp_q[1] * 44130107199L) - ((int128)tmp_q[2] * 22348983976L) + ((int128)tmp_q[3] * 49202499879L) + ((int128)tmp_q[4] * 86873588723L) - ((int128)tmp_q[5] * 923460903L) + ((int128)tmp_q[6] * 42331789584L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

