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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11965777633039403288UL) + ((((uint64_t)op[1] * 17836010218736789271UL) + ((uint64_t)op[2] * 13754557045491267891UL) + ((uint64_t)op[3] * 8707548857154336511UL) + ((uint64_t)op[4] * 15990346062517146501UL) + ((uint64_t)op[5] * 15101668354132704613UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 15101668354132704613UL) + ((uint64_t)op[1] * 11965777633039403288UL) + ((((uint64_t)op[2] * 17836010218736789271UL) + ((uint64_t)op[3] * 13754557045491267891UL) + ((uint64_t)op[4] * 8707548857154336511UL) + ((uint64_t)op[5] * 15990346062517146501UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 15990346062517146501UL) + ((uint64_t)op[1] * 15101668354132704613UL) + ((uint64_t)op[2] * 11965777633039403288UL) + ((((uint64_t)op[3] * 17836010218736789271UL) + ((uint64_t)op[4] * 13754557045491267891UL) + ((uint64_t)op[5] * 8707548857154336511UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 8707548857154336511UL) + ((uint64_t)op[1] * 15990346062517146501UL) + ((uint64_t)op[2] * 15101668354132704613UL) + ((uint64_t)op[3] * 11965777633039403288UL) + ((((uint64_t)op[4] * 17836010218736789271UL) + ((uint64_t)op[5] * 13754557045491267891UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 13754557045491267891UL) + ((uint64_t)op[1] * 8707548857154336511UL) + ((uint64_t)op[2] * 15990346062517146501UL) + ((uint64_t)op[3] * 15101668354132704613UL) + ((uint64_t)op[4] * 11965777633039403288UL) + ((uint64_t)op[5] * 16614542508791264581UL);
	tmp_q[5] = ((uint64_t)op[0] * 17836010218736789271UL) + ((uint64_t)op[1] * 13754557045491267891UL) + ((uint64_t)op[2] * 8707548857154336511UL) + ((uint64_t)op[3] * 15990346062517146501UL) + ((uint64_t)op[4] * 15101668354132704613UL) + ((uint64_t)op[5] * 11965777633039403288UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1213876468L) + ((((int128)tmp_q[1] * 3272172087L) + ((int128)tmp_q[2] * 953332757L) + ((int128)tmp_q[3] * 211944855L) - ((int128)tmp_q[4] * 2479206513L) - ((int128)tmp_q[5] * 215692219L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 215692219L) + ((int128)tmp_q[1] * 1213876468L) + ((((int128)tmp_q[2] * 3272172087L) + ((int128)tmp_q[3] * 953332757L) + ((int128)tmp_q[4] * 211944855L) - ((int128)tmp_q[5] * 2479206513L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 2479206513L) - ((int128)tmp_q[1] * 215692219L) + ((int128)tmp_q[2] * 1213876468L) + ((((int128)tmp_q[3] * 3272172087L) + ((int128)tmp_q[4] * 953332757L) + ((int128)tmp_q[5] * 211944855L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 211944855L) - ((int128)tmp_q[1] * 2479206513L) - ((int128)tmp_q[2] * 215692219L) + ((int128)tmp_q[3] * 1213876468L) + ((((int128)tmp_q[4] * 3272172087L) + ((int128)tmp_q[5] * 953332757L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 953332757L) + ((int128)tmp_q[1] * 211944855L) - ((int128)tmp_q[2] * 2479206513L) - ((int128)tmp_q[3] * 215692219L) + ((int128)tmp_q[4] * 1213876468L) + ((int128)tmp_q[5] * 9816516261L);
	tmp_zero[5] = ((int128)tmp_q[0] * 3272172087L) + ((int128)tmp_q[1] * 953332757L) + ((int128)tmp_q[2] * 211944855L) - ((int128)tmp_q[3] * 2479206513L) - ((int128)tmp_q[4] * 215692219L) + ((int128)tmp_q[5] * 1213876468L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

