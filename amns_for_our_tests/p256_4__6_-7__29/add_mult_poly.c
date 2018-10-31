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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6486183303198996051UL) + ((((uint64_t)op[1] * 12637401507274043734UL) + ((uint64_t)op[2] * 7503288225922872534UL) + ((uint64_t)op[3] * 16470230636052561753UL) + ((uint64_t)op[4] * 13042143521881234050UL) + ((uint64_t)op[5] * 11690515774286714693UL)) * 18446744073709551609);
	tmp_q[1] = ((uint64_t)op[0] * 11690515774286714693UL) + ((uint64_t)op[1] * 6486183303198996051UL) + ((((uint64_t)op[2] * 12637401507274043734UL) + ((uint64_t)op[3] * 7503288225922872534UL) + ((uint64_t)op[4] * 16470230636052561753UL) + ((uint64_t)op[5] * 13042143521881234050UL)) * 18446744073709551609);
	tmp_q[2] = ((uint64_t)op[0] * 13042143521881234050UL) + ((uint64_t)op[1] * 11690515774286714693UL) + ((uint64_t)op[2] * 6486183303198996051UL) + ((((uint64_t)op[3] * 12637401507274043734UL) + ((uint64_t)op[4] * 7503288225922872534UL) + ((uint64_t)op[5] * 16470230636052561753UL)) * 18446744073709551609);
	tmp_q[3] = ((uint64_t)op[0] * 16470230636052561753UL) + ((uint64_t)op[1] * 13042143521881234050UL) + ((uint64_t)op[2] * 11690515774286714693UL) + ((uint64_t)op[3] * 6486183303198996051UL) + ((((uint64_t)op[4] * 12637401507274043734UL) + ((uint64_t)op[5] * 7503288225922872534UL)) * 18446744073709551609);
	tmp_q[4] = ((uint64_t)op[0] * 7503288225922872534UL) + ((uint64_t)op[1] * 16470230636052561753UL) + ((uint64_t)op[2] * 13042143521881234050UL) + ((uint64_t)op[3] * 11690515774286714693UL) + ((uint64_t)op[4] * 6486183303198996051UL) + ((uint64_t)op[5] * 3771909817629451942UL);
	tmp_q[5] = ((uint64_t)op[0] * 12637401507274043734UL) + ((uint64_t)op[1] * 7503288225922872534UL) + ((uint64_t)op[2] * 16470230636052561753UL) + ((uint64_t)op[3] * 13042143521881234050UL) + ((uint64_t)op[4] * 11690515774286714693UL) + ((uint64_t)op[5] * 6486183303198996051UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2655785884728L) - ((((int128)tmp_q[1] * 1545157667331L) + ((int128)tmp_q[2] * 1195258760035L) - ((int128)tmp_q[3] * 665497536280L) - ((int128)tmp_q[4] * 912115140910L) + ((int128)tmp_q[5] * 2170207079745L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 2170207079745L) - ((int128)tmp_q[1] * 2655785884728L) - ((((int128)tmp_q[2] * 1545157667331L) + ((int128)tmp_q[3] * 1195258760035L) - ((int128)tmp_q[4] * 665497536280L) - ((int128)tmp_q[5] * 912115140910L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 912115140910L) + ((int128)tmp_q[1] * 2170207079745L) - ((int128)tmp_q[2] * 2655785884728L) - ((((int128)tmp_q[3] * 1545157667331L) + ((int128)tmp_q[4] * 1195258760035L) - ((int128)tmp_q[5] * 665497536280L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 665497536280L) - ((int128)tmp_q[1] * 912115140910L) + ((int128)tmp_q[2] * 2170207079745L) - ((int128)tmp_q[3] * 2655785884728L) - ((((int128)tmp_q[4] * 1545157667331L) + ((int128)tmp_q[5] * 1195258760035L)) * 7);
	tmp_zero[4] = ((int128)tmp_q[0] * 1195258760035L) - ((int128)tmp_q[1] * 665497536280L) - ((int128)tmp_q[2] * 912115140910L) + ((int128)tmp_q[3] * 2170207079745L) - ((int128)tmp_q[4] * 2655785884728L) - ((int128)tmp_q[5] * 10816103671317L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1545157667331L) + ((int128)tmp_q[1] * 1195258760035L) - ((int128)tmp_q[2] * 665497536280L) - ((int128)tmp_q[3] * 912115140910L) + ((int128)tmp_q[4] * 2170207079745L) - ((int128)tmp_q[5] * 2655785884728L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

