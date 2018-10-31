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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 7);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 14);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 14);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 643565794296212732UL) + ((((uint64_t)op[1] * 14664827049846755143UL) + ((uint64_t)op[2] * 13921845439709352192UL) + ((uint64_t)op[3] * 4992570488823606054UL) + ((uint64_t)op[4] * 5288508215308981662UL) + ((uint64_t)op[5] * 14902779321630543702UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 14902779321630543702UL) + ((uint64_t)op[1] * 643565794296212732UL) + ((((uint64_t)op[2] * 14664827049846755143UL) + ((uint64_t)op[3] * 13921845439709352192UL) + ((uint64_t)op[4] * 4992570488823606054UL) + ((uint64_t)op[5] * 5288508215308981662UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 5288508215308981662UL) + ((uint64_t)op[1] * 14902779321630543702UL) + ((uint64_t)op[2] * 643565794296212732UL) + ((((uint64_t)op[3] * 14664827049846755143UL) + ((uint64_t)op[4] * 13921845439709352192UL) + ((uint64_t)op[5] * 4992570488823606054UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 4992570488823606054UL) + ((uint64_t)op[1] * 5288508215308981662UL) + ((uint64_t)op[2] * 14902779321630543702UL) + ((uint64_t)op[3] * 643565794296212732UL) + ((((uint64_t)op[4] * 14664827049846755143UL) + ((uint64_t)op[5] * 13921845439709352192UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 13921845439709352192UL) + ((uint64_t)op[1] * 4992570488823606054UL) + ((uint64_t)op[2] * 5288508215308981662UL) + ((uint64_t)op[3] * 14902779321630543702UL) + ((uint64_t)op[4] * 643565794296212732UL) + ((uint64_t)op[5] * 10420068980379527921UL);
	tmp_q[5] = ((uint64_t)op[0] * 14664827049846755143UL) + ((uint64_t)op[1] * 13921845439709352192UL) + ((uint64_t)op[2] * 4992570488823606054UL) + ((uint64_t)op[3] * 5288508215308981662UL) + ((uint64_t)op[4] * 14902779321630543702UL) + ((uint64_t)op[5] * 643565794296212732UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 5838000745008L) + ((((int128)tmp_q[1] * 853138512018L) - ((int128)tmp_q[2] * 736580815914L) + ((int128)tmp_q[3] * 550272985290L) - ((int128)tmp_q[4] * 4032830075028L) + ((int128)tmp_q[5] * 788487777067L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 788487777067L) - ((int128)tmp_q[1] * 5838000745008L) + ((((int128)tmp_q[2] * 853138512018L) - ((int128)tmp_q[3] * 736580815914L) + ((int128)tmp_q[4] * 550272985290L) - ((int128)tmp_q[5] * 4032830075028L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 4032830075028L) + ((int128)tmp_q[1] * 788487777067L) - ((int128)tmp_q[2] * 5838000745008L) + ((((int128)tmp_q[3] * 853138512018L) - ((int128)tmp_q[4] * 736580815914L) + ((int128)tmp_q[5] * 550272985290L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 550272985290L) - ((int128)tmp_q[1] * 4032830075028L) + ((int128)tmp_q[2] * 788487777067L) - ((int128)tmp_q[3] * 5838000745008L) + ((((int128)tmp_q[4] * 853138512018L) - ((int128)tmp_q[5] * 736580815914L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 736580815914L) + ((int128)tmp_q[1] * 550272985290L) - ((int128)tmp_q[2] * 4032830075028L) + ((int128)tmp_q[3] * 788487777067L) - ((int128)tmp_q[4] * 5838000745008L) + ((int128)tmp_q[5] * 5971969584126L);
	tmp_zero[5] = ((int128)tmp_q[0] * 853138512018L) - ((int128)tmp_q[1] * 736580815914L) + ((int128)tmp_q[2] * 550272985290L) - ((int128)tmp_q[3] * 4032830075028L) + ((int128)tmp_q[4] * 788487777067L) - ((int128)tmp_q[5] * 5838000745008L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

