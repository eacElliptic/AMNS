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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1553847072079644732UL) + ((((uint64_t)op[1] * 17061411889543850760UL) + ((uint64_t)op[2] * 16400321845983074204UL) + ((uint64_t)op[3] * 11049370438583318402UL) + ((uint64_t)op[4] * 16937877437919231729UL) + ((uint64_t)op[5] * 15609565805321449832UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 15609565805321449832UL) + ((uint64_t)op[1] * 1553847072079644732UL) + ((((uint64_t)op[2] * 17061411889543850760UL) + ((uint64_t)op[3] * 16400321845983074204UL) + ((uint64_t)op[4] * 11049370438583318402UL) + ((uint64_t)op[5] * 16937877437919231729UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 16937877437919231729UL) + ((uint64_t)op[1] * 15609565805321449832UL) + ((uint64_t)op[2] * 1553847072079644732UL) + ((((uint64_t)op[3] * 17061411889543850760UL) + ((uint64_t)op[4] * 16400321845983074204UL) + ((uint64_t)op[5] * 11049370438583318402UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 11049370438583318402UL) + ((uint64_t)op[1] * 16937877437919231729UL) + ((uint64_t)op[2] * 15609565805321449832UL) + ((uint64_t)op[3] * 1553847072079644732UL) + ((((uint64_t)op[4] * 17061411889543850760UL) + ((uint64_t)op[5] * 16400321845983074204UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 16400321845983074204UL) + ((uint64_t)op[1] * 11049370438583318402UL) + ((uint64_t)op[2] * 16937877437919231729UL) + ((uint64_t)op[3] * 15609565805321449832UL) + ((uint64_t)op[4] * 1553847072079644732UL) + ((uint64_t)op[5] * 11520083152881047336UL);
	tmp_q[5] = ((uint64_t)op[0] * 17061411889543850760UL) + ((uint64_t)op[1] * 16400321845983074204UL) + ((uint64_t)op[2] * 11049370438583318402UL) + ((uint64_t)op[3] * 16937877437919231729UL) + ((uint64_t)op[4] * 15609565805321449832UL) + ((uint64_t)op[5] * 1553847072079644732UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 881612004104L) + ((((int128)tmp_q[1] * 1541592625402L) + ((int128)tmp_q[2] * 1267751108339L) - ((int128)tmp_q[3] * 4838544509256L) - ((int128)tmp_q[4] * 4317651125604L) - ((int128)tmp_q[5] * 2399154689856L)) * 5);
	tmp_zero[1] = -((int128)tmp_q[0] * 2399154689856L) + ((int128)tmp_q[1] * 881612004104L) + ((((int128)tmp_q[2] * 1541592625402L) + ((int128)tmp_q[3] * 1267751108339L) - ((int128)tmp_q[4] * 4838544509256L) - ((int128)tmp_q[5] * 4317651125604L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 4317651125604L) - ((int128)tmp_q[1] * 2399154689856L) + ((int128)tmp_q[2] * 881612004104L) + ((((int128)tmp_q[3] * 1541592625402L) + ((int128)tmp_q[4] * 1267751108339L) - ((int128)tmp_q[5] * 4838544509256L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 4838544509256L) - ((int128)tmp_q[1] * 4317651125604L) - ((int128)tmp_q[2] * 2399154689856L) + ((int128)tmp_q[3] * 881612004104L) + ((((int128)tmp_q[4] * 1541592625402L) + ((int128)tmp_q[5] * 1267751108339L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 1267751108339L) - ((int128)tmp_q[1] * 4838544509256L) - ((int128)tmp_q[2] * 4317651125604L) - ((int128)tmp_q[3] * 2399154689856L) + ((int128)tmp_q[4] * 881612004104L) + ((int128)tmp_q[5] * 7707963127010L);
	tmp_zero[5] = ((int128)tmp_q[0] * 1541592625402L) + ((int128)tmp_q[1] * 1267751108339L) - ((int128)tmp_q[2] * 4838544509256L) - ((int128)tmp_q[3] * 4317651125604L) - ((int128)tmp_q[4] * 2399154689856L) + ((int128)tmp_q[5] * 881612004104L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

