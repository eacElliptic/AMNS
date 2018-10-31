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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9113999272779712352UL) + ((((uint64_t)op[1] * 16837793937263908757UL) + ((uint64_t)op[2] * 11952934458646483339UL) + ((uint64_t)op[3] * 15260958045742789735UL) + ((uint64_t)op[4] * 15244315437944617672UL) + ((uint64_t)op[5] * 14232556373308306632UL) + ((uint64_t)op[6] * 2587268184077596548UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 2587268184077596548UL) + ((uint64_t)op[1] * 9113999272779712352UL) + ((((uint64_t)op[2] * 16837793937263908757UL) + ((uint64_t)op[3] * 11952934458646483339UL) + ((uint64_t)op[4] * 15260958045742789735UL) + ((uint64_t)op[5] * 15244315437944617672UL) + ((uint64_t)op[6] * 14232556373308306632UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 14232556373308306632UL) + ((uint64_t)op[1] * 2587268184077596548UL) + ((uint64_t)op[2] * 9113999272779712352UL) + ((((uint64_t)op[3] * 16837793937263908757UL) + ((uint64_t)op[4] * 11952934458646483339UL) + ((uint64_t)op[5] * 15260958045742789735UL) + ((uint64_t)op[6] * 15244315437944617672UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 15244315437944617672UL) + ((uint64_t)op[1] * 14232556373308306632UL) + ((uint64_t)op[2] * 2587268184077596548UL) + ((uint64_t)op[3] * 9113999272779712352UL) + ((((uint64_t)op[4] * 16837793937263908757UL) + ((uint64_t)op[5] * 11952934458646483339UL) + ((uint64_t)op[6] * 15260958045742789735UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 15260958045742789735UL) + ((uint64_t)op[1] * 15244315437944617672UL) + ((uint64_t)op[2] * 14232556373308306632UL) + ((uint64_t)op[3] * 2587268184077596548UL) + ((uint64_t)op[4] * 9113999272779712352UL) + ((((uint64_t)op[5] * 16837793937263908757UL) + ((uint64_t)op[6] * 11952934458646483339UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 11952934458646483339UL) + ((uint64_t)op[1] * 15260958045742789735UL) + ((uint64_t)op[2] * 15244315437944617672UL) + ((uint64_t)op[3] * 14232556373308306632UL) + ((uint64_t)op[4] * 2587268184077596548UL) + ((uint64_t)op[5] * 9113999272779712352UL) + ((uint64_t)op[6] * 4826850409336928577UL);
	tmp_q[6] = ((uint64_t)op[0] * 16837793937263908757UL) + ((uint64_t)op[1] * 11952934458646483339UL) + ((uint64_t)op[2] * 15260958045742789735UL) + ((uint64_t)op[3] * 15244315437944617672UL) + ((uint64_t)op[4] * 14232556373308306632UL) + ((uint64_t)op[5] * 2587268184077596548UL) + ((uint64_t)op[6] * 9113999272779712352UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 50712183238L) - ((((int128)tmp_q[1] * 13550987351L) - ((int128)tmp_q[2] * 18723632979L) + ((int128)tmp_q[3] * 10652273764L) - ((int128)tmp_q[4] * 44120674887L) + ((int128)tmp_q[5] * 1225298277L) + ((int128)tmp_q[6] * 43330407997L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 43330407997L) - ((int128)tmp_q[1] * 50712183238L) - ((((int128)tmp_q[2] * 13550987351L) - ((int128)tmp_q[3] * 18723632979L) + ((int128)tmp_q[4] * 10652273764L) - ((int128)tmp_q[5] * 44120674887L) + ((int128)tmp_q[6] * 1225298277L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 1225298277L) + ((int128)tmp_q[1] * 43330407997L) - ((int128)tmp_q[2] * 50712183238L) - ((((int128)tmp_q[3] * 13550987351L) - ((int128)tmp_q[4] * 18723632979L) + ((int128)tmp_q[5] * 10652273764L) - ((int128)tmp_q[6] * 44120674887L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 44120674887L) + ((int128)tmp_q[1] * 1225298277L) + ((int128)tmp_q[2] * 43330407997L) - ((int128)tmp_q[3] * 50712183238L) - ((((int128)tmp_q[4] * 13550987351L) - ((int128)tmp_q[5] * 18723632979L) + ((int128)tmp_q[6] * 10652273764L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 10652273764L) - ((int128)tmp_q[1] * 44120674887L) + ((int128)tmp_q[2] * 1225298277L) + ((int128)tmp_q[3] * 43330407997L) - ((int128)tmp_q[4] * 50712183238L) - ((((int128)tmp_q[5] * 13550987351L) - ((int128)tmp_q[6] * 18723632979L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 18723632979L) + ((int128)tmp_q[1] * 10652273764L) - ((int128)tmp_q[2] * 44120674887L) + ((int128)tmp_q[3] * 1225298277L) + ((int128)tmp_q[4] * 43330407997L) - ((int128)tmp_q[5] * 50712183238L) - ((int128)tmp_q[6] * 40652962053L);
	tmp_zero[6] = ((int128)tmp_q[0] * 13550987351L) - ((int128)tmp_q[1] * 18723632979L) + ((int128)tmp_q[2] * 10652273764L) - ((int128)tmp_q[3] * 44120674887L) + ((int128)tmp_q[4] * 1225298277L) + ((int128)tmp_q[5] * 43330407997L) - ((int128)tmp_q[6] * 50712183238L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

