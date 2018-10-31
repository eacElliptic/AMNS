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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 2096064140284147677UL) + ((((uint64_t)op[1] * 12913000543279940263UL) + ((uint64_t)op[2] * 4663226907849868756UL) + ((uint64_t)op[3] * 6081226914545767641UL) + ((uint64_t)op[4] * 192360136093759907UL) + ((uint64_t)op[5] * 764922671307422184UL) + ((uint64_t)op[6] * 9613073742325900153UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 9613073742325900153UL) + ((uint64_t)op[1] * 2096064140284147677UL) + ((((uint64_t)op[2] * 12913000543279940263UL) + ((uint64_t)op[3] * 4663226907849868756UL) + ((uint64_t)op[4] * 6081226914545767641UL) + ((uint64_t)op[5] * 192360136093759907UL) + ((uint64_t)op[6] * 764922671307422184UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 764922671307422184UL) + ((uint64_t)op[1] * 9613073742325900153UL) + ((uint64_t)op[2] * 2096064140284147677UL) + ((((uint64_t)op[3] * 12913000543279940263UL) + ((uint64_t)op[4] * 4663226907849868756UL) + ((uint64_t)op[5] * 6081226914545767641UL) + ((uint64_t)op[6] * 192360136093759907UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 192360136093759907UL) + ((uint64_t)op[1] * 764922671307422184UL) + ((uint64_t)op[2] * 9613073742325900153UL) + ((uint64_t)op[3] * 2096064140284147677UL) + ((((uint64_t)op[4] * 12913000543279940263UL) + ((uint64_t)op[5] * 4663226907849868756UL) + ((uint64_t)op[6] * 6081226914545767641UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 6081226914545767641UL) + ((uint64_t)op[1] * 192360136093759907UL) + ((uint64_t)op[2] * 764922671307422184UL) + ((uint64_t)op[3] * 9613073742325900153UL) + ((uint64_t)op[4] * 2096064140284147677UL) + ((((uint64_t)op[5] * 12913000543279940263UL) + ((uint64_t)op[6] * 4663226907849868756UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 4663226907849868756UL) + ((uint64_t)op[1] * 6081226914545767641UL) + ((uint64_t)op[2] * 192360136093759907UL) + ((uint64_t)op[3] * 764922671307422184UL) + ((uint64_t)op[4] * 9613073742325900153UL) + ((uint64_t)op[5] * 2096064140284147677UL) + ((uint64_t)op[6] * 16604027508121375377UL);
	tmp_q[6] = ((uint64_t)op[0] * 12913000543279940263UL) + ((uint64_t)op[1] * 4663226907849868756UL) + ((uint64_t)op[2] * 6081226914545767641UL) + ((uint64_t)op[3] * 192360136093759907UL) + ((uint64_t)op[4] * 764922671307422184UL) + ((uint64_t)op[5] * 9613073742325900153UL) + ((uint64_t)op[6] * 2096064140284147677UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 51220205991L) + ((((int128)tmp_q[1] * 7157688197L) - ((int128)tmp_q[2] * 24760737370L) - ((int128)tmp_q[3] * 7505695756L) + ((int128)tmp_q[4] * 56540977868L) + ((int128)tmp_q[5] * 51624536138L) + ((int128)tmp_q[6] * 5453096623L)) * 7);
	tmp_zero[1] = ((int128)tmp_q[0] * 5453096623L) - ((int128)tmp_q[1] * 51220205991L) + ((((int128)tmp_q[2] * 7157688197L) - ((int128)tmp_q[3] * 24760737370L) - ((int128)tmp_q[4] * 7505695756L) + ((int128)tmp_q[5] * 56540977868L) + ((int128)tmp_q[6] * 51624536138L)) * 7);
	tmp_zero[2] = ((int128)tmp_q[0] * 51624536138L) + ((int128)tmp_q[1] * 5453096623L) - ((int128)tmp_q[2] * 51220205991L) + ((((int128)tmp_q[3] * 7157688197L) - ((int128)tmp_q[4] * 24760737370L) - ((int128)tmp_q[5] * 7505695756L) + ((int128)tmp_q[6] * 56540977868L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 56540977868L) + ((int128)tmp_q[1] * 51624536138L) + ((int128)tmp_q[2] * 5453096623L) - ((int128)tmp_q[3] * 51220205991L) + ((((int128)tmp_q[4] * 7157688197L) - ((int128)tmp_q[5] * 24760737370L) - ((int128)tmp_q[6] * 7505695756L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 7505695756L) + ((int128)tmp_q[1] * 56540977868L) + ((int128)tmp_q[2] * 51624536138L) + ((int128)tmp_q[3] * 5453096623L) - ((int128)tmp_q[4] * 51220205991L) + ((((int128)tmp_q[5] * 7157688197L) - ((int128)tmp_q[6] * 24760737370L)) * 7);
	tmp_zero[5] = -((int128)tmp_q[0] * 24760737370L) - ((int128)tmp_q[1] * 7505695756L) + ((int128)tmp_q[2] * 56540977868L) + ((int128)tmp_q[3] * 51624536138L) + ((int128)tmp_q[4] * 5453096623L) - ((int128)tmp_q[5] * 51220205991L) + ((int128)tmp_q[6] * 50103817379L);
	tmp_zero[6] = ((int128)tmp_q[0] * 7157688197L) - ((int128)tmp_q[1] * 24760737370L) - ((int128)tmp_q[2] * 7505695756L) + ((int128)tmp_q[3] * 56540977868L) + ((int128)tmp_q[4] * 51624536138L) + ((int128)tmp_q[5] * 5453096623L) - ((int128)tmp_q[6] * 51220205991L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

