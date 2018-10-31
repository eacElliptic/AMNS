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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[4] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16913659306981172769UL) + ((((uint64_t)op[1] * 6959560249532003418UL) + ((uint64_t)op[2] * 13460051344501404069UL) + ((uint64_t)op[3] * 15347771185691530733UL) + ((uint64_t)op[4] * 4399184555524707231UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 4399184555524707231UL) + ((uint64_t)op[1] * 16913659306981172769UL) + ((((uint64_t)op[2] * 6959560249532003418UL) + ((uint64_t)op[3] * 13460051344501404069UL) + ((uint64_t)op[4] * 15347771185691530733UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 15347771185691530733UL) + ((uint64_t)op[1] * 4399184555524707231UL) + ((uint64_t)op[2] * 16913659306981172769UL) + ((((uint64_t)op[3] * 6959560249532003418UL) + ((uint64_t)op[4] * 13460051344501404069UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 13460051344501404069UL) + ((uint64_t)op[1] * 15347771185691530733UL) + ((uint64_t)op[2] * 4399184555524707231UL) + ((uint64_t)op[3] * 16913659306981172769UL) + ((uint64_t)op[4] * 13919120499064006836UL);
	tmp_q[4] = ((uint64_t)op[0] * 6959560249532003418UL) + ((uint64_t)op[1] * 13460051344501404069UL) + ((uint64_t)op[2] * 15347771185691530733UL) + ((uint64_t)op[3] * 4399184555524707231UL) + ((uint64_t)op[4] * 16913659306981172769UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 44648520865L) + ((((int128)tmp_q[1] * 121412020815L) - ((int128)tmp_q[2] * 99127117784L) - ((int128)tmp_q[3] * 19129189304L) + ((int128)tmp_q[4] * 201560480555L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 201560480555L) + ((int128)tmp_q[1] * 44648520865L) + ((((int128)tmp_q[2] * 121412020815L) - ((int128)tmp_q[3] * 99127117784L) - ((int128)tmp_q[4] * 19129189304L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 19129189304L) + ((int128)tmp_q[1] * 201560480555L) + ((int128)tmp_q[2] * 44648520865L) + ((((int128)tmp_q[3] * 121412020815L) - ((int128)tmp_q[4] * 99127117784L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 99127117784L) - ((int128)tmp_q[1] * 19129189304L) + ((int128)tmp_q[2] * 201560480555L) + ((int128)tmp_q[3] * 44648520865L) + ((int128)tmp_q[4] * 242824041630L);
	tmp_zero[4] = ((int128)tmp_q[0] * 121412020815L) - ((int128)tmp_q[1] * 99127117784L) - ((int128)tmp_q[2] * 19129189304L) + ((int128)tmp_q[3] * 201560480555L) + ((int128)tmp_q[4] * 44648520865L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

