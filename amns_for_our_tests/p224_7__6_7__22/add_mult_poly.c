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
	tmp_q[0] = ((uint64_t)op[0] * 5039662094746711911UL) + ((((uint64_t)op[1] * 7460634658334075804UL) + ((uint64_t)op[2] * 16547613509314173759UL) + ((uint64_t)op[3] * 15993061458864218383UL) + ((uint64_t)op[4] * 5072311544865901722UL) + ((uint64_t)op[5] * 3430335575574834456UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 3430335575574834456UL) + ((uint64_t)op[1] * 5039662094746711911UL) + ((((uint64_t)op[2] * 7460634658334075804UL) + ((uint64_t)op[3] * 16547613509314173759UL) + ((uint64_t)op[4] * 15993061458864218383UL) + ((uint64_t)op[5] * 5072311544865901722UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 5072311544865901722UL) + ((uint64_t)op[1] * 3430335575574834456UL) + ((uint64_t)op[2] * 5039662094746711911UL) + ((((uint64_t)op[3] * 7460634658334075804UL) + ((uint64_t)op[4] * 16547613509314173759UL) + ((uint64_t)op[5] * 15993061458864218383UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 15993061458864218383UL) + ((uint64_t)op[1] * 5072311544865901722UL) + ((uint64_t)op[2] * 3430335575574834456UL) + ((uint64_t)op[3] * 5039662094746711911UL) + ((((uint64_t)op[4] * 7460634658334075804UL) + ((uint64_t)op[5] * 16547613509314173759UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 16547613509314173759UL) + ((uint64_t)op[1] * 15993061458864218383UL) + ((uint64_t)op[2] * 5072311544865901722UL) + ((uint64_t)op[3] * 3430335575574834456UL) + ((uint64_t)op[4] * 5039662094746711911UL) + ((uint64_t)op[5] * 15330954460919427396UL);
	tmp_q[5] = ((uint64_t)op[0] * 7460634658334075804UL) + ((uint64_t)op[1] * 16547613509314173759UL) + ((uint64_t)op[2] * 15993061458864218383UL) + ((uint64_t)op[3] * 5072311544865901722UL) + ((uint64_t)op[4] * 3430335575574834456UL) + ((uint64_t)op[5] * 5039662094746711911UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 112301451722L) + ((-((int128)tmp_q[1] * 53701337570L) - ((int128)tmp_q[2] * 16200516473L) + ((int128)tmp_q[3] * 100134730766L) - ((int128)tmp_q[4] * 8802360783L) - ((int128)tmp_q[5] * 45931514067L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 45931514067L) - ((int128)tmp_q[1] * 112301451722L) + ((-((int128)tmp_q[2] * 53701337570L) - ((int128)tmp_q[3] * 16200516473L) + ((int128)tmp_q[4] * 100134730766L) - ((int128)tmp_q[5] * 8802360783L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 8802360783L) - ((int128)tmp_q[1] * 45931514067L) - ((int128)tmp_q[2] * 112301451722L) + ((-((int128)tmp_q[3] * 53701337570L) - ((int128)tmp_q[4] * 16200516473L) + ((int128)tmp_q[5] * 100134730766L)) * 7);
	tmp_zero[3] = ((int128)tmp_q[0] * 100134730766L) - ((int128)tmp_q[1] * 8802360783L) - ((int128)tmp_q[2] * 45931514067L) - ((int128)tmp_q[3] * 112301451722L) + ((-((int128)tmp_q[4] * 53701337570L) - ((int128)tmp_q[5] * 16200516473L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 16200516473L) + ((int128)tmp_q[1] * 100134730766L) - ((int128)tmp_q[2] * 8802360783L) - ((int128)tmp_q[3] * 45931514067L) - ((int128)tmp_q[4] * 112301451722L) - ((int128)tmp_q[5] * 375909362990L);
	tmp_zero[5] = -((int128)tmp_q[0] * 53701337570L) - ((int128)tmp_q[1] * 16200516473L) + ((int128)tmp_q[2] * 100134730766L) - ((int128)tmp_q[3] * 8802360783L) - ((int128)tmp_q[4] * 45931514067L) - ((int128)tmp_q[5] * 112301451722L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

