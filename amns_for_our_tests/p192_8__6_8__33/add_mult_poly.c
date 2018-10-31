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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[5]) << 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 4);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[5] * pa[4]) << 4);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[5] * pa[5]) << 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9899910269527336937UL) + ((((uint64_t)op[1] * 6310895555580096006UL) + ((uint64_t)op[2] * 7565442754308123259UL) + ((uint64_t)op[3] * 7212585998473973119UL) + ((uint64_t)op[4] * 3871058431287207574UL) + ((uint64_t)op[5] * 13488846359585776261UL)) * 8);
	tmp_q[1] = ((uint64_t)op[0] * 13488846359585776261UL) + ((uint64_t)op[1] * 9899910269527336937UL) + ((((uint64_t)op[2] * 6310895555580096006UL) + ((uint64_t)op[3] * 7565442754308123259UL) + ((uint64_t)op[4] * 7212585998473973119UL) + ((uint64_t)op[5] * 3871058431287207574UL)) * 8);
	tmp_q[2] = ((uint64_t)op[0] * 3871058431287207574UL) + ((uint64_t)op[1] * 13488846359585776261UL) + ((uint64_t)op[2] * 9899910269527336937UL) + ((((uint64_t)op[3] * 6310895555580096006UL) + ((uint64_t)op[4] * 7565442754308123259UL) + ((uint64_t)op[5] * 7212585998473973119UL)) * 8);
	tmp_q[3] = ((uint64_t)op[0] * 7212585998473973119UL) + ((uint64_t)op[1] * 3871058431287207574UL) + ((uint64_t)op[2] * 13488846359585776261UL) + ((uint64_t)op[3] * 9899910269527336937UL) + ((((uint64_t)op[4] * 6310895555580096006UL) + ((uint64_t)op[5] * 7565442754308123259UL)) * 8);
	tmp_q[4] = ((uint64_t)op[0] * 7565442754308123259UL) + ((uint64_t)op[1] * 7212585998473973119UL) + ((uint64_t)op[2] * 3871058431287207574UL) + ((uint64_t)op[3] * 13488846359585776261UL) + ((uint64_t)op[4] * 9899910269527336937UL) + ((uint64_t)op[5] * 13593676297221664816UL);
	tmp_q[5] = ((uint64_t)op[0] * 6310895555580096006UL) + ((uint64_t)op[1] * 7565442754308123259UL) + ((uint64_t)op[2] * 7212585998473973119UL) + ((uint64_t)op[3] * 3871058431287207574UL) + ((uint64_t)op[4] * 13488846359585776261UL) + ((uint64_t)op[5] * 9899910269527336937UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 40640918125489L) + ((-((int128)tmp_q[1] * 43604604364462L) + ((int128)tmp_q[2] * 90587110346346L) + ((int128)tmp_q[3] * 40521391834112L) + ((int128)tmp_q[4] * 35988281077805L) - ((int128)tmp_q[5] * 9359544825939L)) * 8);
	tmp_zero[1] = -((int128)tmp_q[0] * 9359544825939L) - ((int128)tmp_q[1] * 40640918125489L) + ((-((int128)tmp_q[2] * 43604604364462L) + ((int128)tmp_q[3] * 90587110346346L) + ((int128)tmp_q[4] * 40521391834112L) + ((int128)tmp_q[5] * 35988281077805L)) * 8);
	tmp_zero[2] = ((int128)tmp_q[0] * 35988281077805L) - ((int128)tmp_q[1] * 9359544825939L) - ((int128)tmp_q[2] * 40640918125489L) + ((-((int128)tmp_q[3] * 43604604364462L) + ((int128)tmp_q[4] * 90587110346346L) + ((int128)tmp_q[5] * 40521391834112L)) * 8);
	tmp_zero[3] = ((int128)tmp_q[0] * 40521391834112L) + ((int128)tmp_q[1] * 35988281077805L) - ((int128)tmp_q[2] * 9359544825939L) - ((int128)tmp_q[3] * 40640918125489L) + ((-((int128)tmp_q[4] * 43604604364462L) + ((int128)tmp_q[5] * 90587110346346L)) * 8);
	tmp_zero[4] = ((int128)tmp_q[0] * 90587110346346L) + ((int128)tmp_q[1] * 40521391834112L) + ((int128)tmp_q[2] * 35988281077805L) - ((int128)tmp_q[3] * 9359544825939L) - ((int128)tmp_q[4] * 40640918125489L) - ((int128)tmp_q[5] * 348836834915696L);
	tmp_zero[5] = -((int128)tmp_q[0] * 43604604364462L) + ((int128)tmp_q[1] * 90587110346346L) + ((int128)tmp_q[2] * 40521391834112L) + ((int128)tmp_q[3] * 35988281077805L) - ((int128)tmp_q[4] * 9359544825939L) - ((int128)tmp_q[5] * 40640918125489L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

