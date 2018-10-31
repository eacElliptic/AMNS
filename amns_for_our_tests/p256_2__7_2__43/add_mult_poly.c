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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17273777949762454223UL) + ((((uint64_t)op[1] * 13060820781749187367UL) + ((uint64_t)op[2] * 15949311739685682825UL) + ((uint64_t)op[3] * 7477028583557955120UL) + ((uint64_t)op[4] * 1206417849174560101UL) + ((uint64_t)op[5] * 11114798450518643656UL) + ((uint64_t)op[6] * 8336784512415969697UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 8336784512415969697UL) + ((uint64_t)op[1] * 17273777949762454223UL) + ((((uint64_t)op[2] * 13060820781749187367UL) + ((uint64_t)op[3] * 15949311739685682825UL) + ((uint64_t)op[4] * 7477028583557955120UL) + ((uint64_t)op[5] * 1206417849174560101UL) + ((uint64_t)op[6] * 11114798450518643656UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 11114798450518643656UL) + ((uint64_t)op[1] * 8336784512415969697UL) + ((uint64_t)op[2] * 17273777949762454223UL) + ((((uint64_t)op[3] * 13060820781749187367UL) + ((uint64_t)op[4] * 15949311739685682825UL) + ((uint64_t)op[5] * 7477028583557955120UL) + ((uint64_t)op[6] * 1206417849174560101UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 1206417849174560101UL) + ((uint64_t)op[1] * 11114798450518643656UL) + ((uint64_t)op[2] * 8336784512415969697UL) + ((uint64_t)op[3] * 17273777949762454223UL) + ((((uint64_t)op[4] * 13060820781749187367UL) + ((uint64_t)op[5] * 15949311739685682825UL) + ((uint64_t)op[6] * 7477028583557955120UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 7477028583557955120UL) + ((uint64_t)op[1] * 1206417849174560101UL) + ((uint64_t)op[2] * 11114798450518643656UL) + ((uint64_t)op[3] * 8336784512415969697UL) + ((uint64_t)op[4] * 17273777949762454223UL) + ((((uint64_t)op[5] * 13060820781749187367UL) + ((uint64_t)op[6] * 15949311739685682825UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 15949311739685682825UL) + ((uint64_t)op[1] * 7477028583557955120UL) + ((uint64_t)op[2] * 1206417849174560101UL) + ((uint64_t)op[3] * 11114798450518643656UL) + ((uint64_t)op[4] * 8336784512415969697UL) + ((uint64_t)op[5] * 17273777949762454223UL) + ((uint64_t)op[6] * 7674897489788823118UL);
	tmp_q[6] = ((uint64_t)op[0] * 13060820781749187367UL) + ((uint64_t)op[1] * 15949311739685682825UL) + ((uint64_t)op[2] * 7477028583557955120UL) + ((uint64_t)op[3] * 1206417849174560101UL) + ((uint64_t)op[4] * 11114798450518643656UL) + ((uint64_t)op[5] * 8336784512415969697UL) + ((uint64_t)op[6] * 17273777949762454223UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 31716489759L) + ((-((int128)tmp_q[1] * 61677798547L) - ((int128)tmp_q[2] * 17101986957L) + ((int128)tmp_q[3] * 9741410091L) + ((int128)tmp_q[4] * 44261496386L) + ((int128)tmp_q[5] * 17613149605L) + ((int128)tmp_q[6] * 48258533097L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 48258533097L) - ((int128)tmp_q[1] * 31716489759L) + ((-((int128)tmp_q[2] * 61677798547L) - ((int128)tmp_q[3] * 17101986957L) + ((int128)tmp_q[4] * 9741410091L) + ((int128)tmp_q[5] * 44261496386L) + ((int128)tmp_q[6] * 17613149605L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 17613149605L) + ((int128)tmp_q[1] * 48258533097L) - ((int128)tmp_q[2] * 31716489759L) + ((-((int128)tmp_q[3] * 61677798547L) - ((int128)tmp_q[4] * 17101986957L) + ((int128)tmp_q[5] * 9741410091L) + ((int128)tmp_q[6] * 44261496386L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 44261496386L) + ((int128)tmp_q[1] * 17613149605L) + ((int128)tmp_q[2] * 48258533097L) - ((int128)tmp_q[3] * 31716489759L) + ((-((int128)tmp_q[4] * 61677798547L) - ((int128)tmp_q[5] * 17101986957L) + ((int128)tmp_q[6] * 9741410091L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 9741410091L) + ((int128)tmp_q[1] * 44261496386L) + ((int128)tmp_q[2] * 17613149605L) + ((int128)tmp_q[3] * 48258533097L) - ((int128)tmp_q[4] * 31716489759L) + ((-((int128)tmp_q[5] * 61677798547L) - ((int128)tmp_q[6] * 17101986957L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 17101986957L) + ((int128)tmp_q[1] * 9741410091L) + ((int128)tmp_q[2] * 44261496386L) + ((int128)tmp_q[3] * 17613149605L) + ((int128)tmp_q[4] * 48258533097L) - ((int128)tmp_q[5] * 31716489759L) - ((int128)tmp_q[6] * 123355597094L);
	tmp_zero[6] = -((int128)tmp_q[0] * 61677798547L) - ((int128)tmp_q[1] * 17101986957L) + ((int128)tmp_q[2] * 9741410091L) + ((int128)tmp_q[3] * 44261496386L) + ((int128)tmp_q[4] * 17613149605L) + ((int128)tmp_q[5] * 48258533097L) - ((int128)tmp_q[6] * 31716489759L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

