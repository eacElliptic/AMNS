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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13888963903426036915UL) + ((((uint64_t)op[1] * 8868248862779443635UL) + ((uint64_t)op[2] * 18034514773342196496UL) + ((uint64_t)op[3] * 8327602946410918602UL) + ((uint64_t)op[4] * 4442795477752404017UL) + ((uint64_t)op[5] * 11790003922702346988UL) + ((uint64_t)op[6] * 12052033435639267976UL) + ((uint64_t)op[7] * 9119436979492299092UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 9119436979492299092UL) + ((uint64_t)op[1] * 13888963903426036915UL) + ((((uint64_t)op[2] * 8868248862779443635UL) + ((uint64_t)op[3] * 18034514773342196496UL) + ((uint64_t)op[4] * 8327602946410918602UL) + ((uint64_t)op[5] * 4442795477752404017UL) + ((uint64_t)op[6] * 11790003922702346988UL) + ((uint64_t)op[7] * 12052033435639267976UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 12052033435639267976UL) + ((uint64_t)op[1] * 9119436979492299092UL) + ((uint64_t)op[2] * 13888963903426036915UL) + ((((uint64_t)op[3] * 8868248862779443635UL) + ((uint64_t)op[4] * 18034514773342196496UL) + ((uint64_t)op[5] * 8327602946410918602UL) + ((uint64_t)op[6] * 4442795477752404017UL) + ((uint64_t)op[7] * 11790003922702346988UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 11790003922702346988UL) + ((uint64_t)op[1] * 12052033435639267976UL) + ((uint64_t)op[2] * 9119436979492299092UL) + ((uint64_t)op[3] * 13888963903426036915UL) + ((((uint64_t)op[4] * 8868248862779443635UL) + ((uint64_t)op[5] * 18034514773342196496UL) + ((uint64_t)op[6] * 8327602946410918602UL) + ((uint64_t)op[7] * 4442795477752404017UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 4442795477752404017UL) + ((uint64_t)op[1] * 11790003922702346988UL) + ((uint64_t)op[2] * 12052033435639267976UL) + ((uint64_t)op[3] * 9119436979492299092UL) + ((uint64_t)op[4] * 13888963903426036915UL) + ((((uint64_t)op[5] * 8868248862779443635UL) + ((uint64_t)op[6] * 18034514773342196496UL) + ((uint64_t)op[7] * 8327602946410918602UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 8327602946410918602UL) + ((uint64_t)op[1] * 4442795477752404017UL) + ((uint64_t)op[2] * 11790003922702346988UL) + ((uint64_t)op[3] * 12052033435639267976UL) + ((uint64_t)op[4] * 9119436979492299092UL) + ((uint64_t)op[5] * 13888963903426036915UL) + ((((uint64_t)op[6] * 8868248862779443635UL) + ((uint64_t)op[7] * 18034514773342196496UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 18034514773342196496UL) + ((uint64_t)op[1] * 8327602946410918602UL) + ((uint64_t)op[2] * 4442795477752404017UL) + ((uint64_t)op[3] * 11790003922702346988UL) + ((uint64_t)op[4] * 12052033435639267976UL) + ((uint64_t)op[5] * 9119436979492299092UL) + ((uint64_t)op[6] * 13888963903426036915UL) + ((uint64_t)op[7] * 17736497725558887270UL);
	tmp_q[7] = ((uint64_t)op[0] * 8868248862779443635UL) + ((uint64_t)op[1] * 18034514773342196496UL) + ((uint64_t)op[2] * 8327602946410918602UL) + ((uint64_t)op[3] * 4442795477752404017UL) + ((uint64_t)op[4] * 11790003922702346988UL) + ((uint64_t)op[5] * 12052033435639267976UL) + ((uint64_t)op[6] * 9119436979492299092UL) + ((uint64_t)op[7] * 13888963903426036915UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 87532609288475L) + ((((int128)tmp_q[1] * 125025646985761L) + ((int128)tmp_q[2] * 70643106331314L) - ((int128)tmp_q[3] * 109647580087582L) + ((int128)tmp_q[4] * 103547046001959L) - ((int128)tmp_q[5] * 58720097870744L) - ((int128)tmp_q[6] * 139349311673588L) + ((int128)tmp_q[7] * 17304350099652L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 17304350099652L) + ((int128)tmp_q[1] * 87532609288475L) + ((((int128)tmp_q[2] * 125025646985761L) + ((int128)tmp_q[3] * 70643106331314L) - ((int128)tmp_q[4] * 109647580087582L) + ((int128)tmp_q[5] * 103547046001959L) - ((int128)tmp_q[6] * 58720097870744L) - ((int128)tmp_q[7] * 139349311673588L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 139349311673588L) + ((int128)tmp_q[1] * 17304350099652L) + ((int128)tmp_q[2] * 87532609288475L) + ((((int128)tmp_q[3] * 125025646985761L) + ((int128)tmp_q[4] * 70643106331314L) - ((int128)tmp_q[5] * 109647580087582L) + ((int128)tmp_q[6] * 103547046001959L) - ((int128)tmp_q[7] * 58720097870744L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 58720097870744L) - ((int128)tmp_q[1] * 139349311673588L) + ((int128)tmp_q[2] * 17304350099652L) + ((int128)tmp_q[3] * 87532609288475L) + ((((int128)tmp_q[4] * 125025646985761L) + ((int128)tmp_q[5] * 70643106331314L) - ((int128)tmp_q[6] * 109647580087582L) + ((int128)tmp_q[7] * 103547046001959L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 103547046001959L) - ((int128)tmp_q[1] * 58720097870744L) - ((int128)tmp_q[2] * 139349311673588L) + ((int128)tmp_q[3] * 17304350099652L) + ((int128)tmp_q[4] * 87532609288475L) + ((((int128)tmp_q[5] * 125025646985761L) + ((int128)tmp_q[6] * 70643106331314L) - ((int128)tmp_q[7] * 109647580087582L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 109647580087582L) + ((int128)tmp_q[1] * 103547046001959L) - ((int128)tmp_q[2] * 58720097870744L) - ((int128)tmp_q[3] * 139349311673588L) + ((int128)tmp_q[4] * 17304350099652L) + ((int128)tmp_q[5] * 87532609288475L) + ((((int128)tmp_q[6] * 125025646985761L) + ((int128)tmp_q[7] * 70643106331314L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 70643106331314L) - ((int128)tmp_q[1] * 109647580087582L) + ((int128)tmp_q[2] * 103547046001959L) - ((int128)tmp_q[3] * 58720097870744L) - ((int128)tmp_q[4] * 139349311673588L) + ((int128)tmp_q[5] * 17304350099652L) + ((int128)tmp_q[6] * 87532609288475L) + ((int128)tmp_q[7] * 250051293971522L);
	tmp_zero[7] = ((int128)tmp_q[0] * 125025646985761L) + ((int128)tmp_q[1] * 70643106331314L) - ((int128)tmp_q[2] * 109647580087582L) + ((int128)tmp_q[3] * 103547046001959L) - ((int128)tmp_q[4] * 58720097870744L) - ((int128)tmp_q[5] * 139349311673588L) + ((int128)tmp_q[6] * 17304350099652L) + ((int128)tmp_q[7] * 87532609288475L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

