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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) << 1);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) << 2);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) << 1);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11789567066766345231UL) + ((((uint64_t)op[1] * 2028871269219441641UL) + ((uint64_t)op[2] * 15338898177291282674UL) + ((uint64_t)op[3] * 12439336914687665895UL) + ((uint64_t)op[4] * 18011779295393838587UL) + ((uint64_t)op[5] * 1500806634566570764UL) + ((uint64_t)op[6] * 8853821436126993525UL) + ((uint64_t)op[7] * 11361451194444780761UL) + ((uint64_t)op[8] * 5102579330155346788UL) + ((uint64_t)op[9] * 8305064368330360518UL) + ((uint64_t)op[10] * 14344461159847157371UL) + ((uint64_t)op[11] * 17578588421717713409UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 17578588421717713409UL) + ((uint64_t)op[1] * 11789567066766345231UL) + ((((uint64_t)op[2] * 2028871269219441641UL) + ((uint64_t)op[3] * 15338898177291282674UL) + ((uint64_t)op[4] * 12439336914687665895UL) + ((uint64_t)op[5] * 18011779295393838587UL) + ((uint64_t)op[6] * 1500806634566570764UL) + ((uint64_t)op[7] * 8853821436126993525UL) + ((uint64_t)op[8] * 11361451194444780761UL) + ((uint64_t)op[9] * 5102579330155346788UL) + ((uint64_t)op[10] * 8305064368330360518UL) + ((uint64_t)op[11] * 14344461159847157371UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 14344461159847157371UL) + ((uint64_t)op[1] * 17578588421717713409UL) + ((uint64_t)op[2] * 11789567066766345231UL) + ((((uint64_t)op[3] * 2028871269219441641UL) + ((uint64_t)op[4] * 15338898177291282674UL) + ((uint64_t)op[5] * 12439336914687665895UL) + ((uint64_t)op[6] * 18011779295393838587UL) + ((uint64_t)op[7] * 1500806634566570764UL) + ((uint64_t)op[8] * 8853821436126993525UL) + ((uint64_t)op[9] * 11361451194444780761UL) + ((uint64_t)op[10] * 5102579330155346788UL) + ((uint64_t)op[11] * 8305064368330360518UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 8305064368330360518UL) + ((uint64_t)op[1] * 14344461159847157371UL) + ((uint64_t)op[2] * 17578588421717713409UL) + ((uint64_t)op[3] * 11789567066766345231UL) + ((((uint64_t)op[4] * 2028871269219441641UL) + ((uint64_t)op[5] * 15338898177291282674UL) + ((uint64_t)op[6] * 12439336914687665895UL) + ((uint64_t)op[7] * 18011779295393838587UL) + ((uint64_t)op[8] * 1500806634566570764UL) + ((uint64_t)op[9] * 8853821436126993525UL) + ((uint64_t)op[10] * 11361451194444780761UL) + ((uint64_t)op[11] * 5102579330155346788UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 5102579330155346788UL) + ((uint64_t)op[1] * 8305064368330360518UL) + ((uint64_t)op[2] * 14344461159847157371UL) + ((uint64_t)op[3] * 17578588421717713409UL) + ((uint64_t)op[4] * 11789567066766345231UL) + ((((uint64_t)op[5] * 2028871269219441641UL) + ((uint64_t)op[6] * 15338898177291282674UL) + ((uint64_t)op[7] * 12439336914687665895UL) + ((uint64_t)op[8] * 18011779295393838587UL) + ((uint64_t)op[9] * 1500806634566570764UL) + ((uint64_t)op[10] * 8853821436126993525UL) + ((uint64_t)op[11] * 11361451194444780761UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 11361451194444780761UL) + ((uint64_t)op[1] * 5102579330155346788UL) + ((uint64_t)op[2] * 8305064368330360518UL) + ((uint64_t)op[3] * 14344461159847157371UL) + ((uint64_t)op[4] * 17578588421717713409UL) + ((uint64_t)op[5] * 11789567066766345231UL) + ((((uint64_t)op[6] * 2028871269219441641UL) + ((uint64_t)op[7] * 15338898177291282674UL) + ((uint64_t)op[8] * 12439336914687665895UL) + ((uint64_t)op[9] * 18011779295393838587UL) + ((uint64_t)op[10] * 1500806634566570764UL) + ((uint64_t)op[11] * 8853821436126993525UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 8853821436126993525UL) + ((uint64_t)op[1] * 11361451194444780761UL) + ((uint64_t)op[2] * 5102579330155346788UL) + ((uint64_t)op[3] * 8305064368330360518UL) + ((uint64_t)op[4] * 14344461159847157371UL) + ((uint64_t)op[5] * 17578588421717713409UL) + ((uint64_t)op[6] * 11789567066766345231UL) + ((((uint64_t)op[7] * 2028871269219441641UL) + ((uint64_t)op[8] * 15338898177291282674UL) + ((uint64_t)op[9] * 12439336914687665895UL) + ((uint64_t)op[10] * 18011779295393838587UL) + ((uint64_t)op[11] * 1500806634566570764UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 1500806634566570764UL) + ((uint64_t)op[1] * 8853821436126993525UL) + ((uint64_t)op[2] * 11361451194444780761UL) + ((uint64_t)op[3] * 5102579330155346788UL) + ((uint64_t)op[4] * 8305064368330360518UL) + ((uint64_t)op[5] * 14344461159847157371UL) + ((uint64_t)op[6] * 17578588421717713409UL) + ((uint64_t)op[7] * 11789567066766345231UL) + ((((uint64_t)op[8] * 2028871269219441641UL) + ((uint64_t)op[9] * 15338898177291282674UL) + ((uint64_t)op[10] * 12439336914687665895UL) + ((uint64_t)op[11] * 18011779295393838587UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 18011779295393838587UL) + ((uint64_t)op[1] * 1500806634566570764UL) + ((uint64_t)op[2] * 8853821436126993525UL) + ((uint64_t)op[3] * 11361451194444780761UL) + ((uint64_t)op[4] * 5102579330155346788UL) + ((uint64_t)op[5] * 8305064368330360518UL) + ((uint64_t)op[6] * 14344461159847157371UL) + ((uint64_t)op[7] * 17578588421717713409UL) + ((uint64_t)op[8] * 11789567066766345231UL) + ((((uint64_t)op[9] * 2028871269219441641UL) + ((uint64_t)op[10] * 15338898177291282674UL) + ((uint64_t)op[11] * 12439336914687665895UL)) * 18446744073709551614);
	tmp_q[9] = ((uint64_t)op[0] * 12439336914687665895UL) + ((uint64_t)op[1] * 18011779295393838587UL) + ((uint64_t)op[2] * 1500806634566570764UL) + ((uint64_t)op[3] * 8853821436126993525UL) + ((uint64_t)op[4] * 11361451194444780761UL) + ((uint64_t)op[5] * 5102579330155346788UL) + ((uint64_t)op[6] * 8305064368330360518UL) + ((uint64_t)op[7] * 14344461159847157371UL) + ((uint64_t)op[8] * 17578588421717713409UL) + ((uint64_t)op[9] * 11789567066766345231UL) + ((((uint64_t)op[10] * 2028871269219441641UL) + ((uint64_t)op[11] * 15338898177291282674UL)) * 18446744073709551614);
	tmp_q[10] = ((uint64_t)op[0] * 15338898177291282674UL) + ((uint64_t)op[1] * 12439336914687665895UL) + ((uint64_t)op[2] * 18011779295393838587UL) + ((uint64_t)op[3] * 1500806634566570764UL) + ((uint64_t)op[4] * 8853821436126993525UL) + ((uint64_t)op[5] * 11361451194444780761UL) + ((uint64_t)op[6] * 5102579330155346788UL) + ((uint64_t)op[7] * 8305064368330360518UL) + ((uint64_t)op[8] * 14344461159847157371UL) + ((uint64_t)op[9] * 17578588421717713409UL) + ((uint64_t)op[10] * 11789567066766345231UL) + ((uint64_t)op[11] * 14389001535270668334UL);
	tmp_q[11] = ((uint64_t)op[0] * 2028871269219441641UL) + ((uint64_t)op[1] * 15338898177291282674UL) + ((uint64_t)op[2] * 12439336914687665895UL) + ((uint64_t)op[3] * 18011779295393838587UL) + ((uint64_t)op[4] * 1500806634566570764UL) + ((uint64_t)op[5] * 8853821436126993525UL) + ((uint64_t)op[6] * 11361451194444780761UL) + ((uint64_t)op[7] * 5102579330155346788UL) + ((uint64_t)op[8] * 8305064368330360518UL) + ((uint64_t)op[9] * 14344461159847157371UL) + ((uint64_t)op[10] * 17578588421717713409UL) + ((uint64_t)op[11] * 11789567066766345231UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4844439148275L) - ((((int128)tmp_q[1] * 337957930922L) + ((int128)tmp_q[2] * 5400063526923L) - ((int128)tmp_q[3] * 79282358532L) + ((int128)tmp_q[4] * 642884628788L) - ((int128)tmp_q[5] * 1430294983036L) + ((int128)tmp_q[6] * 5468272947938L) + ((int128)tmp_q[7] * 3887917449271L) + ((int128)tmp_q[8] * 1431129859387L) + ((int128)tmp_q[9] * 5041118435097L) - ((int128)tmp_q[10] * 2471683051640L) - ((int128)tmp_q[11] * 3824142592011L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 3824142592011L) - ((int128)tmp_q[1] * 4844439148275L) - ((((int128)tmp_q[2] * 337957930922L) + ((int128)tmp_q[3] * 5400063526923L) - ((int128)tmp_q[4] * 79282358532L) + ((int128)tmp_q[5] * 642884628788L) - ((int128)tmp_q[6] * 1430294983036L) + ((int128)tmp_q[7] * 5468272947938L) + ((int128)tmp_q[8] * 3887917449271L) + ((int128)tmp_q[9] * 1431129859387L) + ((int128)tmp_q[10] * 5041118435097L) - ((int128)tmp_q[11] * 2471683051640L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 2471683051640L) - ((int128)tmp_q[1] * 3824142592011L) - ((int128)tmp_q[2] * 4844439148275L) - ((((int128)tmp_q[3] * 337957930922L) + ((int128)tmp_q[4] * 5400063526923L) - ((int128)tmp_q[5] * 79282358532L) + ((int128)tmp_q[6] * 642884628788L) - ((int128)tmp_q[7] * 1430294983036L) + ((int128)tmp_q[8] * 5468272947938L) + ((int128)tmp_q[9] * 3887917449271L) + ((int128)tmp_q[10] * 1431129859387L) + ((int128)tmp_q[11] * 5041118435097L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 5041118435097L) - ((int128)tmp_q[1] * 2471683051640L) - ((int128)tmp_q[2] * 3824142592011L) - ((int128)tmp_q[3] * 4844439148275L) - ((((int128)tmp_q[4] * 337957930922L) + ((int128)tmp_q[5] * 5400063526923L) - ((int128)tmp_q[6] * 79282358532L) + ((int128)tmp_q[7] * 642884628788L) - ((int128)tmp_q[8] * 1430294983036L) + ((int128)tmp_q[9] * 5468272947938L) + ((int128)tmp_q[10] * 3887917449271L) + ((int128)tmp_q[11] * 1431129859387L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 1431129859387L) + ((int128)tmp_q[1] * 5041118435097L) - ((int128)tmp_q[2] * 2471683051640L) - ((int128)tmp_q[3] * 3824142592011L) - ((int128)tmp_q[4] * 4844439148275L) - ((((int128)tmp_q[5] * 337957930922L) + ((int128)tmp_q[6] * 5400063526923L) - ((int128)tmp_q[7] * 79282358532L) + ((int128)tmp_q[8] * 642884628788L) - ((int128)tmp_q[9] * 1430294983036L) + ((int128)tmp_q[10] * 5468272947938L) + ((int128)tmp_q[11] * 3887917449271L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 3887917449271L) + ((int128)tmp_q[1] * 1431129859387L) + ((int128)tmp_q[2] * 5041118435097L) - ((int128)tmp_q[3] * 2471683051640L) - ((int128)tmp_q[4] * 3824142592011L) - ((int128)tmp_q[5] * 4844439148275L) - ((((int128)tmp_q[6] * 337957930922L) + ((int128)tmp_q[7] * 5400063526923L) - ((int128)tmp_q[8] * 79282358532L) + ((int128)tmp_q[9] * 642884628788L) - ((int128)tmp_q[10] * 1430294983036L) + ((int128)tmp_q[11] * 5468272947938L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 5468272947938L) + ((int128)tmp_q[1] * 3887917449271L) + ((int128)tmp_q[2] * 1431129859387L) + ((int128)tmp_q[3] * 5041118435097L) - ((int128)tmp_q[4] * 2471683051640L) - ((int128)tmp_q[5] * 3824142592011L) - ((int128)tmp_q[6] * 4844439148275L) - ((((int128)tmp_q[7] * 337957930922L) + ((int128)tmp_q[8] * 5400063526923L) - ((int128)tmp_q[9] * 79282358532L) + ((int128)tmp_q[10] * 642884628788L) - ((int128)tmp_q[11] * 1430294983036L)) * 2);
	tmp_zero[7] = -((int128)tmp_q[0] * 1430294983036L) + ((int128)tmp_q[1] * 5468272947938L) + ((int128)tmp_q[2] * 3887917449271L) + ((int128)tmp_q[3] * 1431129859387L) + ((int128)tmp_q[4] * 5041118435097L) - ((int128)tmp_q[5] * 2471683051640L) - ((int128)tmp_q[6] * 3824142592011L) - ((int128)tmp_q[7] * 4844439148275L) - ((((int128)tmp_q[8] * 337957930922L) + ((int128)tmp_q[9] * 5400063526923L) - ((int128)tmp_q[10] * 79282358532L) + ((int128)tmp_q[11] * 642884628788L)) * 2);
	tmp_zero[8] = ((int128)tmp_q[0] * 642884628788L) - ((int128)tmp_q[1] * 1430294983036L) + ((int128)tmp_q[2] * 5468272947938L) + ((int128)tmp_q[3] * 3887917449271L) + ((int128)tmp_q[4] * 1431129859387L) + ((int128)tmp_q[5] * 5041118435097L) - ((int128)tmp_q[6] * 2471683051640L) - ((int128)tmp_q[7] * 3824142592011L) - ((int128)tmp_q[8] * 4844439148275L) - ((((int128)tmp_q[9] * 337957930922L) + ((int128)tmp_q[10] * 5400063526923L) - ((int128)tmp_q[11] * 79282358532L)) * 2);
	tmp_zero[9] = -((int128)tmp_q[0] * 79282358532L) + ((int128)tmp_q[1] * 642884628788L) - ((int128)tmp_q[2] * 1430294983036L) + ((int128)tmp_q[3] * 5468272947938L) + ((int128)tmp_q[4] * 3887917449271L) + ((int128)tmp_q[5] * 1431129859387L) + ((int128)tmp_q[6] * 5041118435097L) - ((int128)tmp_q[7] * 2471683051640L) - ((int128)tmp_q[8] * 3824142592011L) - ((int128)tmp_q[9] * 4844439148275L) - ((((int128)tmp_q[10] * 337957930922L) + ((int128)tmp_q[11] * 5400063526923L)) * 2);
	tmp_zero[10] = ((int128)tmp_q[0] * 5400063526923L) - ((int128)tmp_q[1] * 79282358532L) + ((int128)tmp_q[2] * 642884628788L) - ((int128)tmp_q[3] * 1430294983036L) + ((int128)tmp_q[4] * 5468272947938L) + ((int128)tmp_q[5] * 3887917449271L) + ((int128)tmp_q[6] * 1431129859387L) + ((int128)tmp_q[7] * 5041118435097L) - ((int128)tmp_q[8] * 2471683051640L) - ((int128)tmp_q[9] * 3824142592011L) - ((int128)tmp_q[10] * 4844439148275L) - ((int128)tmp_q[11] * 675915861844L);
	tmp_zero[11] = ((int128)tmp_q[0] * 337957930922L) + ((int128)tmp_q[1] * 5400063526923L) - ((int128)tmp_q[2] * 79282358532L) + ((int128)tmp_q[3] * 642884628788L) - ((int128)tmp_q[4] * 1430294983036L) + ((int128)tmp_q[5] * 5468272947938L) + ((int128)tmp_q[6] * 3887917449271L) + ((int128)tmp_q[7] * 1431129859387L) + ((int128)tmp_q[8] * 5041118435097L) - ((int128)tmp_q[9] * 2471683051640L) - ((int128)tmp_q[10] * 3824142592011L) - ((int128)tmp_q[11] * 4844439148275L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
	rop[8] = (op[8] + tmp_zero[8]) >> WORD_SIZE;
	rop[9] = (op[9] + tmp_zero[9]) >> WORD_SIZE;
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
	rop[11] = (op[11] + tmp_zero[11]) >> WORD_SIZE;
}

