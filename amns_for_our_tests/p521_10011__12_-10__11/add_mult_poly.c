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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 10);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 10);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 10);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 10);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 10);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 10);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 10);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 10);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 10);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 10);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 10);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 20);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 20);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 20);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 10);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 20);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 10);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 20);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 10);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12885181610280865701UL) + ((((uint64_t)op[1] * 7907227396645182103UL) + ((uint64_t)op[2] * 12894134397415638745UL) + ((uint64_t)op[3] * 13403095976817539657UL) + ((uint64_t)op[4] * 1038689232975079401UL) + ((uint64_t)op[5] * 3949237454037429749UL) + ((uint64_t)op[6] * 17414591964777584421UL) + ((uint64_t)op[7] * 14978318491058378933UL) + ((uint64_t)op[8] * 9004146829020004597UL) + ((uint64_t)op[9] * 4165204667067858452UL) + ((uint64_t)op[10] * 13433053810240343853UL) + ((uint64_t)op[11] * 4038616898171338885UL)) * 18446744073709551606);
	tmp_q[1] = ((uint64_t)op[0] * 4038616898171338885UL) + ((uint64_t)op[1] * 12885181610280865701UL) + ((((uint64_t)op[2] * 7907227396645182103UL) + ((uint64_t)op[3] * 12894134397415638745UL) + ((uint64_t)op[4] * 13403095976817539657UL) + ((uint64_t)op[5] * 1038689232975079401UL) + ((uint64_t)op[6] * 3949237454037429749UL) + ((uint64_t)op[7] * 17414591964777584421UL) + ((uint64_t)op[8] * 14978318491058378933UL) + ((uint64_t)op[9] * 9004146829020004597UL) + ((uint64_t)op[10] * 4165204667067858452UL) + ((uint64_t)op[11] * 13433053810240343853UL)) * 18446744073709551606);
	tmp_q[2] = ((uint64_t)op[0] * 13433053810240343853UL) + ((uint64_t)op[1] * 4038616898171338885UL) + ((uint64_t)op[2] * 12885181610280865701UL) + ((((uint64_t)op[3] * 7907227396645182103UL) + ((uint64_t)op[4] * 12894134397415638745UL) + ((uint64_t)op[5] * 13403095976817539657UL) + ((uint64_t)op[6] * 1038689232975079401UL) + ((uint64_t)op[7] * 3949237454037429749UL) + ((uint64_t)op[8] * 17414591964777584421UL) + ((uint64_t)op[9] * 14978318491058378933UL) + ((uint64_t)op[10] * 9004146829020004597UL) + ((uint64_t)op[11] * 4165204667067858452UL)) * 18446744073709551606);
	tmp_q[3] = ((uint64_t)op[0] * 4165204667067858452UL) + ((uint64_t)op[1] * 13433053810240343853UL) + ((uint64_t)op[2] * 4038616898171338885UL) + ((uint64_t)op[3] * 12885181610280865701UL) + ((((uint64_t)op[4] * 7907227396645182103UL) + ((uint64_t)op[5] * 12894134397415638745UL) + ((uint64_t)op[6] * 13403095976817539657UL) + ((uint64_t)op[7] * 1038689232975079401UL) + ((uint64_t)op[8] * 3949237454037429749UL) + ((uint64_t)op[9] * 17414591964777584421UL) + ((uint64_t)op[10] * 14978318491058378933UL) + ((uint64_t)op[11] * 9004146829020004597UL)) * 18446744073709551606);
	tmp_q[4] = ((uint64_t)op[0] * 9004146829020004597UL) + ((uint64_t)op[1] * 4165204667067858452UL) + ((uint64_t)op[2] * 13433053810240343853UL) + ((uint64_t)op[3] * 4038616898171338885UL) + ((uint64_t)op[4] * 12885181610280865701UL) + ((((uint64_t)op[5] * 7907227396645182103UL) + ((uint64_t)op[6] * 12894134397415638745UL) + ((uint64_t)op[7] * 13403095976817539657UL) + ((uint64_t)op[8] * 1038689232975079401UL) + ((uint64_t)op[9] * 3949237454037429749UL) + ((uint64_t)op[10] * 17414591964777584421UL) + ((uint64_t)op[11] * 14978318491058378933UL)) * 18446744073709551606);
	tmp_q[5] = ((uint64_t)op[0] * 14978318491058378933UL) + ((uint64_t)op[1] * 9004146829020004597UL) + ((uint64_t)op[2] * 4165204667067858452UL) + ((uint64_t)op[3] * 13433053810240343853UL) + ((uint64_t)op[4] * 4038616898171338885UL) + ((uint64_t)op[5] * 12885181610280865701UL) + ((((uint64_t)op[6] * 7907227396645182103UL) + ((uint64_t)op[7] * 12894134397415638745UL) + ((uint64_t)op[8] * 13403095976817539657UL) + ((uint64_t)op[9] * 1038689232975079401UL) + ((uint64_t)op[10] * 3949237454037429749UL) + ((uint64_t)op[11] * 17414591964777584421UL)) * 18446744073709551606);
	tmp_q[6] = ((uint64_t)op[0] * 17414591964777584421UL) + ((uint64_t)op[1] * 14978318491058378933UL) + ((uint64_t)op[2] * 9004146829020004597UL) + ((uint64_t)op[3] * 4165204667067858452UL) + ((uint64_t)op[4] * 13433053810240343853UL) + ((uint64_t)op[5] * 4038616898171338885UL) + ((uint64_t)op[6] * 12885181610280865701UL) + ((((uint64_t)op[7] * 7907227396645182103UL) + ((uint64_t)op[8] * 12894134397415638745UL) + ((uint64_t)op[9] * 13403095976817539657UL) + ((uint64_t)op[10] * 1038689232975079401UL) + ((uint64_t)op[11] * 3949237454037429749UL)) * 18446744073709551606);
	tmp_q[7] = ((uint64_t)op[0] * 3949237454037429749UL) + ((uint64_t)op[1] * 17414591964777584421UL) + ((uint64_t)op[2] * 14978318491058378933UL) + ((uint64_t)op[3] * 9004146829020004597UL) + ((uint64_t)op[4] * 4165204667067858452UL) + ((uint64_t)op[5] * 13433053810240343853UL) + ((uint64_t)op[6] * 4038616898171338885UL) + ((uint64_t)op[7] * 12885181610280865701UL) + ((((uint64_t)op[8] * 7907227396645182103UL) + ((uint64_t)op[9] * 12894134397415638745UL) + ((uint64_t)op[10] * 13403095976817539657UL) + ((uint64_t)op[11] * 1038689232975079401UL)) * 18446744073709551606);
	tmp_q[8] = ((uint64_t)op[0] * 1038689232975079401UL) + ((uint64_t)op[1] * 3949237454037429749UL) + ((uint64_t)op[2] * 17414591964777584421UL) + ((uint64_t)op[3] * 14978318491058378933UL) + ((uint64_t)op[4] * 9004146829020004597UL) + ((uint64_t)op[5] * 4165204667067858452UL) + ((uint64_t)op[6] * 13433053810240343853UL) + ((uint64_t)op[7] * 4038616898171338885UL) + ((uint64_t)op[8] * 12885181610280865701UL) + ((((uint64_t)op[9] * 7907227396645182103UL) + ((uint64_t)op[10] * 12894134397415638745UL) + ((uint64_t)op[11] * 13403095976817539657UL)) * 18446744073709551606);
	tmp_q[9] = ((uint64_t)op[0] * 13403095976817539657UL) + ((uint64_t)op[1] * 1038689232975079401UL) + ((uint64_t)op[2] * 3949237454037429749UL) + ((uint64_t)op[3] * 17414591964777584421UL) + ((uint64_t)op[4] * 14978318491058378933UL) + ((uint64_t)op[5] * 9004146829020004597UL) + ((uint64_t)op[6] * 4165204667067858452UL) + ((uint64_t)op[7] * 13433053810240343853UL) + ((uint64_t)op[8] * 4038616898171338885UL) + ((uint64_t)op[9] * 12885181610280865701UL) + ((((uint64_t)op[10] * 7907227396645182103UL) + ((uint64_t)op[11] * 12894134397415638745UL)) * 18446744073709551606);
	tmp_q[10] = ((uint64_t)op[0] * 12894134397415638745UL) + ((uint64_t)op[1] * 13403095976817539657UL) + ((uint64_t)op[2] * 1038689232975079401UL) + ((uint64_t)op[3] * 3949237454037429749UL) + ((uint64_t)op[4] * 17414591964777584421UL) + ((uint64_t)op[5] * 14978318491058378933UL) + ((uint64_t)op[6] * 9004146829020004597UL) + ((uint64_t)op[7] * 4165204667067858452UL) + ((uint64_t)op[8] * 13433053810240343853UL) + ((uint64_t)op[9] * 4038616898171338885UL) + ((uint64_t)op[10] * 12885181610280865701UL) + ((uint64_t)op[11] * 13161446402095937050UL);
	tmp_q[11] = ((uint64_t)op[0] * 7907227396645182103UL) + ((uint64_t)op[1] * 12894134397415638745UL) + ((uint64_t)op[2] * 13403095976817539657UL) + ((uint64_t)op[3] * 1038689232975079401UL) + ((uint64_t)op[4] * 3949237454037429749UL) + ((uint64_t)op[5] * 17414591964777584421UL) + ((uint64_t)op[6] * 14978318491058378933UL) + ((uint64_t)op[7] * 9004146829020004597UL) + ((uint64_t)op[8] * 4165204667067858452UL) + ((uint64_t)op[9] * 13433053810240343853UL) + ((uint64_t)op[10] * 4038616898171338885UL) + ((uint64_t)op[11] * 12885181610280865701UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4246452819567L) - ((-((int128)tmp_q[1] * 2390158565510L) - ((int128)tmp_q[2] * 528575677410L) - ((int128)tmp_q[3] * 228933567964L) - ((int128)tmp_q[4] * 1435629064887L) + ((int128)tmp_q[5] * 1618988219213L) - ((int128)tmp_q[6] * 1781836344293L) - ((int128)tmp_q[7] * 2565792588581L) + ((int128)tmp_q[8] * 5404511451548L) - ((int128)tmp_q[9] * 5601342120431L) + ((int128)tmp_q[10] * 3921299628540L) - ((int128)tmp_q[11] * 951071846257L)) * 10);
	tmp_zero[1] = -((int128)tmp_q[0] * 951071846257L) + ((int128)tmp_q[1] * 4246452819567L) - ((-((int128)tmp_q[2] * 2390158565510L) - ((int128)tmp_q[3] * 528575677410L) - ((int128)tmp_q[4] * 228933567964L) - ((int128)tmp_q[5] * 1435629064887L) + ((int128)tmp_q[6] * 1618988219213L) - ((int128)tmp_q[7] * 1781836344293L) - ((int128)tmp_q[8] * 2565792588581L) + ((int128)tmp_q[9] * 5404511451548L) - ((int128)tmp_q[10] * 5601342120431L) + ((int128)tmp_q[11] * 3921299628540L)) * 10);
	tmp_zero[2] = ((int128)tmp_q[0] * 3921299628540L) - ((int128)tmp_q[1] * 951071846257L) + ((int128)tmp_q[2] * 4246452819567L) - ((-((int128)tmp_q[3] * 2390158565510L) - ((int128)tmp_q[4] * 528575677410L) - ((int128)tmp_q[5] * 228933567964L) - ((int128)tmp_q[6] * 1435629064887L) + ((int128)tmp_q[7] * 1618988219213L) - ((int128)tmp_q[8] * 1781836344293L) - ((int128)tmp_q[9] * 2565792588581L) + ((int128)tmp_q[10] * 5404511451548L) - ((int128)tmp_q[11] * 5601342120431L)) * 10);
	tmp_zero[3] = -((int128)tmp_q[0] * 5601342120431L) + ((int128)tmp_q[1] * 3921299628540L) - ((int128)tmp_q[2] * 951071846257L) + ((int128)tmp_q[3] * 4246452819567L) - ((-((int128)tmp_q[4] * 2390158565510L) - ((int128)tmp_q[5] * 528575677410L) - ((int128)tmp_q[6] * 228933567964L) - ((int128)tmp_q[7] * 1435629064887L) + ((int128)tmp_q[8] * 1618988219213L) - ((int128)tmp_q[9] * 1781836344293L) - ((int128)tmp_q[10] * 2565792588581L) + ((int128)tmp_q[11] * 5404511451548L)) * 10);
	tmp_zero[4] = ((int128)tmp_q[0] * 5404511451548L) - ((int128)tmp_q[1] * 5601342120431L) + ((int128)tmp_q[2] * 3921299628540L) - ((int128)tmp_q[3] * 951071846257L) + ((int128)tmp_q[4] * 4246452819567L) - ((-((int128)tmp_q[5] * 2390158565510L) - ((int128)tmp_q[6] * 528575677410L) - ((int128)tmp_q[7] * 228933567964L) - ((int128)tmp_q[8] * 1435629064887L) + ((int128)tmp_q[9] * 1618988219213L) - ((int128)tmp_q[10] * 1781836344293L) - ((int128)tmp_q[11] * 2565792588581L)) * 10);
	tmp_zero[5] = -((int128)tmp_q[0] * 2565792588581L) + ((int128)tmp_q[1] * 5404511451548L) - ((int128)tmp_q[2] * 5601342120431L) + ((int128)tmp_q[3] * 3921299628540L) - ((int128)tmp_q[4] * 951071846257L) + ((int128)tmp_q[5] * 4246452819567L) - ((-((int128)tmp_q[6] * 2390158565510L) - ((int128)tmp_q[7] * 528575677410L) - ((int128)tmp_q[8] * 228933567964L) - ((int128)tmp_q[9] * 1435629064887L) + ((int128)tmp_q[10] * 1618988219213L) - ((int128)tmp_q[11] * 1781836344293L)) * 10);
	tmp_zero[6] = -((int128)tmp_q[0] * 1781836344293L) - ((int128)tmp_q[1] * 2565792588581L) + ((int128)tmp_q[2] * 5404511451548L) - ((int128)tmp_q[3] * 5601342120431L) + ((int128)tmp_q[4] * 3921299628540L) - ((int128)tmp_q[5] * 951071846257L) + ((int128)tmp_q[6] * 4246452819567L) - ((-((int128)tmp_q[7] * 2390158565510L) - ((int128)tmp_q[8] * 528575677410L) - ((int128)tmp_q[9] * 228933567964L) - ((int128)tmp_q[10] * 1435629064887L) + ((int128)tmp_q[11] * 1618988219213L)) * 10);
	tmp_zero[7] = ((int128)tmp_q[0] * 1618988219213L) - ((int128)tmp_q[1] * 1781836344293L) - ((int128)tmp_q[2] * 2565792588581L) + ((int128)tmp_q[3] * 5404511451548L) - ((int128)tmp_q[4] * 5601342120431L) + ((int128)tmp_q[5] * 3921299628540L) - ((int128)tmp_q[6] * 951071846257L) + ((int128)tmp_q[7] * 4246452819567L) - ((-((int128)tmp_q[8] * 2390158565510L) - ((int128)tmp_q[9] * 528575677410L) - ((int128)tmp_q[10] * 228933567964L) - ((int128)tmp_q[11] * 1435629064887L)) * 10);
	tmp_zero[8] = -((int128)tmp_q[0] * 1435629064887L) + ((int128)tmp_q[1] * 1618988219213L) - ((int128)tmp_q[2] * 1781836344293L) - ((int128)tmp_q[3] * 2565792588581L) + ((int128)tmp_q[4] * 5404511451548L) - ((int128)tmp_q[5] * 5601342120431L) + ((int128)tmp_q[6] * 3921299628540L) - ((int128)tmp_q[7] * 951071846257L) + ((int128)tmp_q[8] * 4246452819567L) - ((-((int128)tmp_q[9] * 2390158565510L) - ((int128)tmp_q[10] * 528575677410L) - ((int128)tmp_q[11] * 228933567964L)) * 10);
	tmp_zero[9] = -((int128)tmp_q[0] * 228933567964L) - ((int128)tmp_q[1] * 1435629064887L) + ((int128)tmp_q[2] * 1618988219213L) - ((int128)tmp_q[3] * 1781836344293L) - ((int128)tmp_q[4] * 2565792588581L) + ((int128)tmp_q[5] * 5404511451548L) - ((int128)tmp_q[6] * 5601342120431L) + ((int128)tmp_q[7] * 3921299628540L) - ((int128)tmp_q[8] * 951071846257L) + ((int128)tmp_q[9] * 4246452819567L) - ((-((int128)tmp_q[10] * 2390158565510L) - ((int128)tmp_q[11] * 528575677410L)) * 10);
	tmp_zero[10] = -((int128)tmp_q[0] * 528575677410L) - ((int128)tmp_q[1] * 228933567964L) - ((int128)tmp_q[2] * 1435629064887L) + ((int128)tmp_q[3] * 1618988219213L) - ((int128)tmp_q[4] * 1781836344293L) - ((int128)tmp_q[5] * 2565792588581L) + ((int128)tmp_q[6] * 5404511451548L) - ((int128)tmp_q[7] * 5601342120431L) + ((int128)tmp_q[8] * 3921299628540L) - ((int128)tmp_q[9] * 951071846257L) + ((int128)tmp_q[10] * 4246452819567L) + ((int128)tmp_q[11] * 23901585655100L);
	tmp_zero[11] = -((int128)tmp_q[0] * 2390158565510L) - ((int128)tmp_q[1] * 528575677410L) - ((int128)tmp_q[2] * 228933567964L) - ((int128)tmp_q[3] * 1435629064887L) + ((int128)tmp_q[4] * 1618988219213L) - ((int128)tmp_q[5] * 1781836344293L) - ((int128)tmp_q[6] * 2565792588581L) + ((int128)tmp_q[7] * 5404511451548L) - ((int128)tmp_q[8] * 5601342120431L) + ((int128)tmp_q[9] * 3921299628540L) - ((int128)tmp_q[10] * 951071846257L) + ((int128)tmp_q[11] * 4246452819567L);

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

