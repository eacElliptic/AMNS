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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[8] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 18164864355073873467UL) + ((((uint64_t)op[1] * 10164694433116726296UL) + ((uint64_t)op[2] * 10609789699974493468UL) + ((uint64_t)op[3] * 446823990102881178UL) + ((uint64_t)op[4] * 1707947933092213125UL) + ((uint64_t)op[5] * 12422810313576280597UL) + ((uint64_t)op[6] * 12426480616074726530UL) + ((uint64_t)op[7] * 11761751278650308237UL) + ((uint64_t)op[8] * 10189203805452990300UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 10189203805452990300UL) + ((uint64_t)op[1] * 18164864355073873467UL) + ((((uint64_t)op[2] * 10164694433116726296UL) + ((uint64_t)op[3] * 10609789699974493468UL) + ((uint64_t)op[4] * 446823990102881178UL) + ((uint64_t)op[5] * 1707947933092213125UL) + ((uint64_t)op[6] * 12422810313576280597UL) + ((uint64_t)op[7] * 12426480616074726530UL) + ((uint64_t)op[8] * 11761751278650308237UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 11761751278650308237UL) + ((uint64_t)op[1] * 10189203805452990300UL) + ((uint64_t)op[2] * 18164864355073873467UL) + ((((uint64_t)op[3] * 10164694433116726296UL) + ((uint64_t)op[4] * 10609789699974493468UL) + ((uint64_t)op[5] * 446823990102881178UL) + ((uint64_t)op[6] * 1707947933092213125UL) + ((uint64_t)op[7] * 12422810313576280597UL) + ((uint64_t)op[8] * 12426480616074726530UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 12426480616074726530UL) + ((uint64_t)op[1] * 11761751278650308237UL) + ((uint64_t)op[2] * 10189203805452990300UL) + ((uint64_t)op[3] * 18164864355073873467UL) + ((((uint64_t)op[4] * 10164694433116726296UL) + ((uint64_t)op[5] * 10609789699974493468UL) + ((uint64_t)op[6] * 446823990102881178UL) + ((uint64_t)op[7] * 1707947933092213125UL) + ((uint64_t)op[8] * 12422810313576280597UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 12422810313576280597UL) + ((uint64_t)op[1] * 12426480616074726530UL) + ((uint64_t)op[2] * 11761751278650308237UL) + ((uint64_t)op[3] * 10189203805452990300UL) + ((uint64_t)op[4] * 18164864355073873467UL) + ((((uint64_t)op[5] * 10164694433116726296UL) + ((uint64_t)op[6] * 10609789699974493468UL) + ((uint64_t)op[7] * 446823990102881178UL) + ((uint64_t)op[8] * 1707947933092213125UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 1707947933092213125UL) + ((uint64_t)op[1] * 12422810313576280597UL) + ((uint64_t)op[2] * 12426480616074726530UL) + ((uint64_t)op[3] * 11761751278650308237UL) + ((uint64_t)op[4] * 10189203805452990300UL) + ((uint64_t)op[5] * 18164864355073873467UL) + ((((uint64_t)op[6] * 10164694433116726296UL) + ((uint64_t)op[7] * 10609789699974493468UL) + ((uint64_t)op[8] * 446823990102881178UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 446823990102881178UL) + ((uint64_t)op[1] * 1707947933092213125UL) + ((uint64_t)op[2] * 12422810313576280597UL) + ((uint64_t)op[3] * 12426480616074726530UL) + ((uint64_t)op[4] * 11761751278650308237UL) + ((uint64_t)op[5] * 10189203805452990300UL) + ((uint64_t)op[6] * 18164864355073873467UL) + ((((uint64_t)op[7] * 10164694433116726296UL) + ((uint64_t)op[8] * 10609789699974493468UL)) * 2);
	tmp_q[7] = ((uint64_t)op[0] * 10609789699974493468UL) + ((uint64_t)op[1] * 446823990102881178UL) + ((uint64_t)op[2] * 1707947933092213125UL) + ((uint64_t)op[3] * 12422810313576280597UL) + ((uint64_t)op[4] * 12426480616074726530UL) + ((uint64_t)op[5] * 11761751278650308237UL) + ((uint64_t)op[6] * 10189203805452990300UL) + ((uint64_t)op[7] * 18164864355073873467UL) + ((uint64_t)op[8] * 1882644792523900976UL);
	tmp_q[8] = ((uint64_t)op[0] * 10164694433116726296UL) + ((uint64_t)op[1] * 10609789699974493468UL) + ((uint64_t)op[2] * 446823990102881178UL) + ((uint64_t)op[3] * 1707947933092213125UL) + ((uint64_t)op[4] * 12422810313576280597UL) + ((uint64_t)op[5] * 12426480616074726530UL) + ((uint64_t)op[6] * 11761751278650308237UL) + ((uint64_t)op[7] * 10189203805452990300UL) + ((uint64_t)op[8] * 18164864355073873467UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 2523267277865L) + ((((int128)tmp_q[1] * 2148794128443L) + ((int128)tmp_q[2] * 1826938335416L) + ((int128)tmp_q[3] * 1360710426803L) - ((int128)tmp_q[4] * 4353136553459L) - ((int128)tmp_q[5] * 804467892746L) - ((int128)tmp_q[6] * 144834102050L) - ((int128)tmp_q[7] * 3425557318563L) + ((int128)tmp_q[8] * 2388994726578L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 2388994726578L) - ((int128)tmp_q[1] * 2523267277865L) + ((((int128)tmp_q[2] * 2148794128443L) + ((int128)tmp_q[3] * 1826938335416L) + ((int128)tmp_q[4] * 1360710426803L) - ((int128)tmp_q[5] * 4353136553459L) - ((int128)tmp_q[6] * 804467892746L) - ((int128)tmp_q[7] * 144834102050L) - ((int128)tmp_q[8] * 3425557318563L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 3425557318563L) + ((int128)tmp_q[1] * 2388994726578L) - ((int128)tmp_q[2] * 2523267277865L) + ((((int128)tmp_q[3] * 2148794128443L) + ((int128)tmp_q[4] * 1826938335416L) + ((int128)tmp_q[5] * 1360710426803L) - ((int128)tmp_q[6] * 4353136553459L) - ((int128)tmp_q[7] * 804467892746L) - ((int128)tmp_q[8] * 144834102050L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 144834102050L) - ((int128)tmp_q[1] * 3425557318563L) + ((int128)tmp_q[2] * 2388994726578L) - ((int128)tmp_q[3] * 2523267277865L) + ((((int128)tmp_q[4] * 2148794128443L) + ((int128)tmp_q[5] * 1826938335416L) + ((int128)tmp_q[6] * 1360710426803L) - ((int128)tmp_q[7] * 4353136553459L) - ((int128)tmp_q[8] * 804467892746L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 804467892746L) - ((int128)tmp_q[1] * 144834102050L) - ((int128)tmp_q[2] * 3425557318563L) + ((int128)tmp_q[3] * 2388994726578L) - ((int128)tmp_q[4] * 2523267277865L) + ((((int128)tmp_q[5] * 2148794128443L) + ((int128)tmp_q[6] * 1826938335416L) + ((int128)tmp_q[7] * 1360710426803L) - ((int128)tmp_q[8] * 4353136553459L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 4353136553459L) - ((int128)tmp_q[1] * 804467892746L) - ((int128)tmp_q[2] * 144834102050L) - ((int128)tmp_q[3] * 3425557318563L) + ((int128)tmp_q[4] * 2388994726578L) - ((int128)tmp_q[5] * 2523267277865L) + ((((int128)tmp_q[6] * 2148794128443L) + ((int128)tmp_q[7] * 1826938335416L) + ((int128)tmp_q[8] * 1360710426803L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 1360710426803L) - ((int128)tmp_q[1] * 4353136553459L) - ((int128)tmp_q[2] * 804467892746L) - ((int128)tmp_q[3] * 144834102050L) - ((int128)tmp_q[4] * 3425557318563L) + ((int128)tmp_q[5] * 2388994726578L) - ((int128)tmp_q[6] * 2523267277865L) + ((((int128)tmp_q[7] * 2148794128443L) + ((int128)tmp_q[8] * 1826938335416L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 1826938335416L) + ((int128)tmp_q[1] * 1360710426803L) - ((int128)tmp_q[2] * 4353136553459L) - ((int128)tmp_q[3] * 804467892746L) - ((int128)tmp_q[4] * 144834102050L) - ((int128)tmp_q[5] * 3425557318563L) + ((int128)tmp_q[6] * 2388994726578L) - ((int128)tmp_q[7] * 2523267277865L) + ((int128)tmp_q[8] * 4297588256886L);
	tmp_zero[8] = ((int128)tmp_q[0] * 2148794128443L) + ((int128)tmp_q[1] * 1826938335416L) + ((int128)tmp_q[2] * 1360710426803L) - ((int128)tmp_q[3] * 4353136553459L) - ((int128)tmp_q[4] * 804467892746L) - ((int128)tmp_q[5] * 144834102050L) - ((int128)tmp_q[6] * 3425557318563L) + ((int128)tmp_q[7] * 2388994726578L) - ((int128)tmp_q[8] * 2523267277865L);

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
}

