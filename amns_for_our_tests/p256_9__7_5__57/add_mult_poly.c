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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 6070273197117422846UL) + ((((uint64_t)op[1] * 10702859500784045706UL) + ((uint64_t)op[2] * 12970089250144367838UL) + ((uint64_t)op[3] * 5743662953381719076UL) + ((uint64_t)op[4] * 9613304865984797390UL) + ((uint64_t)op[5] * 8105231720640713465UL) + ((uint64_t)op[6] * 12066746116129163176UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 12066746116129163176UL) + ((uint64_t)op[1] * 6070273197117422846UL) + ((((uint64_t)op[2] * 10702859500784045706UL) + ((uint64_t)op[3] * 12970089250144367838UL) + ((uint64_t)op[4] * 5743662953381719076UL) + ((uint64_t)op[5] * 9613304865984797390UL) + ((uint64_t)op[6] * 8105231720640713465UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 8105231720640713465UL) + ((uint64_t)op[1] * 12066746116129163176UL) + ((uint64_t)op[2] * 6070273197117422846UL) + ((((uint64_t)op[3] * 10702859500784045706UL) + ((uint64_t)op[4] * 12970089250144367838UL) + ((uint64_t)op[5] * 5743662953381719076UL) + ((uint64_t)op[6] * 9613304865984797390UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 9613304865984797390UL) + ((uint64_t)op[1] * 8105231720640713465UL) + ((uint64_t)op[2] * 12066746116129163176UL) + ((uint64_t)op[3] * 6070273197117422846UL) + ((((uint64_t)op[4] * 10702859500784045706UL) + ((uint64_t)op[5] * 12970089250144367838UL) + ((uint64_t)op[6] * 5743662953381719076UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 5743662953381719076UL) + ((uint64_t)op[1] * 9613304865984797390UL) + ((uint64_t)op[2] * 8105231720640713465UL) + ((uint64_t)op[3] * 12066746116129163176UL) + ((uint64_t)op[4] * 6070273197117422846UL) + ((((uint64_t)op[5] * 10702859500784045706UL) + ((uint64_t)op[6] * 12970089250144367838UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 12970089250144367838UL) + ((uint64_t)op[1] * 5743662953381719076UL) + ((uint64_t)op[2] * 9613304865984797390UL) + ((uint64_t)op[3] * 8105231720640713465UL) + ((uint64_t)op[4] * 12066746116129163176UL) + ((uint64_t)op[5] * 6070273197117422846UL) + ((uint64_t)op[6] * 16620809356501125298UL);
	tmp_q[6] = ((uint64_t)op[0] * 10702859500784045706UL) + ((uint64_t)op[1] * 12970089250144367838UL) + ((uint64_t)op[2] * 5743662953381719076UL) + ((uint64_t)op[3] * 9613304865984797390UL) + ((uint64_t)op[4] * 8105231720640713465UL) + ((uint64_t)op[5] * 12066746116129163176UL) + ((uint64_t)op[6] * 6070273197117422846UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 34474813128L) + ((((int128)tmp_q[1] * 56339285762L) + ((int128)tmp_q[2] * 21012271307L) + ((int128)tmp_q[3] * 38184859812L) + ((int128)tmp_q[4] * 49128047966L) - ((int128)tmp_q[5] * 3257067862L) + ((int128)tmp_q[6] * 35243808746L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 35243808746L) + ((int128)tmp_q[1] * 34474813128L) + ((((int128)tmp_q[2] * 56339285762L) + ((int128)tmp_q[3] * 21012271307L) + ((int128)tmp_q[4] * 38184859812L) + ((int128)tmp_q[5] * 49128047966L) - ((int128)tmp_q[6] * 3257067862L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 3257067862L) + ((int128)tmp_q[1] * 35243808746L) + ((int128)tmp_q[2] * 34474813128L) + ((((int128)tmp_q[3] * 56339285762L) + ((int128)tmp_q[4] * 21012271307L) + ((int128)tmp_q[5] * 38184859812L) + ((int128)tmp_q[6] * 49128047966L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 49128047966L) - ((int128)tmp_q[1] * 3257067862L) + ((int128)tmp_q[2] * 35243808746L) + ((int128)tmp_q[3] * 34474813128L) + ((((int128)tmp_q[4] * 56339285762L) + ((int128)tmp_q[5] * 21012271307L) + ((int128)tmp_q[6] * 38184859812L)) * 5);
	tmp_zero[4] = ((int128)tmp_q[0] * 38184859812L) + ((int128)tmp_q[1] * 49128047966L) - ((int128)tmp_q[2] * 3257067862L) + ((int128)tmp_q[3] * 35243808746L) + ((int128)tmp_q[4] * 34474813128L) + ((((int128)tmp_q[5] * 56339285762L) + ((int128)tmp_q[6] * 21012271307L)) * 5);
	tmp_zero[5] = ((int128)tmp_q[0] * 21012271307L) + ((int128)tmp_q[1] * 38184859812L) + ((int128)tmp_q[2] * 49128047966L) - ((int128)tmp_q[3] * 3257067862L) + ((int128)tmp_q[4] * 35243808746L) + ((int128)tmp_q[5] * 34474813128L) + ((int128)tmp_q[6] * 281696428810L);
	tmp_zero[6] = ((int128)tmp_q[0] * 56339285762L) + ((int128)tmp_q[1] * 21012271307L) + ((int128)tmp_q[2] * 38184859812L) + ((int128)tmp_q[3] * 49128047966L) - ((int128)tmp_q[4] * 3257067862L) + ((int128)tmp_q[5] * 35243808746L) + ((int128)tmp_q[6] * 34474813128L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

