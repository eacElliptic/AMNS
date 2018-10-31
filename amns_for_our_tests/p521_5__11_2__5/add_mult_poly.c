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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) << 1);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) << 2);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) << 1);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 16154456111883867305UL) + ((((uint64_t)op[1] * 8869479548860326006UL) + ((uint64_t)op[2] * 5119788002511585082UL) + ((uint64_t)op[3] * 7346190131884729535UL) + ((uint64_t)op[4] * 7534914907632187535UL) + ((uint64_t)op[5] * 16979964985818080258UL) + ((uint64_t)op[6] * 12948633822410685167UL) + ((uint64_t)op[7] * 12760992487808256779UL) + ((uint64_t)op[8] * 1336353340633487607UL) + ((uint64_t)op[9] * 8560930429279022193UL) + ((uint64_t)op[10] * 12315049366561132540UL)) * 2);
	tmp_q[1] = ((uint64_t)op[0] * 12315049366561132540UL) + ((uint64_t)op[1] * 16154456111883867305UL) + ((((uint64_t)op[2] * 8869479548860326006UL) + ((uint64_t)op[3] * 5119788002511585082UL) + ((uint64_t)op[4] * 7346190131884729535UL) + ((uint64_t)op[5] * 7534914907632187535UL) + ((uint64_t)op[6] * 16979964985818080258UL) + ((uint64_t)op[7] * 12948633822410685167UL) + ((uint64_t)op[8] * 12760992487808256779UL) + ((uint64_t)op[9] * 1336353340633487607UL) + ((uint64_t)op[10] * 8560930429279022193UL)) * 2);
	tmp_q[2] = ((uint64_t)op[0] * 8560930429279022193UL) + ((uint64_t)op[1] * 12315049366561132540UL) + ((uint64_t)op[2] * 16154456111883867305UL) + ((((uint64_t)op[3] * 8869479548860326006UL) + ((uint64_t)op[4] * 5119788002511585082UL) + ((uint64_t)op[5] * 7346190131884729535UL) + ((uint64_t)op[6] * 7534914907632187535UL) + ((uint64_t)op[7] * 16979964985818080258UL) + ((uint64_t)op[8] * 12948633822410685167UL) + ((uint64_t)op[9] * 12760992487808256779UL) + ((uint64_t)op[10] * 1336353340633487607UL)) * 2);
	tmp_q[3] = ((uint64_t)op[0] * 1336353340633487607UL) + ((uint64_t)op[1] * 8560930429279022193UL) + ((uint64_t)op[2] * 12315049366561132540UL) + ((uint64_t)op[3] * 16154456111883867305UL) + ((((uint64_t)op[4] * 8869479548860326006UL) + ((uint64_t)op[5] * 5119788002511585082UL) + ((uint64_t)op[6] * 7346190131884729535UL) + ((uint64_t)op[7] * 7534914907632187535UL) + ((uint64_t)op[8] * 16979964985818080258UL) + ((uint64_t)op[9] * 12948633822410685167UL) + ((uint64_t)op[10] * 12760992487808256779UL)) * 2);
	tmp_q[4] = ((uint64_t)op[0] * 12760992487808256779UL) + ((uint64_t)op[1] * 1336353340633487607UL) + ((uint64_t)op[2] * 8560930429279022193UL) + ((uint64_t)op[3] * 12315049366561132540UL) + ((uint64_t)op[4] * 16154456111883867305UL) + ((((uint64_t)op[5] * 8869479548860326006UL) + ((uint64_t)op[6] * 5119788002511585082UL) + ((uint64_t)op[7] * 7346190131884729535UL) + ((uint64_t)op[8] * 7534914907632187535UL) + ((uint64_t)op[9] * 16979964985818080258UL) + ((uint64_t)op[10] * 12948633822410685167UL)) * 2);
	tmp_q[5] = ((uint64_t)op[0] * 12948633822410685167UL) + ((uint64_t)op[1] * 12760992487808256779UL) + ((uint64_t)op[2] * 1336353340633487607UL) + ((uint64_t)op[3] * 8560930429279022193UL) + ((uint64_t)op[4] * 12315049366561132540UL) + ((uint64_t)op[5] * 16154456111883867305UL) + ((((uint64_t)op[6] * 8869479548860326006UL) + ((uint64_t)op[7] * 5119788002511585082UL) + ((uint64_t)op[8] * 7346190131884729535UL) + ((uint64_t)op[9] * 7534914907632187535UL) + ((uint64_t)op[10] * 16979964985818080258UL)) * 2);
	tmp_q[6] = ((uint64_t)op[0] * 16979964985818080258UL) + ((uint64_t)op[1] * 12948633822410685167UL) + ((uint64_t)op[2] * 12760992487808256779UL) + ((uint64_t)op[3] * 1336353340633487607UL) + ((uint64_t)op[4] * 8560930429279022193UL) + ((uint64_t)op[5] * 12315049366561132540UL) + ((uint64_t)op[6] * 16154456111883867305UL) + ((((uint64_t)op[7] * 8869479548860326006UL) + ((uint64_t)op[8] * 5119788002511585082UL) + ((uint64_t)op[9] * 7346190131884729535UL) + ((uint64_t)op[10] * 7534914907632187535UL)) * 2);
	tmp_q[7] = ((uint64_t)op[0] * 7534914907632187535UL) + ((uint64_t)op[1] * 16979964985818080258UL) + ((uint64_t)op[2] * 12948633822410685167UL) + ((uint64_t)op[3] * 12760992487808256779UL) + ((uint64_t)op[4] * 1336353340633487607UL) + ((uint64_t)op[5] * 8560930429279022193UL) + ((uint64_t)op[6] * 12315049366561132540UL) + ((uint64_t)op[7] * 16154456111883867305UL) + ((((uint64_t)op[8] * 8869479548860326006UL) + ((uint64_t)op[9] * 5119788002511585082UL) + ((uint64_t)op[10] * 7346190131884729535UL)) * 2);
	tmp_q[8] = ((uint64_t)op[0] * 7346190131884729535UL) + ((uint64_t)op[1] * 7534914907632187535UL) + ((uint64_t)op[2] * 16979964985818080258UL) + ((uint64_t)op[3] * 12948633822410685167UL) + ((uint64_t)op[4] * 12760992487808256779UL) + ((uint64_t)op[5] * 1336353340633487607UL) + ((uint64_t)op[6] * 8560930429279022193UL) + ((uint64_t)op[7] * 12315049366561132540UL) + ((uint64_t)op[8] * 16154456111883867305UL) + ((((uint64_t)op[9] * 8869479548860326006UL) + ((uint64_t)op[10] * 5119788002511585082UL)) * 2);
	tmp_q[9] = ((uint64_t)op[0] * 5119788002511585082UL) + ((uint64_t)op[1] * 7346190131884729535UL) + ((uint64_t)op[2] * 7534914907632187535UL) + ((uint64_t)op[3] * 16979964985818080258UL) + ((uint64_t)op[4] * 12948633822410685167UL) + ((uint64_t)op[5] * 12760992487808256779UL) + ((uint64_t)op[6] * 1336353340633487607UL) + ((uint64_t)op[7] * 8560930429279022193UL) + ((uint64_t)op[8] * 12315049366561132540UL) + ((uint64_t)op[9] * 16154456111883867305UL) + ((uint64_t)op[10] * 17738959097720652012UL);
	tmp_q[10] = ((uint64_t)op[0] * 8869479548860326006UL) + ((uint64_t)op[1] * 5119788002511585082UL) + ((uint64_t)op[2] * 7346190131884729535UL) + ((uint64_t)op[3] * 7534914907632187535UL) + ((uint64_t)op[4] * 16979964985818080258UL) + ((uint64_t)op[5] * 12948633822410685167UL) + ((uint64_t)op[6] * 12760992487808256779UL) + ((uint64_t)op[7] * 1336353340633487607UL) + ((uint64_t)op[8] * 8560930429279022193UL) + ((uint64_t)op[9] * 12315049366561132540UL) + ((uint64_t)op[10] * 16154456111883867305UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 4988315393927L) + ((((int128)tmp_q[1] * 8983251195634L) - ((int128)tmp_q[2] * 78450187215630L) - ((int128)tmp_q[3] * 58388975464995L) - ((int128)tmp_q[4] * 84836772043500L) - ((int128)tmp_q[5] * 60566330150110L) + ((int128)tmp_q[6] * 36050225585023L) - ((int128)tmp_q[7] * 31575692979972L) + ((int128)tmp_q[8] * 26452640604723L) + ((int128)tmp_q[9] * 80674299699329L) + ((int128)tmp_q[10] * 6647200163020L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 6647200163020L) + ((int128)tmp_q[1] * 4988315393927L) + ((((int128)tmp_q[2] * 8983251195634L) - ((int128)tmp_q[3] * 78450187215630L) - ((int128)tmp_q[4] * 58388975464995L) - ((int128)tmp_q[5] * 84836772043500L) - ((int128)tmp_q[6] * 60566330150110L) + ((int128)tmp_q[7] * 36050225585023L) - ((int128)tmp_q[8] * 31575692979972L) + ((int128)tmp_q[9] * 26452640604723L) + ((int128)tmp_q[10] * 80674299699329L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 80674299699329L) + ((int128)tmp_q[1] * 6647200163020L) + ((int128)tmp_q[2] * 4988315393927L) + ((((int128)tmp_q[3] * 8983251195634L) - ((int128)tmp_q[4] * 78450187215630L) - ((int128)tmp_q[5] * 58388975464995L) - ((int128)tmp_q[6] * 84836772043500L) - ((int128)tmp_q[7] * 60566330150110L) + ((int128)tmp_q[8] * 36050225585023L) - ((int128)tmp_q[9] * 31575692979972L) + ((int128)tmp_q[10] * 26452640604723L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 26452640604723L) + ((int128)tmp_q[1] * 80674299699329L) + ((int128)tmp_q[2] * 6647200163020L) + ((int128)tmp_q[3] * 4988315393927L) + ((((int128)tmp_q[4] * 8983251195634L) - ((int128)tmp_q[5] * 78450187215630L) - ((int128)tmp_q[6] * 58388975464995L) - ((int128)tmp_q[7] * 84836772043500L) - ((int128)tmp_q[8] * 60566330150110L) + ((int128)tmp_q[9] * 36050225585023L) - ((int128)tmp_q[10] * 31575692979972L)) * 2);
	tmp_zero[4] = -((int128)tmp_q[0] * 31575692979972L) + ((int128)tmp_q[1] * 26452640604723L) + ((int128)tmp_q[2] * 80674299699329L) + ((int128)tmp_q[3] * 6647200163020L) + ((int128)tmp_q[4] * 4988315393927L) + ((((int128)tmp_q[5] * 8983251195634L) - ((int128)tmp_q[6] * 78450187215630L) - ((int128)tmp_q[7] * 58388975464995L) - ((int128)tmp_q[8] * 84836772043500L) - ((int128)tmp_q[9] * 60566330150110L) + ((int128)tmp_q[10] * 36050225585023L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 36050225585023L) - ((int128)tmp_q[1] * 31575692979972L) + ((int128)tmp_q[2] * 26452640604723L) + ((int128)tmp_q[3] * 80674299699329L) + ((int128)tmp_q[4] * 6647200163020L) + ((int128)tmp_q[5] * 4988315393927L) + ((((int128)tmp_q[6] * 8983251195634L) - ((int128)tmp_q[7] * 78450187215630L) - ((int128)tmp_q[8] * 58388975464995L) - ((int128)tmp_q[9] * 84836772043500L) - ((int128)tmp_q[10] * 60566330150110L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 60566330150110L) + ((int128)tmp_q[1] * 36050225585023L) - ((int128)tmp_q[2] * 31575692979972L) + ((int128)tmp_q[3] * 26452640604723L) + ((int128)tmp_q[4] * 80674299699329L) + ((int128)tmp_q[5] * 6647200163020L) + ((int128)tmp_q[6] * 4988315393927L) + ((((int128)tmp_q[7] * 8983251195634L) - ((int128)tmp_q[8] * 78450187215630L) - ((int128)tmp_q[9] * 58388975464995L) - ((int128)tmp_q[10] * 84836772043500L)) * 2);
	tmp_zero[7] = -((int128)tmp_q[0] * 84836772043500L) - ((int128)tmp_q[1] * 60566330150110L) + ((int128)tmp_q[2] * 36050225585023L) - ((int128)tmp_q[3] * 31575692979972L) + ((int128)tmp_q[4] * 26452640604723L) + ((int128)tmp_q[5] * 80674299699329L) + ((int128)tmp_q[6] * 6647200163020L) + ((int128)tmp_q[7] * 4988315393927L) + ((((int128)tmp_q[8] * 8983251195634L) - ((int128)tmp_q[9] * 78450187215630L) - ((int128)tmp_q[10] * 58388975464995L)) * 2);
	tmp_zero[8] = -((int128)tmp_q[0] * 58388975464995L) - ((int128)tmp_q[1] * 84836772043500L) - ((int128)tmp_q[2] * 60566330150110L) + ((int128)tmp_q[3] * 36050225585023L) - ((int128)tmp_q[4] * 31575692979972L) + ((int128)tmp_q[5] * 26452640604723L) + ((int128)tmp_q[6] * 80674299699329L) + ((int128)tmp_q[7] * 6647200163020L) + ((int128)tmp_q[8] * 4988315393927L) + ((((int128)tmp_q[9] * 8983251195634L) - ((int128)tmp_q[10] * 78450187215630L)) * 2);
	tmp_zero[9] = -((int128)tmp_q[0] * 78450187215630L) - ((int128)tmp_q[1] * 58388975464995L) - ((int128)tmp_q[2] * 84836772043500L) - ((int128)tmp_q[3] * 60566330150110L) + ((int128)tmp_q[4] * 36050225585023L) - ((int128)tmp_q[5] * 31575692979972L) + ((int128)tmp_q[6] * 26452640604723L) + ((int128)tmp_q[7] * 80674299699329L) + ((int128)tmp_q[8] * 6647200163020L) + ((int128)tmp_q[9] * 4988315393927L) + ((int128)tmp_q[10] * 17966502391268L);
	tmp_zero[10] = ((int128)tmp_q[0] * 8983251195634L) - ((int128)tmp_q[1] * 78450187215630L) - ((int128)tmp_q[2] * 58388975464995L) - ((int128)tmp_q[3] * 84836772043500L) - ((int128)tmp_q[4] * 60566330150110L) + ((int128)tmp_q[5] * 36050225585023L) - ((int128)tmp_q[6] * 31575692979972L) + ((int128)tmp_q[7] * 26452640604723L) + ((int128)tmp_q[8] * 80674299699329L) + ((int128)tmp_q[9] * 6647200163020L) + ((int128)tmp_q[10] * 4988315393927L);

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
}

