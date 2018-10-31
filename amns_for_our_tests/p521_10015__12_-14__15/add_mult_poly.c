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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[11] + (int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2] + (int128)pa[11] * pb[1]) * 14);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[11] + (int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3] + (int128)pa[11] * pb[2]) * 14);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[11] + (int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4] + (int128)pa[11] * pb[3]) * 14);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[11] + (int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5] + (int128)pa[11] * pb[4]) * 14);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[11] + (int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6] + (int128)pa[11] * pb[5]) * 14);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[11] + (int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7] + (int128)pa[11] * pb[6]) * 14);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[11] + (int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8] + (int128)pa[11] * pb[7]) * 14);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[11] + (int128)pa[9] * pb[10] + (int128)pa[10] * pb[9] + (int128)pa[11] * pb[8]) * 14);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[11] + (int128)pa[10] * pb[10] + (int128)pa[11] * pb[9]) * 14);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] - (((int128)pa[10] * pb[11] + (int128)pa[11] * pb[10]) * 14);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0] - (((int128)pa[11] * pb[11]) * 14);
	tmp_prod_result[11] = (int128)pa[0] * pb[11] + (int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1] + (int128)pa[11] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2] + (int128)pa[11] * pa[1]) << 1) + (int128)pa[6] * pa[6]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3] + (int128)pa[11] * pa[2]) * 28);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4] + (int128)pa[11] * pa[3]) << 1) + (int128)pa[7] * pa[7]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5] + (int128)pa[11] * pa[4]) * 28);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6] + (int128)pa[11] * pa[5]) << 1) + (int128)pa[8] * pa[8]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7] + (int128)pa[11] * pa[6]) * 28);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[10] * pa[8] + (int128)pa[11] * pa[7]) << 1) + (int128)pa[9] * pa[9]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[10] * pa[9] + (int128)pa[11] * pa[8]) * 28);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((((int128)pa[11] * pa[9]) << 1) + (int128)pa[10] * pa[10]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) - (((int128)pa[11] * pa[10]) * 28);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5] - (((int128)pa[11] * pa[11]) * 14);
	tmp_prod_result[11] = (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1] + (int128)pa[11] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10772573315394169047UL) + ((((uint64_t)op[1] * 14590325118758063687UL) + ((uint64_t)op[2] * 2933481310287690214UL) + ((uint64_t)op[3] * 6478547575863544192UL) + ((uint64_t)op[4] * 15427621959397184703UL) + ((uint64_t)op[5] * 2472737004397238086UL) + ((uint64_t)op[6] * 13353477544330242624UL) + ((uint64_t)op[7] * 16751716314857504272UL) + ((uint64_t)op[8] * 6704711288411505111UL) + ((uint64_t)op[9] * 960559382600009376UL) + ((uint64_t)op[10] * 13336336383532263213UL) + ((uint64_t)op[11] * 10630860590770986248UL)) * 18446744073709551602);
	tmp_q[1] = ((uint64_t)op[0] * 10630860590770986248UL) + ((uint64_t)op[1] * 10772573315394169047UL) + ((((uint64_t)op[2] * 14590325118758063687UL) + ((uint64_t)op[3] * 2933481310287690214UL) + ((uint64_t)op[4] * 6478547575863544192UL) + ((uint64_t)op[5] * 15427621959397184703UL) + ((uint64_t)op[6] * 2472737004397238086UL) + ((uint64_t)op[7] * 13353477544330242624UL) + ((uint64_t)op[8] * 16751716314857504272UL) + ((uint64_t)op[9] * 6704711288411505111UL) + ((uint64_t)op[10] * 960559382600009376UL) + ((uint64_t)op[11] * 13336336383532263213UL)) * 18446744073709551602);
	tmp_q[2] = ((uint64_t)op[0] * 13336336383532263213UL) + ((uint64_t)op[1] * 10630860590770986248UL) + ((uint64_t)op[2] * 10772573315394169047UL) + ((((uint64_t)op[3] * 14590325118758063687UL) + ((uint64_t)op[4] * 2933481310287690214UL) + ((uint64_t)op[5] * 6478547575863544192UL) + ((uint64_t)op[6] * 15427621959397184703UL) + ((uint64_t)op[7] * 2472737004397238086UL) + ((uint64_t)op[8] * 13353477544330242624UL) + ((uint64_t)op[9] * 16751716314857504272UL) + ((uint64_t)op[10] * 6704711288411505111UL) + ((uint64_t)op[11] * 960559382600009376UL)) * 18446744073709551602);
	tmp_q[3] = ((uint64_t)op[0] * 960559382600009376UL) + ((uint64_t)op[1] * 13336336383532263213UL) + ((uint64_t)op[2] * 10630860590770986248UL) + ((uint64_t)op[3] * 10772573315394169047UL) + ((((uint64_t)op[4] * 14590325118758063687UL) + ((uint64_t)op[5] * 2933481310287690214UL) + ((uint64_t)op[6] * 6478547575863544192UL) + ((uint64_t)op[7] * 15427621959397184703UL) + ((uint64_t)op[8] * 2472737004397238086UL) + ((uint64_t)op[9] * 13353477544330242624UL) + ((uint64_t)op[10] * 16751716314857504272UL) + ((uint64_t)op[11] * 6704711288411505111UL)) * 18446744073709551602);
	tmp_q[4] = ((uint64_t)op[0] * 6704711288411505111UL) + ((uint64_t)op[1] * 960559382600009376UL) + ((uint64_t)op[2] * 13336336383532263213UL) + ((uint64_t)op[3] * 10630860590770986248UL) + ((uint64_t)op[4] * 10772573315394169047UL) + ((((uint64_t)op[5] * 14590325118758063687UL) + ((uint64_t)op[6] * 2933481310287690214UL) + ((uint64_t)op[7] * 6478547575863544192UL) + ((uint64_t)op[8] * 15427621959397184703UL) + ((uint64_t)op[9] * 2472737004397238086UL) + ((uint64_t)op[10] * 13353477544330242624UL) + ((uint64_t)op[11] * 16751716314857504272UL)) * 18446744073709551602);
	tmp_q[5] = ((uint64_t)op[0] * 16751716314857504272UL) + ((uint64_t)op[1] * 6704711288411505111UL) + ((uint64_t)op[2] * 960559382600009376UL) + ((uint64_t)op[3] * 13336336383532263213UL) + ((uint64_t)op[4] * 10630860590770986248UL) + ((uint64_t)op[5] * 10772573315394169047UL) + ((((uint64_t)op[6] * 14590325118758063687UL) + ((uint64_t)op[7] * 2933481310287690214UL) + ((uint64_t)op[8] * 6478547575863544192UL) + ((uint64_t)op[9] * 15427621959397184703UL) + ((uint64_t)op[10] * 2472737004397238086UL) + ((uint64_t)op[11] * 13353477544330242624UL)) * 18446744073709551602);
	tmp_q[6] = ((uint64_t)op[0] * 13353477544330242624UL) + ((uint64_t)op[1] * 16751716314857504272UL) + ((uint64_t)op[2] * 6704711288411505111UL) + ((uint64_t)op[3] * 960559382600009376UL) + ((uint64_t)op[4] * 13336336383532263213UL) + ((uint64_t)op[5] * 10630860590770986248UL) + ((uint64_t)op[6] * 10772573315394169047UL) + ((((uint64_t)op[7] * 14590325118758063687UL) + ((uint64_t)op[8] * 2933481310287690214UL) + ((uint64_t)op[9] * 6478547575863544192UL) + ((uint64_t)op[10] * 15427621959397184703UL) + ((uint64_t)op[11] * 2472737004397238086UL)) * 18446744073709551602);
	tmp_q[7] = ((uint64_t)op[0] * 2472737004397238086UL) + ((uint64_t)op[1] * 13353477544330242624UL) + ((uint64_t)op[2] * 16751716314857504272UL) + ((uint64_t)op[3] * 6704711288411505111UL) + ((uint64_t)op[4] * 960559382600009376UL) + ((uint64_t)op[5] * 13336336383532263213UL) + ((uint64_t)op[6] * 10630860590770986248UL) + ((uint64_t)op[7] * 10772573315394169047UL) + ((((uint64_t)op[8] * 14590325118758063687UL) + ((uint64_t)op[9] * 2933481310287690214UL) + ((uint64_t)op[10] * 6478547575863544192UL) + ((uint64_t)op[11] * 15427621959397184703UL)) * 18446744073709551602);
	tmp_q[8] = ((uint64_t)op[0] * 15427621959397184703UL) + ((uint64_t)op[1] * 2472737004397238086UL) + ((uint64_t)op[2] * 13353477544330242624UL) + ((uint64_t)op[3] * 16751716314857504272UL) + ((uint64_t)op[4] * 6704711288411505111UL) + ((uint64_t)op[5] * 960559382600009376UL) + ((uint64_t)op[6] * 13336336383532263213UL) + ((uint64_t)op[7] * 10630860590770986248UL) + ((uint64_t)op[8] * 10772573315394169047UL) + ((((uint64_t)op[9] * 14590325118758063687UL) + ((uint64_t)op[10] * 2933481310287690214UL) + ((uint64_t)op[11] * 6478547575863544192UL)) * 18446744073709551602);
	tmp_q[9] = ((uint64_t)op[0] * 6478547575863544192UL) + ((uint64_t)op[1] * 15427621959397184703UL) + ((uint64_t)op[2] * 2472737004397238086UL) + ((uint64_t)op[3] * 13353477544330242624UL) + ((uint64_t)op[4] * 16751716314857504272UL) + ((uint64_t)op[5] * 6704711288411505111UL) + ((uint64_t)op[6] * 960559382600009376UL) + ((uint64_t)op[7] * 13336336383532263213UL) + ((uint64_t)op[8] * 10630860590770986248UL) + ((uint64_t)op[9] * 10772573315394169047UL) + ((((uint64_t)op[10] * 14590325118758063687UL) + ((uint64_t)op[11] * 2933481310287690214UL)) * 18446744073709551602);
	tmp_q[10] = ((uint64_t)op[0] * 2933481310287690214UL) + ((uint64_t)op[1] * 6478547575863544192UL) + ((uint64_t)op[2] * 15427621959397184703UL) + ((uint64_t)op[3] * 2472737004397238086UL) + ((uint64_t)op[4] * 13353477544330242624UL) + ((uint64_t)op[5] * 16751716314857504272UL) + ((uint64_t)op[6] * 6704711288411505111UL) + ((uint64_t)op[7] * 960559382600009376UL) + ((uint64_t)op[8] * 13336336383532263213UL) + ((uint64_t)op[9] * 10630860590770986248UL) + ((uint64_t)op[10] * 10772573315394169047UL) + ((uint64_t)op[11] * 17096377221901727774UL);
	tmp_q[11] = ((uint64_t)op[0] * 14590325118758063687UL) + ((uint64_t)op[1] * 2933481310287690214UL) + ((uint64_t)op[2] * 6478547575863544192UL) + ((uint64_t)op[3] * 15427621959397184703UL) + ((uint64_t)op[4] * 2472737004397238086UL) + ((uint64_t)op[5] * 13353477544330242624UL) + ((uint64_t)op[6] * 16751716314857504272UL) + ((uint64_t)op[7] * 6704711288411505111UL) + ((uint64_t)op[8] * 960559382600009376UL) + ((uint64_t)op[9] * 13336336383532263213UL) + ((uint64_t)op[10] * 10630860590770986248UL) + ((uint64_t)op[11] * 10772573315394169047UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 3041132233341L) - ((-((int128)tmp_q[1] * 4956728313805L) - ((int128)tmp_q[2] * 415353690310L) + ((int128)tmp_q[3] * 514907176984L) + ((int128)tmp_q[4] * 1794096871536L) - ((int128)tmp_q[5] * 913823888042L) - ((int128)tmp_q[6] * 2072316227513L) - ((int128)tmp_q[7] * 3857335043068L) - ((int128)tmp_q[8] * 594632556850L) - ((int128)tmp_q[9] * 236413275618L) - ((int128)tmp_q[10] * 6710765005821L) - ((int128)tmp_q[11] * 6435916234364L)) * 14);
	tmp_zero[1] = -((int128)tmp_q[0] * 6435916234364L) + ((int128)tmp_q[1] * 3041132233341L) - ((-((int128)tmp_q[2] * 4956728313805L) - ((int128)tmp_q[3] * 415353690310L) + ((int128)tmp_q[4] * 514907176984L) + ((int128)tmp_q[5] * 1794096871536L) - ((int128)tmp_q[6] * 913823888042L) - ((int128)tmp_q[7] * 2072316227513L) - ((int128)tmp_q[8] * 3857335043068L) - ((int128)tmp_q[9] * 594632556850L) - ((int128)tmp_q[10] * 236413275618L) - ((int128)tmp_q[11] * 6710765005821L)) * 14);
	tmp_zero[2] = -((int128)tmp_q[0] * 6710765005821L) - ((int128)tmp_q[1] * 6435916234364L) + ((int128)tmp_q[2] * 3041132233341L) - ((-((int128)tmp_q[3] * 4956728313805L) - ((int128)tmp_q[4] * 415353690310L) + ((int128)tmp_q[5] * 514907176984L) + ((int128)tmp_q[6] * 1794096871536L) - ((int128)tmp_q[7] * 913823888042L) - ((int128)tmp_q[8] * 2072316227513L) - ((int128)tmp_q[9] * 3857335043068L) - ((int128)tmp_q[10] * 594632556850L) - ((int128)tmp_q[11] * 236413275618L)) * 14);
	tmp_zero[3] = -((int128)tmp_q[0] * 236413275618L) - ((int128)tmp_q[1] * 6710765005821L) - ((int128)tmp_q[2] * 6435916234364L) + ((int128)tmp_q[3] * 3041132233341L) - ((-((int128)tmp_q[4] * 4956728313805L) - ((int128)tmp_q[5] * 415353690310L) + ((int128)tmp_q[6] * 514907176984L) + ((int128)tmp_q[7] * 1794096871536L) - ((int128)tmp_q[8] * 913823888042L) - ((int128)tmp_q[9] * 2072316227513L) - ((int128)tmp_q[10] * 3857335043068L) - ((int128)tmp_q[11] * 594632556850L)) * 14);
	tmp_zero[4] = -((int128)tmp_q[0] * 594632556850L) - ((int128)tmp_q[1] * 236413275618L) - ((int128)tmp_q[2] * 6710765005821L) - ((int128)tmp_q[3] * 6435916234364L) + ((int128)tmp_q[4] * 3041132233341L) - ((-((int128)tmp_q[5] * 4956728313805L) - ((int128)tmp_q[6] * 415353690310L) + ((int128)tmp_q[7] * 514907176984L) + ((int128)tmp_q[8] * 1794096871536L) - ((int128)tmp_q[9] * 913823888042L) - ((int128)tmp_q[10] * 2072316227513L) - ((int128)tmp_q[11] * 3857335043068L)) * 14);
	tmp_zero[5] = -((int128)tmp_q[0] * 3857335043068L) - ((int128)tmp_q[1] * 594632556850L) - ((int128)tmp_q[2] * 236413275618L) - ((int128)tmp_q[3] * 6710765005821L) - ((int128)tmp_q[4] * 6435916234364L) + ((int128)tmp_q[5] * 3041132233341L) - ((-((int128)tmp_q[6] * 4956728313805L) - ((int128)tmp_q[7] * 415353690310L) + ((int128)tmp_q[8] * 514907176984L) + ((int128)tmp_q[9] * 1794096871536L) - ((int128)tmp_q[10] * 913823888042L) - ((int128)tmp_q[11] * 2072316227513L)) * 14);
	tmp_zero[6] = -((int128)tmp_q[0] * 2072316227513L) - ((int128)tmp_q[1] * 3857335043068L) - ((int128)tmp_q[2] * 594632556850L) - ((int128)tmp_q[3] * 236413275618L) - ((int128)tmp_q[4] * 6710765005821L) - ((int128)tmp_q[5] * 6435916234364L) + ((int128)tmp_q[6] * 3041132233341L) - ((-((int128)tmp_q[7] * 4956728313805L) - ((int128)tmp_q[8] * 415353690310L) + ((int128)tmp_q[9] * 514907176984L) + ((int128)tmp_q[10] * 1794096871536L) - ((int128)tmp_q[11] * 913823888042L)) * 14);
	tmp_zero[7] = -((int128)tmp_q[0] * 913823888042L) - ((int128)tmp_q[1] * 2072316227513L) - ((int128)tmp_q[2] * 3857335043068L) - ((int128)tmp_q[3] * 594632556850L) - ((int128)tmp_q[4] * 236413275618L) - ((int128)tmp_q[5] * 6710765005821L) - ((int128)tmp_q[6] * 6435916234364L) + ((int128)tmp_q[7] * 3041132233341L) - ((-((int128)tmp_q[8] * 4956728313805L) - ((int128)tmp_q[9] * 415353690310L) + ((int128)tmp_q[10] * 514907176984L) + ((int128)tmp_q[11] * 1794096871536L)) * 14);
	tmp_zero[8] = ((int128)tmp_q[0] * 1794096871536L) - ((int128)tmp_q[1] * 913823888042L) - ((int128)tmp_q[2] * 2072316227513L) - ((int128)tmp_q[3] * 3857335043068L) - ((int128)tmp_q[4] * 594632556850L) - ((int128)tmp_q[5] * 236413275618L) - ((int128)tmp_q[6] * 6710765005821L) - ((int128)tmp_q[7] * 6435916234364L) + ((int128)tmp_q[8] * 3041132233341L) - ((-((int128)tmp_q[9] * 4956728313805L) - ((int128)tmp_q[10] * 415353690310L) + ((int128)tmp_q[11] * 514907176984L)) * 14);
	tmp_zero[9] = ((int128)tmp_q[0] * 514907176984L) + ((int128)tmp_q[1] * 1794096871536L) - ((int128)tmp_q[2] * 913823888042L) - ((int128)tmp_q[3] * 2072316227513L) - ((int128)tmp_q[4] * 3857335043068L) - ((int128)tmp_q[5] * 594632556850L) - ((int128)tmp_q[6] * 236413275618L) - ((int128)tmp_q[7] * 6710765005821L) - ((int128)tmp_q[8] * 6435916234364L) + ((int128)tmp_q[9] * 3041132233341L) - ((-((int128)tmp_q[10] * 4956728313805L) - ((int128)tmp_q[11] * 415353690310L)) * 14);
	tmp_zero[10] = -((int128)tmp_q[0] * 415353690310L) + ((int128)tmp_q[1] * 514907176984L) + ((int128)tmp_q[2] * 1794096871536L) - ((int128)tmp_q[3] * 913823888042L) - ((int128)tmp_q[4] * 2072316227513L) - ((int128)tmp_q[5] * 3857335043068L) - ((int128)tmp_q[6] * 594632556850L) - ((int128)tmp_q[7] * 236413275618L) - ((int128)tmp_q[8] * 6710765005821L) - ((int128)tmp_q[9] * 6435916234364L) + ((int128)tmp_q[10] * 3041132233341L) + ((int128)tmp_q[11] * 69394196393270L);
	tmp_zero[11] = -((int128)tmp_q[0] * 4956728313805L) - ((int128)tmp_q[1] * 415353690310L) + ((int128)tmp_q[2] * 514907176984L) + ((int128)tmp_q[3] * 1794096871536L) - ((int128)tmp_q[4] * 913823888042L) - ((int128)tmp_q[5] * 2072316227513L) - ((int128)tmp_q[6] * 3857335043068L) - ((int128)tmp_q[7] * 594632556850L) - ((int128)tmp_q[8] * 236413275618L) - ((int128)tmp_q[9] * 6710765005821L) - ((int128)tmp_q[10] * 6435916234364L) + ((int128)tmp_q[11] * 3041132233341L);

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

