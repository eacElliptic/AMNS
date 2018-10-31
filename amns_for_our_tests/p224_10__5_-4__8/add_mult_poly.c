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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[4] + (int128)pa[4] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1]) << 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[4] * pa[2]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[4] * pa[3]) << 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 313946906954887821UL) + ((((uint64_t)op[1] * 11009783368483206963UL) + ((uint64_t)op[2] * 16303965199559541060UL) + ((uint64_t)op[3] * 9851549103627281917UL) + ((uint64_t)op[4] * 11580205357454119140UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 11580205357454119140UL) + ((uint64_t)op[1] * 313946906954887821UL) + ((((uint64_t)op[2] * 11009783368483206963UL) + ((uint64_t)op[3] * 16303965199559541060UL) + ((uint64_t)op[4] * 9851549103627281917UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 9851549103627281917UL) + ((uint64_t)op[1] * 11580205357454119140UL) + ((uint64_t)op[2] * 313946906954887821UL) + ((((uint64_t)op[3] * 11009783368483206963UL) + ((uint64_t)op[4] * 16303965199559541060UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 16303965199559541060UL) + ((uint64_t)op[1] * 9851549103627281917UL) + ((uint64_t)op[2] * 11580205357454119140UL) + ((uint64_t)op[3] * 313946906954887821UL) + ((uint64_t)op[4] * 11301098747195826996UL);
	tmp_q[4] = ((uint64_t)op[0] * 11009783368483206963UL) + ((uint64_t)op[1] * 16303965199559541060UL) + ((uint64_t)op[2] * 9851549103627281917UL) + ((uint64_t)op[3] * 11580205357454119140UL) + ((uint64_t)op[4] * 313946906954887821UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 4498584260341L) - ((-((int128)tmp_q[1] * 17249999610090L) - ((int128)tmp_q[2] * 10504655518080L) - ((int128)tmp_q[3] * 1484052917211L) - ((int128)tmp_q[4] * 8282568629432L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 8282568629432L) - ((int128)tmp_q[1] * 4498584260341L) - ((-((int128)tmp_q[2] * 17249999610090L) - ((int128)tmp_q[3] * 10504655518080L) - ((int128)tmp_q[4] * 1484052917211L)) * 4);
	tmp_zero[2] = -((int128)tmp_q[0] * 1484052917211L) - ((int128)tmp_q[1] * 8282568629432L) - ((int128)tmp_q[2] * 4498584260341L) - ((-((int128)tmp_q[3] * 17249999610090L) - ((int128)tmp_q[4] * 10504655518080L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 10504655518080L) - ((int128)tmp_q[1] * 1484052917211L) - ((int128)tmp_q[2] * 8282568629432L) - ((int128)tmp_q[3] * 4498584260341L) + ((int128)tmp_q[4] * 68999998440360L);
	tmp_zero[4] = -((int128)tmp_q[0] * 17249999610090L) - ((int128)tmp_q[1] * 10504655518080L) - ((int128)tmp_q[2] * 1484052917211L) - ((int128)tmp_q[3] * 8282568629432L) - ((int128)tmp_q[4] * 4498584260341L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
}

