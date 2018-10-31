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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[7]) * 5);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 10);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 10);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((int128)pa[7] * pa[6]) * 10);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[7] * pa[7]) * 5);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 4329446415213474737UL) + ((((uint64_t)op[1] * 7792448788223422054UL) + ((uint64_t)op[2] * 7426712856774093435UL) + ((uint64_t)op[3] * 11433975928305562951UL) + ((uint64_t)op[4] * 8613079278975451153UL) + ((uint64_t)op[5] * 6217949247535709052UL) + ((uint64_t)op[6] * 558223588085038633UL) + ((uint64_t)op[7] * 6208043528424576184UL)) * 5);
	tmp_q[1] = ((uint64_t)op[0] * 6208043528424576184UL) + ((uint64_t)op[1] * 4329446415213474737UL) + ((((uint64_t)op[2] * 7792448788223422054UL) + ((uint64_t)op[3] * 7426712856774093435UL) + ((uint64_t)op[4] * 11433975928305562951UL) + ((uint64_t)op[5] * 8613079278975451153UL) + ((uint64_t)op[6] * 6217949247535709052UL) + ((uint64_t)op[7] * 558223588085038633UL)) * 5);
	tmp_q[2] = ((uint64_t)op[0] * 558223588085038633UL) + ((uint64_t)op[1] * 6208043528424576184UL) + ((uint64_t)op[2] * 4329446415213474737UL) + ((((uint64_t)op[3] * 7792448788223422054UL) + ((uint64_t)op[4] * 7426712856774093435UL) + ((uint64_t)op[5] * 11433975928305562951UL) + ((uint64_t)op[6] * 8613079278975451153UL) + ((uint64_t)op[7] * 6217949247535709052UL)) * 5);
	tmp_q[3] = ((uint64_t)op[0] * 6217949247535709052UL) + ((uint64_t)op[1] * 558223588085038633UL) + ((uint64_t)op[2] * 6208043528424576184UL) + ((uint64_t)op[3] * 4329446415213474737UL) + ((((uint64_t)op[4] * 7792448788223422054UL) + ((uint64_t)op[5] * 7426712856774093435UL) + ((uint64_t)op[6] * 11433975928305562951UL) + ((uint64_t)op[7] * 8613079278975451153UL)) * 5);
	tmp_q[4] = ((uint64_t)op[0] * 8613079278975451153UL) + ((uint64_t)op[1] * 6217949247535709052UL) + ((uint64_t)op[2] * 558223588085038633UL) + ((uint64_t)op[3] * 6208043528424576184UL) + ((uint64_t)op[4] * 4329446415213474737UL) + ((((uint64_t)op[5] * 7792448788223422054UL) + ((uint64_t)op[6] * 7426712856774093435UL) + ((uint64_t)op[7] * 11433975928305562951UL)) * 5);
	tmp_q[5] = ((uint64_t)op[0] * 11433975928305562951UL) + ((uint64_t)op[1] * 8613079278975451153UL) + ((uint64_t)op[2] * 6217949247535709052UL) + ((uint64_t)op[3] * 558223588085038633UL) + ((uint64_t)op[4] * 6208043528424576184UL) + ((uint64_t)op[5] * 4329446415213474737UL) + ((((uint64_t)op[6] * 7792448788223422054UL) + ((uint64_t)op[7] * 7426712856774093435UL)) * 5);
	tmp_q[6] = ((uint64_t)op[0] * 7426712856774093435UL) + ((uint64_t)op[1] * 11433975928305562951UL) + ((uint64_t)op[2] * 8613079278975451153UL) + ((uint64_t)op[3] * 6217949247535709052UL) + ((uint64_t)op[4] * 558223588085038633UL) + ((uint64_t)op[5] * 6208043528424576184UL) + ((uint64_t)op[6] * 4329446415213474737UL) + ((uint64_t)op[7] * 2068755793698007038UL);
	tmp_q[7] = ((uint64_t)op[0] * 7792448788223422054UL) + ((uint64_t)op[1] * 7426712856774093435UL) + ((uint64_t)op[2] * 11433975928305562951UL) + ((uint64_t)op[3] * 8613079278975451153UL) + ((uint64_t)op[4] * 6217949247535709052UL) + ((uint64_t)op[5] * 558223588085038633UL) + ((uint64_t)op[6] * 6208043528424576184UL) + ((uint64_t)op[7] * 4329446415213474737UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 83707349482823L) + ((-((int128)tmp_q[1] * 110795289161852L) - ((int128)tmp_q[2] * 30413452675667L) - ((int128)tmp_q[3] * 90868050250066L) - ((int128)tmp_q[4] * 26597126366533L) + ((int128)tmp_q[5] * 28324570043613L) - ((int128)tmp_q[6] * 43557981862783L) + ((int128)tmp_q[7] * 76240746690520L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 76240746690520L) - ((int128)tmp_q[1] * 83707349482823L) + ((-((int128)tmp_q[2] * 110795289161852L) - ((int128)tmp_q[3] * 30413452675667L) - ((int128)tmp_q[4] * 90868050250066L) - ((int128)tmp_q[5] * 26597126366533L) + ((int128)tmp_q[6] * 28324570043613L) - ((int128)tmp_q[7] * 43557981862783L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 43557981862783L) + ((int128)tmp_q[1] * 76240746690520L) - ((int128)tmp_q[2] * 83707349482823L) + ((-((int128)tmp_q[3] * 110795289161852L) - ((int128)tmp_q[4] * 30413452675667L) - ((int128)tmp_q[5] * 90868050250066L) - ((int128)tmp_q[6] * 26597126366533L) + ((int128)tmp_q[7] * 28324570043613L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 28324570043613L) - ((int128)tmp_q[1] * 43557981862783L) + ((int128)tmp_q[2] * 76240746690520L) - ((int128)tmp_q[3] * 83707349482823L) + ((-((int128)tmp_q[4] * 110795289161852L) - ((int128)tmp_q[5] * 30413452675667L) - ((int128)tmp_q[6] * 90868050250066L) - ((int128)tmp_q[7] * 26597126366533L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 26597126366533L) + ((int128)tmp_q[1] * 28324570043613L) - ((int128)tmp_q[2] * 43557981862783L) + ((int128)tmp_q[3] * 76240746690520L) - ((int128)tmp_q[4] * 83707349482823L) + ((-((int128)tmp_q[5] * 110795289161852L) - ((int128)tmp_q[6] * 30413452675667L) - ((int128)tmp_q[7] * 90868050250066L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 90868050250066L) - ((int128)tmp_q[1] * 26597126366533L) + ((int128)tmp_q[2] * 28324570043613L) - ((int128)tmp_q[3] * 43557981862783L) + ((int128)tmp_q[4] * 76240746690520L) - ((int128)tmp_q[5] * 83707349482823L) + ((-((int128)tmp_q[6] * 110795289161852L) - ((int128)tmp_q[7] * 30413452675667L)) * 5);
	tmp_zero[6] = -((int128)tmp_q[0] * 30413452675667L) - ((int128)tmp_q[1] * 90868050250066L) - ((int128)tmp_q[2] * 26597126366533L) + ((int128)tmp_q[3] * 28324570043613L) - ((int128)tmp_q[4] * 43557981862783L) + ((int128)tmp_q[5] * 76240746690520L) - ((int128)tmp_q[6] * 83707349482823L) - ((int128)tmp_q[7] * 553976445809260L);
	tmp_zero[7] = -((int128)tmp_q[0] * 110795289161852L) - ((int128)tmp_q[1] * 30413452675667L) - ((int128)tmp_q[2] * 90868050250066L) - ((int128)tmp_q[3] * 26597126366533L) + ((int128)tmp_q[4] * 28324570043613L) - ((int128)tmp_q[5] * 43557981862783L) + ((int128)tmp_q[6] * 76240746690520L) - ((int128)tmp_q[7] * 83707349482823L);

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

