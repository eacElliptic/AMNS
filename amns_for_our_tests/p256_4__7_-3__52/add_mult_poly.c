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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1134860617448312340UL) + ((((uint64_t)op[1] * 5988568102294077950UL) + ((uint64_t)op[2] * 9976287446033776195UL) + ((uint64_t)op[3] * 12763119056777018575UL) + ((uint64_t)op[4] * 2525263425217575679UL) + ((uint64_t)op[5] * 6829217900803033760UL) + ((uint64_t)op[6] * 7454264475460903730UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 7454264475460903730UL) + ((uint64_t)op[1] * 1134860617448312340UL) + ((((uint64_t)op[2] * 5988568102294077950UL) + ((uint64_t)op[3] * 9976287446033776195UL) + ((uint64_t)op[4] * 12763119056777018575UL) + ((uint64_t)op[5] * 2525263425217575679UL) + ((uint64_t)op[6] * 6829217900803033760UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 6829217900803033760UL) + ((uint64_t)op[1] * 7454264475460903730UL) + ((uint64_t)op[2] * 1134860617448312340UL) + ((((uint64_t)op[3] * 5988568102294077950UL) + ((uint64_t)op[4] * 9976287446033776195UL) + ((uint64_t)op[5] * 12763119056777018575UL) + ((uint64_t)op[6] * 2525263425217575679UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 2525263425217575679UL) + ((uint64_t)op[1] * 6829217900803033760UL) + ((uint64_t)op[2] * 7454264475460903730UL) + ((uint64_t)op[3] * 1134860617448312340UL) + ((((uint64_t)op[4] * 5988568102294077950UL) + ((uint64_t)op[5] * 9976287446033776195UL) + ((uint64_t)op[6] * 12763119056777018575UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 12763119056777018575UL) + ((uint64_t)op[1] * 2525263425217575679UL) + ((uint64_t)op[2] * 6829217900803033760UL) + ((uint64_t)op[3] * 7454264475460903730UL) + ((uint64_t)op[4] * 1134860617448312340UL) + ((((uint64_t)op[5] * 5988568102294077950UL) + ((uint64_t)op[6] * 9976287446033776195UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 9976287446033776195UL) + ((uint64_t)op[1] * 12763119056777018575UL) + ((uint64_t)op[2] * 2525263425217575679UL) + ((uint64_t)op[3] * 6829217900803033760UL) + ((uint64_t)op[4] * 7454264475460903730UL) + ((uint64_t)op[5] * 1134860617448312340UL) + ((uint64_t)op[6] * 481039766827317766UL);
	tmp_q[6] = ((uint64_t)op[0] * 5988568102294077950UL) + ((uint64_t)op[1] * 9976287446033776195UL) + ((uint64_t)op[2] * 12763119056777018575UL) + ((uint64_t)op[3] * 2525263425217575679UL) + ((uint64_t)op[4] * 6829217900803033760UL) + ((uint64_t)op[5] * 7454264475460903730UL) + ((uint64_t)op[6] * 1134860617448312340UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 3795927721L) - ((((int128)tmp_q[1] * 78469127775L) + ((int128)tmp_q[2] * 13384671986L) + ((int128)tmp_q[3] * 63089044029L) + ((int128)tmp_q[4] * 27139852433L) - ((int128)tmp_q[5] * 3322943987L) + ((int128)tmp_q[6] * 38152387748L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 38152387748L) - ((int128)tmp_q[1] * 3795927721L) - ((((int128)tmp_q[2] * 78469127775L) + ((int128)tmp_q[3] * 13384671986L) + ((int128)tmp_q[4] * 63089044029L) + ((int128)tmp_q[5] * 27139852433L) - ((int128)tmp_q[6] * 3322943987L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 3322943987L) + ((int128)tmp_q[1] * 38152387748L) - ((int128)tmp_q[2] * 3795927721L) - ((((int128)tmp_q[3] * 78469127775L) + ((int128)tmp_q[4] * 13384671986L) + ((int128)tmp_q[5] * 63089044029L) + ((int128)tmp_q[6] * 27139852433L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 27139852433L) - ((int128)tmp_q[1] * 3322943987L) + ((int128)tmp_q[2] * 38152387748L) - ((int128)tmp_q[3] * 3795927721L) - ((((int128)tmp_q[4] * 78469127775L) + ((int128)tmp_q[5] * 13384671986L) + ((int128)tmp_q[6] * 63089044029L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 63089044029L) + ((int128)tmp_q[1] * 27139852433L) - ((int128)tmp_q[2] * 3322943987L) + ((int128)tmp_q[3] * 38152387748L) - ((int128)tmp_q[4] * 3795927721L) - ((((int128)tmp_q[5] * 78469127775L) + ((int128)tmp_q[6] * 13384671986L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 13384671986L) + ((int128)tmp_q[1] * 63089044029L) + ((int128)tmp_q[2] * 27139852433L) - ((int128)tmp_q[3] * 3322943987L) + ((int128)tmp_q[4] * 38152387748L) - ((int128)tmp_q[5] * 3795927721L) - ((int128)tmp_q[6] * 235407383325L);
	tmp_zero[6] = ((int128)tmp_q[0] * 78469127775L) + ((int128)tmp_q[1] * 13384671986L) + ((int128)tmp_q[2] * 63089044029L) + ((int128)tmp_q[3] * 27139852433L) - ((int128)tmp_q[4] * 3322943987L) + ((int128)tmp_q[5] * 38152387748L) - ((int128)tmp_q[6] * 3795927721L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

