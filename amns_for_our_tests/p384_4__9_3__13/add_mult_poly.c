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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[8] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 12218504154860753631UL) + ((((uint64_t)op[1] * 11770817965669754942UL) + ((uint64_t)op[2] * 15323887648432681155UL) + ((uint64_t)op[3] * 13477609993288347499UL) + ((uint64_t)op[4] * 10536615971490207119UL) + ((uint64_t)op[5] * 2285939478845294911UL) + ((uint64_t)op[6] * 17953314295429502776UL) + ((uint64_t)op[7] * 5274015244445443323UL) + ((uint64_t)op[8] * 14027813177956902485UL)) * 3);
	tmp_q[1] = ((uint64_t)op[0] * 14027813177956902485UL) + ((uint64_t)op[1] * 12218504154860753631UL) + ((((uint64_t)op[2] * 11770817965669754942UL) + ((uint64_t)op[3] * 15323887648432681155UL) + ((uint64_t)op[4] * 13477609993288347499UL) + ((uint64_t)op[5] * 10536615971490207119UL) + ((uint64_t)op[6] * 2285939478845294911UL) + ((uint64_t)op[7] * 17953314295429502776UL) + ((uint64_t)op[8] * 5274015244445443323UL)) * 3);
	tmp_q[2] = ((uint64_t)op[0] * 5274015244445443323UL) + ((uint64_t)op[1] * 14027813177956902485UL) + ((uint64_t)op[2] * 12218504154860753631UL) + ((((uint64_t)op[3] * 11770817965669754942UL) + ((uint64_t)op[4] * 15323887648432681155UL) + ((uint64_t)op[5] * 13477609993288347499UL) + ((uint64_t)op[6] * 10536615971490207119UL) + ((uint64_t)op[7] * 2285939478845294911UL) + ((uint64_t)op[8] * 17953314295429502776UL)) * 3);
	tmp_q[3] = ((uint64_t)op[0] * 17953314295429502776UL) + ((uint64_t)op[1] * 5274015244445443323UL) + ((uint64_t)op[2] * 14027813177956902485UL) + ((uint64_t)op[3] * 12218504154860753631UL) + ((((uint64_t)op[4] * 11770817965669754942UL) + ((uint64_t)op[5] * 15323887648432681155UL) + ((uint64_t)op[6] * 13477609993288347499UL) + ((uint64_t)op[7] * 10536615971490207119UL) + ((uint64_t)op[8] * 2285939478845294911UL)) * 3);
	tmp_q[4] = ((uint64_t)op[0] * 2285939478845294911UL) + ((uint64_t)op[1] * 17953314295429502776UL) + ((uint64_t)op[2] * 5274015244445443323UL) + ((uint64_t)op[3] * 14027813177956902485UL) + ((uint64_t)op[4] * 12218504154860753631UL) + ((((uint64_t)op[5] * 11770817965669754942UL) + ((uint64_t)op[6] * 15323887648432681155UL) + ((uint64_t)op[7] * 13477609993288347499UL) + ((uint64_t)op[8] * 10536615971490207119UL)) * 3);
	tmp_q[5] = ((uint64_t)op[0] * 10536615971490207119UL) + ((uint64_t)op[1] * 2285939478845294911UL) + ((uint64_t)op[2] * 17953314295429502776UL) + ((uint64_t)op[3] * 5274015244445443323UL) + ((uint64_t)op[4] * 14027813177956902485UL) + ((uint64_t)op[5] * 12218504154860753631UL) + ((((uint64_t)op[6] * 11770817965669754942UL) + ((uint64_t)op[7] * 15323887648432681155UL) + ((uint64_t)op[8] * 13477609993288347499UL)) * 3);
	tmp_q[6] = ((uint64_t)op[0] * 13477609993288347499UL) + ((uint64_t)op[1] * 10536615971490207119UL) + ((uint64_t)op[2] * 2285939478845294911UL) + ((uint64_t)op[3] * 17953314295429502776UL) + ((uint64_t)op[4] * 5274015244445443323UL) + ((uint64_t)op[5] * 14027813177956902485UL) + ((uint64_t)op[6] * 12218504154860753631UL) + ((((uint64_t)op[7] * 11770817965669754942UL) + ((uint64_t)op[8] * 15323887648432681155UL)) * 3);
	tmp_q[7] = ((uint64_t)op[0] * 15323887648432681155UL) + ((uint64_t)op[1] * 13477609993288347499UL) + ((uint64_t)op[2] * 10536615971490207119UL) + ((uint64_t)op[3] * 2285939478845294911UL) + ((uint64_t)op[4] * 17953314295429502776UL) + ((uint64_t)op[5] * 5274015244445443323UL) + ((uint64_t)op[6] * 14027813177956902485UL) + ((uint64_t)op[7] * 12218504154860753631UL) + ((uint64_t)op[8] * 16865709823299713210UL);
	tmp_q[8] = ((uint64_t)op[0] * 11770817965669754942UL) + ((uint64_t)op[1] * 15323887648432681155UL) + ((uint64_t)op[2] * 13477609993288347499UL) + ((uint64_t)op[3] * 10536615971490207119UL) + ((uint64_t)op[4] * 2285939478845294911UL) + ((uint64_t)op[5] * 17953314295429502776UL) + ((uint64_t)op[6] * 5274015244445443323UL) + ((uint64_t)op[7] * 14027813177956902485UL) + ((uint64_t)op[8] * 12218504154860753631UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1940590764383L) + ((-((int128)tmp_q[1] * 425515553317L) - ((int128)tmp_q[2] * 2670424997107L) + ((int128)tmp_q[3] * 4734117166657L) - ((int128)tmp_q[4] * 2834986306724L) + ((int128)tmp_q[5] * 2360485043362L) - ((int128)tmp_q[6] * 2676637074176L) + ((int128)tmp_q[7] * 2145975177802L) - ((int128)tmp_q[8] * 3711865599913L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 3711865599913L) - ((int128)tmp_q[1] * 1940590764383L) + ((-((int128)tmp_q[2] * 425515553317L) - ((int128)tmp_q[3] * 2670424997107L) + ((int128)tmp_q[4] * 4734117166657L) - ((int128)tmp_q[5] * 2834986306724L) + ((int128)tmp_q[6] * 2360485043362L) - ((int128)tmp_q[7] * 2676637074176L) + ((int128)tmp_q[8] * 2145975177802L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 2145975177802L) - ((int128)tmp_q[1] * 3711865599913L) - ((int128)tmp_q[2] * 1940590764383L) + ((-((int128)tmp_q[3] * 425515553317L) - ((int128)tmp_q[4] * 2670424997107L) + ((int128)tmp_q[5] * 4734117166657L) - ((int128)tmp_q[6] * 2834986306724L) + ((int128)tmp_q[7] * 2360485043362L) - ((int128)tmp_q[8] * 2676637074176L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 2676637074176L) + ((int128)tmp_q[1] * 2145975177802L) - ((int128)tmp_q[2] * 3711865599913L) - ((int128)tmp_q[3] * 1940590764383L) + ((-((int128)tmp_q[4] * 425515553317L) - ((int128)tmp_q[5] * 2670424997107L) + ((int128)tmp_q[6] * 4734117166657L) - ((int128)tmp_q[7] * 2834986306724L) + ((int128)tmp_q[8] * 2360485043362L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 2360485043362L) - ((int128)tmp_q[1] * 2676637074176L) + ((int128)tmp_q[2] * 2145975177802L) - ((int128)tmp_q[3] * 3711865599913L) - ((int128)tmp_q[4] * 1940590764383L) + ((-((int128)tmp_q[5] * 425515553317L) - ((int128)tmp_q[6] * 2670424997107L) + ((int128)tmp_q[7] * 4734117166657L) - ((int128)tmp_q[8] * 2834986306724L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 2834986306724L) + ((int128)tmp_q[1] * 2360485043362L) - ((int128)tmp_q[2] * 2676637074176L) + ((int128)tmp_q[3] * 2145975177802L) - ((int128)tmp_q[4] * 3711865599913L) - ((int128)tmp_q[5] * 1940590764383L) + ((-((int128)tmp_q[6] * 425515553317L) - ((int128)tmp_q[7] * 2670424997107L) + ((int128)tmp_q[8] * 4734117166657L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 4734117166657L) - ((int128)tmp_q[1] * 2834986306724L) + ((int128)tmp_q[2] * 2360485043362L) - ((int128)tmp_q[3] * 2676637074176L) + ((int128)tmp_q[4] * 2145975177802L) - ((int128)tmp_q[5] * 3711865599913L) - ((int128)tmp_q[6] * 1940590764383L) + ((-((int128)tmp_q[7] * 425515553317L) - ((int128)tmp_q[8] * 2670424997107L)) * 3);
	tmp_zero[7] = -((int128)tmp_q[0] * 2670424997107L) + ((int128)tmp_q[1] * 4734117166657L) - ((int128)tmp_q[2] * 2834986306724L) + ((int128)tmp_q[3] * 2360485043362L) - ((int128)tmp_q[4] * 2676637074176L) + ((int128)tmp_q[5] * 2145975177802L) - ((int128)tmp_q[6] * 3711865599913L) - ((int128)tmp_q[7] * 1940590764383L) - ((int128)tmp_q[8] * 1276546659951L);
	tmp_zero[8] = -((int128)tmp_q[0] * 425515553317L) - ((int128)tmp_q[1] * 2670424997107L) + ((int128)tmp_q[2] * 4734117166657L) - ((int128)tmp_q[3] * 2834986306724L) + ((int128)tmp_q[4] * 2360485043362L) - ((int128)tmp_q[5] * 2676637074176L) + ((int128)tmp_q[6] * 2145975177802L) - ((int128)tmp_q[7] * 3711865599913L) - ((int128)tmp_q[8] * 1940590764383L);

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

