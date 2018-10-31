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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[9] + (int128)pa[9] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[9]) << 1);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2]) << 2);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4]) << 2);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6]) << 2);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[9] * pa[7]) << 1) + (int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[9] * pa[8]) << 2);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[9] * pa[9]) << 1);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 10773746088520450151UL) + ((((uint64_t)op[1] * 4592429164023824795UL) + ((uint64_t)op[2] * 15721672951029712405UL) + ((uint64_t)op[3] * 16556497540425869029UL) + ((uint64_t)op[4] * 10437050973780786078UL) + ((uint64_t)op[5] * 17998220869296846599UL) + ((uint64_t)op[6] * 6717411263826142866UL) + ((uint64_t)op[7] * 13091019643316963401UL) + ((uint64_t)op[8] * 16725393099077756620UL) + ((uint64_t)op[9] * 9805895641677525383UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 9805895641677525383UL) + ((uint64_t)op[1] * 10773746088520450151UL) + ((((uint64_t)op[2] * 4592429164023824795UL) + ((uint64_t)op[3] * 15721672951029712405UL) + ((uint64_t)op[4] * 16556497540425869029UL) + ((uint64_t)op[5] * 10437050973780786078UL) + ((uint64_t)op[6] * 17998220869296846599UL) + ((uint64_t)op[7] * 6717411263826142866UL) + ((uint64_t)op[8] * 13091019643316963401UL) + ((uint64_t)op[9] * 16725393099077756620UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 16725393099077756620UL) + ((uint64_t)op[1] * 9805895641677525383UL) + ((uint64_t)op[2] * 10773746088520450151UL) + ((((uint64_t)op[3] * 4592429164023824795UL) + ((uint64_t)op[4] * 15721672951029712405UL) + ((uint64_t)op[5] * 16556497540425869029UL) + ((uint64_t)op[6] * 10437050973780786078UL) + ((uint64_t)op[7] * 17998220869296846599UL) + ((uint64_t)op[8] * 6717411263826142866UL) + ((uint64_t)op[9] * 13091019643316963401UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 13091019643316963401UL) + ((uint64_t)op[1] * 16725393099077756620UL) + ((uint64_t)op[2] * 9805895641677525383UL) + ((uint64_t)op[3] * 10773746088520450151UL) + ((((uint64_t)op[4] * 4592429164023824795UL) + ((uint64_t)op[5] * 15721672951029712405UL) + ((uint64_t)op[6] * 16556497540425869029UL) + ((uint64_t)op[7] * 10437050973780786078UL) + ((uint64_t)op[8] * 17998220869296846599UL) + ((uint64_t)op[9] * 6717411263826142866UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 6717411263826142866UL) + ((uint64_t)op[1] * 13091019643316963401UL) + ((uint64_t)op[2] * 16725393099077756620UL) + ((uint64_t)op[3] * 9805895641677525383UL) + ((uint64_t)op[4] * 10773746088520450151UL) + ((((uint64_t)op[5] * 4592429164023824795UL) + ((uint64_t)op[6] * 15721672951029712405UL) + ((uint64_t)op[7] * 16556497540425869029UL) + ((uint64_t)op[8] * 10437050973780786078UL) + ((uint64_t)op[9] * 17998220869296846599UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 17998220869296846599UL) + ((uint64_t)op[1] * 6717411263826142866UL) + ((uint64_t)op[2] * 13091019643316963401UL) + ((uint64_t)op[3] * 16725393099077756620UL) + ((uint64_t)op[4] * 9805895641677525383UL) + ((uint64_t)op[5] * 10773746088520450151UL) + ((((uint64_t)op[6] * 4592429164023824795UL) + ((uint64_t)op[7] * 15721672951029712405UL) + ((uint64_t)op[8] * 16556497540425869029UL) + ((uint64_t)op[9] * 10437050973780786078UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 10437050973780786078UL) + ((uint64_t)op[1] * 17998220869296846599UL) + ((uint64_t)op[2] * 6717411263826142866UL) + ((uint64_t)op[3] * 13091019643316963401UL) + ((uint64_t)op[4] * 16725393099077756620UL) + ((uint64_t)op[5] * 9805895641677525383UL) + ((uint64_t)op[6] * 10773746088520450151UL) + ((((uint64_t)op[7] * 4592429164023824795UL) + ((uint64_t)op[8] * 15721672951029712405UL) + ((uint64_t)op[9] * 16556497540425869029UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 16556497540425869029UL) + ((uint64_t)op[1] * 10437050973780786078UL) + ((uint64_t)op[2] * 17998220869296846599UL) + ((uint64_t)op[3] * 6717411263826142866UL) + ((uint64_t)op[4] * 13091019643316963401UL) + ((uint64_t)op[5] * 16725393099077756620UL) + ((uint64_t)op[6] * 9805895641677525383UL) + ((uint64_t)op[7] * 10773746088520450151UL) + ((((uint64_t)op[8] * 4592429164023824795UL) + ((uint64_t)op[9] * 15721672951029712405UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 15721672951029712405UL) + ((uint64_t)op[1] * 16556497540425869029UL) + ((uint64_t)op[2] * 10437050973780786078UL) + ((uint64_t)op[3] * 17998220869296846599UL) + ((uint64_t)op[4] * 6717411263826142866UL) + ((uint64_t)op[5] * 13091019643316963401UL) + ((uint64_t)op[6] * 16725393099077756620UL) + ((uint64_t)op[7] * 9805895641677525383UL) + ((uint64_t)op[8] * 10773746088520450151UL) + ((uint64_t)op[9] * 9261885745661902026UL);
	tmp_q[9] = ((uint64_t)op[0] * 4592429164023824795UL) + ((uint64_t)op[1] * 15721672951029712405UL) + ((uint64_t)op[2] * 16556497540425869029UL) + ((uint64_t)op[3] * 10437050973780786078UL) + ((uint64_t)op[4] * 17998220869296846599UL) + ((uint64_t)op[5] * 6717411263826142866UL) + ((uint64_t)op[6] * 13091019643316963401UL) + ((uint64_t)op[7] * 16725393099077756620UL) + ((uint64_t)op[8] * 9805895641677525383UL) + ((uint64_t)op[9] * 10773746088520450151UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 1241807495495201L) - ((-((int128)tmp_q[1] * 1203014628693242L) + ((int128)tmp_q[2] * 726719569547476L) + ((int128)tmp_q[3] * 2689739605078305L) - ((int128)tmp_q[4] * 3295676602270602L) - ((int128)tmp_q[5] * 125853016067341L) + ((int128)tmp_q[6] * 3986845734040463L) + ((int128)tmp_q[7] * 656809476034922L) - ((int128)tmp_q[8] * 62671343314165L) + ((int128)tmp_q[9] * 1114881539696239L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 1114881539696239L) + ((int128)tmp_q[1] * 1241807495495201L) - ((-((int128)tmp_q[2] * 1203014628693242L) + ((int128)tmp_q[3] * 726719569547476L) + ((int128)tmp_q[4] * 2689739605078305L) - ((int128)tmp_q[5] * 3295676602270602L) - ((int128)tmp_q[6] * 125853016067341L) + ((int128)tmp_q[7] * 3986845734040463L) + ((int128)tmp_q[8] * 656809476034922L) - ((int128)tmp_q[9] * 62671343314165L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 62671343314165L) + ((int128)tmp_q[1] * 1114881539696239L) + ((int128)tmp_q[2] * 1241807495495201L) - ((-((int128)tmp_q[3] * 1203014628693242L) + ((int128)tmp_q[4] * 726719569547476L) + ((int128)tmp_q[5] * 2689739605078305L) - ((int128)tmp_q[6] * 3295676602270602L) - ((int128)tmp_q[7] * 125853016067341L) + ((int128)tmp_q[8] * 3986845734040463L) + ((int128)tmp_q[9] * 656809476034922L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 656809476034922L) - ((int128)tmp_q[1] * 62671343314165L) + ((int128)tmp_q[2] * 1114881539696239L) + ((int128)tmp_q[3] * 1241807495495201L) - ((-((int128)tmp_q[4] * 1203014628693242L) + ((int128)tmp_q[5] * 726719569547476L) + ((int128)tmp_q[6] * 2689739605078305L) - ((int128)tmp_q[7] * 3295676602270602L) - ((int128)tmp_q[8] * 125853016067341L) + ((int128)tmp_q[9] * 3986845734040463L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 3986845734040463L) + ((int128)tmp_q[1] * 656809476034922L) - ((int128)tmp_q[2] * 62671343314165L) + ((int128)tmp_q[3] * 1114881539696239L) + ((int128)tmp_q[4] * 1241807495495201L) - ((-((int128)tmp_q[5] * 1203014628693242L) + ((int128)tmp_q[6] * 726719569547476L) + ((int128)tmp_q[7] * 2689739605078305L) - ((int128)tmp_q[8] * 3295676602270602L) - ((int128)tmp_q[9] * 125853016067341L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 125853016067341L) + ((int128)tmp_q[1] * 3986845734040463L) + ((int128)tmp_q[2] * 656809476034922L) - ((int128)tmp_q[3] * 62671343314165L) + ((int128)tmp_q[4] * 1114881539696239L) + ((int128)tmp_q[5] * 1241807495495201L) - ((-((int128)tmp_q[6] * 1203014628693242L) + ((int128)tmp_q[7] * 726719569547476L) + ((int128)tmp_q[8] * 2689739605078305L) - ((int128)tmp_q[9] * 3295676602270602L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 3295676602270602L) - ((int128)tmp_q[1] * 125853016067341L) + ((int128)tmp_q[2] * 3986845734040463L) + ((int128)tmp_q[3] * 656809476034922L) - ((int128)tmp_q[4] * 62671343314165L) + ((int128)tmp_q[5] * 1114881539696239L) + ((int128)tmp_q[6] * 1241807495495201L) - ((-((int128)tmp_q[7] * 1203014628693242L) + ((int128)tmp_q[8] * 726719569547476L) + ((int128)tmp_q[9] * 2689739605078305L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 2689739605078305L) - ((int128)tmp_q[1] * 3295676602270602L) - ((int128)tmp_q[2] * 125853016067341L) + ((int128)tmp_q[3] * 3986845734040463L) + ((int128)tmp_q[4] * 656809476034922L) - ((int128)tmp_q[5] * 62671343314165L) + ((int128)tmp_q[6] * 1114881539696239L) + ((int128)tmp_q[7] * 1241807495495201L) - ((-((int128)tmp_q[8] * 1203014628693242L) + ((int128)tmp_q[9] * 726719569547476L)) * 2);
	tmp_zero[8] = ((int128)tmp_q[0] * 726719569547476L) + ((int128)tmp_q[1] * 2689739605078305L) - ((int128)tmp_q[2] * 3295676602270602L) - ((int128)tmp_q[3] * 125853016067341L) + ((int128)tmp_q[4] * 3986845734040463L) + ((int128)tmp_q[5] * 656809476034922L) - ((int128)tmp_q[6] * 62671343314165L) + ((int128)tmp_q[7] * 1114881539696239L) + ((int128)tmp_q[8] * 1241807495495201L) + ((int128)tmp_q[9] * 2406029257386484L);
	tmp_zero[9] = -((int128)tmp_q[0] * 1203014628693242L) + ((int128)tmp_q[1] * 726719569547476L) + ((int128)tmp_q[2] * 2689739605078305L) - ((int128)tmp_q[3] * 3295676602270602L) - ((int128)tmp_q[4] * 125853016067341L) + ((int128)tmp_q[5] * 3986845734040463L) + ((int128)tmp_q[6] * 656809476034922L) - ((int128)tmp_q[7] * 62671343314165L) + ((int128)tmp_q[8] * 1114881539696239L) + ((int128)tmp_q[9] * 1241807495495201L);

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
}

