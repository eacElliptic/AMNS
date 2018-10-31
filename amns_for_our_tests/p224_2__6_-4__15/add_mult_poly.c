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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1]) << 2);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2]) << 2);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3]) << 2);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[5] + (int128)pa[5] * pb[4]) << 2);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[5]) << 2);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1]) << 1) + (int128)pa[3] * pa[3]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2]) << 3);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[5] * pa[3]) << 1) + (int128)pa[4] * pa[4]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[5] * pa[4]) << 3);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[5] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 17666537222290543887UL) + ((((uint64_t)op[1] * 16056000695422086332UL) + ((uint64_t)op[2] * 1779189931874062160UL) + ((uint64_t)op[3] * 162169907545555511UL) + ((uint64_t)op[4] * 12938373705442302403UL) + ((uint64_t)op[5] * 16473793704853350357UL)) * 18446744073709551612);
	tmp_q[1] = ((uint64_t)op[0] * 16473793704853350357UL) + ((uint64_t)op[1] * 17666537222290543887UL) + ((((uint64_t)op[2] * 16056000695422086332UL) + ((uint64_t)op[3] * 1779189931874062160UL) + ((uint64_t)op[4] * 162169907545555511UL) + ((uint64_t)op[5] * 12938373705442302403UL)) * 18446744073709551612);
	tmp_q[2] = ((uint64_t)op[0] * 12938373705442302403UL) + ((uint64_t)op[1] * 16473793704853350357UL) + ((uint64_t)op[2] * 17666537222290543887UL) + ((((uint64_t)op[3] * 16056000695422086332UL) + ((uint64_t)op[4] * 1779189931874062160UL) + ((uint64_t)op[5] * 162169907545555511UL)) * 18446744073709551612);
	tmp_q[3] = ((uint64_t)op[0] * 162169907545555511UL) + ((uint64_t)op[1] * 12938373705442302403UL) + ((uint64_t)op[2] * 16473793704853350357UL) + ((uint64_t)op[3] * 17666537222290543887UL) + ((((uint64_t)op[4] * 16056000695422086332UL) + ((uint64_t)op[5] * 1779189931874062160UL)) * 18446744073709551612);
	tmp_q[4] = ((uint64_t)op[0] * 1779189931874062160UL) + ((uint64_t)op[1] * 162169907545555511UL) + ((uint64_t)op[2] * 12938373705442302403UL) + ((uint64_t)op[3] * 16473793704853350357UL) + ((uint64_t)op[4] * 17666537222290543887UL) + ((uint64_t)op[5] * 9562973513149861136UL);
	tmp_q[5] = ((uint64_t)op[0] * 16056000695422086332UL) + ((uint64_t)op[1] * 1779189931874062160UL) + ((uint64_t)op[2] * 162169907545555511UL) + ((uint64_t)op[3] * 12938373705442302403UL) + ((uint64_t)op[4] * 16473793704853350357UL) + ((uint64_t)op[5] * 17666537222290543887UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 61812291025L) - ((-((int128)tmp_q[1] * 17772690621L) + ((int128)tmp_q[2] * 11958125505L) - ((int128)tmp_q[3] * 17677614546L) + ((int128)tmp_q[4] * 98894483928L) - ((int128)tmp_q[5] * 29271183307L)) * 4);
	tmp_zero[1] = -((int128)tmp_q[0] * 29271183307L) + ((int128)tmp_q[1] * 61812291025L) - ((-((int128)tmp_q[2] * 17772690621L) + ((int128)tmp_q[3] * 11958125505L) - ((int128)tmp_q[4] * 17677614546L) + ((int128)tmp_q[5] * 98894483928L)) * 4);
	tmp_zero[2] = ((int128)tmp_q[0] * 98894483928L) - ((int128)tmp_q[1] * 29271183307L) + ((int128)tmp_q[2] * 61812291025L) - ((-((int128)tmp_q[3] * 17772690621L) + ((int128)tmp_q[4] * 11958125505L) - ((int128)tmp_q[5] * 17677614546L)) * 4);
	tmp_zero[3] = -((int128)tmp_q[0] * 17677614546L) + ((int128)tmp_q[1] * 98894483928L) - ((int128)tmp_q[2] * 29271183307L) + ((int128)tmp_q[3] * 61812291025L) - ((-((int128)tmp_q[4] * 17772690621L) + ((int128)tmp_q[5] * 11958125505L)) * 4);
	tmp_zero[4] = ((int128)tmp_q[0] * 11958125505L) - ((int128)tmp_q[1] * 17677614546L) + ((int128)tmp_q[2] * 98894483928L) - ((int128)tmp_q[3] * 29271183307L) + ((int128)tmp_q[4] * 61812291025L) + ((int128)tmp_q[5] * 71090762484L);
	tmp_zero[5] = -((int128)tmp_q[0] * 17772690621L) + ((int128)tmp_q[1] * 11958125505L) - ((int128)tmp_q[2] * 17677614546L) + ((int128)tmp_q[3] * 98894483928L) - ((int128)tmp_q[4] * 29271183307L) + ((int128)tmp_q[5] * 61812291025L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
}

