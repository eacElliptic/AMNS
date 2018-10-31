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
	tmp_q[0] = ((uint64_t)op[0] * 13726072779734294839UL) + ((((uint64_t)op[1] * 14733020954686674005UL) + ((uint64_t)op[2] * 3498383493595768667UL) + ((uint64_t)op[3] * 8361472621463891746UL) + ((uint64_t)op[4] * 14909877303284415588UL) + ((uint64_t)op[5] * 364714425179575154UL) + ((uint64_t)op[6] * 12615109125814804768UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 12615109125814804768UL) + ((uint64_t)op[1] * 13726072779734294839UL) + ((((uint64_t)op[2] * 14733020954686674005UL) + ((uint64_t)op[3] * 3498383493595768667UL) + ((uint64_t)op[4] * 8361472621463891746UL) + ((uint64_t)op[5] * 14909877303284415588UL) + ((uint64_t)op[6] * 364714425179575154UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 364714425179575154UL) + ((uint64_t)op[1] * 12615109125814804768UL) + ((uint64_t)op[2] * 13726072779734294839UL) + ((((uint64_t)op[3] * 14733020954686674005UL) + ((uint64_t)op[4] * 3498383493595768667UL) + ((uint64_t)op[5] * 8361472621463891746UL) + ((uint64_t)op[6] * 14909877303284415588UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 14909877303284415588UL) + ((uint64_t)op[1] * 364714425179575154UL) + ((uint64_t)op[2] * 12615109125814804768UL) + ((uint64_t)op[3] * 13726072779734294839UL) + ((((uint64_t)op[4] * 14733020954686674005UL) + ((uint64_t)op[5] * 3498383493595768667UL) + ((uint64_t)op[6] * 8361472621463891746UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 8361472621463891746UL) + ((uint64_t)op[1] * 14909877303284415588UL) + ((uint64_t)op[2] * 364714425179575154UL) + ((uint64_t)op[3] * 12615109125814804768UL) + ((uint64_t)op[4] * 13726072779734294839UL) + ((((uint64_t)op[5] * 14733020954686674005UL) + ((uint64_t)op[6] * 3498383493595768667UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 3498383493595768667UL) + ((uint64_t)op[1] * 8361472621463891746UL) + ((uint64_t)op[2] * 14909877303284415588UL) + ((uint64_t)op[3] * 364714425179575154UL) + ((uint64_t)op[4] * 12615109125814804768UL) + ((uint64_t)op[5] * 13726072779734294839UL) + ((uint64_t)op[6] * 11141169357068632833UL);
	tmp_q[6] = ((uint64_t)op[0] * 14733020954686674005UL) + ((uint64_t)op[1] * 3498383493595768667UL) + ((uint64_t)op[2] * 8361472621463891746UL) + ((uint64_t)op[3] * 14909877303284415588UL) + ((uint64_t)op[4] * 364714425179575154UL) + ((uint64_t)op[5] * 12615109125814804768UL) + ((uint64_t)op[6] * 13726072779734294839UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 7220519239L) - ((-((int128)tmp_q[1] * 13907315816L) + ((int128)tmp_q[2] * 46548650931L) - ((int128)tmp_q[3] * 17943812269L) + ((int128)tmp_q[4] * 27354316578L) + ((int128)tmp_q[5] * 7985741139L) - ((int128)tmp_q[6] * 13239542221L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 13239542221L) + ((int128)tmp_q[1] * 7220519239L) - ((-((int128)tmp_q[2] * 13907315816L) + ((int128)tmp_q[3] * 46548650931L) - ((int128)tmp_q[4] * 17943812269L) + ((int128)tmp_q[5] * 27354316578L) + ((int128)tmp_q[6] * 7985741139L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 7985741139L) - ((int128)tmp_q[1] * 13239542221L) + ((int128)tmp_q[2] * 7220519239L) - ((-((int128)tmp_q[3] * 13907315816L) + ((int128)tmp_q[4] * 46548650931L) - ((int128)tmp_q[5] * 17943812269L) + ((int128)tmp_q[6] * 27354316578L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 27354316578L) + ((int128)tmp_q[1] * 7985741139L) - ((int128)tmp_q[2] * 13239542221L) + ((int128)tmp_q[3] * 7220519239L) - ((-((int128)tmp_q[4] * 13907315816L) + ((int128)tmp_q[5] * 46548650931L) - ((int128)tmp_q[6] * 17943812269L)) * 3);
	tmp_zero[4] = -((int128)tmp_q[0] * 17943812269L) + ((int128)tmp_q[1] * 27354316578L) + ((int128)tmp_q[2] * 7985741139L) - ((int128)tmp_q[3] * 13239542221L) + ((int128)tmp_q[4] * 7220519239L) - ((-((int128)tmp_q[5] * 13907315816L) + ((int128)tmp_q[6] * 46548650931L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 46548650931L) - ((int128)tmp_q[1] * 17943812269L) + ((int128)tmp_q[2] * 27354316578L) + ((int128)tmp_q[3] * 7985741139L) - ((int128)tmp_q[4] * 13239542221L) + ((int128)tmp_q[5] * 7220519239L) + ((int128)tmp_q[6] * 41721947448L);
	tmp_zero[6] = -((int128)tmp_q[0] * 13907315816L) + ((int128)tmp_q[1] * 46548650931L) - ((int128)tmp_q[2] * 17943812269L) + ((int128)tmp_q[3] * 27354316578L) + ((int128)tmp_q[4] * 7985741139L) - ((int128)tmp_q[5] * 13239542221L) + ((int128)tmp_q[6] * 7220519239L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

