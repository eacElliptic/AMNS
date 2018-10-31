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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[9] + (int128)pa[9] * pb[8]) * 3);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] - (((int128)pa[9] * pb[9]) * 3);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5]) << 1) + (int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((((int128)pa[9] * pa[7]) << 1) + (int128)pa[8] * pa[8]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[9] * pa[8]) * 6);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] - (((int128)pa[9] * pa[9]) * 3);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 9945979388966912969UL) + ((((uint64_t)op[1] * 17629625551606160677UL) + ((uint64_t)op[2] * 1470842542593020985UL) + ((uint64_t)op[3] * 18202253505623178728UL) + ((uint64_t)op[4] * 3031247507902701173UL) + ((uint64_t)op[5] * 2803312231869001829UL) + ((uint64_t)op[6] * 12948586714636451487UL) + ((uint64_t)op[7] * 8134912577998312117UL) + ((uint64_t)op[8] * 1530701680287441903UL) + ((uint64_t)op[9] * 403101767405704575UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 403101767405704575UL) + ((uint64_t)op[1] * 9945979388966912969UL) + ((((uint64_t)op[2] * 17629625551606160677UL) + ((uint64_t)op[3] * 1470842542593020985UL) + ((uint64_t)op[4] * 18202253505623178728UL) + ((uint64_t)op[5] * 3031247507902701173UL) + ((uint64_t)op[6] * 2803312231869001829UL) + ((uint64_t)op[7] * 12948586714636451487UL) + ((uint64_t)op[8] * 8134912577998312117UL) + ((uint64_t)op[9] * 1530701680287441903UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 1530701680287441903UL) + ((uint64_t)op[1] * 403101767405704575UL) + ((uint64_t)op[2] * 9945979388966912969UL) + ((((uint64_t)op[3] * 17629625551606160677UL) + ((uint64_t)op[4] * 1470842542593020985UL) + ((uint64_t)op[5] * 18202253505623178728UL) + ((uint64_t)op[6] * 3031247507902701173UL) + ((uint64_t)op[7] * 2803312231869001829UL) + ((uint64_t)op[8] * 12948586714636451487UL) + ((uint64_t)op[9] * 8134912577998312117UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 8134912577998312117UL) + ((uint64_t)op[1] * 1530701680287441903UL) + ((uint64_t)op[2] * 403101767405704575UL) + ((uint64_t)op[3] * 9945979388966912969UL) + ((((uint64_t)op[4] * 17629625551606160677UL) + ((uint64_t)op[5] * 1470842542593020985UL) + ((uint64_t)op[6] * 18202253505623178728UL) + ((uint64_t)op[7] * 3031247507902701173UL) + ((uint64_t)op[8] * 2803312231869001829UL) + ((uint64_t)op[9] * 12948586714636451487UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 12948586714636451487UL) + ((uint64_t)op[1] * 8134912577998312117UL) + ((uint64_t)op[2] * 1530701680287441903UL) + ((uint64_t)op[3] * 403101767405704575UL) + ((uint64_t)op[4] * 9945979388966912969UL) + ((((uint64_t)op[5] * 17629625551606160677UL) + ((uint64_t)op[6] * 1470842542593020985UL) + ((uint64_t)op[7] * 18202253505623178728UL) + ((uint64_t)op[8] * 3031247507902701173UL) + ((uint64_t)op[9] * 2803312231869001829UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 2803312231869001829UL) + ((uint64_t)op[1] * 12948586714636451487UL) + ((uint64_t)op[2] * 8134912577998312117UL) + ((uint64_t)op[3] * 1530701680287441903UL) + ((uint64_t)op[4] * 403101767405704575UL) + ((uint64_t)op[5] * 9945979388966912969UL) + ((((uint64_t)op[6] * 17629625551606160677UL) + ((uint64_t)op[7] * 1470842542593020985UL) + ((uint64_t)op[8] * 18202253505623178728UL) + ((uint64_t)op[9] * 3031247507902701173UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 3031247507902701173UL) + ((uint64_t)op[1] * 2803312231869001829UL) + ((uint64_t)op[2] * 12948586714636451487UL) + ((uint64_t)op[3] * 8134912577998312117UL) + ((uint64_t)op[4] * 1530701680287441903UL) + ((uint64_t)op[5] * 403101767405704575UL) + ((uint64_t)op[6] * 9945979388966912969UL) + ((((uint64_t)op[7] * 17629625551606160677UL) + ((uint64_t)op[8] * 1470842542593020985UL) + ((uint64_t)op[9] * 18202253505623178728UL)) * 18446744073709551613);
	tmp_q[7] = ((uint64_t)op[0] * 18202253505623178728UL) + ((uint64_t)op[1] * 3031247507902701173UL) + ((uint64_t)op[2] * 2803312231869001829UL) + ((uint64_t)op[3] * 12948586714636451487UL) + ((uint64_t)op[4] * 8134912577998312117UL) + ((uint64_t)op[5] * 1530701680287441903UL) + ((uint64_t)op[6] * 403101767405704575UL) + ((uint64_t)op[7] * 9945979388966912969UL) + ((((uint64_t)op[8] * 17629625551606160677UL) + ((uint64_t)op[9] * 1470842542593020985UL)) * 18446744073709551613);
	tmp_q[8] = ((uint64_t)op[0] * 1470842542593020985UL) + ((uint64_t)op[1] * 18202253505623178728UL) + ((uint64_t)op[2] * 3031247507902701173UL) + ((uint64_t)op[3] * 2803312231869001829UL) + ((uint64_t)op[4] * 12948586714636451487UL) + ((uint64_t)op[5] * 8134912577998312117UL) + ((uint64_t)op[6] * 1530701680287441903UL) + ((uint64_t)op[7] * 403101767405704575UL) + ((uint64_t)op[8] * 9945979388966912969UL) + ((uint64_t)op[9] * 2451355566310172817UL);
	tmp_q[9] = ((uint64_t)op[0] * 17629625551606160677UL) + ((uint64_t)op[1] * 1470842542593020985UL) + ((uint64_t)op[2] * 18202253505623178728UL) + ((uint64_t)op[3] * 3031247507902701173UL) + ((uint64_t)op[4] * 2803312231869001829UL) + ((uint64_t)op[5] * 12948586714636451487UL) + ((uint64_t)op[6] * 8134912577998312117UL) + ((uint64_t)op[7] * 1530701680287441903UL) + ((uint64_t)op[8] * 403101767405704575UL) + ((uint64_t)op[9] * 9945979388966912969UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 173487364602087L) - ((((int128)tmp_q[1] * 1098654326321703L) - ((int128)tmp_q[2] * 263603796747827L) + ((int128)tmp_q[3] * 1495847326958773L) + ((int128)tmp_q[4] * 1950366959122615L) + ((int128)tmp_q[5] * 599747169748435L) + ((int128)tmp_q[6] * 2240019851360387L) + ((int128)tmp_q[7] * 2792694875843928L) - ((int128)tmp_q[8] * 101400704013961L) - ((int128)tmp_q[9] * 1636883694301409L)) * 3);
	tmp_zero[1] = -((int128)tmp_q[0] * 1636883694301409L) - ((int128)tmp_q[1] * 173487364602087L) - ((((int128)tmp_q[2] * 1098654326321703L) - ((int128)tmp_q[3] * 263603796747827L) + ((int128)tmp_q[4] * 1495847326958773L) + ((int128)tmp_q[5] * 1950366959122615L) + ((int128)tmp_q[6] * 599747169748435L) + ((int128)tmp_q[7] * 2240019851360387L) + ((int128)tmp_q[8] * 2792694875843928L) - ((int128)tmp_q[9] * 101400704013961L)) * 3);
	tmp_zero[2] = -((int128)tmp_q[0] * 101400704013961L) - ((int128)tmp_q[1] * 1636883694301409L) - ((int128)tmp_q[2] * 173487364602087L) - ((((int128)tmp_q[3] * 1098654326321703L) - ((int128)tmp_q[4] * 263603796747827L) + ((int128)tmp_q[5] * 1495847326958773L) + ((int128)tmp_q[6] * 1950366959122615L) + ((int128)tmp_q[7] * 599747169748435L) + ((int128)tmp_q[8] * 2240019851360387L) + ((int128)tmp_q[9] * 2792694875843928L)) * 3);
	tmp_zero[3] = ((int128)tmp_q[0] * 2792694875843928L) - ((int128)tmp_q[1] * 101400704013961L) - ((int128)tmp_q[2] * 1636883694301409L) - ((int128)tmp_q[3] * 173487364602087L) - ((((int128)tmp_q[4] * 1098654326321703L) - ((int128)tmp_q[5] * 263603796747827L) + ((int128)tmp_q[6] * 1495847326958773L) + ((int128)tmp_q[7] * 1950366959122615L) + ((int128)tmp_q[8] * 599747169748435L) + ((int128)tmp_q[9] * 2240019851360387L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 2240019851360387L) + ((int128)tmp_q[1] * 2792694875843928L) - ((int128)tmp_q[2] * 101400704013961L) - ((int128)tmp_q[3] * 1636883694301409L) - ((int128)tmp_q[4] * 173487364602087L) - ((((int128)tmp_q[5] * 1098654326321703L) - ((int128)tmp_q[6] * 263603796747827L) + ((int128)tmp_q[7] * 1495847326958773L) + ((int128)tmp_q[8] * 1950366959122615L) + ((int128)tmp_q[9] * 599747169748435L)) * 3);
	tmp_zero[5] = ((int128)tmp_q[0] * 599747169748435L) + ((int128)tmp_q[1] * 2240019851360387L) + ((int128)tmp_q[2] * 2792694875843928L) - ((int128)tmp_q[3] * 101400704013961L) - ((int128)tmp_q[4] * 1636883694301409L) - ((int128)tmp_q[5] * 173487364602087L) - ((((int128)tmp_q[6] * 1098654326321703L) - ((int128)tmp_q[7] * 263603796747827L) + ((int128)tmp_q[8] * 1495847326958773L) + ((int128)tmp_q[9] * 1950366959122615L)) * 3);
	tmp_zero[6] = ((int128)tmp_q[0] * 1950366959122615L) + ((int128)tmp_q[1] * 599747169748435L) + ((int128)tmp_q[2] * 2240019851360387L) + ((int128)tmp_q[3] * 2792694875843928L) - ((int128)tmp_q[4] * 101400704013961L) - ((int128)tmp_q[5] * 1636883694301409L) - ((int128)tmp_q[6] * 173487364602087L) - ((((int128)tmp_q[7] * 1098654326321703L) - ((int128)tmp_q[8] * 263603796747827L) + ((int128)tmp_q[9] * 1495847326958773L)) * 3);
	tmp_zero[7] = ((int128)tmp_q[0] * 1495847326958773L) + ((int128)tmp_q[1] * 1950366959122615L) + ((int128)tmp_q[2] * 599747169748435L) + ((int128)tmp_q[3] * 2240019851360387L) + ((int128)tmp_q[4] * 2792694875843928L) - ((int128)tmp_q[5] * 101400704013961L) - ((int128)tmp_q[6] * 1636883694301409L) - ((int128)tmp_q[7] * 173487364602087L) - ((((int128)tmp_q[8] * 1098654326321703L) - ((int128)tmp_q[9] * 263603796747827L)) * 3);
	tmp_zero[8] = -((int128)tmp_q[0] * 263603796747827L) + ((int128)tmp_q[1] * 1495847326958773L) + ((int128)tmp_q[2] * 1950366959122615L) + ((int128)tmp_q[3] * 599747169748435L) + ((int128)tmp_q[4] * 2240019851360387L) + ((int128)tmp_q[5] * 2792694875843928L) - ((int128)tmp_q[6] * 101400704013961L) - ((int128)tmp_q[7] * 1636883694301409L) - ((int128)tmp_q[8] * 173487364602087L) - ((int128)tmp_q[9] * 3295962978965109L);
	tmp_zero[9] = ((int128)tmp_q[0] * 1098654326321703L) - ((int128)tmp_q[1] * 263603796747827L) + ((int128)tmp_q[2] * 1495847326958773L) + ((int128)tmp_q[3] * 1950366959122615L) + ((int128)tmp_q[4] * 599747169748435L) + ((int128)tmp_q[5] * 2240019851360387L) + ((int128)tmp_q[6] * 2792694875843928L) - ((int128)tmp_q[7] * 101400704013961L) - ((int128)tmp_q[8] * 1636883694301409L) - ((int128)tmp_q[9] * 173487364602087L);

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

