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
	tmp_q[0] = ((uint64_t)op[0] * 3869359117327632853UL) + ((((uint64_t)op[1] * 1126566556316948740UL) + ((uint64_t)op[2] * 1838084532748814302UL) + ((uint64_t)op[3] * 12611838636449518723UL) + ((uint64_t)op[4] * 5674683486801572332UL) + ((uint64_t)op[5] * 2217132845002094455UL) + ((uint64_t)op[6] * 5309596880153394035UL) + ((uint64_t)op[7] * 12557795114661437905UL) + ((uint64_t)op[8] * 2360773370770457083UL) + ((uint64_t)op[9] * 12530282819695681310UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 12530282819695681310UL) + ((uint64_t)op[1] * 3869359117327632853UL) + ((((uint64_t)op[2] * 1126566556316948740UL) + ((uint64_t)op[3] * 1838084532748814302UL) + ((uint64_t)op[4] * 12611838636449518723UL) + ((uint64_t)op[5] * 5674683486801572332UL) + ((uint64_t)op[6] * 2217132845002094455UL) + ((uint64_t)op[7] * 5309596880153394035UL) + ((uint64_t)op[8] * 12557795114661437905UL) + ((uint64_t)op[9] * 2360773370770457083UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 2360773370770457083UL) + ((uint64_t)op[1] * 12530282819695681310UL) + ((uint64_t)op[2] * 3869359117327632853UL) + ((((uint64_t)op[3] * 1126566556316948740UL) + ((uint64_t)op[4] * 1838084532748814302UL) + ((uint64_t)op[5] * 12611838636449518723UL) + ((uint64_t)op[6] * 5674683486801572332UL) + ((uint64_t)op[7] * 2217132845002094455UL) + ((uint64_t)op[8] * 5309596880153394035UL) + ((uint64_t)op[9] * 12557795114661437905UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 12557795114661437905UL) + ((uint64_t)op[1] * 2360773370770457083UL) + ((uint64_t)op[2] * 12530282819695681310UL) + ((uint64_t)op[3] * 3869359117327632853UL) + ((((uint64_t)op[4] * 1126566556316948740UL) + ((uint64_t)op[5] * 1838084532748814302UL) + ((uint64_t)op[6] * 12611838636449518723UL) + ((uint64_t)op[7] * 5674683486801572332UL) + ((uint64_t)op[8] * 2217132845002094455UL) + ((uint64_t)op[9] * 5309596880153394035UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 5309596880153394035UL) + ((uint64_t)op[1] * 12557795114661437905UL) + ((uint64_t)op[2] * 2360773370770457083UL) + ((uint64_t)op[3] * 12530282819695681310UL) + ((uint64_t)op[4] * 3869359117327632853UL) + ((((uint64_t)op[5] * 1126566556316948740UL) + ((uint64_t)op[6] * 1838084532748814302UL) + ((uint64_t)op[7] * 12611838636449518723UL) + ((uint64_t)op[8] * 5674683486801572332UL) + ((uint64_t)op[9] * 2217132845002094455UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 2217132845002094455UL) + ((uint64_t)op[1] * 5309596880153394035UL) + ((uint64_t)op[2] * 12557795114661437905UL) + ((uint64_t)op[3] * 2360773370770457083UL) + ((uint64_t)op[4] * 12530282819695681310UL) + ((uint64_t)op[5] * 3869359117327632853UL) + ((((uint64_t)op[6] * 1126566556316948740UL) + ((uint64_t)op[7] * 1838084532748814302UL) + ((uint64_t)op[8] * 12611838636449518723UL) + ((uint64_t)op[9] * 5674683486801572332UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 5674683486801572332UL) + ((uint64_t)op[1] * 2217132845002094455UL) + ((uint64_t)op[2] * 5309596880153394035UL) + ((uint64_t)op[3] * 12557795114661437905UL) + ((uint64_t)op[4] * 2360773370770457083UL) + ((uint64_t)op[5] * 12530282819695681310UL) + ((uint64_t)op[6] * 3869359117327632853UL) + ((((uint64_t)op[7] * 1126566556316948740UL) + ((uint64_t)op[8] * 1838084532748814302UL) + ((uint64_t)op[9] * 12611838636449518723UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 12611838636449518723UL) + ((uint64_t)op[1] * 5674683486801572332UL) + ((uint64_t)op[2] * 2217132845002094455UL) + ((uint64_t)op[3] * 5309596880153394035UL) + ((uint64_t)op[4] * 12557795114661437905UL) + ((uint64_t)op[5] * 2360773370770457083UL) + ((uint64_t)op[6] * 12530282819695681310UL) + ((uint64_t)op[7] * 3869359117327632853UL) + ((((uint64_t)op[8] * 1126566556316948740UL) + ((uint64_t)op[9] * 1838084532748814302UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 1838084532748814302UL) + ((uint64_t)op[1] * 12611838636449518723UL) + ((uint64_t)op[2] * 5674683486801572332UL) + ((uint64_t)op[3] * 2217132845002094455UL) + ((uint64_t)op[4] * 5309596880153394035UL) + ((uint64_t)op[5] * 12557795114661437905UL) + ((uint64_t)op[6] * 2360773370770457083UL) + ((uint64_t)op[7] * 12530282819695681310UL) + ((uint64_t)op[8] * 3869359117327632853UL) + ((uint64_t)op[9] * 16193610961075654136UL);
	tmp_q[9] = ((uint64_t)op[0] * 1126566556316948740UL) + ((uint64_t)op[1] * 1838084532748814302UL) + ((uint64_t)op[2] * 12611838636449518723UL) + ((uint64_t)op[3] * 5674683486801572332UL) + ((uint64_t)op[4] * 2217132845002094455UL) + ((uint64_t)op[5] * 5309596880153394035UL) + ((uint64_t)op[6] * 12557795114661437905UL) + ((uint64_t)op[7] * 2360773370770457083UL) + ((uint64_t)op[8] * 12530282819695681310UL) + ((uint64_t)op[9] * 3869359117327632853UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 1068032049586209L) - ((-((int128)tmp_q[1] * 732356205479644L) + ((int128)tmp_q[2] * 1587152653072592L) + ((int128)tmp_q[3] * 663076964773806L) - ((int128)tmp_q[4] * 1066242998266592L) + ((int128)tmp_q[5] * 2189623624226721L) + ((int128)tmp_q[6] * 1441199565655268L) + ((int128)tmp_q[7] * 1478012021826713L) - ((int128)tmp_q[8] * 1902570122274931L) - ((int128)tmp_q[9] * 2861786442825818L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 2861786442825818L) - ((int128)tmp_q[1] * 1068032049586209L) - ((-((int128)tmp_q[2] * 732356205479644L) + ((int128)tmp_q[3] * 1587152653072592L) + ((int128)tmp_q[4] * 663076964773806L) - ((int128)tmp_q[5] * 1066242998266592L) + ((int128)tmp_q[6] * 2189623624226721L) + ((int128)tmp_q[7] * 1441199565655268L) + ((int128)tmp_q[8] * 1478012021826713L) - ((int128)tmp_q[9] * 1902570122274931L)) * 2);
	tmp_zero[2] = -((int128)tmp_q[0] * 1902570122274931L) - ((int128)tmp_q[1] * 2861786442825818L) - ((int128)tmp_q[2] * 1068032049586209L) - ((-((int128)tmp_q[3] * 732356205479644L) + ((int128)tmp_q[4] * 1587152653072592L) + ((int128)tmp_q[5] * 663076964773806L) - ((int128)tmp_q[6] * 1066242998266592L) + ((int128)tmp_q[7] * 2189623624226721L) + ((int128)tmp_q[8] * 1441199565655268L) + ((int128)tmp_q[9] * 1478012021826713L)) * 2);
	tmp_zero[3] = ((int128)tmp_q[0] * 1478012021826713L) - ((int128)tmp_q[1] * 1902570122274931L) - ((int128)tmp_q[2] * 2861786442825818L) - ((int128)tmp_q[3] * 1068032049586209L) - ((-((int128)tmp_q[4] * 732356205479644L) + ((int128)tmp_q[5] * 1587152653072592L) + ((int128)tmp_q[6] * 663076964773806L) - ((int128)tmp_q[7] * 1066242998266592L) + ((int128)tmp_q[8] * 2189623624226721L) + ((int128)tmp_q[9] * 1441199565655268L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 1441199565655268L) + ((int128)tmp_q[1] * 1478012021826713L) - ((int128)tmp_q[2] * 1902570122274931L) - ((int128)tmp_q[3] * 2861786442825818L) - ((int128)tmp_q[4] * 1068032049586209L) - ((-((int128)tmp_q[5] * 732356205479644L) + ((int128)tmp_q[6] * 1587152653072592L) + ((int128)tmp_q[7] * 663076964773806L) - ((int128)tmp_q[8] * 1066242998266592L) + ((int128)tmp_q[9] * 2189623624226721L)) * 2);
	tmp_zero[5] = ((int128)tmp_q[0] * 2189623624226721L) + ((int128)tmp_q[1] * 1441199565655268L) + ((int128)tmp_q[2] * 1478012021826713L) - ((int128)tmp_q[3] * 1902570122274931L) - ((int128)tmp_q[4] * 2861786442825818L) - ((int128)tmp_q[5] * 1068032049586209L) - ((-((int128)tmp_q[6] * 732356205479644L) + ((int128)tmp_q[7] * 1587152653072592L) + ((int128)tmp_q[8] * 663076964773806L) - ((int128)tmp_q[9] * 1066242998266592L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 1066242998266592L) + ((int128)tmp_q[1] * 2189623624226721L) + ((int128)tmp_q[2] * 1441199565655268L) + ((int128)tmp_q[3] * 1478012021826713L) - ((int128)tmp_q[4] * 1902570122274931L) - ((int128)tmp_q[5] * 2861786442825818L) - ((int128)tmp_q[6] * 1068032049586209L) - ((-((int128)tmp_q[7] * 732356205479644L) + ((int128)tmp_q[8] * 1587152653072592L) + ((int128)tmp_q[9] * 663076964773806L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 663076964773806L) - ((int128)tmp_q[1] * 1066242998266592L) + ((int128)tmp_q[2] * 2189623624226721L) + ((int128)tmp_q[3] * 1441199565655268L) + ((int128)tmp_q[4] * 1478012021826713L) - ((int128)tmp_q[5] * 1902570122274931L) - ((int128)tmp_q[6] * 2861786442825818L) - ((int128)tmp_q[7] * 1068032049586209L) - ((-((int128)tmp_q[8] * 732356205479644L) + ((int128)tmp_q[9] * 1587152653072592L)) * 2);
	tmp_zero[8] = ((int128)tmp_q[0] * 1587152653072592L) + ((int128)tmp_q[1] * 663076964773806L) - ((int128)tmp_q[2] * 1066242998266592L) + ((int128)tmp_q[3] * 2189623624226721L) + ((int128)tmp_q[4] * 1441199565655268L) + ((int128)tmp_q[5] * 1478012021826713L) - ((int128)tmp_q[6] * 1902570122274931L) - ((int128)tmp_q[7] * 2861786442825818L) - ((int128)tmp_q[8] * 1068032049586209L) + ((int128)tmp_q[9] * 1464712410959288L);
	tmp_zero[9] = -((int128)tmp_q[0] * 732356205479644L) + ((int128)tmp_q[1] * 1587152653072592L) + ((int128)tmp_q[2] * 663076964773806L) - ((int128)tmp_q[3] * 1066242998266592L) + ((int128)tmp_q[4] * 2189623624226721L) + ((int128)tmp_q[5] * 1441199565655268L) + ((int128)tmp_q[6] * 1478012021826713L) - ((int128)tmp_q[7] * 1902570122274931L) - ((int128)tmp_q[8] * 2861786442825818L) - ((int128)tmp_q[9] * 1068032049586209L);

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

