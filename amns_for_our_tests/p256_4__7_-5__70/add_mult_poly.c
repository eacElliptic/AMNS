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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1]) * 5);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2]) * 5);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3]) * 5);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4]) * 5);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[6] + (int128)pa[6] * pb[5]) * 5);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[6]) * 5);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1]) * 10);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2]) << 1) + (int128)pa[4] * pa[4]) * 5);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3]) * 10);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[6] * pa[4]) << 1) + (int128)pa[5] * pa[5]) * 5);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[6] * pa[5]) * 10);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[6] * pa[6]) * 5);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 11236834908608899785UL) + ((((uint64_t)op[1] * 13929577311704013087UL) + ((uint64_t)op[2] * 1130472805151230593UL) + ((uint64_t)op[3] * 1274812289532549312UL) + ((uint64_t)op[4] * 16006942133368431883UL) + ((uint64_t)op[5] * 2933066524258901547UL) + ((uint64_t)op[6] * 12981935686609147442UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 12981935686609147442UL) + ((uint64_t)op[1] * 11236834908608899785UL) + ((((uint64_t)op[2] * 13929577311704013087UL) + ((uint64_t)op[3] * 1130472805151230593UL) + ((uint64_t)op[4] * 1274812289532549312UL) + ((uint64_t)op[5] * 16006942133368431883UL) + ((uint64_t)op[6] * 2933066524258901547UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 2933066524258901547UL) + ((uint64_t)op[1] * 12981935686609147442UL) + ((uint64_t)op[2] * 11236834908608899785UL) + ((((uint64_t)op[3] * 13929577311704013087UL) + ((uint64_t)op[4] * 1130472805151230593UL) + ((uint64_t)op[5] * 1274812289532549312UL) + ((uint64_t)op[6] * 16006942133368431883UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 16006942133368431883UL) + ((uint64_t)op[1] * 2933066524258901547UL) + ((uint64_t)op[2] * 12981935686609147442UL) + ((uint64_t)op[3] * 11236834908608899785UL) + ((((uint64_t)op[4] * 13929577311704013087UL) + ((uint64_t)op[5] * 1130472805151230593UL) + ((uint64_t)op[6] * 1274812289532549312UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 1274812289532549312UL) + ((uint64_t)op[1] * 16006942133368431883UL) + ((uint64_t)op[2] * 2933066524258901547UL) + ((uint64_t)op[3] * 12981935686609147442UL) + ((uint64_t)op[4] * 11236834908608899785UL) + ((((uint64_t)op[5] * 13929577311704013087UL) + ((uint64_t)op[6] * 1130472805151230593UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 1130472805151230593UL) + ((uint64_t)op[1] * 1274812289532549312UL) + ((uint64_t)op[2] * 16006942133368431883UL) + ((uint64_t)op[3] * 2933066524258901547UL) + ((uint64_t)op[4] * 12981935686609147442UL) + ((uint64_t)op[5] * 11236834908608899785UL) + ((uint64_t)op[6] * 4139089736318141029UL);
	tmp_q[6] = ((uint64_t)op[0] * 13929577311704013087UL) + ((uint64_t)op[1] * 1130472805151230593UL) + ((uint64_t)op[2] * 1274812289532549312UL) + ((uint64_t)op[3] * 16006942133368431883UL) + ((uint64_t)op[4] * 2933066524258901547UL) + ((uint64_t)op[5] * 12981935686609147442UL) + ((uint64_t)op[6] * 11236834908608899785UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 21446049357L) - ((-((int128)tmp_q[1] * 18822999610L) - ((int128)tmp_q[2] * 30285002826L) - ((int128)tmp_q[3] * 40166210216L) + ((int128)tmp_q[4] * 18883318146L) - ((int128)tmp_q[5] * 34554334615L) + ((int128)tmp_q[6] * 40615181235L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 40615181235L) - ((int128)tmp_q[1] * 21446049357L) - ((-((int128)tmp_q[2] * 18822999610L) - ((int128)tmp_q[3] * 30285002826L) - ((int128)tmp_q[4] * 40166210216L) + ((int128)tmp_q[5] * 18883318146L) - ((int128)tmp_q[6] * 34554334615L)) * 5);
	tmp_zero[2] = -((int128)tmp_q[0] * 34554334615L) + ((int128)tmp_q[1] * 40615181235L) - ((int128)tmp_q[2] * 21446049357L) - ((-((int128)tmp_q[3] * 18822999610L) - ((int128)tmp_q[4] * 30285002826L) - ((int128)tmp_q[5] * 40166210216L) + ((int128)tmp_q[6] * 18883318146L)) * 5);
	tmp_zero[3] = ((int128)tmp_q[0] * 18883318146L) - ((int128)tmp_q[1] * 34554334615L) + ((int128)tmp_q[2] * 40615181235L) - ((int128)tmp_q[3] * 21446049357L) - ((-((int128)tmp_q[4] * 18822999610L) - ((int128)tmp_q[5] * 30285002826L) - ((int128)tmp_q[6] * 40166210216L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 40166210216L) + ((int128)tmp_q[1] * 18883318146L) - ((int128)tmp_q[2] * 34554334615L) + ((int128)tmp_q[3] * 40615181235L) - ((int128)tmp_q[4] * 21446049357L) - ((-((int128)tmp_q[5] * 18822999610L) - ((int128)tmp_q[6] * 30285002826L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 30285002826L) - ((int128)tmp_q[1] * 40166210216L) + ((int128)tmp_q[2] * 18883318146L) - ((int128)tmp_q[3] * 34554334615L) + ((int128)tmp_q[4] * 40615181235L) - ((int128)tmp_q[5] * 21446049357L) + ((int128)tmp_q[6] * 94114998050L);
	tmp_zero[6] = -((int128)tmp_q[0] * 18822999610L) - ((int128)tmp_q[1] * 30285002826L) - ((int128)tmp_q[2] * 40166210216L) + ((int128)tmp_q[3] * 18883318146L) - ((int128)tmp_q[4] * 34554334615L) + ((int128)tmp_q[5] * 40615181235L) - ((int128)tmp_q[6] * 21446049357L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

