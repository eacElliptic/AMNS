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
	tmp_q[0] = ((uint64_t)op[0] * 10775693794249288623UL) + ((((uint64_t)op[1] * 2151175711009165943UL) + ((uint64_t)op[2] * 11054907801860823595UL) + ((uint64_t)op[3] * 5729991550902513601UL) + ((uint64_t)op[4] * 12204061623497279750UL) + ((uint64_t)op[5] * 487917035510503881UL) + ((uint64_t)op[6] * 9653425029663906252UL)) * 18446744073709551611);
	tmp_q[1] = ((uint64_t)op[0] * 9653425029663906252UL) + ((uint64_t)op[1] * 10775693794249288623UL) + ((((uint64_t)op[2] * 2151175711009165943UL) + ((uint64_t)op[3] * 11054907801860823595UL) + ((uint64_t)op[4] * 5729991550902513601UL) + ((uint64_t)op[5] * 12204061623497279750UL) + ((uint64_t)op[6] * 487917035510503881UL)) * 18446744073709551611);
	tmp_q[2] = ((uint64_t)op[0] * 487917035510503881UL) + ((uint64_t)op[1] * 9653425029663906252UL) + ((uint64_t)op[2] * 10775693794249288623UL) + ((((uint64_t)op[3] * 2151175711009165943UL) + ((uint64_t)op[4] * 11054907801860823595UL) + ((uint64_t)op[5] * 5729991550902513601UL) + ((uint64_t)op[6] * 12204061623497279750UL)) * 18446744073709551611);
	tmp_q[3] = ((uint64_t)op[0] * 12204061623497279750UL) + ((uint64_t)op[1] * 487917035510503881UL) + ((uint64_t)op[2] * 9653425029663906252UL) + ((uint64_t)op[3] * 10775693794249288623UL) + ((((uint64_t)op[4] * 2151175711009165943UL) + ((uint64_t)op[5] * 11054907801860823595UL) + ((uint64_t)op[6] * 5729991550902513601UL)) * 18446744073709551611);
	tmp_q[4] = ((uint64_t)op[0] * 5729991550902513601UL) + ((uint64_t)op[1] * 12204061623497279750UL) + ((uint64_t)op[2] * 487917035510503881UL) + ((uint64_t)op[3] * 9653425029663906252UL) + ((uint64_t)op[4] * 10775693794249288623UL) + ((((uint64_t)op[5] * 2151175711009165943UL) + ((uint64_t)op[6] * 11054907801860823595UL)) * 18446744073709551611);
	tmp_q[5] = ((uint64_t)op[0] * 11054907801860823595UL) + ((uint64_t)op[1] * 5729991550902513601UL) + ((uint64_t)op[2] * 12204061623497279750UL) + ((uint64_t)op[3] * 487917035510503881UL) + ((uint64_t)op[4] * 9653425029663906252UL) + ((uint64_t)op[5] * 10775693794249288623UL) + ((uint64_t)op[6] * 7690865518663721901UL);
	tmp_q[6] = ((uint64_t)op[0] * 2151175711009165943UL) + ((uint64_t)op[1] * 11054907801860823595UL) + ((uint64_t)op[2] * 5729991550902513601UL) + ((uint64_t)op[3] * 12204061623497279750UL) + ((uint64_t)op[4] * 487917035510503881UL) + ((uint64_t)op[5] * 9653425029663906252UL) + ((uint64_t)op[6] * 10775693794249288623UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 55076905396L) - ((((int128)tmp_q[1] * 31586598468L) - ((int128)tmp_q[2] * 49539541263L) - ((int128)tmp_q[3] * 10006499640L) - ((int128)tmp_q[4] * 76546951654L) + ((int128)tmp_q[5] * 8673003797L) + ((int128)tmp_q[6] * 442142849L)) * 5);
	tmp_zero[1] = ((int128)tmp_q[0] * 442142849L) - ((int128)tmp_q[1] * 55076905396L) - ((((int128)tmp_q[2] * 31586598468L) - ((int128)tmp_q[3] * 49539541263L) - ((int128)tmp_q[4] * 10006499640L) - ((int128)tmp_q[5] * 76546951654L) + ((int128)tmp_q[6] * 8673003797L)) * 5);
	tmp_zero[2] = ((int128)tmp_q[0] * 8673003797L) + ((int128)tmp_q[1] * 442142849L) - ((int128)tmp_q[2] * 55076905396L) - ((((int128)tmp_q[3] * 31586598468L) - ((int128)tmp_q[4] * 49539541263L) - ((int128)tmp_q[5] * 10006499640L) - ((int128)tmp_q[6] * 76546951654L)) * 5);
	tmp_zero[3] = -((int128)tmp_q[0] * 76546951654L) + ((int128)tmp_q[1] * 8673003797L) + ((int128)tmp_q[2] * 442142849L) - ((int128)tmp_q[3] * 55076905396L) - ((((int128)tmp_q[4] * 31586598468L) - ((int128)tmp_q[5] * 49539541263L) - ((int128)tmp_q[6] * 10006499640L)) * 5);
	tmp_zero[4] = -((int128)tmp_q[0] * 10006499640L) - ((int128)tmp_q[1] * 76546951654L) + ((int128)tmp_q[2] * 8673003797L) + ((int128)tmp_q[3] * 442142849L) - ((int128)tmp_q[4] * 55076905396L) - ((((int128)tmp_q[5] * 31586598468L) - ((int128)tmp_q[6] * 49539541263L)) * 5);
	tmp_zero[5] = -((int128)tmp_q[0] * 49539541263L) - ((int128)tmp_q[1] * 10006499640L) - ((int128)tmp_q[2] * 76546951654L) + ((int128)tmp_q[3] * 8673003797L) + ((int128)tmp_q[4] * 442142849L) - ((int128)tmp_q[5] * 55076905396L) - ((int128)tmp_q[6] * 157932992340L);
	tmp_zero[6] = ((int128)tmp_q[0] * 31586598468L) - ((int128)tmp_q[1] * 49539541263L) - ((int128)tmp_q[2] * 10006499640L) - ((int128)tmp_q[3] * 76546951654L) + ((int128)tmp_q[4] * 8673003797L) + ((int128)tmp_q[5] * 442142849L) - ((int128)tmp_q[6] * 55076905396L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
}

