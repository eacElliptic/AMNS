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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] + (((int128)pa[1] * pb[10] + (int128)pa[2] * pb[9] + (int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3] + (int128)pa[9] * pb[2] + (int128)pa[10] * pb[1]) * 7);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] + (((int128)pa[2] * pb[10] + (int128)pa[3] * pb[9] + (int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4] + (int128)pa[9] * pb[3] + (int128)pa[10] * pb[2]) * 7);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] + (((int128)pa[3] * pb[10] + (int128)pa[4] * pb[9] + (int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5] + (int128)pa[9] * pb[4] + (int128)pa[10] * pb[3]) * 7);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] + (((int128)pa[4] * pb[10] + (int128)pa[5] * pb[9] + (int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6] + (int128)pa[9] * pb[5] + (int128)pa[10] * pb[4]) * 7);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] + (((int128)pa[5] * pb[10] + (int128)pa[6] * pb[9] + (int128)pa[7] * pb[8] + (int128)pa[8] * pb[7] + (int128)pa[9] * pb[6] + (int128)pa[10] * pb[5]) * 7);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] + (((int128)pa[6] * pb[10] + (int128)pa[7] * pb[9] + (int128)pa[8] * pb[8] + (int128)pa[9] * pb[7] + (int128)pa[10] * pb[6]) * 7);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] + (((int128)pa[7] * pb[10] + (int128)pa[8] * pb[9] + (int128)pa[9] * pb[8] + (int128)pa[10] * pb[7]) * 7);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] + (((int128)pa[8] * pb[10] + (int128)pa[9] * pb[9] + (int128)pa[10] * pb[8]) * 7);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0] + (((int128)pa[9] * pb[10] + (int128)pa[10] * pb[9]) * 7);
	tmp_prod_result[9] = (int128)pa[0] * pb[9] + (int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1] + (int128)pa[9] * pb[0] + (((int128)pa[10] * pb[10]) * 7);
	tmp_prod_result[10] = (int128)pa[0] * pb[10] + (int128)pa[1] * pb[9] + (int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2] + (int128)pa[9] * pb[1] + (int128)pa[10] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] + (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3] + (int128)pa[9] * pa[2] + (int128)pa[10] * pa[1]) * 14);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) + (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4] + (int128)pa[9] * pa[3] + (int128)pa[10] * pa[2]) << 1) + (int128)pa[6] * pa[6]) * 7);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] + (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5] + (int128)pa[9] * pa[4] + (int128)pa[10] * pa[3]) * 14);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) + (((((int128)pa[8] * pa[6] + (int128)pa[9] * pa[5] + (int128)pa[10] * pa[4]) << 1) + (int128)pa[7] * pa[7]) * 7);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] + (((int128)pa[8] * pa[7] + (int128)pa[9] * pa[6] + (int128)pa[10] * pa[5]) * 14);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) + (((((int128)pa[9] * pa[7] + (int128)pa[10] * pa[6]) << 1) + (int128)pa[8] * pa[8]) * 7);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] + (((int128)pa[9] * pa[8] + (int128)pa[10] * pa[7]) * 14);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) + (((((int128)pa[10] * pa[8]) << 1) + (int128)pa[9] * pa[9]) * 7);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4] + (((int128)pa[10] * pa[9]) * 14);
	tmp_prod_result[9] = (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1] + (int128)pa[9] * pa[0]) << 1) + (((int128)pa[10] * pa[10]) * 7);
	tmp_prod_result[10] = (((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2] + (int128)pa[9] * pa[1] + (int128)pa[10] * pa[0]) << 1) + (int128)pa[5] * pa[5];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 1473647133326713317UL) + ((((uint64_t)op[1] * 4270738813664669008UL) + ((uint64_t)op[2] * 11554351881453359700UL) + ((uint64_t)op[3] * 10731662810726815966UL) + ((uint64_t)op[4] * 14725829993045922307UL) + ((uint64_t)op[5] * 13182533234761706359UL) + ((uint64_t)op[6] * 15285701718870050433UL) + ((uint64_t)op[7] * 13886158810488315615UL) + ((uint64_t)op[8] * 4605647320554802783UL) + ((uint64_t)op[9] * 3637739900359287000UL) + ((uint64_t)op[10] * 10466828628376630569UL)) * 7);
	tmp_q[1] = ((uint64_t)op[0] * 10466828628376630569UL) + ((uint64_t)op[1] * 1473647133326713317UL) + ((((uint64_t)op[2] * 4270738813664669008UL) + ((uint64_t)op[3] * 11554351881453359700UL) + ((uint64_t)op[4] * 10731662810726815966UL) + ((uint64_t)op[5] * 14725829993045922307UL) + ((uint64_t)op[6] * 13182533234761706359UL) + ((uint64_t)op[7] * 15285701718870050433UL) + ((uint64_t)op[8] * 13886158810488315615UL) + ((uint64_t)op[9] * 4605647320554802783UL) + ((uint64_t)op[10] * 3637739900359287000UL)) * 7);
	tmp_q[2] = ((uint64_t)op[0] * 3637739900359287000UL) + ((uint64_t)op[1] * 10466828628376630569UL) + ((uint64_t)op[2] * 1473647133326713317UL) + ((((uint64_t)op[3] * 4270738813664669008UL) + ((uint64_t)op[4] * 11554351881453359700UL) + ((uint64_t)op[5] * 10731662810726815966UL) + ((uint64_t)op[6] * 14725829993045922307UL) + ((uint64_t)op[7] * 13182533234761706359UL) + ((uint64_t)op[8] * 15285701718870050433UL) + ((uint64_t)op[9] * 13886158810488315615UL) + ((uint64_t)op[10] * 4605647320554802783UL)) * 7);
	tmp_q[3] = ((uint64_t)op[0] * 4605647320554802783UL) + ((uint64_t)op[1] * 3637739900359287000UL) + ((uint64_t)op[2] * 10466828628376630569UL) + ((uint64_t)op[3] * 1473647133326713317UL) + ((((uint64_t)op[4] * 4270738813664669008UL) + ((uint64_t)op[5] * 11554351881453359700UL) + ((uint64_t)op[6] * 10731662810726815966UL) + ((uint64_t)op[7] * 14725829993045922307UL) + ((uint64_t)op[8] * 13182533234761706359UL) + ((uint64_t)op[9] * 15285701718870050433UL) + ((uint64_t)op[10] * 13886158810488315615UL)) * 7);
	tmp_q[4] = ((uint64_t)op[0] * 13886158810488315615UL) + ((uint64_t)op[1] * 4605647320554802783UL) + ((uint64_t)op[2] * 3637739900359287000UL) + ((uint64_t)op[3] * 10466828628376630569UL) + ((uint64_t)op[4] * 1473647133326713317UL) + ((((uint64_t)op[5] * 4270738813664669008UL) + ((uint64_t)op[6] * 11554351881453359700UL) + ((uint64_t)op[7] * 10731662810726815966UL) + ((uint64_t)op[8] * 14725829993045922307UL) + ((uint64_t)op[9] * 13182533234761706359UL) + ((uint64_t)op[10] * 15285701718870050433UL)) * 7);
	tmp_q[5] = ((uint64_t)op[0] * 15285701718870050433UL) + ((uint64_t)op[1] * 13886158810488315615UL) + ((uint64_t)op[2] * 4605647320554802783UL) + ((uint64_t)op[3] * 3637739900359287000UL) + ((uint64_t)op[4] * 10466828628376630569UL) + ((uint64_t)op[5] * 1473647133326713317UL) + ((((uint64_t)op[6] * 4270738813664669008UL) + ((uint64_t)op[7] * 11554351881453359700UL) + ((uint64_t)op[8] * 10731662810726815966UL) + ((uint64_t)op[9] * 14725829993045922307UL) + ((uint64_t)op[10] * 13182533234761706359UL)) * 7);
	tmp_q[6] = ((uint64_t)op[0] * 13182533234761706359UL) + ((uint64_t)op[1] * 15285701718870050433UL) + ((uint64_t)op[2] * 13886158810488315615UL) + ((uint64_t)op[3] * 4605647320554802783UL) + ((uint64_t)op[4] * 3637739900359287000UL) + ((uint64_t)op[5] * 10466828628376630569UL) + ((uint64_t)op[6] * 1473647133326713317UL) + ((((uint64_t)op[7] * 4270738813664669008UL) + ((uint64_t)op[8] * 11554351881453359700UL) + ((uint64_t)op[9] * 10731662810726815966UL) + ((uint64_t)op[10] * 14725829993045922307UL)) * 7);
	tmp_q[7] = ((uint64_t)op[0] * 14725829993045922307UL) + ((uint64_t)op[1] * 13182533234761706359UL) + ((uint64_t)op[2] * 15285701718870050433UL) + ((uint64_t)op[3] * 13886158810488315615UL) + ((uint64_t)op[4] * 4605647320554802783UL) + ((uint64_t)op[5] * 3637739900359287000UL) + ((uint64_t)op[6] * 10466828628376630569UL) + ((uint64_t)op[7] * 1473647133326713317UL) + ((((uint64_t)op[8] * 4270738813664669008UL) + ((uint64_t)op[9] * 11554351881453359700UL) + ((uint64_t)op[10] * 10731662810726815966UL)) * 7);
	tmp_q[8] = ((uint64_t)op[0] * 10731662810726815966UL) + ((uint64_t)op[1] * 14725829993045922307UL) + ((uint64_t)op[2] * 13182533234761706359UL) + ((uint64_t)op[3] * 15285701718870050433UL) + ((uint64_t)op[4] * 13886158810488315615UL) + ((uint64_t)op[5] * 4605647320554802783UL) + ((uint64_t)op[6] * 3637739900359287000UL) + ((uint64_t)op[7] * 10466828628376630569UL) + ((uint64_t)op[8] * 1473647133326713317UL) + ((((uint64_t)op[9] * 4270738813664669008UL) + ((uint64_t)op[10] * 11554351881453359700UL)) * 7);
	tmp_q[9] = ((uint64_t)op[0] * 11554351881453359700UL) + ((uint64_t)op[1] * 10731662810726815966UL) + ((uint64_t)op[2] * 14725829993045922307UL) + ((uint64_t)op[3] * 13182533234761706359UL) + ((uint64_t)op[4] * 15285701718870050433UL) + ((uint64_t)op[5] * 13886158810488315615UL) + ((uint64_t)op[6] * 4605647320554802783UL) + ((uint64_t)op[7] * 3637739900359287000UL) + ((uint64_t)op[8] * 10466828628376630569UL) + ((uint64_t)op[9] * 1473647133326713317UL) + ((uint64_t)op[10] * 11448427621943131440UL);
	tmp_q[10] = ((uint64_t)op[0] * 4270738813664669008UL) + ((uint64_t)op[1] * 11554351881453359700UL) + ((uint64_t)op[2] * 10731662810726815966UL) + ((uint64_t)op[3] * 14725829993045922307UL) + ((uint64_t)op[4] * 13182533234761706359UL) + ((uint64_t)op[5] * 15285701718870050433UL) + ((uint64_t)op[6] * 13886158810488315615UL) + ((uint64_t)op[7] * 4605647320554802783UL) + ((uint64_t)op[8] * 3637739900359287000UL) + ((uint64_t)op[9] * 10466828628376630569UL) + ((uint64_t)op[10] * 1473647133326713317UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 91166946050654L) + ((((int128)tmp_q[1] * 51578632348061L) - ((int128)tmp_q[2] * 67991227840900L) - ((int128)tmp_q[3] * 62826941623363L) + ((int128)tmp_q[4] * 67510802185672L) + ((int128)tmp_q[5] * 52828854998924L) + ((int128)tmp_q[6] * 63540702042889L) - ((int128)tmp_q[7] * 18873667409176L) - ((int128)tmp_q[8] * 67417275944514L) - ((int128)tmp_q[9] * 64993861568874L) - ((int128)tmp_q[10] * 12260796939648L)) * 7);
	tmp_zero[1] = -((int128)tmp_q[0] * 12260796939648L) + ((int128)tmp_q[1] * 91166946050654L) + ((((int128)tmp_q[2] * 51578632348061L) - ((int128)tmp_q[3] * 67991227840900L) - ((int128)tmp_q[4] * 62826941623363L) + ((int128)tmp_q[5] * 67510802185672L) + ((int128)tmp_q[6] * 52828854998924L) + ((int128)tmp_q[7] * 63540702042889L) - ((int128)tmp_q[8] * 18873667409176L) - ((int128)tmp_q[9] * 67417275944514L) - ((int128)tmp_q[10] * 64993861568874L)) * 7);
	tmp_zero[2] = -((int128)tmp_q[0] * 64993861568874L) - ((int128)tmp_q[1] * 12260796939648L) + ((int128)tmp_q[2] * 91166946050654L) + ((((int128)tmp_q[3] * 51578632348061L) - ((int128)tmp_q[4] * 67991227840900L) - ((int128)tmp_q[5] * 62826941623363L) + ((int128)tmp_q[6] * 67510802185672L) + ((int128)tmp_q[7] * 52828854998924L) + ((int128)tmp_q[8] * 63540702042889L) - ((int128)tmp_q[9] * 18873667409176L) - ((int128)tmp_q[10] * 67417275944514L)) * 7);
	tmp_zero[3] = -((int128)tmp_q[0] * 67417275944514L) - ((int128)tmp_q[1] * 64993861568874L) - ((int128)tmp_q[2] * 12260796939648L) + ((int128)tmp_q[3] * 91166946050654L) + ((((int128)tmp_q[4] * 51578632348061L) - ((int128)tmp_q[5] * 67991227840900L) - ((int128)tmp_q[6] * 62826941623363L) + ((int128)tmp_q[7] * 67510802185672L) + ((int128)tmp_q[8] * 52828854998924L) + ((int128)tmp_q[9] * 63540702042889L) - ((int128)tmp_q[10] * 18873667409176L)) * 7);
	tmp_zero[4] = -((int128)tmp_q[0] * 18873667409176L) - ((int128)tmp_q[1] * 67417275944514L) - ((int128)tmp_q[2] * 64993861568874L) - ((int128)tmp_q[3] * 12260796939648L) + ((int128)tmp_q[4] * 91166946050654L) + ((((int128)tmp_q[5] * 51578632348061L) - ((int128)tmp_q[6] * 67991227840900L) - ((int128)tmp_q[7] * 62826941623363L) + ((int128)tmp_q[8] * 67510802185672L) + ((int128)tmp_q[9] * 52828854998924L) + ((int128)tmp_q[10] * 63540702042889L)) * 7);
	tmp_zero[5] = ((int128)tmp_q[0] * 63540702042889L) - ((int128)tmp_q[1] * 18873667409176L) - ((int128)tmp_q[2] * 67417275944514L) - ((int128)tmp_q[3] * 64993861568874L) - ((int128)tmp_q[4] * 12260796939648L) + ((int128)tmp_q[5] * 91166946050654L) + ((((int128)tmp_q[6] * 51578632348061L) - ((int128)tmp_q[7] * 67991227840900L) - ((int128)tmp_q[8] * 62826941623363L) + ((int128)tmp_q[9] * 67510802185672L) + ((int128)tmp_q[10] * 52828854998924L)) * 7);
	tmp_zero[6] = ((int128)tmp_q[0] * 52828854998924L) + ((int128)tmp_q[1] * 63540702042889L) - ((int128)tmp_q[2] * 18873667409176L) - ((int128)tmp_q[3] * 67417275944514L) - ((int128)tmp_q[4] * 64993861568874L) - ((int128)tmp_q[5] * 12260796939648L) + ((int128)tmp_q[6] * 91166946050654L) + ((((int128)tmp_q[7] * 51578632348061L) - ((int128)tmp_q[8] * 67991227840900L) - ((int128)tmp_q[9] * 62826941623363L) + ((int128)tmp_q[10] * 67510802185672L)) * 7);
	tmp_zero[7] = ((int128)tmp_q[0] * 67510802185672L) + ((int128)tmp_q[1] * 52828854998924L) + ((int128)tmp_q[2] * 63540702042889L) - ((int128)tmp_q[3] * 18873667409176L) - ((int128)tmp_q[4] * 67417275944514L) - ((int128)tmp_q[5] * 64993861568874L) - ((int128)tmp_q[6] * 12260796939648L) + ((int128)tmp_q[7] * 91166946050654L) + ((((int128)tmp_q[8] * 51578632348061L) - ((int128)tmp_q[9] * 67991227840900L) - ((int128)tmp_q[10] * 62826941623363L)) * 7);
	tmp_zero[8] = -((int128)tmp_q[0] * 62826941623363L) + ((int128)tmp_q[1] * 67510802185672L) + ((int128)tmp_q[2] * 52828854998924L) + ((int128)tmp_q[3] * 63540702042889L) - ((int128)tmp_q[4] * 18873667409176L) - ((int128)tmp_q[5] * 67417275944514L) - ((int128)tmp_q[6] * 64993861568874L) - ((int128)tmp_q[7] * 12260796939648L) + ((int128)tmp_q[8] * 91166946050654L) + ((((int128)tmp_q[9] * 51578632348061L) - ((int128)tmp_q[10] * 67991227840900L)) * 7);
	tmp_zero[9] = -((int128)tmp_q[0] * 67991227840900L) - ((int128)tmp_q[1] * 62826941623363L) + ((int128)tmp_q[2] * 67510802185672L) + ((int128)tmp_q[3] * 52828854998924L) + ((int128)tmp_q[4] * 63540702042889L) - ((int128)tmp_q[5] * 18873667409176L) - ((int128)tmp_q[6] * 67417275944514L) - ((int128)tmp_q[7] * 64993861568874L) - ((int128)tmp_q[8] * 12260796939648L) + ((int128)tmp_q[9] * 91166946050654L) + ((int128)tmp_q[10] * 361050426436427L);
	tmp_zero[10] = ((int128)tmp_q[0] * 51578632348061L) - ((int128)tmp_q[1] * 67991227840900L) - ((int128)tmp_q[2] * 62826941623363L) + ((int128)tmp_q[3] * 67510802185672L) + ((int128)tmp_q[4] * 52828854998924L) + ((int128)tmp_q[5] * 63540702042889L) - ((int128)tmp_q[6] * 18873667409176L) - ((int128)tmp_q[7] * 67417275944514L) - ((int128)tmp_q[8] * 64993861568874L) - ((int128)tmp_q[9] * 12260796939648L) + ((int128)tmp_q[10] * 91166946050654L);

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
	rop[10] = (op[10] + tmp_zero[10]) >> WORD_SIZE;
}

