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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 6);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 6);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 6);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 6);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 6);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 6);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 6);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 6);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 12);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 6);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 12);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 6);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 12);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 6);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 13935517013557060919UL) + ((((uint64_t)op[1] * 3717670403085400850UL) + ((uint64_t)op[2] * 13073438454849504884UL) + ((uint64_t)op[3] * 13748234964203316216UL) + ((uint64_t)op[4] * 995649586476698732UL) + ((uint64_t)op[5] * 15822679364996510287UL) + ((uint64_t)op[6] * 6376517619820191277UL) + ((uint64_t)op[7] * 13854284533629248036UL)) * 18446744073709551610);
	tmp_q[1] = ((uint64_t)op[0] * 13854284533629248036UL) + ((uint64_t)op[1] * 13935517013557060919UL) + ((((uint64_t)op[2] * 3717670403085400850UL) + ((uint64_t)op[3] * 13073438454849504884UL) + ((uint64_t)op[4] * 13748234964203316216UL) + ((uint64_t)op[5] * 995649586476698732UL) + ((uint64_t)op[6] * 15822679364996510287UL) + ((uint64_t)op[7] * 6376517619820191277UL)) * 18446744073709551610);
	tmp_q[2] = ((uint64_t)op[0] * 6376517619820191277UL) + ((uint64_t)op[1] * 13854284533629248036UL) + ((uint64_t)op[2] * 13935517013557060919UL) + ((((uint64_t)op[3] * 3717670403085400850UL) + ((uint64_t)op[4] * 13073438454849504884UL) + ((uint64_t)op[5] * 13748234964203316216UL) + ((uint64_t)op[6] * 995649586476698732UL) + ((uint64_t)op[7] * 15822679364996510287UL)) * 18446744073709551610);
	tmp_q[3] = ((uint64_t)op[0] * 15822679364996510287UL) + ((uint64_t)op[1] * 6376517619820191277UL) + ((uint64_t)op[2] * 13854284533629248036UL) + ((uint64_t)op[3] * 13935517013557060919UL) + ((((uint64_t)op[4] * 3717670403085400850UL) + ((uint64_t)op[5] * 13073438454849504884UL) + ((uint64_t)op[6] * 13748234964203316216UL) + ((uint64_t)op[7] * 995649586476698732UL)) * 18446744073709551610);
	tmp_q[4] = ((uint64_t)op[0] * 995649586476698732UL) + ((uint64_t)op[1] * 15822679364996510287UL) + ((uint64_t)op[2] * 6376517619820191277UL) + ((uint64_t)op[3] * 13854284533629248036UL) + ((uint64_t)op[4] * 13935517013557060919UL) + ((((uint64_t)op[5] * 3717670403085400850UL) + ((uint64_t)op[6] * 13073438454849504884UL) + ((uint64_t)op[7] * 13748234964203316216UL)) * 18446744073709551610);
	tmp_q[5] = ((uint64_t)op[0] * 13748234964203316216UL) + ((uint64_t)op[1] * 995649586476698732UL) + ((uint64_t)op[2] * 15822679364996510287UL) + ((uint64_t)op[3] * 6376517619820191277UL) + ((uint64_t)op[4] * 13854284533629248036UL) + ((uint64_t)op[5] * 13935517013557060919UL) + ((((uint64_t)op[6] * 3717670403085400850UL) + ((uint64_t)op[7] * 13073438454849504884UL)) * 18446744073709551610);
	tmp_q[6] = ((uint64_t)op[0] * 13073438454849504884UL) + ((uint64_t)op[1] * 13748234964203316216UL) + ((uint64_t)op[2] * 995649586476698732UL) + ((uint64_t)op[3] * 15822679364996510287UL) + ((uint64_t)op[4] * 6376517619820191277UL) + ((uint64_t)op[5] * 13854284533629248036UL) + ((uint64_t)op[6] * 13935517013557060919UL) + ((uint64_t)op[7] * 14587465728906698132UL);
	tmp_q[7] = ((uint64_t)op[0] * 3717670403085400850UL) + ((uint64_t)op[1] * 13073438454849504884UL) + ((uint64_t)op[2] * 13748234964203316216UL) + ((uint64_t)op[3] * 995649586476698732UL) + ((uint64_t)op[4] * 15822679364996510287UL) + ((uint64_t)op[5] * 6376517619820191277UL) + ((uint64_t)op[6] * 13854284533629248036UL) + ((uint64_t)op[7] * 13935517013557060919UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 123601098257811L) - ((-((int128)tmp_q[1] * 77340052583213L) - ((int128)tmp_q[2] * 74738944617972L) + ((int128)tmp_q[3] * 101082891351778L) + ((int128)tmp_q[4] * 132582827041893L) - ((int128)tmp_q[5] * 107750945548335L) - ((int128)tmp_q[6] * 14083057406505L) + ((int128)tmp_q[7] * 57910798242646L)) * 6);
	tmp_zero[1] = ((int128)tmp_q[0] * 57910798242646L) - ((int128)tmp_q[1] * 123601098257811L) - ((-((int128)tmp_q[2] * 77340052583213L) - ((int128)tmp_q[3] * 74738944617972L) + ((int128)tmp_q[4] * 101082891351778L) + ((int128)tmp_q[5] * 132582827041893L) - ((int128)tmp_q[6] * 107750945548335L) - ((int128)tmp_q[7] * 14083057406505L)) * 6);
	tmp_zero[2] = -((int128)tmp_q[0] * 14083057406505L) + ((int128)tmp_q[1] * 57910798242646L) - ((int128)tmp_q[2] * 123601098257811L) - ((-((int128)tmp_q[3] * 77340052583213L) - ((int128)tmp_q[4] * 74738944617972L) + ((int128)tmp_q[5] * 101082891351778L) + ((int128)tmp_q[6] * 132582827041893L) - ((int128)tmp_q[7] * 107750945548335L)) * 6);
	tmp_zero[3] = -((int128)tmp_q[0] * 107750945548335L) - ((int128)tmp_q[1] * 14083057406505L) + ((int128)tmp_q[2] * 57910798242646L) - ((int128)tmp_q[3] * 123601098257811L) - ((-((int128)tmp_q[4] * 77340052583213L) - ((int128)tmp_q[5] * 74738944617972L) + ((int128)tmp_q[6] * 101082891351778L) + ((int128)tmp_q[7] * 132582827041893L)) * 6);
	tmp_zero[4] = ((int128)tmp_q[0] * 132582827041893L) - ((int128)tmp_q[1] * 107750945548335L) - ((int128)tmp_q[2] * 14083057406505L) + ((int128)tmp_q[3] * 57910798242646L) - ((int128)tmp_q[4] * 123601098257811L) - ((-((int128)tmp_q[5] * 77340052583213L) - ((int128)tmp_q[6] * 74738944617972L) + ((int128)tmp_q[7] * 101082891351778L)) * 6);
	tmp_zero[5] = ((int128)tmp_q[0] * 101082891351778L) + ((int128)tmp_q[1] * 132582827041893L) - ((int128)tmp_q[2] * 107750945548335L) - ((int128)tmp_q[3] * 14083057406505L) + ((int128)tmp_q[4] * 57910798242646L) - ((int128)tmp_q[5] * 123601098257811L) - ((-((int128)tmp_q[6] * 77340052583213L) - ((int128)tmp_q[7] * 74738944617972L)) * 6);
	tmp_zero[6] = -((int128)tmp_q[0] * 74738944617972L) + ((int128)tmp_q[1] * 101082891351778L) + ((int128)tmp_q[2] * 132582827041893L) - ((int128)tmp_q[3] * 107750945548335L) - ((int128)tmp_q[4] * 14083057406505L) + ((int128)tmp_q[5] * 57910798242646L) - ((int128)tmp_q[6] * 123601098257811L) + ((int128)tmp_q[7] * 464040315499278L);
	tmp_zero[7] = -((int128)tmp_q[0] * 77340052583213L) - ((int128)tmp_q[1] * 74738944617972L) + ((int128)tmp_q[2] * 101082891351778L) + ((int128)tmp_q[3] * 132582827041893L) - ((int128)tmp_q[4] * 107750945548335L) - ((int128)tmp_q[5] * 14083057406505L) + ((int128)tmp_q[6] * 57910798242646L) - ((int128)tmp_q[7] * 123601098257811L);

	//~ computation of : (op + tmp_zero)/mont_phi
	rop[0] = (op[0] + tmp_zero[0]) >> WORD_SIZE;
	rop[1] = (op[1] + tmp_zero[1]) >> WORD_SIZE;
	rop[2] = (op[2] + tmp_zero[2]) >> WORD_SIZE;
	rop[3] = (op[3] + tmp_zero[3]) >> WORD_SIZE;
	rop[4] = (op[4] + tmp_zero[4]) >> WORD_SIZE;
	rop[5] = (op[5] + tmp_zero[5]) >> WORD_SIZE;
	rop[6] = (op[6] + tmp_zero[6]) >> WORD_SIZE;
	rop[7] = (op[7] + tmp_zero[7]) >> WORD_SIZE;
}

