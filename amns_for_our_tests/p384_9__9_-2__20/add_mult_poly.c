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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[8] + (int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2] + (int128)pa[8] * pb[1]) << 1);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[8] + (int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3] + (int128)pa[8] * pb[2]) << 1);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[8] + (int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4] + (int128)pa[8] * pb[3]) << 1);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[8] + (int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5] + (int128)pa[8] * pb[4]) << 1);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[8] + (int128)pa[6] * pb[7] + (int128)pa[7] * pb[6] + (int128)pa[8] * pb[5]) << 1);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[8] + (int128)pa[7] * pb[7] + (int128)pa[8] * pb[6]) << 1);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[8] + (int128)pa[8] * pb[7]) << 1);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0] - (((int128)pa[8] * pb[8]) << 1);
	tmp_prod_result[8] = (int128)pa[0] * pb[8] + (int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1] + (int128)pa[8] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2] + (int128)pa[8] * pa[1]) << 2);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3] + (int128)pa[8] * pa[2]) << 1) + (int128)pa[5] * pa[5]) << 1);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4] + (int128)pa[8] * pa[3]) << 2);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((((int128)pa[7] * pa[5] + (int128)pa[8] * pa[4]) << 1) + (int128)pa[6] * pa[6]) << 1);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((int128)pa[7] * pa[6] + (int128)pa[8] * pa[5]) << 2);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((((int128)pa[8] * pa[6]) << 1) + (int128)pa[7] * pa[7]) << 1);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[8] * pa[7]) << 2);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1) - (((int128)pa[8] * pa[8]) << 1);
	tmp_prod_result[8] = (((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1] + (int128)pa[8] * pa[0]) << 1) + (int128)pa[4] * pa[4];

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3822963429905716575UL) + ((((uint64_t)op[1] * 3644959619242900462UL) + ((uint64_t)op[2] * 7797900946705423642UL) + ((uint64_t)op[3] * 4608394859859829353UL) + ((uint64_t)op[4] * 5384868672479376308UL) + ((uint64_t)op[5] * 13772592247788727259UL) + ((uint64_t)op[6] * 12159577304765828173UL) + ((uint64_t)op[7] * 15017299836659572829UL) + ((uint64_t)op[8] * 15509550373127096701UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 15509550373127096701UL) + ((uint64_t)op[1] * 3822963429905716575UL) + ((((uint64_t)op[2] * 3644959619242900462UL) + ((uint64_t)op[3] * 7797900946705423642UL) + ((uint64_t)op[4] * 4608394859859829353UL) + ((uint64_t)op[5] * 5384868672479376308UL) + ((uint64_t)op[6] * 13772592247788727259UL) + ((uint64_t)op[7] * 12159577304765828173UL) + ((uint64_t)op[8] * 15017299836659572829UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 15017299836659572829UL) + ((uint64_t)op[1] * 15509550373127096701UL) + ((uint64_t)op[2] * 3822963429905716575UL) + ((((uint64_t)op[3] * 3644959619242900462UL) + ((uint64_t)op[4] * 7797900946705423642UL) + ((uint64_t)op[5] * 4608394859859829353UL) + ((uint64_t)op[6] * 5384868672479376308UL) + ((uint64_t)op[7] * 13772592247788727259UL) + ((uint64_t)op[8] * 12159577304765828173UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 12159577304765828173UL) + ((uint64_t)op[1] * 15017299836659572829UL) + ((uint64_t)op[2] * 15509550373127096701UL) + ((uint64_t)op[3] * 3822963429905716575UL) + ((((uint64_t)op[4] * 3644959619242900462UL) + ((uint64_t)op[5] * 7797900946705423642UL) + ((uint64_t)op[6] * 4608394859859829353UL) + ((uint64_t)op[7] * 5384868672479376308UL) + ((uint64_t)op[8] * 13772592247788727259UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 13772592247788727259UL) + ((uint64_t)op[1] * 12159577304765828173UL) + ((uint64_t)op[2] * 15017299836659572829UL) + ((uint64_t)op[3] * 15509550373127096701UL) + ((uint64_t)op[4] * 3822963429905716575UL) + ((((uint64_t)op[5] * 3644959619242900462UL) + ((uint64_t)op[6] * 7797900946705423642UL) + ((uint64_t)op[7] * 4608394859859829353UL) + ((uint64_t)op[8] * 5384868672479376308UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 5384868672479376308UL) + ((uint64_t)op[1] * 13772592247788727259UL) + ((uint64_t)op[2] * 12159577304765828173UL) + ((uint64_t)op[3] * 15017299836659572829UL) + ((uint64_t)op[4] * 15509550373127096701UL) + ((uint64_t)op[5] * 3822963429905716575UL) + ((((uint64_t)op[6] * 3644959619242900462UL) + ((uint64_t)op[7] * 7797900946705423642UL) + ((uint64_t)op[8] * 4608394859859829353UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 4608394859859829353UL) + ((uint64_t)op[1] * 5384868672479376308UL) + ((uint64_t)op[2] * 13772592247788727259UL) + ((uint64_t)op[3] * 12159577304765828173UL) + ((uint64_t)op[4] * 15017299836659572829UL) + ((uint64_t)op[5] * 15509550373127096701UL) + ((uint64_t)op[6] * 3822963429905716575UL) + ((((uint64_t)op[7] * 3644959619242900462UL) + ((uint64_t)op[8] * 7797900946705423642UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 7797900946705423642UL) + ((uint64_t)op[1] * 4608394859859829353UL) + ((uint64_t)op[2] * 5384868672479376308UL) + ((uint64_t)op[3] * 13772592247788727259UL) + ((uint64_t)op[4] * 12159577304765828173UL) + ((uint64_t)op[5] * 15017299836659572829UL) + ((uint64_t)op[6] * 15509550373127096701UL) + ((uint64_t)op[7] * 3822963429905716575UL) + ((uint64_t)op[8] * 11156824835223750692UL);
	tmp_q[8] = ((uint64_t)op[0] * 3644959619242900462UL) + ((uint64_t)op[1] * 7797900946705423642UL) + ((uint64_t)op[2] * 4608394859859829353UL) + ((uint64_t)op[3] * 5384868672479376308UL) + ((uint64_t)op[4] * 13772592247788727259UL) + ((uint64_t)op[5] * 12159577304765828173UL) + ((uint64_t)op[6] * 15017299836659572829UL) + ((uint64_t)op[7] * 15509550373127096701UL) + ((uint64_t)op[8] * 3822963429905716575UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 929026096663L) - ((-((int128)tmp_q[1] * 619939290435L) + ((int128)tmp_q[2] * 4158725733956L) - ((int128)tmp_q[3] * 3091511603724L) - ((int128)tmp_q[4] * 861035561379L) + ((int128)tmp_q[5] * 278838302280L) - ((int128)tmp_q[6] * 1799914358798L) + ((int128)tmp_q[7] * 963825465472L) - ((int128)tmp_q[8] * 3868204846797L)) * 2);
	tmp_zero[1] = -((int128)tmp_q[0] * 3868204846797L) - ((int128)tmp_q[1] * 929026096663L) - ((-((int128)tmp_q[2] * 619939290435L) + ((int128)tmp_q[3] * 4158725733956L) - ((int128)tmp_q[4] * 3091511603724L) - ((int128)tmp_q[5] * 861035561379L) + ((int128)tmp_q[6] * 278838302280L) - ((int128)tmp_q[7] * 1799914358798L) + ((int128)tmp_q[8] * 963825465472L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 963825465472L) - ((int128)tmp_q[1] * 3868204846797L) - ((int128)tmp_q[2] * 929026096663L) - ((-((int128)tmp_q[3] * 619939290435L) + ((int128)tmp_q[4] * 4158725733956L) - ((int128)tmp_q[5] * 3091511603724L) - ((int128)tmp_q[6] * 861035561379L) + ((int128)tmp_q[7] * 278838302280L) - ((int128)tmp_q[8] * 1799914358798L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 1799914358798L) + ((int128)tmp_q[1] * 963825465472L) - ((int128)tmp_q[2] * 3868204846797L) - ((int128)tmp_q[3] * 929026096663L) - ((-((int128)tmp_q[4] * 619939290435L) + ((int128)tmp_q[5] * 4158725733956L) - ((int128)tmp_q[6] * 3091511603724L) - ((int128)tmp_q[7] * 861035561379L) + ((int128)tmp_q[8] * 278838302280L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 278838302280L) - ((int128)tmp_q[1] * 1799914358798L) + ((int128)tmp_q[2] * 963825465472L) - ((int128)tmp_q[3] * 3868204846797L) - ((int128)tmp_q[4] * 929026096663L) - ((-((int128)tmp_q[5] * 619939290435L) + ((int128)tmp_q[6] * 4158725733956L) - ((int128)tmp_q[7] * 3091511603724L) - ((int128)tmp_q[8] * 861035561379L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 861035561379L) + ((int128)tmp_q[1] * 278838302280L) - ((int128)tmp_q[2] * 1799914358798L) + ((int128)tmp_q[3] * 963825465472L) - ((int128)tmp_q[4] * 3868204846797L) - ((int128)tmp_q[5] * 929026096663L) - ((-((int128)tmp_q[6] * 619939290435L) + ((int128)tmp_q[7] * 4158725733956L) - ((int128)tmp_q[8] * 3091511603724L)) * 2);
	tmp_zero[6] = -((int128)tmp_q[0] * 3091511603724L) - ((int128)tmp_q[1] * 861035561379L) + ((int128)tmp_q[2] * 278838302280L) - ((int128)tmp_q[3] * 1799914358798L) + ((int128)tmp_q[4] * 963825465472L) - ((int128)tmp_q[5] * 3868204846797L) - ((int128)tmp_q[6] * 929026096663L) - ((-((int128)tmp_q[7] * 619939290435L) + ((int128)tmp_q[8] * 4158725733956L)) * 2);
	tmp_zero[7] = ((int128)tmp_q[0] * 4158725733956L) - ((int128)tmp_q[1] * 3091511603724L) - ((int128)tmp_q[2] * 861035561379L) + ((int128)tmp_q[3] * 278838302280L) - ((int128)tmp_q[4] * 1799914358798L) + ((int128)tmp_q[5] * 963825465472L) - ((int128)tmp_q[6] * 3868204846797L) - ((int128)tmp_q[7] * 929026096663L) + ((int128)tmp_q[8] * 1239878580870L);
	tmp_zero[8] = -((int128)tmp_q[0] * 619939290435L) + ((int128)tmp_q[1] * 4158725733956L) - ((int128)tmp_q[2] * 3091511603724L) - ((int128)tmp_q[3] * 861035561379L) + ((int128)tmp_q[4] * 278838302280L) - ((int128)tmp_q[5] * 1799914358798L) + ((int128)tmp_q[6] * 963825465472L) - ((int128)tmp_q[7] * 3868204846797L) - ((int128)tmp_q[8] * 929026096663L);

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

