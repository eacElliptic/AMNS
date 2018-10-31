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
	tmp_q[0] = ((uint64_t)op[0] * 15401694195148755465UL) + ((((uint64_t)op[1] * 12809721689538425942UL) + ((uint64_t)op[2] * 1165109566965045184UL) + ((uint64_t)op[3] * 8959605271747099477UL) + ((uint64_t)op[4] * 9143314483462275409UL) + ((uint64_t)op[5] * 13306355277525666595UL) + ((uint64_t)op[6] * 17187645934115330343UL) + ((uint64_t)op[7] * 6142898538869318662UL) + ((uint64_t)op[8] * 12542054912784686891UL) + ((uint64_t)op[9] * 18105265991440669982UL)) * 18446744073709551614);
	tmp_q[1] = ((uint64_t)op[0] * 18105265991440669982UL) + ((uint64_t)op[1] * 15401694195148755465UL) + ((((uint64_t)op[2] * 12809721689538425942UL) + ((uint64_t)op[3] * 1165109566965045184UL) + ((uint64_t)op[4] * 8959605271747099477UL) + ((uint64_t)op[5] * 9143314483462275409UL) + ((uint64_t)op[6] * 13306355277525666595UL) + ((uint64_t)op[7] * 17187645934115330343UL) + ((uint64_t)op[8] * 6142898538869318662UL) + ((uint64_t)op[9] * 12542054912784686891UL)) * 18446744073709551614);
	tmp_q[2] = ((uint64_t)op[0] * 12542054912784686891UL) + ((uint64_t)op[1] * 18105265991440669982UL) + ((uint64_t)op[2] * 15401694195148755465UL) + ((((uint64_t)op[3] * 12809721689538425942UL) + ((uint64_t)op[4] * 1165109566965045184UL) + ((uint64_t)op[5] * 8959605271747099477UL) + ((uint64_t)op[6] * 9143314483462275409UL) + ((uint64_t)op[7] * 13306355277525666595UL) + ((uint64_t)op[8] * 17187645934115330343UL) + ((uint64_t)op[9] * 6142898538869318662UL)) * 18446744073709551614);
	tmp_q[3] = ((uint64_t)op[0] * 6142898538869318662UL) + ((uint64_t)op[1] * 12542054912784686891UL) + ((uint64_t)op[2] * 18105265991440669982UL) + ((uint64_t)op[3] * 15401694195148755465UL) + ((((uint64_t)op[4] * 12809721689538425942UL) + ((uint64_t)op[5] * 1165109566965045184UL) + ((uint64_t)op[6] * 8959605271747099477UL) + ((uint64_t)op[7] * 9143314483462275409UL) + ((uint64_t)op[8] * 13306355277525666595UL) + ((uint64_t)op[9] * 17187645934115330343UL)) * 18446744073709551614);
	tmp_q[4] = ((uint64_t)op[0] * 17187645934115330343UL) + ((uint64_t)op[1] * 6142898538869318662UL) + ((uint64_t)op[2] * 12542054912784686891UL) + ((uint64_t)op[3] * 18105265991440669982UL) + ((uint64_t)op[4] * 15401694195148755465UL) + ((((uint64_t)op[5] * 12809721689538425942UL) + ((uint64_t)op[6] * 1165109566965045184UL) + ((uint64_t)op[7] * 8959605271747099477UL) + ((uint64_t)op[8] * 9143314483462275409UL) + ((uint64_t)op[9] * 13306355277525666595UL)) * 18446744073709551614);
	tmp_q[5] = ((uint64_t)op[0] * 13306355277525666595UL) + ((uint64_t)op[1] * 17187645934115330343UL) + ((uint64_t)op[2] * 6142898538869318662UL) + ((uint64_t)op[3] * 12542054912784686891UL) + ((uint64_t)op[4] * 18105265991440669982UL) + ((uint64_t)op[5] * 15401694195148755465UL) + ((((uint64_t)op[6] * 12809721689538425942UL) + ((uint64_t)op[7] * 1165109566965045184UL) + ((uint64_t)op[8] * 8959605271747099477UL) + ((uint64_t)op[9] * 9143314483462275409UL)) * 18446744073709551614);
	tmp_q[6] = ((uint64_t)op[0] * 9143314483462275409UL) + ((uint64_t)op[1] * 13306355277525666595UL) + ((uint64_t)op[2] * 17187645934115330343UL) + ((uint64_t)op[3] * 6142898538869318662UL) + ((uint64_t)op[4] * 12542054912784686891UL) + ((uint64_t)op[5] * 18105265991440669982UL) + ((uint64_t)op[6] * 15401694195148755465UL) + ((((uint64_t)op[7] * 12809721689538425942UL) + ((uint64_t)op[8] * 1165109566965045184UL) + ((uint64_t)op[9] * 8959605271747099477UL)) * 18446744073709551614);
	tmp_q[7] = ((uint64_t)op[0] * 8959605271747099477UL) + ((uint64_t)op[1] * 9143314483462275409UL) + ((uint64_t)op[2] * 13306355277525666595UL) + ((uint64_t)op[3] * 17187645934115330343UL) + ((uint64_t)op[4] * 6142898538869318662UL) + ((uint64_t)op[5] * 12542054912784686891UL) + ((uint64_t)op[6] * 18105265991440669982UL) + ((uint64_t)op[7] * 15401694195148755465UL) + ((((uint64_t)op[8] * 12809721689538425942UL) + ((uint64_t)op[9] * 1165109566965045184UL)) * 18446744073709551614);
	tmp_q[8] = ((uint64_t)op[0] * 1165109566965045184UL) + ((uint64_t)op[1] * 8959605271747099477UL) + ((uint64_t)op[2] * 9143314483462275409UL) + ((uint64_t)op[3] * 13306355277525666595UL) + ((uint64_t)op[4] * 17187645934115330343UL) + ((uint64_t)op[5] * 6142898538869318662UL) + ((uint64_t)op[6] * 12542054912784686891UL) + ((uint64_t)op[7] * 18105265991440669982UL) + ((uint64_t)op[8] * 15401694195148755465UL) + ((uint64_t)op[9] * 11274044768342251348UL);
	tmp_q[9] = ((uint64_t)op[0] * 12809721689538425942UL) + ((uint64_t)op[1] * 1165109566965045184UL) + ((uint64_t)op[2] * 8959605271747099477UL) + ((uint64_t)op[3] * 9143314483462275409UL) + ((uint64_t)op[4] * 13306355277525666595UL) + ((uint64_t)op[5] * 17187645934115330343UL) + ((uint64_t)op[6] * 6142898538869318662UL) + ((uint64_t)op[7] * 12542054912784686891UL) + ((uint64_t)op[8] * 18105265991440669982UL) + ((uint64_t)op[9] * 15401694195148755465UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = -((int128)tmp_q[0] * 694342057148613L) - ((-((int128)tmp_q[1] * 1926707237202823L) - ((int128)tmp_q[2] * 741020355415697L) - ((int128)tmp_q[3] * 1337485991562969L) + ((int128)tmp_q[4] * 1093783566106096L) - ((int128)tmp_q[5] * 1814159096194085L) + ((int128)tmp_q[6] * 795401488945070L) - ((int128)tmp_q[7] * 4832493406492634L) + ((int128)tmp_q[8] * 992155818376525L) + ((int128)tmp_q[9] * 1202875037979088L)) * 2);
	tmp_zero[1] = ((int128)tmp_q[0] * 1202875037979088L) - ((int128)tmp_q[1] * 694342057148613L) - ((-((int128)tmp_q[2] * 1926707237202823L) - ((int128)tmp_q[3] * 741020355415697L) - ((int128)tmp_q[4] * 1337485991562969L) + ((int128)tmp_q[5] * 1093783566106096L) - ((int128)tmp_q[6] * 1814159096194085L) + ((int128)tmp_q[7] * 795401488945070L) - ((int128)tmp_q[8] * 4832493406492634L) + ((int128)tmp_q[9] * 992155818376525L)) * 2);
	tmp_zero[2] = ((int128)tmp_q[0] * 992155818376525L) + ((int128)tmp_q[1] * 1202875037979088L) - ((int128)tmp_q[2] * 694342057148613L) - ((-((int128)tmp_q[3] * 1926707237202823L) - ((int128)tmp_q[4] * 741020355415697L) - ((int128)tmp_q[5] * 1337485991562969L) + ((int128)tmp_q[6] * 1093783566106096L) - ((int128)tmp_q[7] * 1814159096194085L) + ((int128)tmp_q[8] * 795401488945070L) - ((int128)tmp_q[9] * 4832493406492634L)) * 2);
	tmp_zero[3] = -((int128)tmp_q[0] * 4832493406492634L) + ((int128)tmp_q[1] * 992155818376525L) + ((int128)tmp_q[2] * 1202875037979088L) - ((int128)tmp_q[3] * 694342057148613L) - ((-((int128)tmp_q[4] * 1926707237202823L) - ((int128)tmp_q[5] * 741020355415697L) - ((int128)tmp_q[6] * 1337485991562969L) + ((int128)tmp_q[7] * 1093783566106096L) - ((int128)tmp_q[8] * 1814159096194085L) + ((int128)tmp_q[9] * 795401488945070L)) * 2);
	tmp_zero[4] = ((int128)tmp_q[0] * 795401488945070L) - ((int128)tmp_q[1] * 4832493406492634L) + ((int128)tmp_q[2] * 992155818376525L) + ((int128)tmp_q[3] * 1202875037979088L) - ((int128)tmp_q[4] * 694342057148613L) - ((-((int128)tmp_q[5] * 1926707237202823L) - ((int128)tmp_q[6] * 741020355415697L) - ((int128)tmp_q[7] * 1337485991562969L) + ((int128)tmp_q[8] * 1093783566106096L) - ((int128)tmp_q[9] * 1814159096194085L)) * 2);
	tmp_zero[5] = -((int128)tmp_q[0] * 1814159096194085L) + ((int128)tmp_q[1] * 795401488945070L) - ((int128)tmp_q[2] * 4832493406492634L) + ((int128)tmp_q[3] * 992155818376525L) + ((int128)tmp_q[4] * 1202875037979088L) - ((int128)tmp_q[5] * 694342057148613L) - ((-((int128)tmp_q[6] * 1926707237202823L) - ((int128)tmp_q[7] * 741020355415697L) - ((int128)tmp_q[8] * 1337485991562969L) + ((int128)tmp_q[9] * 1093783566106096L)) * 2);
	tmp_zero[6] = ((int128)tmp_q[0] * 1093783566106096L) - ((int128)tmp_q[1] * 1814159096194085L) + ((int128)tmp_q[2] * 795401488945070L) - ((int128)tmp_q[3] * 4832493406492634L) + ((int128)tmp_q[4] * 992155818376525L) + ((int128)tmp_q[5] * 1202875037979088L) - ((int128)tmp_q[6] * 694342057148613L) - ((-((int128)tmp_q[7] * 1926707237202823L) - ((int128)tmp_q[8] * 741020355415697L) - ((int128)tmp_q[9] * 1337485991562969L)) * 2);
	tmp_zero[7] = -((int128)tmp_q[0] * 1337485991562969L) + ((int128)tmp_q[1] * 1093783566106096L) - ((int128)tmp_q[2] * 1814159096194085L) + ((int128)tmp_q[3] * 795401488945070L) - ((int128)tmp_q[4] * 4832493406492634L) + ((int128)tmp_q[5] * 992155818376525L) + ((int128)tmp_q[6] * 1202875037979088L) - ((int128)tmp_q[7] * 694342057148613L) - ((-((int128)tmp_q[8] * 1926707237202823L) - ((int128)tmp_q[9] * 741020355415697L)) * 2);
	tmp_zero[8] = -((int128)tmp_q[0] * 741020355415697L) - ((int128)tmp_q[1] * 1337485991562969L) + ((int128)tmp_q[2] * 1093783566106096L) - ((int128)tmp_q[3] * 1814159096194085L) + ((int128)tmp_q[4] * 795401488945070L) - ((int128)tmp_q[5] * 4832493406492634L) + ((int128)tmp_q[6] * 992155818376525L) + ((int128)tmp_q[7] * 1202875037979088L) - ((int128)tmp_q[8] * 694342057148613L) + ((int128)tmp_q[9] * 3853414474405646L);
	tmp_zero[9] = -((int128)tmp_q[0] * 1926707237202823L) - ((int128)tmp_q[1] * 741020355415697L) - ((int128)tmp_q[2] * 1337485991562969L) + ((int128)tmp_q[3] * 1093783566106096L) - ((int128)tmp_q[4] * 1814159096194085L) + ((int128)tmp_q[5] * 795401488945070L) - ((int128)tmp_q[6] * 4832493406492634L) + ((int128)tmp_q[7] * 992155818376525L) + ((int128)tmp_q[8] * 1202875037979088L) - ((int128)tmp_q[9] * 694342057148613L);

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

