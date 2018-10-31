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

	tmp_prod_result[0] = (int128)pa[0] * pb[0] - (((int128)pa[1] * pb[7] + (int128)pa[2] * pb[6] + (int128)pa[3] * pb[5] + (int128)pa[4] * pb[4] + (int128)pa[5] * pb[3] + (int128)pa[6] * pb[2] + (int128)pa[7] * pb[1]) * 3);
	tmp_prod_result[1] = (int128)pa[0] * pb[1] + (int128)pa[1] * pb[0] - (((int128)pa[2] * pb[7] + (int128)pa[3] * pb[6] + (int128)pa[4] * pb[5] + (int128)pa[5] * pb[4] + (int128)pa[6] * pb[3] + (int128)pa[7] * pb[2]) * 3);
	tmp_prod_result[2] = (int128)pa[0] * pb[2] + (int128)pa[1] * pb[1] + (int128)pa[2] * pb[0] - (((int128)pa[3] * pb[7] + (int128)pa[4] * pb[6] + (int128)pa[5] * pb[5] + (int128)pa[6] * pb[4] + (int128)pa[7] * pb[3]) * 3);
	tmp_prod_result[3] = (int128)pa[0] * pb[3] + (int128)pa[1] * pb[2] + (int128)pa[2] * pb[1] + (int128)pa[3] * pb[0] - (((int128)pa[4] * pb[7] + (int128)pa[5] * pb[6] + (int128)pa[6] * pb[5] + (int128)pa[7] * pb[4]) * 3);
	tmp_prod_result[4] = (int128)pa[0] * pb[4] + (int128)pa[1] * pb[3] + (int128)pa[2] * pb[2] + (int128)pa[3] * pb[1] + (int128)pa[4] * pb[0] - (((int128)pa[5] * pb[7] + (int128)pa[6] * pb[6] + (int128)pa[7] * pb[5]) * 3);
	tmp_prod_result[5] = (int128)pa[0] * pb[5] + (int128)pa[1] * pb[4] + (int128)pa[2] * pb[3] + (int128)pa[3] * pb[2] + (int128)pa[4] * pb[1] + (int128)pa[5] * pb[0] - (((int128)pa[6] * pb[7] + (int128)pa[7] * pb[6]) * 3);
	tmp_prod_result[6] = (int128)pa[0] * pb[6] + (int128)pa[1] * pb[5] + (int128)pa[2] * pb[4] + (int128)pa[3] * pb[3] + (int128)pa[4] * pb[2] + (int128)pa[5] * pb[1] + (int128)pa[6] * pb[0] - (((int128)pa[7] * pb[7]) * 3);
	tmp_prod_result[7] = (int128)pa[0] * pb[7] + (int128)pa[1] * pb[6] + (int128)pa[2] * pb[5] + (int128)pa[3] * pb[4] + (int128)pa[4] * pb[3] + (int128)pa[5] * pb[2] + (int128)pa[6] * pb[1] + (int128)pa[7] * pb[0];

	internal_reduction(rop, tmp_prod_result);
}

//~ Computes pa(X)^2 mod(X^n - c)
void square_mod_poly(int64_t *rop, int64_t *pa){

	int128 tmp_prod_result[NB_COEFF];

	tmp_prod_result[0] = (int128)pa[0] * pa[0] - (((((int128)pa[5] * pa[3] + (int128)pa[6] * pa[2] + (int128)pa[7] * pa[1]) << 1) + (int128)pa[4] * pa[4]) * 3);
	tmp_prod_result[1] = (((int128)pa[1] * pa[0]) << 1) - (((int128)pa[5] * pa[4] + (int128)pa[6] * pa[3] + (int128)pa[7] * pa[2]) * 6);
	tmp_prod_result[2] = (((int128)pa[2] * pa[0]) << 1) + (int128)pa[1] * pa[1] - (((((int128)pa[6] * pa[4] + (int128)pa[7] * pa[3]) << 1) + (int128)pa[5] * pa[5]) * 3);
	tmp_prod_result[3] = (((int128)pa[2] * pa[1] + (int128)pa[3] * pa[0]) << 1) - (((int128)pa[6] * pa[5] + (int128)pa[7] * pa[4]) * 6);
	tmp_prod_result[4] = (((int128)pa[3] * pa[1] + (int128)pa[4] * pa[0]) << 1) + (int128)pa[2] * pa[2] - (((((int128)pa[7] * pa[5]) << 1) + (int128)pa[6] * pa[6]) * 3);
	tmp_prod_result[5] = (((int128)pa[3] * pa[2] + (int128)pa[4] * pa[1] + (int128)pa[5] * pa[0]) << 1) - (((int128)pa[7] * pa[6]) * 6);
	tmp_prod_result[6] = (((int128)pa[4] * pa[2] + (int128)pa[5] * pa[1] + (int128)pa[6] * pa[0]) << 1) + (int128)pa[3] * pa[3] - (((int128)pa[7] * pa[7]) * 3);
	tmp_prod_result[7] = (((int128)pa[4] * pa[3] + (int128)pa[5] * pa[2] + (int128)pa[6] * pa[1] + (int128)pa[7] * pa[0]) << 1);

	internal_reduction(rop, tmp_prod_result);
}

//~ performs the internal reduction on 'op' and puts the result in 'rop'
//~ IMPORTANT : We take 'mont_phi = 1 << WORD_SIZE', so operations modulo mont_phi are automatically done using the appropriate variable type.
void internal_reduction(int64_t *rop, int128 *op){

	uint64_t tmp_q[NB_COEFF];
	int128 tmp_zero[NB_COEFF];

	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)
	tmp_q[0] = ((uint64_t)op[0] * 3167336677092148830UL) + ((((uint64_t)op[1] * 17205756469588458269UL) + ((uint64_t)op[2] * 5131112348993135182UL) + ((uint64_t)op[3] * 10270546136756576274UL) + ((uint64_t)op[4] * 14495888836521561415UL) + ((uint64_t)op[5] * 3479328476312549859UL) + ((uint64_t)op[6] * 17675662722558573465UL) + ((uint64_t)op[7] * 110619624264431647UL)) * 18446744073709551613);
	tmp_q[1] = ((uint64_t)op[0] * 110619624264431647UL) + ((uint64_t)op[1] * 3167336677092148830UL) + ((((uint64_t)op[2] * 17205756469588458269UL) + ((uint64_t)op[3] * 5131112348993135182UL) + ((uint64_t)op[4] * 10270546136756576274UL) + ((uint64_t)op[5] * 14495888836521561415UL) + ((uint64_t)op[6] * 3479328476312549859UL) + ((uint64_t)op[7] * 17675662722558573465UL)) * 18446744073709551613);
	tmp_q[2] = ((uint64_t)op[0] * 17675662722558573465UL) + ((uint64_t)op[1] * 110619624264431647UL) + ((uint64_t)op[2] * 3167336677092148830UL) + ((((uint64_t)op[3] * 17205756469588458269UL) + ((uint64_t)op[4] * 5131112348993135182UL) + ((uint64_t)op[5] * 10270546136756576274UL) + ((uint64_t)op[6] * 14495888836521561415UL) + ((uint64_t)op[7] * 3479328476312549859UL)) * 18446744073709551613);
	tmp_q[3] = ((uint64_t)op[0] * 3479328476312549859UL) + ((uint64_t)op[1] * 17675662722558573465UL) + ((uint64_t)op[2] * 110619624264431647UL) + ((uint64_t)op[3] * 3167336677092148830UL) + ((((uint64_t)op[4] * 17205756469588458269UL) + ((uint64_t)op[5] * 5131112348993135182UL) + ((uint64_t)op[6] * 10270546136756576274UL) + ((uint64_t)op[7] * 14495888836521561415UL)) * 18446744073709551613);
	tmp_q[4] = ((uint64_t)op[0] * 14495888836521561415UL) + ((uint64_t)op[1] * 3479328476312549859UL) + ((uint64_t)op[2] * 17675662722558573465UL) + ((uint64_t)op[3] * 110619624264431647UL) + ((uint64_t)op[4] * 3167336677092148830UL) + ((((uint64_t)op[5] * 17205756469588458269UL) + ((uint64_t)op[6] * 5131112348993135182UL) + ((uint64_t)op[7] * 10270546136756576274UL)) * 18446744073709551613);
	tmp_q[5] = ((uint64_t)op[0] * 10270546136756576274UL) + ((uint64_t)op[1] * 14495888836521561415UL) + ((uint64_t)op[2] * 3479328476312549859UL) + ((uint64_t)op[3] * 17675662722558573465UL) + ((uint64_t)op[4] * 110619624264431647UL) + ((uint64_t)op[5] * 3167336677092148830UL) + ((((uint64_t)op[6] * 17205756469588458269UL) + ((uint64_t)op[7] * 5131112348993135182UL)) * 18446744073709551613);
	tmp_q[6] = ((uint64_t)op[0] * 5131112348993135182UL) + ((uint64_t)op[1] * 10270546136756576274UL) + ((uint64_t)op[2] * 14495888836521561415UL) + ((uint64_t)op[3] * 3479328476312549859UL) + ((uint64_t)op[4] * 17675662722558573465UL) + ((uint64_t)op[5] * 110619624264431647UL) + ((uint64_t)op[6] * 3167336677092148830UL) + ((uint64_t)op[7] * 3722962812363280041UL);
	tmp_q[7] = ((uint64_t)op[0] * 17205756469588458269UL) + ((uint64_t)op[1] * 5131112348993135182UL) + ((uint64_t)op[2] * 10270546136756576274UL) + ((uint64_t)op[3] * 14495888836521561415UL) + ((uint64_t)op[4] * 3479328476312549859UL) + ((uint64_t)op[5] * 17675662722558573465UL) + ((uint64_t)op[6] * 110619624264431647UL) + ((uint64_t)op[7] * 3167336677092148830UL);

	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)
	tmp_zero[0] = ((int128)tmp_q[0] * 47621296962404L) - ((-((int128)tmp_q[1] * 154916624563747L) - ((int128)tmp_q[2] * 179779124883783L) - ((int128)tmp_q[3] * 9011440055382L) + ((int128)tmp_q[4] * 48616209220755L) - ((int128)tmp_q[5] * 12946466528100L) + ((int128)tmp_q[6] * 47385891966744L) + ((int128)tmp_q[7] * 113910925639684L)) * 3);
	tmp_zero[1] = ((int128)tmp_q[0] * 113910925639684L) + ((int128)tmp_q[1] * 47621296962404L) - ((-((int128)tmp_q[2] * 154916624563747L) - ((int128)tmp_q[3] * 179779124883783L) - ((int128)tmp_q[4] * 9011440055382L) + ((int128)tmp_q[5] * 48616209220755L) - ((int128)tmp_q[6] * 12946466528100L) + ((int128)tmp_q[7] * 47385891966744L)) * 3);
	tmp_zero[2] = ((int128)tmp_q[0] * 47385891966744L) + ((int128)tmp_q[1] * 113910925639684L) + ((int128)tmp_q[2] * 47621296962404L) - ((-((int128)tmp_q[3] * 154916624563747L) - ((int128)tmp_q[4] * 179779124883783L) - ((int128)tmp_q[5] * 9011440055382L) + ((int128)tmp_q[6] * 48616209220755L) - ((int128)tmp_q[7] * 12946466528100L)) * 3);
	tmp_zero[3] = -((int128)tmp_q[0] * 12946466528100L) + ((int128)tmp_q[1] * 47385891966744L) + ((int128)tmp_q[2] * 113910925639684L) + ((int128)tmp_q[3] * 47621296962404L) - ((-((int128)tmp_q[4] * 154916624563747L) - ((int128)tmp_q[5] * 179779124883783L) - ((int128)tmp_q[6] * 9011440055382L) + ((int128)tmp_q[7] * 48616209220755L)) * 3);
	tmp_zero[4] = ((int128)tmp_q[0] * 48616209220755L) - ((int128)tmp_q[1] * 12946466528100L) + ((int128)tmp_q[2] * 47385891966744L) + ((int128)tmp_q[3] * 113910925639684L) + ((int128)tmp_q[4] * 47621296962404L) - ((-((int128)tmp_q[5] * 154916624563747L) - ((int128)tmp_q[6] * 179779124883783L) - ((int128)tmp_q[7] * 9011440055382L)) * 3);
	tmp_zero[5] = -((int128)tmp_q[0] * 9011440055382L) + ((int128)tmp_q[1] * 48616209220755L) - ((int128)tmp_q[2] * 12946466528100L) + ((int128)tmp_q[3] * 47385891966744L) + ((int128)tmp_q[4] * 113910925639684L) + ((int128)tmp_q[5] * 47621296962404L) - ((-((int128)tmp_q[6] * 154916624563747L) - ((int128)tmp_q[7] * 179779124883783L)) * 3);
	tmp_zero[6] = -((int128)tmp_q[0] * 179779124883783L) - ((int128)tmp_q[1] * 9011440055382L) + ((int128)tmp_q[2] * 48616209220755L) - ((int128)tmp_q[3] * 12946466528100L) + ((int128)tmp_q[4] * 47385891966744L) + ((int128)tmp_q[5] * 113910925639684L) + ((int128)tmp_q[6] * 47621296962404L) + ((int128)tmp_q[7] * 464749873691241L);
	tmp_zero[7] = -((int128)tmp_q[0] * 154916624563747L) - ((int128)tmp_q[1] * 179779124883783L) - ((int128)tmp_q[2] * 9011440055382L) + ((int128)tmp_q[3] * 48616209220755L) - ((int128)tmp_q[4] * 12946466528100L) + ((int128)tmp_q[5] * 47385891966744L) + ((int128)tmp_q[6] * 113910925639684L) + ((int128)tmp_q[7] * 47621296962404L);

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

