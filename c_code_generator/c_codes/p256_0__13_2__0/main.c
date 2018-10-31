#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <gmp.h>
#include <openssl/bn.h>

#include "structs_data.h"
#include "add_mult_poly.c"
#include "useful_functs.c"
#include "amns_init.c"

#define BILLION 1000000000L


//~ Compilation command : gcc -Wall -O3 main.c -o main -lgmp -lcrypto
//~ Execution command : ./main

//~ Important : polynomials representations form is P(X) = a0 + ... + an.X^n = (a0, ..., an).


int main(void){

	int i, nbiter;
	mpz_t A, B, C, E;
	mpz_inits (A, B, C, E, NULL);

	BN_CTX *ctx = BN_CTX_new();
	BN_CTX_start(ctx);

	BIGNUM *opA = BN_CTX_get(ctx);
	BIGNUM *opB = BN_CTX_get(ctx);
	BIGNUM *opP = BN_CTX_get(ctx);
	BIGNUM *opAA = BN_CTX_get(ctx);
	BIGNUM *opBB = BN_CTX_get(ctx);
	BN_MONT_CTX *mont_ctx = BN_MONT_CTX_new();


	int pa[NB_COEFF];
	int pb[NB_COEFF];

	struct timespec start1, end1;
	struct timespec start2, end2;
	struct timespec start3, end3;
	struct timespec start4, end4;
	struct timespec start5, end5;
	struct timespec start6, end6;
	uint64_t diff1, diff2, diff3, diff4, diff5, diff6;

	unsigned long seed = time(NULL);
	gmp_randstate_t r;
	gmp_randinit_default(r);
	gmp_randseed_ui(r, seed);


	init_data();


	mpz_urandomm(A, r, modul_p);
	mpz_urandomm(B, r, modul_p);
	mpz_set(E, A);

	BN_dec2bn(&opA, mpz_get_str (NULL, 10, A));
	BN_dec2bn(&opB, mpz_get_str (NULL, 10, B));
	BN_dec2bn(&opP, mpz_get_str (NULL, 10, modul_p));
	BN_copy(opAA, opA);
	BN_copy(opBB, opB);
	BN_MONT_CTX_set(mont_ctx, opP, ctx);


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start3);
	from_int_to_amns(pa, A);
	from_int_to_amns(pb, B);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end3);
	diff3 = BILLION * (end3.tv_sec - start3.tv_sec) + (end3.tv_nsec - start3.tv_nsec);

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start5);
	BN_to_montgomery(opA, opA, mont_ctx, ctx);
	BN_to_montgomery(opB, opB, mont_ctx, ctx);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end5);
	diff5 = BILLION * (end5.tv_sec - start5.tv_sec) + (end5.tv_nsec - start5.tv_nsec);


	nbiter = 1 << 25;


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start1);
	for (i=0; i<nbiter; i++) {
		mpz_mul (A, A, B);
		mpz_mod (A, A, modul_p);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end1);
	diff1 = BILLION * (end1.tv_sec - start1.tv_sec) + (end1.tv_nsec - start1.tv_nsec);


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start2);
	for (i=0; i<nbiter; i++) {
		mult_mod_poly(pa, pa, pb);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end2);
	diff2 = BILLION * (end2.tv_sec - start2.tv_sec) + (end2.tv_nsec - start2.tv_nsec);


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start4);
	for (i=0; i<nbiter; i++) {
		BN_mod_mul_montgomery(opA, opA, opB, mont_ctx, ctx);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end4);
	diff4 = BILLION * (end4.tv_sec - start4.tv_sec) + (end4.tv_nsec - start4.tv_nsec);


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start6);
	for (i=0; i<nbiter; i++) {
		BN_mod_mul(opAA, opAA, opBB, opP, ctx);
	}
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end6);
	diff6 = BILLION * (end6.tv_sec - start6.tv_sec) + (end6.tv_nsec - start6.tv_nsec);


	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start3);
	from_amns_to_int(C, pa);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end3);
	diff3 += BILLION * (end3.tv_sec - start3.tv_sec) + (end3.tv_nsec - start3.tv_nsec);

	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start5);
	BN_from_montgomery(opA, opA, mont_ctx, ctx);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end5);
	diff5 += BILLION * (end5.tv_sec - start5.tv_sec) + (end5.tv_nsec - start5.tv_nsec);


	printf("%lu %d %lu %lu %lu %lu\n\n", mpz_sizeinbase (modul_p, 2), NB_COEFF, (diff2+diff3), (diff4+diff5), diff6, diff1);

	printf("nbiter = %d\n\n", nbiter);
	gmp_printf("A       : %Zd\n", E);
	gmp_printf("B       : %Zd\n\n", B);
	gmp_printf("r_gmp   : %Zd\n", A);
	gmp_printf("r_amns  : %Zd\n", C);
	gmp_printf("r_ssld  : %s\n", BN_bn2dec(opAA));
	gmp_printf("r_sslm  : %s\n", BN_bn2dec(opA));
	printf("\ntime using gmp			= %lu nanoseconds\n", diff1);
	printf("\ntime using openssl		= %lu nanoseconds\n", diff6);
	printf("\ntime using amns prod		= %lu nanoseconds\n", diff2);
	printf("total time using amns prod	= %lu nanoseconds\n", (diff2+diff3));
	printf("\ntime using openssl mont		= %lu nanoseconds\n", diff4);
	printf("total time using openssl mont	= %lu nanoseconds\n", (diff4+diff5));


	mpz_clears (A, B, C, E, NULL);
	gmp_randclear(r);

	BN_MONT_CTX_free(mont_ctx);
	BN_CTX_end(ctx);
	BN_CTX_free(ctx);

	free_data();
	return 0;
}

