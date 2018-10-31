
def build_main_file(dir_path, small_int):
	with open(dir_path+"/main.c", "w") as f:
		
		f.write("#include <stdlib.h>\n")
		f.write("#include <stdio.h>\n")
		f.write("#include <stdint.h>\n")
		f.write("#include <time.h>\n")
		f.write("#include <gmp.h>\n")
		f.write("#include <openssl/bn.h>\n\n")
		
		f.write("#include \"structs_data.h\"\n")
		f.write("#include \"add_mult_poly.c\"\n")
		f.write("#include \"useful_functs.c\"\n")
		f.write("#include \"amns_init.c\"\n\n")
		f.write("#define BILLION 1000000000L\n\n\n")
		
		f.write("//~ Compilation command : gcc -Wall -O3 main.c -o main -lgmp -lcrypto\n")
		f.write("//~ Execution command : ./main\n\n")
		
		f.write("//~ Important : polynomials representations form is P(X) = a0 + ... + an.X^n = (a0, ..., an).\n\n\n")
		
		f.write("int main(void){\n\n")
		f.write("	int i, nbiter;\n")
		f.write("	mpz_t A, B, C, E;\n")
		f.write("	mpz_inits (A, B, C, E, NULL);\n\n")
		
		f.write("	BN_CTX *ctx = BN_CTX_new();\n")
		f.write("	BN_CTX_start(ctx);\n\n")
		f.write("	BIGNUM *opA = BN_CTX_get(ctx);\n")
		f.write("	BIGNUM *opB = BN_CTX_get(ctx);\n")
		f.write("	BIGNUM *opP = BN_CTX_get(ctx);\n")
		f.write("	BIGNUM *opAA = BN_CTX_get(ctx);\n")
		f.write("	BIGNUM *opBB = BN_CTX_get(ctx);\n")
		f.write("	BN_MONT_CTX *mont_ctx = BN_MONT_CTX_new();\n\n\n")
		
		f.write("	" + small_int + " pa[NB_COEFF];\n")
		f.write("	" + small_int + " pb[NB_COEFF];\n\n")
		
		f.write("	struct timespec start1, end1;\n")
		f.write("	struct timespec start2, end2;\n")
		f.write("	struct timespec start3, end3;\n")
		f.write("	struct timespec start4, end4;\n")
		f.write("	struct timespec start5, end5;\n")
		f.write("	struct timespec start6, end6;\n")
		f.write("	uint64_t diff1, diff2, diff3, diff4, diff5, diff6;\n\n")
		
		f.write("	unsigned long seed = time(NULL);\n")
		f.write("	gmp_randstate_t r;\n")
		f.write("	gmp_randinit_default(r);\n")
		f.write("	gmp_randseed_ui(r, seed);\n\n\n")
		
		f.write("	init_data();\n\n\n")
		
		f.write("	mpz_urandomm(A, r, modul_p);\n")
		f.write("	mpz_urandomm(B, r, modul_p);\n")
		f.write("	mpz_set(E, A);\n\n")
		
		f.write("	BN_dec2bn(&opA, mpz_get_str (NULL, 10, A));\n")
		f.write("	BN_dec2bn(&opB, mpz_get_str (NULL, 10, B));\n")
		f.write("	BN_dec2bn(&opP, mpz_get_str (NULL, 10, modul_p));\n")
		f.write("	BN_copy(opAA, opA);\n")
		f.write("	BN_copy(opBB, opB);\n")
		f.write("	BN_MONT_CTX_set(mont_ctx, opP, ctx);\n\n\n")
		
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start3);\n")
		f.write("	from_int_to_amns(pa, A);\n")
		f.write("	from_int_to_amns(pb, B);\n")
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end3);\n")
		f.write("	diff3 = BILLION * (end3.tv_sec - start3.tv_sec) + (end3.tv_nsec - start3.tv_nsec);\n\n")
		
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start5);\n")
		f.write("	BN_to_montgomery(opA, opA, mont_ctx, ctx);\n")
		f.write("	BN_to_montgomery(opB, opB, mont_ctx, ctx);\n")
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end5);\n")
		f.write("	diff5 = BILLION * (end5.tv_sec - start5.tv_sec) + (end5.tv_nsec - start5.tv_nsec);\n\n\n")

		f.write("	nbiter = 1 << 25;\n\n\n")

		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start1);\n")
		f.write("	for (i=0; i<nbiter; i++) {\n")
		f.write("		mpz_mul (A, A, B);\n")
		f.write("		mpz_mod (A, A, modul_p);\n")
		f.write("	}\n")
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end1);\n")
		f.write("	diff1 = BILLION * (end1.tv_sec - start1.tv_sec) + (end1.tv_nsec - start1.tv_nsec);\n\n\n")
		
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start2);\n")
		f.write("	for (i=0; i<nbiter; i++) {\n")
		f.write("		mult_mod_poly(pa, pa, pb);\n")
		f.write("	}\n")
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end2);\n")
		f.write("	diff2 = BILLION * (end2.tv_sec - start2.tv_sec) + (end2.tv_nsec - start2.tv_nsec);\n\n\n")
		
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start4);\n")
		f.write("	for (i=0; i<nbiter; i++) {\n")
		f.write("		BN_mod_mul_montgomery(opA, opA, opB, mont_ctx, ctx);\n")
		f.write("	}\n")
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end4);\n")
		f.write("	diff4 = BILLION * (end4.tv_sec - start4.tv_sec) + (end4.tv_nsec - start4.tv_nsec);\n\n\n")
		
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start6);\n")
		f.write("	for (i=0; i<nbiter; i++) {\n")
		f.write("		BN_mod_mul(opAA, opAA, opBB, opP, ctx);\n")
		f.write("	}\n")
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end6);\n")
		f.write("	diff6 = BILLION * (end6.tv_sec - start6.tv_sec) + (end6.tv_nsec - start6.tv_nsec);\n\n\n")
		
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start3);\n")
		f.write("	from_amns_to_int(C, pa);\n")
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end3);\n")
		f.write("	diff3 += BILLION * (end3.tv_sec - start3.tv_sec) + (end3.tv_nsec - start3.tv_nsec);\n\n")
		
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &start5);\n")
		f.write("	BN_from_montgomery(opA, opA, mont_ctx, ctx);\n")
		f.write("	clock_gettime(CLOCK_PROCESS_CPUTIME_ID , &end5);\n")
		f.write("	diff5 += BILLION * (end5.tv_sec - start5.tv_sec) + (end5.tv_nsec - start5.tv_nsec);\n\n\n")
		
		
		f.write("	printf(\"%lu %d %lu %lu %lu %lu\\n\\n\", mpz_sizeinbase (modul_p, 2), NB_COEFF, (diff2+diff3), (diff4+diff5), diff6, diff1);\n\n")
		
		f.write("	printf(\"nbiter = %d\\n\\n\", nbiter);\n")
		f.write("	gmp_printf(\"A       : %Zd\\n\", E);\n")
		f.write("	gmp_printf(\"B       : %Zd\\n\\n\", B);\n")
		f.write("	gmp_printf(\"r_gmp   : %Zd\\n\", A);\n")
		f.write("	gmp_printf(\"r_amns  : %Zd\\n\", C);\n")
		f.write("	gmp_printf(\"r_ssld  : %s\\n\", BN_bn2dec(opAA));\n")
		f.write("	gmp_printf(\"r_sslm  : %s\\n\", BN_bn2dec(opA));\n")
		f.write("	printf(\"\\ntime using gmp			= %lu nanoseconds\\n\", diff1);\n")
		f.write("	printf(\"\\ntime using openssl		= %lu nanoseconds\\n\", diff6);\n")
		f.write("	printf(\"\\ntime using amns prod		= %lu nanoseconds\\n\", diff2);\n")
		f.write("	printf(\"total time using amns prod	= %lu nanoseconds\\n\", (diff2+diff3));\n")
		f.write("	printf(\"\\ntime using openssl mont		= %lu nanoseconds\\n\", diff4);\n")
		f.write("	printf(\"total time using openssl mont	= %lu nanoseconds\\n\", (diff4+diff5));\n\n\n")
		
		f.write("	mpz_clears (A, B, C, E, NULL);\n")
		f.write("	gmp_randclear(r);\n\n")
		f.write("	BN_MONT_CTX_free(mont_ctx);\n")
		f.write("	BN_CTX_end(ctx);\n")
		f.write("	BN_CTX_free(ctx);\n\n")
		f.write("	free_data();\n")
		f.write("	return 0;\n")
		f.write("}\n\n")













