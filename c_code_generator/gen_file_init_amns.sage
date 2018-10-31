
def build_amns_init_h_file(dir_path):
	with open(dir_path+"/amns_init.h", "w") as f:
		f.write("#ifndef AMNS_INIT\n")
		f.write("#define AMNS_INIT\n\n\n")
		f.write("void init_data();\n\n")
		f.write("void free_data();\n\n")
		f.write("void compute_rho_pows();\n\n")
		f.write("#endif\n\n")



def build_amns_init_c_file(dir_path, p, gmm, lambd, lambd_mod_phi, small_int):
	with open(dir_path+"/amns_init.c", "w") as f:
		
		f.write("#include \"amns_init.h\"\n\n\n")
		
		f.write("void init_data(){\n")
		f.write("	int i;\n")
		f.write("	for(i=0; i<POLY_DEG; i++)\n")
		f.write("		mpz_init (gama_pow[i]);\n\n")
		f.write("	mpz_init (modul_p);\n\n\n")
		f.write("	mpz_set_str (modul_p, \"" + str(p) + "\", 10);\n\n")
		f.write("	mpz_set_str (gama_pow[0], \"" + str(gmm) + "\", 10);\n")
		f.write("	for(i=1; i<POLY_DEG; i++){\n")
		f.write("		mpz_mul (gama_pow[i], gama_pow[i-1], gama_pow[0]);\n")
		f.write("		mpz_mod (gama_pow[i], gama_pow[i], modul_p);\n")
		f.write("	}\n\n")
		f.write("	//~ Note : lambda = " + str(lambd) + " (see : modular multiplication in 'add_mult_poly.c').\n\n")
		f.write("	amns_rho = 1L << RHO_LOG2;\n\n")
		f.write("	//~ IMPORTANT : initialisations above must be done before those below.\n")
		f.write("	compute_rho_pows();\n")
		f.write("}\n\n\n")
		
		f.write("//~ computes representatives of powers of 'phi' in the AMNS.\n")
		f.write("void compute_rho_pows(){\n")
		f.write("	int i, l;\n")
		f.write("	" + small_int + " tmp_rho[NB_COEFF];\n\n")
		f.write("	for(i=0; i<NB_COEFF; i++)\n")
		f.write("		tmp_rho[i] = rho_rep[i];\n\n")
		f.write("	//~ computation of a representative of 'rho'\n")
		f.write("	from_mont_domain(rho_rep, rho_rep);\n\n")
		f.write("	//~ computation of representatives of (rho)^i (for i=2,3,...)\n")
		f.write("	l = NB_COEFF - 2;\n")
		f.write("	if (l > 0){\n")
		f.write("		mult_mod_poly(RHO_POWS[0], rho_rep, tmp_rho);\n")
		f.write("		for(i=1; i<l; i++)\n")
		f.write("			mult_mod_poly(RHO_POWS[i], RHO_POWS[i-1], tmp_rho);\n")
		f.write("	}\n")
		f.write("}\n\n\n")
		
		f.write("void free_data(){\n")
		f.write("	int i;\n")
		f.write("	for(i=0; i<POLY_DEG; i++)\n")
		f.write("		mpz_clear (gama_pow[i]);\n")
		f.write("\n")
		f.write("	mpz_clear (modul_p);\n")
		f.write("}\n\n")









