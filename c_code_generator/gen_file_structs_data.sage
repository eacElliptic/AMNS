
def build_structs_data_file(dir_path, word_size, poly_deg, phi__rho_rep, rho_log2, nb_max_add, big_int_name, big_int, unsigned_big_int, small_int, unsigned_small_int):
	with open(dir_path+"/structs_data.h", "w") as f:
		
		f.write("#ifndef STRUCTS_DATA\n")
		f.write("#define STRUCTS_DATA\n\n\n")
		
		f.write("//~ IMPORTANT : We take 'phi = 1 << WORD_SIZE'\n")
		f.write("#define WORD_SIZE " + str(word_size) + "\n")
		f.write("#define POLY_DEG " + str(poly_deg) + "\n")
		f.write("#define NB_COEFF " + str(poly_deg+1) + "\n")
		f.write("#define NB_ADD_MAX " + str(nb_max_add) + "\n\n")
		
		f.write("#define RHO_LOG2 " + str(rho_log2) + "\n")
		f.write("//~ We will take : rho = 1 << RHO_LOG2.\n\n\n")
		
		f.write("typedef " + big_int_name + " " + big_int + ";\n")
		f.write("typedef unsigned " + big_int_name + " " + unsigned_big_int + ";\n\n\n")
		
		f.write("static " + unsigned_small_int + " amns_rho;\n\n")
		
		f.write("//~ will contain a representative of 'rho' in the amns\n")
		f.write("//~ important : this initial value is that of 'rho*phi'; correct value will be put during initialisation step\n")
		f.write("static " + small_int + " rho_rep[NB_COEFF] = {" + str(phi__rho_rep)[1:-1] + "};\n\n")
		
		f.write("//~ representatives of (RHO)^i (for i=2,4,...) in the amns (for convertions : int to amns)\n")
		f.write("static " + small_int + " RHO_POWS[(NB_COEFF - 2)][NB_COEFF];\n\n")
		
		f.write("static mpz_t modul_p;\n")
		f.write("static mpz_t gama_pow[POLY_DEG];\n\n")
		
		f.write("#endif\n\n")
	





