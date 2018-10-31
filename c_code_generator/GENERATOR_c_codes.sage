from os import makedirs as create_dir
from shutil import rmtree as delete_code
load("c_codes_gen_utils.sage")
load("gen_file_main.sage")
load("gen_file_init_amns.sage")
load("gen_file_structs_data.sage")
load("gen_file_add_mult_poly.sage")
load("gen_file_useful_functs.sage")



#~ assumes that polynomial representation form is P(X) = a_0 + ... + a_n.X^n = [a_0, ..., a_n].
#~ assumes 'amns_data' is a list containning (in this order) : [nb_max_add, n, lambd, rho_log2, gmm, ri, list(phi__rho_rep), list(neg_inv_ri)]
#~ assumes 'target_archi_info' is a list containning (in this order) : [word_size, small_int_name, unsigned_small_int_name, big_int_name, big_int_new_name, unsigned_big_int_name]
#~ WARNING : 'amns_data' is supposed generated for 'target_archi_info'
def build_amns_c_codes(p, amns_data, target_archi_info, p_num=0, amns_num=0):
	
	word_size = target_archi_info[0]
	small_int = target_archi_info[1]
	unsigned_small_int = target_archi_info[2]
	big_int_name = target_archi_info[3]
	big_int = target_archi_info[4]
	unsigned_big_int = target_archi_info[5]
	
	mont_phi = 1 << word_size
	lambd = amns_data[2]
	lambd_mod_phi = lambd % mont_phi
	nb_max_add = amns_data[0]
	n = amns_data[1]
	rho_log2 = amns_data[3]
	gmm = amns_data[4]
	red_int_pol = amns_data[5]
	phi__rho_rep = amns_data[6]
	neg_inv_ri = amns_data[7]
	
	
	dir_name = 'p' + str(p.nbits()) + '_'  + str(p_num) + '__' + str(n) + '_' + str(lambd) + '__' + str(amns_num)
	dir_path = "c_codes/" + dir_name
	try:
		create_dir(dir_path)
	except OSError:  # if this directory already exist
		delete_code(dir_path)
		create_dir(dir_path)
	
	
	build_main_file(dir_path, small_int)
	
	build_structs_data_file(dir_path, word_size, (n-1), phi__rho_rep, rho_log2, nb_max_add, big_int_name, big_int, unsigned_big_int, small_int, unsigned_small_int)

	build_amns_init_h_file(dir_path)
	build_amns_init_c_file(dir_path, p, gmm, lambd, lambd_mod_phi, small_int)

	build_add_mult_poly_h_file(dir_path, small_int, big_int)
	build_add_mult_poly_c_file(dir_path, n, mont_phi, lambd, small_int, unsigned_small_int, big_int, lambd_mod_phi, red_int_pol, neg_inv_ri)

	build_useful_functs_h_file(dir_path, small_int, unsigned_small_int, big_int)
	build_useful_functs_c_file(dir_path, small_int, big_int, unsigned_small_int)
	
	return;






