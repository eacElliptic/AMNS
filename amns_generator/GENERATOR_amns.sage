from multiprocessing import Process, Queue
load("amns_gen_utils.sage")

#~ ----------------------------------------------------------------------------------------------

#~ representation is of the form : a0 + a1.X + ... + a_{n-1}.X^{n-1}
#~ important : 'base' is supposed to be a power of two 
def compute_rep_in_base(val, base):
	val = Integer(val)
	if val==0 :
		return [0]
	rep = []
	mask = base - 1
	block_size = mask.nbits()
	while val != 0 :
		rep.append(val & mask)
		val >>= block_size
	return rep

#~ assumes that little-endian representation is used 
def compute_val_from_rep(rep, base):
	l = len(rep)
	if l==0 :
		return 0
	val = rep[0]
	for i in range(1, l):
		val += rep[i]*(base**i)
	return val


#~ return 'True' if n = +/- 2^t or n = +/-2^a +/- 2^b , otherwise return 'False'
def has_good_shape(n):
	if n == 0 :
		return False
	abv = Integer(abs(n))
	nb_setBit = abv.popcount()
	if nb_setBit <= 2 :
		#~ i.e : abs(n) = 2^a or 2^a + 2^b
		return True
	tmp = abv >> abv.trailing_zero_bits()
	almost2pow = ((tmp & (tmp + 1)) == 0)
	if almost2pow :
		#~ i.e : abs(n) = 2^a - 2^b, for a > b
		#~ note : '2^a - 2^b' is always in the form '11..1100...000'
		return True
	#~ at this level, 'n' has no special shape for us
	return False

#~ ----------------------------------------------------------------------------------------------

def compute_neg_inv_ri_poly(n, phi, ri_poly, ext_poly):
	R.<x> = QQ[]
	P = ZZ.quo(phi)
	PP.<y> = P[]

	e = R(ext_poly) 
	m = R(ri_poly)
	imy = PP(m.inverse_mod(e))

	return (-imy)

#~ returns a representative of 'op/phi'
def amns_red_int(op, ext_pol, ri_poly, neg_inv_ri, R, PP, phi):
	q = (PP(op)*neg_inv_ri).mod(PP(ext_pol))
	r0 = op + (R(q)*ri_poly).mod(ext_pol)
	return (r0/phi)

#~ returns a representative of '(op1*op2)/phi'
def amns_mont_mult(op1, op2, ext_pol, ri_poly, neg_inv_ri, R, PP, phi):
	#~ op1 = R(op1)
	#~ op2 = R(op2)
	#~ ext_pol = R(ext_pol)
	#~ ri_poly = R(ri_poly)
	#~ neg_inv_ri = PP(neg_inv_ri)
	
	c = (op1*op2).mod(ext_pol)
	q = (PP(c)*neg_inv_ri).mod(PP(ext_pol))
	r0 = c + (R(q)*ri_poly).mod(ext_pol)
	
	return (r0/phi)

#~ returns a representative of 'phi*val'
#~ we use here method 2 of conversion to AMNS
def compute_rep_in_amns(val, n, p, gmm, phi, ext_pol, ri_poly, neg_inv_ri, rho):
	F = GF(p)
	R.<x> = QQ[]
	P = ZZ.quo(phi); PP.<y> = P[]
	
	tmp_rep = Integer(F(val * phi.powermod(n, p)))
	rep = amns_red_int(R(tmp_rep), ext_pol, ri_poly, neg_inv_ri, R, PP, phi)
	for i in xrange(n-2):
		rep = amns_red_int(rep, ext_pol, ri_poly, neg_inv_ri, R, PP, phi)
	
	if F(rep(gmm)) != F(val*phi) :
		print("ERROR : Bad conversion !!!")
	
	if rep.norm(infinity) >= rho :
		print("ERROR : Element outside the AMNS after conversion !!!")
	
	return rep


#~ ----------------------------------------------------------------------------------------------

def nth_root_computation(F, n, lambd, find_all_nthroot, queue):
	res = []
	try:
		res = F(lambd).nth_root(n, all=find_all_nthroot)
		if not find_all_nthroot :
			res = [res]
	except ValueError, msg:
		#~ print msg
		res = []
	finally:
		queue.put(res)

#~ tries to find a nth_root of 'F(lambd)' within 'max_duration' seconds
def find_nth_root_with_timeout(F, n, lambd, nth_root_max_duration_checks, find_all_nthroot):
	res = []
	queue = Queue() # used to get the result
	proc = Process(target=nth_root_computation, args=(F, n, lambd, find_all_nthroot, queue)) # creation of a process calling longfunction with the specified arguments
	proc.start() # lauching the processus on another thread
	try:
		res = queue.get(timeout=nth_root_max_duration_checks) # getting the resultat under 'max_duration' seconds or stop
		proc.join() # proper delete if the computation has take less than timeout seconds
	except Exception, msg:
		proc.terminate() # kill the process
		#~ print ("Timed out!")
	return res

#~ checks if 'lambd' produces AMNS and returns it, if so. 
def check_lambda(word_size, F, p, n, phi, lambd, nth_root_max_duration_checks, find_all_nthroot) :
	
	gmms = find_nth_root_with_timeout(F, n, lambd, nth_root_max_duration_checks, find_all_nthroot)
	
	amns_list = []
	for gmm in gmms :
		gmm = Integer(gmm)
		if (gmm == 1):
			continue; # not a useful nth-root
		
		if (lambd%2) == 0 :
			ri = find_valid_vect_even_lambda(p, n, gmm)
		else :
			ri = find_valid_vect_odd_lambda(p, n, gmm)
		
		if ri != [] : # it should always be the case
			(rho_log2, nb_max_add) = compute_rho_min_log2_with_nb_add_max(word_size, phi, p, n, lambd, ri)
			if nb_max_add != -1 :
				rho = 1 << rho_log2
				R.<x> = QQ[]
				ri = list(ri)
				ri_poly = R(ri)
				ext_pol = R(x**n - lambd)
				neg_inv_ri = compute_neg_inv_ri_poly(n, phi, ri_poly, ext_pol)
				phi__rho_rep = compute_rep_in_amns(rho, n, p, gmm, phi, ext_pol, ri_poly, neg_inv_ri, rho)
				
				amns_list.append([nb_max_add, n, lambd, rho_log2, gmm, ri, list(phi__rho_rep), list(neg_inv_ri)])
	return amns_list

# checks all integers (with good shape) in {-lamb_max, ..., lamb_max}\{0}
def build_amns_candidates_for_n(word_size, F, p, n, phi, lamb_max, nth_root_max_duration_checks, find_all_nthroot):
	amns_cands = []
	
	for c in range(1, lamb_max+1):
		if not has_good_shape(c):
			continue
		tmp1 = check_lambda(word_size, F, p, n, phi,  c, nth_root_max_duration_checks, find_all_nthroot)
		tmp2 = check_lambda(word_size, F, p, n, phi, -c, nth_root_max_duration_checks, find_all_nthroot)
		if tmp1 != [] :
			amns_cands.append(tmp1)
		if tmp2 != [] :
			amns_cands.append(tmp2)
	
	print ("n = " + str(n) + " done.")
	
	return amns_cands

def generate_amns_candidates_with_n_min_max(word_size, prime_module, n_min, n_max, lamb_max, nth_root_max_duration_checks, find_all_nthroot):
	mont_phi = 1 << word_size  # will be used for multiplications in AMNS
	F = GF(prime_module)
	amns_list = []
	for n in range(n_min, n_max + 1):
		cand = build_amns_candidates_for_n(word_size, F, prime_module, n, mont_phi, lamb_max, nth_root_max_duration_checks, find_all_nthroot)
		if cand != [] :
			amns_list.append(cand)
	return flatten(amns_list, max_level=2)

