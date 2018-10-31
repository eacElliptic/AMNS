#~ return the infinite norm of 'vect'
def get_inf_norm(vect):
	tmp = [abs(a) for a in vect]
	return max(tmp)

#~ return an element of 'm' with the smallest infinite norm
def get_smallest_vect(m): 
	sml_vect = m[0]
	sv_norm = get_inf_norm(m[0])
	for v in m[1:] :
		tmp = get_inf_norm(v)
		if sv_norm > tmp :
			sv_norm = tmp
			sml_vect = v
	return list(sml_vect)

#~ -------------- Case where : E(X) = X^n - lambda, with abs(lambda) even ------------------------------
#~ note : each line corresponds to a poly of the form : a0 + a1.X + ... + a_{n-1}.X^{n-1}
def build_lattice_base_even_lambda(p, n, gmm):
	b = []
	l = [p] + [0]*(n-1)
	b.append(l)
	m = identity_matrix(n)
	for i in range(1,n):
		l = [(-gmm.powermod(i, p))%p] + list(m[i][1:])
		b.append(l)
	bb = matrix(b)
	return bb.LLL(delta=1, algorithm='NTL:LLL')

#~ an element 'pol' is valid if 'pol[0] % 2 == 1' 
def get_valid_elements_list_even_lambda(m):
	res = []
	for l in m:
		if l[0]%2 == 1 : # it must be odd for an odd resultant
			res.append(list(l))
	return res

#~ return a small valid vect
def find_valid_vect_even_lambda(p, n, gmm):
	m = build_lattice_base_even_lambda(p, n, gmm)
	vld_elmts = get_valid_elements_list_even_lambda(m)
	if vld_elmts == [] :
		print("WARNING : even lambda, no valid poly found; this shouldn't happen !!!")
	return get_smallest_vect(vld_elmts)

#~ -------------- Case where : E(X) = X^n - lambda, with abs(lambda) odd ------------------------------
#~ note : each line corresponds to a poly of the form : a0 + a1.X + ... + a_{n-1}.X^{n-1}
#~ note : 'p' is supposed odd
def build_lattice_base_odd_lambda(p, n, gmm):
	b = []
	l = [p] + [0]*(n-1)
	b.append(l)
	m = identity_matrix(n)
	for i in range(1,n):
		t = (-gmm.powermod(i, p))%p
		if t%2 == 1 :
			t += p
		l = [t] + list(m[i][1:])
		b.append(l)
	bb = matrix(b)
	return bb.LLL(delta=1, algorithm='NTL:LLL')

# assumes : k lower than 1<<n
def build_weight_vect(k, n):
	vt = list(bin(k)[2:])
	vt = [int(x) for x in vt]
	vt = [0]*(n-len(vt)) + vt
	return vt

#~ here, 'vect' is good if ... (see the article) 
def check_vect_odd_lambda(vect, e2, R):
	nb_odds = 0
	vect2 = []
	for x in vect :
		x2 = x % 2
		vect2.append(x2)
		if x2 == 1:
			nb_odds += 1
	if (nb_odds % 2) == 0 :
		#~ no need to continue if there is an even number of odd elements
		return False
	m = R(vect2)
	return (e2.gcd(m) == 1)

#~ return a small valid vect
def find_valid_vect_odd_lambda(p, n, gmm):
	R.<x> = GF(2)[]
	e2 = R(x**n - 1)
	vld_vects = []
	mm = build_lattice_base_odd_lambda(p, n, gmm)
	vmax = 1 << n
	for k in range(1, vmax): 
		k_vt = build_weight_vect(k, n)
		vsum = 0*mm[0] # to obtain the good type of vector
		for i in range(n):
			vsum += k_vt[i]*mm[i]
		if check_vect_odd_lambda(vsum, e2, R) :
			vld_vects.append(list(vsum))
	if vld_vects != [] :
		return get_smallest_vect(vld_vects)
	print("WARNING : odd lambda, no valid poly found; this shouldn't happen !!!")
	return []


#~ ---------------------------------------------------------------------------------

#~ assumes : E(X) = X^n - lambd
#~ note : we assume that elements coeffs can be negatifs in the AMNS, so one bit will be allocated for that.
def compute_rho_min_log2_with_nb_add_max(word_size, phi, p, n, lambd, ri_poly):
	
	nm_inf = get_inf_norm(ri_poly)
	rho = 2 * abs(lambd) * n * nm_inf  # codition for the montgomery-like reduction to be correct
	
	pw1 = ceil(log(rho, 2))
	pw2 = ceil(ceil(log(p, 2))/n) - 1 # we must be sure that the AMNS can represent all elements modulo p, -1 because elements signed.
	
	
	pw  = max(pw1, pw2)
	rho = 1 << pw                  # we want 'rho' to be a power of two
	
	nb_max_add = -1
	ok = True
	while ok :
		nb_max_add += 1
		if (1 << (word_size - 1)) <= ((nb_max_add+1)*(rho-1)) :
			ok = False # here, 'rho' (and/or 'nb_max_add') is too big for 'word_size', (one bit is used for sign)
		
		if phi < (2 * abs(lambd) * n * rho * (nb_max_add+1)**2) :
			ok = False # here, 'rho' (and/or 'nb_max_add') is too big for the montgomery-like reduction to be correct
		
		if (1 << (2*word_size - 1)) <= (n * abs(lambd) * ((rho-1)**2) * (nb_max_add+1)**2) :
			ok = False # here, 'rho' (and/or 'nb_max_add') is also too big for 'word_size' (because multiplication result coeffs will be too big)
	
	return (pw, (nb_max_add-1))




