def get_sign(k):
	if k < 0 :
		return '-'
	else :
		return '+'


def is_power_of_two(k) :
	if k <= 0 :
		return False
	return ((k & (k - 1)) == 0)
#~ -----------------------------------------------------------------------------------------------

def build_square_code(n, c, big_int):
	c_sign = get_sign(c)
	abs_c = abs(c)
	abs_c_is_2pow = is_power_of_two(abs_c)
	tmp_part1 = ['']*n
	tmp_part3 = ['']*n
	tmp_part2 = ['']*(n-1)
	i=0
	while i < n :
		tmp_part1[i] = 'tmp_prod_result['+str(i)+'] = '
		i=i+1
	i=0
	while i < n :
		j=0
		while j < i :
			pos = i+j
			if pos < n :
				tmp_part3[pos] += ' (' + big_int + ')pa['+ str(i) + '] * pa['+ str(j) + '] +'
			else:
				tmp_part2[pos%n] += ' (' + big_int + ')pa['+ str(i) + '] * pa['+ str(j) + '] +'
			j=j+1
		i=i+1
	for i in range(n):
		if tmp_part3[i] != '' :
			tmp_part3[i] = '((' + tmp_part3[i][1:-2] + ') << 1)'
	for i in range(n, 2*n-1):
		if tmp_part2[i%n] != '' :
			if i%2 == 0 :
				tmp_part2[i%n] = '((' + tmp_part2[i%n][1:-2] + ') << 1)'
			else :
				tmp_part2[i%n] = tmp_part2[i%n][1:-2]
	for i in range(n):
		pos = i+i
		if pos < n :
			if tmp_part3[pos] != '' :
				tmp_part3[pos] += ' + (' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
			else :
				tmp_part3[pos] += '(' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
		else:
			if tmp_part2[pos%n] != '' :
				tmp_part2[pos%n] += ' + (' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
			else :
				tmp_part2[pos%n] += '(' + big_int + ')pa['+ str(i) + '] * pa['+ str(i) + ']'
	for i in range(n, 2*n-1):
		if i%2 == 0 :
			if abs_c == 1 :
				tmp_part2[i%n] = ' ' + c_sign + ' (' + tmp_part2[i%n] +')'
			elif abs_c_is_2pow :
				tmp_part2[i%n] = ' ' + c_sign + ' ((' + tmp_part2[i%n] +') << ' + str(int(log(abs_c, 2))) + ')'
			else :
				tmp_part2[i%n] = ' ' + c_sign + ' ((' + tmp_part2[i%n] +') * '  + str(abs_c) + ')'
		else :
			if abs_c_is_2pow :
				tmp_part2[i%n] = ' ' + c_sign + ' ((' + tmp_part2[i%n] +') << ' + str(int(log(abs_c, 2))+1) + ')'
			else :
				tmp_part2[i%n] = ' ' + c_sign + ' ((' + tmp_part2[i%n] +') * '  + str(abs_c*2) + ')'
	result = ''
	for i in range(n-1):
		result += "	" + tmp_part1[i]+ tmp_part3[i] + tmp_part2[i] + ';\n'
	result += "	" + tmp_part1[n-1]+ tmp_part3[n-1] + ';\n'
	return result


def build_prod_code(n, c, big_int):
	c_sign = get_sign(c)
	abs_c = abs(c)
	abs_c_is_2pow = is_power_of_two(abs_c)
	tmp_part1 = ['']*n
	tmp_part2 = ['']*(n-1)
	i=0
	while i < n :
		tmp_part1[i] = 'tmp_prod_result['+str(i)+'] ='
		i=i+1
	i=0
	while i < n :
		j=0
		while j < n :
			pos = i+j
			if pos < n :
				tmp_part1[pos]+= ' (' + big_int + ')pa['+ str(i) + '] * pb['+ str(j) + '] +'
			else:
				tmp_part2[pos%n]+= ' (' + big_int + ')pa['+ str(i) + '] * pb['+ str(j) + '] +'
			j=j+1
		i=i+1
	i=0
	while i < (n-1) :
		tmp_part1[i] = tmp_part1[i][:-2]
		if abs_c == 1 :
			tmp_part2[i] = ' ' + c_sign + ' (' + tmp_part2[i][1:-2] +')'
		elif abs_c_is_2pow :
			tmp_part2[i] = ' ' + c_sign + ' ((' + tmp_part2[i][1:-2] +') << ' + str(int(log(abs_c, 2))) + ')'
		else :
			tmp_part2[i] = ' ' + c_sign + ' ((' + tmp_part2[i][1:-2] +') * ' + str(abs_c) + ')'
		i += 1
	tmp_part1[n-1] = tmp_part1[n-1][:-2]
	i=0
	result = ''
	while i < (n-1) :
		result += "	" + tmp_part1[i] + tmp_part2[i] + ';\n'
		i += 1
	result +=  "	" + tmp_part1[n-1] + ';\n'
	return result

#~ -----------------------------------------------------------------------------------------------

def build_red_int_code(unsigned_small_int, big_int, n, mont_phi, c, lambda_mod_phi, red_int_coeff, neg_inv_ri_rep_coeff):
	
	result = "	//~ computation of : op*neg_inv_ri_rep_coeff mod((X^n - c), mont_phi)\n"
	result += build_redInt_tmpQ_prod_code(n, mont_phi, lambda_mod_phi, neg_inv_ri_rep_coeff, unsigned_small_int)
	
	result += "\n	//~ computation of : tmp_q*red_int_coeff mod(X^n - c)\n"
	result += build_redInt_tmpZero_prod_code(n, c, red_int_coeff, big_int)
	
	result += "\n	//~ computation of : (op + tmp_zero)/mont_phi\n"
	result += build_redInt_lastStep_prod_code(n)
	
	return result


#~ important : 'lambda_mod_phi' is supposed > 0 and elements of 'neg_inv_ri_rep_coeff' are >= 0
def build_redInt_tmpQ_prod_code(n, mont_phi, lambda_mod_phi, neg_inv_ri_rep_coeff, unsigned_small_int):
	tmp_part1 = ['']*n
	tmp_part2 = ['']*(n-1)
	tmp_partStart = ['']*n
	
	#~ some adjustments, if needed
	neg_inv_ri_rep_coeff += [0] * (n - len(neg_inv_ri_rep_coeff))
	
	for i in xrange(n):
		tmp_partStart[i] = 'tmp_q['+str(i)+'] = '
	
	#~ note : elements of 'neg_inv_ri_rep_coeff' are >= 0
	for i in xrange(n):
		for j in xrange(n):
			
			if neg_inv_ri_rep_coeff[j] == 0 :
				continue
			elif neg_inv_ri_rep_coeff[j] == 1 :
				tmp = ' + (' + unsigned_small_int + ')op['+ str(i) + ']'
			else :
				tmp = ' + ((' + unsigned_small_int + ')op['+ str(i) + '] * ' + str(neg_inv_ri_rep_coeff[j]) + 'UL)'
			
			pos = i+j
			if pos < n :
				tmp_part1[pos] += tmp
			else:
				tmp_part2[pos%n] += tmp
	
	#~ some adjustments
	for i in xrange(n):
		if tmp_part1[i] == '':
			continue
		else :
			tmp_part1[i] = tmp_part1[i][3:]
	
	# note : lambda_mod_phi > 0
	for i in xrange(n-2):
		if (tmp_part2[i] == '') or (lambda_mod_phi == 1) :
			continue
		else :
			tmp_part2[i] = tmp_part2[i][3:]
			tmp_part2[i] = ' + ((' + tmp_part2[i] +') * ' + str(lambda_mod_phi) + ')'
	if (tmp_part2[n-2] == '') or (lambda_mod_phi == 1) :
		pass
	else :
		tred = (neg_inv_ri_rep_coeff[n-1]*lambda_mod_phi) % mont_phi
		tmp_part2[n-2] = ' + ((' + unsigned_small_int + ')op['+ str(n-1) + '] * ' + str(tred) + 'UL)'
		
	result = ''
	for i in xrange(n-1):
		if tmp_part1[i] == '' and tmp_part2[i] == '':
			continue
		result += "	" + tmp_partStart[i] + tmp_part1[i] + tmp_part2[i] + ';\n'
	if tmp_part1[n-1] == '':
		pass
	else :
		result +=  "	" + tmp_partStart[n-1] + tmp_part1[n-1] + ';\n'
	
	return result

def build_redInt_tmpZero_prod_code(n, c, red_int_coeff, big_int):
	tmp_part1 = ['']*n
	tmp_part2 = ['']*(n-1)
	tmp_partStart = ['']*n
	
	#~ some adjustments, if needed
	red_int_coeff += [0] * (n - len(red_int_coeff))
	
	abs_c = abs(c)
	c_sign = get_sign(c)
	
	for i in xrange(n):
		tmp_partStart[i] = 'tmp_zero['+str(i)+'] = '
	
	for i in xrange(n):
		for j in xrange(n):
			
			if red_int_coeff[j] == 0 :
				continue
			if abs(red_int_coeff[j]) == 1 :
				tmp = ' ' + get_sign(red_int_coeff[j]) + ' (' + big_int + ')tmp_q['+ str(i) + ']'
			else : # note : phi=(1<<word_size) >= 2*n*abs_c*rho and rho > inf_norm(red_int_coeff), so overflow shouldn't happen on target architecture.
				tmp = ' ' + get_sign(red_int_coeff[j]) + ' ((' + big_int + ')tmp_q['+ str(i) + '] * ' + str(abs(red_int_coeff[j])) + 'L)'
			
			pos = i+j
			if pos < n :
				tmp_part1[pos] += tmp
			else:
				tmp_part2[pos%n] += tmp
	
	#~ some adjustments
	for i in xrange(n):
		if tmp_part1[i] == '':
			continue
		elif tmp_part1[i][1] == '+':
			tmp_part1[i] = tmp_part1[i][3:]
		else :
			tmp_part1[i] = '-' + tmp_part1[i][3:]
	
	for i in xrange(n-2):
		if tmp_part2[i] == '':
			continue
		else :
			if tmp_part2[i][1] == '+':
				tmp_part2[i] = tmp_part2[i][3:]
			else :
				tmp_part2[i] = '-' + tmp_part2[i][3:]
			
			if abs_c == 1 :  
				tmp_part2[i] = ' ' + c_sign + ' (' + tmp_part2[i] + ')'
			else :
				tmp_part2[i] = ' ' + c_sign + ' ((' + tmp_part2[i] +') * ' + str(abs_c) + ')'
	
	if tmp_part2[n-2] == '' :
		pass
	else : # note : phi=(1<<word_size) >= 2*n*abs_c*rho and rho > inf_norm(red_int_coeff), so overflow shouldn't happen on target architecture.
		if get_sign(red_int_coeff[n-1]) == c_sign :
			tmp_part2[n-2] = ' + ((' + big_int + ')tmp_q['+ str(n-1) + '] * ' + str(abs(red_int_coeff[n-1])*abs_c) + 'L)'
		else :
			tmp_part2[n-2] = ' - ((' + big_int + ')tmp_q['+ str(n-1) + '] * ' + str(abs(red_int_coeff[n-1])*abs_c) + 'L)' # idem
	
	result = ''
	for i in xrange(n-1):
		if tmp_part1[i] == '' and tmp_part2[i] == '':
			continue
		result += "	" + tmp_partStart[i] + tmp_part1[i] + tmp_part2[i] + ';\n'
	if tmp_part1[n-1] == '':
		pass
	else :
		result +=  "	" + tmp_partStart[n-1] + tmp_part1[n-1] + ';\n'
	
	return result

def build_redInt_lastStep_prod_code(n):
	
	result = ''
	for i in xrange(n):
		tmp = "	rop[" + str(i) + "] = (op[" + str(i) + "] + tmp_zero[" + str(i) + "]) >> WORD_SIZE;\n"
		result += tmp
	
	return result




