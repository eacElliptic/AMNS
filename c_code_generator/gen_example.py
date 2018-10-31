load("GENERATOR_c_codes.sage")

#~ -------------------------------- for 32-bits computer -------------------------------------------------------------------

word_size = 32

p = 113556790634660157267263905667270152901197793648196242794250351611737421945657

amns_data = [0, 13, 2, 25, 27909538719875102511551125443225519740089046541151275361691586791170147696467,
				[-49379,336486,-484377,82209,-328187,100938,225053,456205,76663,-559880,-231035,-340281,235204],
				[-353569,-1551726,-75969,-442411,585906,-733358,-205736,-1223599,-861081,-806387,-312442,257680,117074],
				[1049726001,3481490052,789633811,2090964479,32802490,1825549316,3461746829,2777209316,130940849,1318620905,1138226165,3298942106,281380021]]


#~ -------------------------------- for 64-bits computer -----------------------------------------------------------------

word_size = 64

p = 113556790634660157267263905667270152901197793648196242794250351611737421945657

amns_data = [4, 5, 2, 55, 21138829655901097943222712296378255075934887833046131852578355841495498583960,
				[1626219128380213,-982952150664854,-1663301451987807,-1345411614681303,292601798855450],
				[-5889130220341370,-4048668363131559,-2345436401269775,-821621600159061,-1839930390478987],
				[14356718804419984039,1794657206921851374,6508683138508426415,5592215177109344109,10537774877680181295]]


#~ -------------------------------- C codes generation according to the value of word_size ----------------------------------

#~ WARNING : 'amns_data' is supposed generated with the value of 'word_size'; otherwise, computations will be incorrect. 

if word_size == 64 :
	
	target_archi_info = [64, "int64_t", "uint64_t", "__int128", "int128", "uint128"]
	
	build_amns_c_codes(p, amns_data, target_archi_info)
elif word_size == 32 :
	
	target_archi_info = [32, "int", "uint", "long long", "llong", "ullong"]
	
	build_amns_c_codes(p, amns_data, target_archi_info)
else :
	print "Unknown architecture; please specify 'target_archi_info'." 

