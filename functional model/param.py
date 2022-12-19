
import numpy as np
import numpy.random as rnd

posOLD=np.load("data/pos1.npy").tolist()

# move dlcs more rostrally to enable contact from tSts
for i in range(214,370):
	posOLD[i]=posOLD[i]-900


# NU = not useful to obtain sync to swim transitions

class SimulationParams(object):
	def __init__(self):
		# ===== NUMBER OF CELLS AND CELL TYPES SPECIFICATION (UNIVERSAL ORDERING) - FOR ANDREA =====
		self.cell_types = ["rb","dlc","ain","cin","din","mn","dla","xin","tst","tin","tsp","mhr","dinr","ecin"] # cell types
		self.tot_num_x_type = [0, 214, 156, 194, 542, 326, 496, 58, 60, 130, 40, 80, 40, 0, 0] # total number of cells in the whole connectome 
		self.num_types = [self.tot_num_x_type[i]+sum(self.tot_num_x_type[:i]) for i in xrange(len(self.tot_num_x_type))] # indexes of the cells per type
		self.halves = [w+(u-w)/2 for (u,w) in zip(self.num_types[1:],self.num_types[0:-1])] # half number of cells
		self.vect_index = [range(i,j) for (i,j) in zip(self.num_types[0:-1],self.num_types[1:])] # vector of indexes 
		self.left_index = [range(i,j) for (i,j) in zip(self.num_types[0:-1],self.halves)] # left indexes
		self.right_index = [range(i,j) for (i,j) in zip(self.halves,self.num_types[1:])] # right indexes
		self.colors = ["#ffd232","#ff0000","#4646b4","#00aadc","#96501e","#00963c","#ffaa8c","#696969","#ffa500","#FF00FF","#ffa500","#0000ff","k","#cc99ff"]

		self.tot_num_x_type_old = [0, 126, 104, 136, 384, 236, 338, 58, 60, 130, 40, 80, 40, 0, 0] # total number of cells in the whole connectome 
		self.num_types_old = [self.tot_num_x_type_old[i]+sum(self.tot_num_x_type_old[:i]) for i in xrange(len(self.tot_num_x_type_old))] # indexes of the cells per type
		self.halves_old = [w+(u-w)/2 for (u,w) in zip(self.num_types_old[1:],self.num_types_old[0:-1])] # half number of cells
		self.vect_index_old = [range(i,j) for (i,j) in zip(self.num_types_old[0:-1],self.num_types_old[1:])] # vector of indexes 
		self.left_index_old = [range(i,j) for (i,j) in zip(self.num_types_old[0:-1],self.halves_old)] # left indexes
		self.right_index_old = [range(i,j) for (i,j) in zip(self.halves_old,self.num_types_old[1:])] # right indexes

		# ===== NUMERICAL INTEGRATION PARAMETERS - FOR ANDREA =====
		self.varDt = False # variable time step integration
		self.atol = 1e-5 # absolute tolerance
		self.rtol = 1e-5 # relative tolerance
		self.dt = 0.025 # time step integration when varDt=False
		self.RecAll = 0 # this value can be 0="no variable recorded" - 1="record only time-voltage" - 2:"record all synaptic currents"
		self.time_end = 7000 # total time of simulation

		# ===== PARAM FOR INJECTING CURRENT TO SELECTED NEURONS =====
		self.pos_active = 30 # position of the RBs to activate
		self.n_active = 2 # number of active RBs - i.e. number of RBs receiving a 150pA and 5ms long pulse of current
		self.delay = 50.0 # ms from the start of current injection
		self.duration = 5.0 # ms of inj duration
		self.amplitude_mean = 0.3 # in nA
		self.amplitude_std = 0.0 # in nA

		# ===== SYNAPTIC-ANATOMICAL PARAMETERS =====
		self.sa_prop = 0.02 # variability in the cell RC position 
		self.fixed_delay = 1.0 # fixed synaptic delay 
		self.var_delay = 0.0035 # variability in the synaptic delay to multiply by the RC pre-post synaptic cell RC distance
		self.std_on = 1
		self.pos = posOLD # list of RC positions

		tot_num_x_type = self.tot_num_x_type
		num_half = np.divide(tot_num_x_type[1:],2)
		# add xIN pos
		for i in xrange(num_half[7]):
			self.pos.append(-130.0+20*i)
		for i in xrange(num_half[7]):
			self.pos.append(-130.0+20*i)
		# tSt pos
		for i in xrange(num_half[8]):	
			self.pos.append(-50.0*rnd.rand())
		for i in xrange(num_half[8]):
			self.pos.append(-50.0*rnd.rand())
		# tIN pos
		[self.pos.append(i) for i in np.linspace(150,330,num_half[9])]
		[self.pos.append(i) for i in np.linspace(150,330,num_half[9])]
		# tSp pos
		for i in xrange(num_half[10]):
			self.pos.append(10.0*rnd.rand())
		for i in xrange(num_half[10]):
			self.pos.append(10.0*rnd.rand())
		# mhr pos
		[self.pos.append(i) for i in np.linspace(400,550,num_half[11])]
		[self.pos.append(i) for i in np.linspace(400,550,num_half[11])]
		# dINr pos
		[self.pos.append(i) for i in np.linspace(300,500,num_half[12])]
		[self.pos.append(i) for i in np.linspace(300,500,num_half[12])]
		# ecIN pos
		if num_half[13]-num_half[12]>0:
			[self.pos.append(pos[i]) for i in range(366,750)]

		# ===== SYNAPTIC PARAMETERS =====
		self.v=[0,0.2,0.4,0.6,0.8] # xIN recurrent connection synaptic variability multiplicative constants
		self.alpha = 0.993 # xIN->xIN depression rate
		self.beta1 = 0.995 # dIN->dIN depression rate
		self.beta2 = 0.98 # xIN->CPG depression rate

		# synaptic strengths in this format : "pre_post": [syn_type, average strength, std, delay]: delay is "distance" if the cells have a defined distance dependent delay or a number corresponding to the prescribed delay if cell positions are not known
		molt=1e-3
		dli_str=0.04*molt
		mhr_inh = 3.0*molt 
		ain_inh=0.0*molt
		xin_cpg_ampa = 2*0.7*molt 
		xin_cpg_nmda = 2*0.35*molt 
		dli_str2=5*molt
		dli_std2=2*molt
		cin_inh=0.435*molt  
		ampa_str=0.593*molt  
		self.w = { 
			# pre rb
			"rb -> dla": (["ampa", 4.0*molt, 4.0*molt, "distance"],None),  
			"rb -> dlc": (["ampa", 4.0*molt, 4.0*molt, "distance"],None), 

			# pre din
			"din -> din": (["ampa", ampa_str, 0.05*ampa_str, "distance"], ["nmda_sat", 0.2*molt, 0.05*0.15*molt, "distance"]),  
			"din -> cin": (["ampa", ampa_str, 0.05*ampa_str, "distance"],None),  
			"din -> dlc": (["ampa", ampa_str, 0.05*ampa_str, "distance"],None),
			"din -> mn": (["ampa", ampa_str, 0.05*ampa_str, "distance"],None),
			"din -> dla": (["ampa", ampa_str, 0.05*ampa_str, "distance"],None),
			"din -> ain": (["ampa", 0.1*molt, 0.005*molt, 3.0],None),

			# pre ain
			"ain -> dlc": (["inh_ain", ain_inh, 0.05*ain_inh, "distance"],None),
			"ain -> ain": (["inh_ain", ain_inh, 0.05*ain_inh, "distance"],None),
			"ain -> cin": (["inh_ain", ain_inh, 0.05*ain_inh, "distance"],None),
			"ain -> din": (["inh_ain", ain_inh, 0.05*ain_inh, "distance"],None),
			"ain -> mn": (["inh_ain", ain_inh, 0.05*ain_inh, "distance"],None),
			"ain -> dla": (["inh_ain", ain_inh, 0.05*ain_inh, "distance"],None),

			# pre cin
			"cin -> dlc": (["inh", cin_inh, cin_inh*0.05, "distance"],None),
			"cin -> ain": (["inh", cin_inh, cin_inh*0.05, "distance"],None),
			"cin -> cin": (["inh", cin_inh, cin_inh*0.05, "distance"],None),
			"cin -> din": (["inh", cin_inh, cin_inh*0.05, "distance"],None), 
			"cin -> mn": (["inh", cin_inh, cin_inh*0.05, "distance"],None),
			"cin -> dla": (["inh", cin_inh, cin_inh*0.05, "distance"],None),

			# pre mn
			"mn -> ain": (["ampa", ampa_str, ampa_str*0.05, "distance"],None),
			"mn -> cin": (["ampa", ampa_str, ampa_str*0.05, "distance"],None),
			"mn -> din": (["ampa", ampa_str, ampa_str*0.05, "distance"],None),
			"mn -> mn": (["ampa", ampa_str, ampa_str*0.05, "distance"],None),

			# pre tst 
			"tsp -> mhr": (["ampa", 5.0*molt, 0.05*5.0*molt,"distance"],["nmda", 1.0*molt, 0.05*1.0*molt,"distance"]),  

			# pre mhr
			"mhr -> din": (["inh_mhr", mhr_inh, mhr_inh*0.05, "distance"],None), 
			"mhr -> ain": (["inh_mhr", mhr_inh, mhr_inh*0.05, "distance"],None),
			"mhr -> cin": (["inh_mhr", mhr_inh, mhr_inh*0.05, "distance"],None),
			"mhr -> mn": (["inh_mhr", mhr_inh, mhr_inh*0.05, "distance"],None),

			# pre dlc
			"dlc -> xin": (["ampa", dli_str2+dli_str2*0.4, dli_std2, 1.0], ["nmda", dli_str2+dli_str2*0.4, dli_std2, 1.0]),
			"dla -> xin": (["ampa", dli_str2-dli_str2*0.2, dli_std2, 1.0], ["nmda", dli_str2-dli_str2*0.2, dli_std2, 1.0]),

			# pre xin
			"xin -> xin": (["ampa_std", 6.5*molt, 0.0, 1.0], ["nmda_std", 1.4*molt, 0.0, 1.0]), 
			"xin -> din": (["ampa_std", xin_cpg_ampa, 0.0, 1.0], ["nmda_std", xin_cpg_nmda, 0.0, 1.0]),
			"xin -> cin": (["ampa_std", xin_cpg_ampa, 0.0, 1.0], ["nmda_std", xin_cpg_nmda, 0.0, 1.0]),
			"xin -> ain": (["ampa_std", xin_cpg_ampa, 0.0, 1.0], ["nmda_std", xin_cpg_nmda, 0.0, 1.0]),
			"xin -> mn": (["ampa_std", xin_cpg_ampa, 0.0, 1.0], ["nmda_std", xin_cpg_nmda, 0.0, 1.0]),
			
			# pre tst 
			"tst -> dlc": (["ampa", 8.0*molt, 4.0*molt, "distance"],None), 
			"tst -> tin": (["ampa",1.0225*molt,0.56375*0.05*molt,"distance"],["nmda",0.6925*molt,0.54625*0.05*molt,"distance"]),  

			# pre tin
			"tin -> din": (["ampa",0.34790508*molt,0.05*molt,"distance"],["nmda",0.31814063*molt,0.05*0.23*molt,"distance"]),  
			"tin -> xin": (["ampa", dli_str2+dli_str2*0.4, dli_std2, 1.0], ["nmda", dli_str2+dli_str2*0.4, dli_std2, 1.0]),
			}  


def create_params():
	params = SimulationParams()
	return params






