
from abc import ABCMeta
import numpy as np
from neuron import h

molt_gmax=1e-9/1e-5
molt_cap=1

def connection(self,dest=None,InputSynName=None,w=0.0,delay=0.0,gmax=0.0,seed=None):
	if InputSynName=='ampa':
		if hasattr(dest,"syn_ampa")==0:
			dest.syn_ampa = h.syn_ampa(dest.soma(0.5))
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_ampa,self.theta,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['i_ampa']=h.Vector() 
			dest.record['i_ampa'].record(dest.syn_ampa._ref_i,sec=self.soma)
	elif InputSynName=='ampa_dale':
		if hasattr(dest,"syn_ampa_dale")==0:
			dest.syn_ampa_dale = h.syn_ampa_dale(dest.soma(0.5))
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_ampa_dale,self.theta,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['i_ampa_dale']=h.Vector() 
			dest.record['i_ampa_dale'].record(dest.syn_ampa_dale._ref_i,sec=self.soma)
	elif InputSynName=='nmda':
		if hasattr(dest,"syn_nmda")==0:
			dest.syn_nmda = h.syn_nmda(dest.soma(0.5))
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_nmda,self.theta,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['i_nmda']=h.Vector() 
			dest.record['i_nmda'].record(dest.syn_nmda._ref_i,sec=self.soma)
	elif InputSynName=='inh':
		if hasattr(dest,"syn_inh")==0:
			dest.syn_inh = h.syn_inh(dest.soma(0.5))
			#dest.syn_inh.tau_c = np.abs(np.random.normal(4.0,1.0))
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_inh,self.theta,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['i_inh']=h.Vector()
			dest.record['i_inh'].record(dest.syn_inh._ref_i,sec=self.soma)
	elif InputSynName=='inh_mhr':
		if hasattr(dest,"syn_inh_mhr")==0:
			dest.syn_inh_mhr = h.syn_inh_mhr(dest.soma(0.5))
			#dest.syn_inh_mhr.tau_c = 15.0
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_inh_mhr,self.theta,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['i_inh']=h.Vector()
			dest.record['i_inh'].record(dest.syn_inh._ref_i,sec=self.soma)
	elif InputSynName=='inh_ain':
		if hasattr(dest,"syn_inh_mhr")==0:
			dest.syn_inh_mhr = h.syn_inh_mhr(dest.soma(0.5))
			dest.syn_inh_mhr.tau_c = 80.0
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_inh_mhr,self.theta,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['i_inh']=h.Vector()
			dest.record['i_inh'].record(dest.syn_inh._ref_i,sec=self.soma)
	elif InputSynName=='gap':
		dest.gap_list.append(h.Gap(dest.soma(0.5)))
		dest.gap_list[-1].gmax=gmax
		self.connlist.append(h.setpointer(self.soma(0.5)._ref_v,'vgap',dest.gap_list[-1]))
		if dest.RecAll==2:
			dest.record['i_gap']=h.Vector() 
			dest.record['i_gap'].record(dest.gap_list[-1]._ref_i,sec=self.soma)
	elif InputSynName=='ampa_var':
		if hasattr(dest,"syn_ampa_var")==0:
			dest.syn_ampa_var = h.syn_ampa_var(dest.soma(0.5))
			dest.rnd_ampa=h.Random(seed)
			dest.rnd_ampa.uniform(0,1)
			dest.syn_ampa_var.noiseFromRandom(dest.rnd_ampa)
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_ampa_var,0.0,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['IAMPA_VAR'].record(dest.syn_ampa_var._ref_i,sec=self.soma)
	elif InputSynName=='nmda_var':
		if hasattr(dest,"syn_nmda_var")==0:
			dest.syn_nmda_var = h.syn_nmda_var(dest.soma(0.5))
			dest.rnd_nmda=h.Random(seed)
			dest.rnd_nmda.uniform(0,1)
			dest.syn_nmda_var.noiseFromRandom(dest.rnd_nmda)
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_nmda_var,0.0,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['INMDA_VAR'].record(dest.syn_nmda_var._ref_i,sec=self.soma)
	elif InputSynName=='nmda_sat':
		if hasattr(dest,"syn_nmda_sat")==0:
			dest.syn_nmda_sat = h.syn_nmda_sat(dest.soma(0.5))
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_nmda_sat,self.theta,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['i_nmda_sat']=h.Vector() 
			dest.record['i_nmda_sat'].record(dest.syn_nmda_sat._ref_i,sec=self.soma)
	elif InputSynName=='ampa_sat':
		if hasattr(dest,"syn_ampa_sat")==0:
			dest.syn_ampa_sat = h.syn_ampa_sat(dest.soma(0.5))
			print dest.syn_ampa_sat.ampa_sat
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_ampa_sat,self.theta,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['i_ampa_sat']=h.Vector() 
			dest.record['i_ampa_sat'].record(dest.syn_ampa_sat._ref_i,sec=self.soma)
	elif InputSynName=='inh_sat':
		if hasattr(dest,"syn_inh_sat")==0:
			dest.syn_inh_sat = h.syn_inh_sat(dest.soma(0.5))
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_inh_sat,self.theta,delay,w,sec=self.soma))
		if dest.RecAll==2:
			dest.record['i_inh_sat']=h.Vector() 
			dest.record['i_inh_sat'].record(dest.syn_inh_sat._ref_i,sec=self.soma)
	elif InputSynName=='nmda_std':
		if hasattr(dest,"syn_nmda_std")==0:
			dest.syn_nmda_std = h.syn_nmda_std(dest.soma(0.5))
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_nmda_std,0.0,delay,w,sec=self.soma))
		if self.RecAll==2:
			dest.record['INMDA_STD'].record(dest.syn_nmda_std._ref_i,sec=self.soma)
	elif InputSynName=='ampa_std':
		if hasattr(dest,"syn_ampa_std")==0:
			dest.syn_ampa_std = h.syn_ampa_std(dest.soma(0.5))
		self.connlist.append(h.NetCon(self.soma(0.5)._ref_v,dest.syn_ampa_std,0.0,delay,w,sec=self.soma))
		if self.RecAll==2:
			dest.record['IAMPA_STD'].record(dest.syn_ampa_std._ref_i,sec=self.soma)
	else:
		print "connection " + InputSynName + " unknown!!\n"

class Cell(object):
	__metaclass__ = ABCMeta
	
	def __init__(self, cap=0.4*molt_cap, RecAll=0, varDt=False, atol=1e-5, rtol=1e-5, theta=0.0, channels={}):
		'''
			This class defines a single compartment cell standard class with the following parameters
			- cap: capacitance
			- RecAll: 0 save nothing, 1 saves time and voltage, 2 saves also synaptic currents
			- varDt: True set the use of local time step integrator, False uses predefined step dt
			- atol, rtol: in case varDt=1 sets the tolerances of the CVode numerical integration
			- theta: threshold for spike transmission
			- channels: active channels
		'''
	
		# position and side
		self.pos = 0.0
		self.side = 0 # 0=unassigned, 1=left, 2=right
		self.index = None

		# param
		self.RecAll = RecAll
		self.theta = theta
	
		# soma anatomy
		self.soma=h.Section(cell=self)
		self.soma.nseg=1
		self.soma.L=100.0 # um (?)
		self.surface_area=1e3 # cm^2
		self.soma.diam= self.surface_area/(np.pi*self.soma.L) #(self.surface_area*1e8)/(np.pi*self.soma.L)
		self.soma.cm=cap

		# soma ionic channels
		for ch in channels.viewkeys():
			self.soma.insert(ch)
		
		# set current clamp and voltage clamp 
		self.CC = h.IClamp(self.soma(0.5))
		self.CC2 = h.IClamp(self.soma(0.5))
		self.VC = h.VClamp(self.soma(0.5))

		# record spike times
		self.record = {}
		self.threshold = 0.0
		self.nc_spike = h.NetCon(self.soma(0.5)._ref_v,None,self.threshold,0.0,0.0,sec=self.soma)
		self.record['spk'] = h.Vector()
		self.nc_spike.record(self.record['spk'])

		# vector of connections
		self.connlist = []
		self.gap_list=[]

		if varDt:
			# variable time step integrator
			self.Hines = h.CVode()
			self.Hines.active(1)
			self.Hines.use_local_dt(1)
			self.Hines.atol(atol)
			self.Hines.rtol(rtol)

		if RecAll>0:
			self.record['t']=h.Vector()
			self.record['vm']=h.Vector()
			if varDt:
				self.Hines.record(self.soma(0.5)._ref_v,self.record['vm'],self.record['t'],sec=self.soma)
				if RecAll==2:
					for item in channels.viewitems():
						self.record["i_"+item[0]]=h.Vector()
						self.record["i_"+item[0]].record(eval("self.soma(0.5)._ref_"+item[1])) 
			else: 
				self.record['t'].record(h._ref_t)
				self.record['vm'].record(self.soma(0.5)._ref_v)
				if RecAll==2:
					for item in channels.viewitems():
						self.record["i_"+item[0]]=h.Vector()
						self.record["i_"+item[0]].record(eval("self.soma(0.5)._ref_"+item[1])) 

	def connect(self,dest,InputSynName=None,w=0.0,delay=0.0,gmax=0.0,seed=None):
		connection(self,dest=dest,InputSynName=InputSynName,w=w,delay=delay,gmax=gmax,seed=seed)

	# destroy
	def destroy(self):
		del self.record
		#del self.soma
		
class MHRcell(Cell):
	def __init__(self,RecAll=1,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.color = "#a9a9a9"
		self.whatami = "mhr" 
		self.impulses = []

		# channel properties and cell init
		self.channels = {"pas":"i_pas","MN_na":"ina_MN_na","MN_kFast":"ik_MN_kFast","MN_kSlow":"ik_MN_kSlow"}

		# init cell
		Cell.__init__(self, cap=.4*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		self.soma.g_pas = 3.816*molt_gmax
		self.soma.e_pas = -60.0 

		self.soma.gmax_MN_na = 420.0*molt_gmax 
		self.soma.gmax_MN_kFast = 70.0*molt_gmax
		self.soma.gmax_MN_kSlow = 10.0*molt_gmax

		d = 10 # to obtain ~-40/-35mV voltage treshold characteristic of mhr
		# Na constants
		self.soma.alpha_A_m_MN_na = 13.26 +3 # to obtain multiple firing at any level of current injected and avoid the fixed depolarized state
		self.soma.alpha_B_m_MN_na = 0.0
		self.soma.alpha_C_m_MN_na = 3.0
		self.soma.alpha_D_m_MN_na = -3.01-d
		self.soma.alpha_E_m_MN_na = -12.56
	    	
		self.soma.beta_A_m_MN_na = 5.73 
		self.soma.beta_B_m_MN_na = 0.0
		self.soma.beta_C_m_MN_na = 1.0
		self.soma.beta_D_m_MN_na = 6.01-d
		self.soma.beta_E_m_MN_na = 9.69

		self.soma.alpha_A_h_MN_na = 0.06 
		self.soma.alpha_B_h_MN_na = 0.0
		self.soma.alpha_C_h_MN_na = 0.0
		self.soma.alpha_D_h_MN_na = 19.88-d
		self.soma.alpha_E_h_MN_na = 26.0
	    	
		self.soma.beta_A_h_MN_na = 4.08 
		self.soma.beta_B_h_MN_na = 0.0
		self.soma.beta_C_h_MN_na = 0.001
		self.soma.beta_D_h_MN_na = -8.09-d
		self.soma.beta_E_h_MN_na = -10.21 

		# Kf constants
		self.soma.alpha_A_MN_kFast = 3.1  # lower the frequency
		self.soma.alpha_B_MN_kFast = 0.0
		self.soma.alpha_C_MN_kFast = 1.0
		self.soma.alpha_D_MN_kFast = -32.5-d
		self.soma.alpha_E_MN_kFast = -9.3
	    	
		self.soma.beta_A_MN_kFast = 1.1 -0.9 # to obtain multiple firing at any level of current injected and avoid the fixed depolarized state &  lower the frequency
		self.soma.beta_B_MN_kFast = 0.0
		self.soma.beta_C_MN_kFast = 2.0 
		self.soma.beta_D_MN_kFast = 3.98-d
		self.soma.beta_E_MN_kFast = 16.19

		# Ks constants
		self.soma.alpha_A_MN_kSlow = 4.0 
		self.soma.alpha_B_MN_kSlow = 0.0
		self.soma.alpha_C_MN_kSlow = 1.0 
		self.soma.alpha_D_MN_kSlow = -53.0-d
		self.soma.alpha_E_MN_kSlow = -7.74
	    	
		self.soma.beta_A_MN_kSlow = 0.01
		self.soma.beta_B_MN_kSlow = 0.0
		self.soma.beta_C_MN_kSlow = 1.0
		self.soma.beta_D_MN_kSlow = 47.0-d
		self.soma.beta_E_MN_kSlow = 6.1



class DLAcell(Cell):
	def __init__(self,RecAll=1,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="dla"
		self.color="#ffaa8c"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","MN_na":"ina_MN_na","MN_kFast":"ik_MN_kFast","MN_kSlow":"ik_MN_kSlow"}

		# init cell		
		Cell.__init__(self, cap=.4*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# Na constants
		self.soma.cm = .4*molt_cap
		self.soma.g_pas = 0.6964*molt_gmax
		self.soma.e_pas = -63.0
		
		self.soma.gmax_MN_na = 150.0*molt_gmax 
		self.soma.gmax_MN_kFast = 70.0*molt_gmax
		self.soma.gmax_MN_kSlow = 5.0*molt_gmax

		# Na constants
		self.soma.alpha_A_m_MN_na = 13.26
		self.soma.alpha_B_m_MN_na = 0.0
		self.soma.alpha_C_m_MN_na = 1.2
		self.soma.alpha_D_m_MN_na = -9.01
		self.soma.alpha_E_m_MN_na = -12.56
	    	
		self.soma.beta_A_m_MN_na = 5.73 
		self.soma.beta_B_m_MN_na = 0.0
		self.soma.beta_C_m_MN_na = 1.0
		self.soma.beta_D_m_MN_na = 1.01
		self.soma.beta_E_m_MN_na = 9.69

		self.soma.alpha_A_h_MN_na = 0.04
		self.soma.alpha_B_h_MN_na = 0.0
		self.soma.alpha_C_h_MN_na = 0.0
		self.soma.alpha_D_h_MN_na = 14.88
		self.soma.alpha_E_h_MN_na = 26.0
	    	
		self.soma.beta_A_h_MN_na = 2.04
		self.soma.beta_B_h_MN_na = 0.0
		self.soma.beta_C_h_MN_na = 0.001
		self.soma.beta_D_h_MN_na = -13.09
		self.soma.beta_E_h_MN_na = -10.21 

		# Kf constants
		self.soma.alpha_A_MN_kFast = 3.1
		self.soma.alpha_B_MN_kFast = 0.0
		self.soma.alpha_C_MN_kFast = 1.0
		self.soma.alpha_D_MN_kFast = -37.5
		self.soma.alpha_E_MN_kFast = -9.3
	    	
		self.soma.beta_A_MN_kFast = 1.1 
		self.soma.beta_B_MN_kFast = 0.0
		self.soma.beta_C_MN_kFast = 2.0
		self.soma.beta_D_MN_kFast = -1.02
		self.soma.beta_E_MN_kFast = 16.19

		# Ks constants
		self.soma.alpha_A_MN_kSlow = 4.0
		self.soma.alpha_B_MN_kSlow = 0.0
		self.soma.alpha_C_MN_kSlow = 1.0 
		self.soma.alpha_D_MN_kSlow = -58.0
		self.soma.alpha_E_MN_kSlow = -7.74
	    	
		self.soma.beta_A_MN_kSlow = 0.01
		self.soma.beta_B_MN_kSlow = 0.0
		self.soma.beta_C_MN_kSlow = 1.0
		self.soma.beta_D_MN_kSlow = 42.0
		self.soma.beta_E_MN_kSlow = 6.1


class DLCcell(Cell):
	def __init__(self,RecAll=1,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="dlc"
		self.color="#ff0000"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","MN_na":"ina_MN_na","MN_kFast":"ik_MN_kFast","MN_kSlow":"ik_MN_kSlow"}

		# init cell		
		Cell.__init__(self, cap=.4*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# Na constants
		self.soma.cm = .4*molt_cap
		self.soma.g_pas = 2.3364*molt_gmax
		self.soma.e_pas = -66.0 
		
		self.soma.gmax_MN_na = 420.0*molt_gmax 
		self.soma.gmax_MN_kFast = 70.0*molt_gmax
		self.soma.gmax_MN_kSlow = 10.0*molt_gmax

		# Na constants
		self.soma.alpha_A_m_MN_na = 13.26
		self.soma.alpha_B_m_MN_na = 0.0
		self.soma.alpha_C_m_MN_na = 3.0
		self.soma.alpha_D_m_MN_na = -3.01
		self.soma.alpha_E_m_MN_na = -12.56
	    	
		self.soma.beta_A_m_MN_na = 5.73 
		self.soma.beta_B_m_MN_na = 0.0
		self.soma.beta_C_m_MN_na = 1.0
		self.soma.beta_D_m_MN_na = 6.01
		self.soma.beta_E_m_MN_na = 9.69

		self.soma.alpha_A_h_MN_na = 0.06 
		self.soma.alpha_B_h_MN_na = 0.0
		self.soma.alpha_C_h_MN_na = 0.0
		self.soma.alpha_D_h_MN_na = 19.88
		self.soma.alpha_E_h_MN_na = 26.0
	    	
		self.soma.beta_A_h_MN_na = 4.08
		self.soma.beta_B_h_MN_na = 0.0
		self.soma.beta_C_h_MN_na = 0.001
		self.soma.beta_D_h_MN_na = -8.09
		self.soma.beta_E_h_MN_na = -10.21 

		# Kf constants
		self.soma.alpha_A_MN_kFast = 3.1
		self.soma.alpha_B_MN_kFast = 0.0
		self.soma.alpha_C_MN_kFast = 1.0
		self.soma.alpha_D_MN_kFast = -32.5
		self.soma.alpha_E_MN_kFast = -9.3
	    	
		self.soma.beta_A_MN_kFast = 1.1 
		self.soma.beta_B_MN_kFast = 0.0
		self.soma.beta_C_MN_kFast = 2.0
		self.soma.beta_D_MN_kFast = 3.98
		self.soma.beta_E_MN_kFast = 16.19

		# Ks constants
		self.soma.alpha_A_MN_kSlow = 4.0
		self.soma.alpha_B_MN_kSlow = 0.0
		self.soma.alpha_C_MN_kSlow = 1.0 
		self.soma.alpha_D_MN_kSlow = -53.0
		self.soma.alpha_E_MN_kSlow = -7.74
	    	
		self.soma.beta_A_MN_kSlow = 0.01
		self.soma.beta_B_MN_kSlow = 0.0
		self.soma.beta_C_MN_kSlow = 1.0
		self.soma.beta_D_MN_kSlow = 47.0
		self.soma.beta_E_MN_kSlow = 6.1


class aINcell(Cell):
	def __init__(self,RecAll=1,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="aIN"
		self.color="#4646b4"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","MN_na":"ina_MN_na","MN_kFast":"ik_MN_kFast","MN_kSlow":"ik_MN_kSlow"}

		# init cell
		Cell.__init__(self, cap=.4*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# Na constants
		self.soma.cm = .4*molt_cap
		self.soma.g_pas = 1.3514*molt_gmax
		self.soma.e_pas = -54.0 
		
		self.soma.gmax_MN_na = 150.0*molt_gmax 
		self.soma.gmax_MN_kFast = 15.0*molt_gmax
		self.soma.gmax_MN_kSlow = 2.5*molt_gmax

		# Na constants
		self.soma.alpha_A_m_MN_na = 8.67
		self.soma.alpha_B_m_MN_na = 0.0
		self.soma.alpha_C_m_MN_na = 0.5
		self.soma.alpha_D_m_MN_na = -13.01
		self.soma.alpha_E_m_MN_na = -18.56
	    	
		self.soma.beta_A_m_MN_na = 5.73 
		self.soma.beta_B_m_MN_na = 0.0
		self.soma.beta_C_m_MN_na = 1.0
		self.soma.beta_D_m_MN_na = -2.99
		self.soma.beta_E_m_MN_na = 9.69

		self.soma.alpha_A_h_MN_na = 0.04 
		self.soma.alpha_B_h_MN_na = 0.0
		self.soma.alpha_C_h_MN_na = 0.0
		self.soma.alpha_D_h_MN_na = 15.8
		self.soma.alpha_E_h_MN_na = 26.0
	    	
		self.soma.beta_A_h_MN_na = 4.08
		self.soma.beta_B_h_MN_na = 0.0
		self.soma.beta_C_h_MN_na = 0.001
		self.soma.beta_D_h_MN_na = -19.09
		self.soma.beta_E_h_MN_na = -10.21 

		# Kf constants
		self.soma.alpha_A_MN_kFast = 3.1
		self.soma.alpha_B_MN_kFast = 0.0
		self.soma.alpha_C_MN_kFast = 1.0
		self.soma.alpha_D_MN_kFast = -35.5
		self.soma.alpha_E_MN_kFast = -9.3
	    	
		self.soma.beta_A_MN_kFast = 1.1 
		self.soma.beta_B_MN_kFast = 0.0
		self.soma.beta_C_MN_kFast = 1.0
		self.soma.beta_D_MN_kFast = 0.98
		self.soma.beta_E_MN_kFast = 16.19

		# Ks constants
		self.soma.alpha_A_MN_kSlow = 0.2
		self.soma.alpha_B_MN_kSlow = 0.0
		self.soma.alpha_C_MN_kSlow = 1.0 
		self.soma.alpha_D_MN_kSlow = -10.96
		self.soma.alpha_E_MN_kSlow = -7.74
	    	
		self.soma.beta_A_MN_kSlow = 0.05
		self.soma.beta_B_MN_kSlow = 0.0
		self.soma.beta_C_MN_kSlow = 1.0
		self.soma.beta_D_MN_kSlow = -22.07
		self.soma.beta_E_MN_kSlow = 6.1


class DALECell(Cell):
	def __init__(self,RecAll=1,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="dale"
		self.color="#4646b4"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","MN_na":"ina_MN_na","MN_kFast":"ik_MN_kFast","MN_kSlow":"ik_MN_kSlow"}

		# init cell
		Cell.__init__(self, cap=.4*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# Na constants
		self.soma.cm = .4*molt_cap
		self.soma.g_pas = 1.3514*molt_gmax
		self.soma.e_pas = -54.0 
		
		self.soma.gmax_MN_na = 150.0*molt_gmax 
		self.soma.gmax_MN_kFast = 15.0*molt_gmax
		self.soma.gmax_MN_kSlow = 2.5*molt_gmax

		# Na constants
		self.soma.alpha_A_m_MN_na = 8.67
		self.soma.alpha_B_m_MN_na = 0.0
		self.soma.alpha_C_m_MN_na = 1.0
		self.soma.alpha_D_m_MN_na = -1.01
		self.soma.alpha_E_m_MN_na = -12.56
	    	
		self.soma.beta_A_m_MN_na = 3.82 
		self.soma.beta_B_m_MN_na = 0.0
		self.soma.beta_C_m_MN_na = 1.0
		self.soma.beta_D_m_MN_na = 9.01
		self.soma.beta_E_m_MN_na = 9.69

		self.soma.alpha_A_h_MN_na = 0.08
		self.soma.alpha_B_h_MN_na = 0.0
		self.soma.alpha_C_h_MN_na = 0.0
		self.soma.alpha_D_h_MN_na = 38.88
		self.soma.alpha_E_h_MN_na = 26.0
	    	
		self.soma.beta_A_h_MN_na = 4.08
		self.soma.beta_B_h_MN_na = 0.0
		self.soma.beta_C_h_MN_na = 1.0
		self.soma.beta_D_h_MN_na = -5.09
		self.soma.beta_E_h_MN_na = -10.21 

		# Kf constants
		self.soma.alpha_A_MN_kFast = 3.1
		self.soma.alpha_B_MN_kFast = 0.0
		self.soma.alpha_C_MN_kFast = 1.0
		self.soma.alpha_D_MN_kFast = -29.5
		self.soma.alpha_E_MN_kFast = -23.3
	    	
		self.soma.beta_A_MN_kFast = 0.44 
		self.soma.beta_B_MN_kFast = 0.0
		self.soma.beta_C_MN_kFast = 1.0
		self.soma.beta_D_MN_kFast = 6.98
		self.soma.beta_E_MN_kFast = 16.19

		# Ks constants
		self.soma.alpha_A_MN_kSlow = 0.16
		self.soma.alpha_B_MN_kSlow = 0.0
		self.soma.alpha_C_MN_kSlow = 1.0
		self.soma.alpha_D_MN_kSlow = -4.69
		self.soma.alpha_E_MN_kSlow = -7.74
	    	
		self.soma.beta_A_MN_kSlow = 0.04
		self.soma.beta_B_MN_kSlow = 0.0
		self.soma.beta_C_MN_kSlow = 1.0
		self.soma.beta_D_MN_kSlow = -16.07
		self.soma.beta_E_MN_kSlow = 6.1

class dINcell(Cell):
	def __init__(self,RecAll=2,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="dIN"
		self.color="#96501e"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","dIN_na":"ina_dIN_na","dIN_kFast":"ik_dIN_kFast","dIN_kSlow":"ik_dIN_kSlow","dIN_ca":"ica_dIN_ca"}

		# init cell		
		Cell.__init__(self, cap=1.0*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# parameters for dINs
		self.parameters = {"erev_lk": -52.0, 
		                   "gmax_lk": 1.405*molt_gmax,
		                   "erev_na": 50.0,
		                   "gmax_na": 240.5*molt_gmax, # Down to 210.5*molt_gmax reduces mid-cycle dINs
		                   "erev_k": -81.5,
		                   "gmax_kf": 12.0*molt_gmax,
		                   "gmax_ks": 9.6*molt_gmax,
		                   "T_ca": 300.0,
		                   "Sin_ca": 100e-9,
		                   "Sout_ca": 10e-6,
		                   "perm_ca": 1.425e-10/1e-5} # this gives the same results as Bob's dINs

		# ion channels
		self.soma.cm=1.0*molt_cap
		self.soma.g_pas=self.parameters["gmax_lk"]
		self.soma.e_pas=self.parameters["erev_lk"]
		self.soma.gmax_dIN_na=self.parameters["gmax_na"]
		self.soma.gmax_dIN_kFast=self.parameters["gmax_kf"]
		self.soma.gmax_dIN_kSlow=self.parameters["gmax_ks"]
		self.soma.T_dIN_ca=self.parameters["T_ca"]
		self.soma.perm_dIN_ca=self.parameters["perm_ca"]
		self.soma.Sin_dIN_ca=self.parameters["Sin_ca"]
		self.soma_Sout_dIN_ca=self.parameters["Sout_ca"]
		self.soma.ena=self.parameters["erev_na"]
		self.soma.ek=self.parameters["erev_k"]

		#self.soma.alpha_A_m_dIN_na = 8.67 -1
		#self.soma.alpha_E_m_dIN_na = -12.56 -1
		#self.soma.alpha_A_h_dIN_na = 0.08 -0.03

class dINcell_Li2014(Cell):
	def __init__(self,RecAll=2,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="dIN"
		self.color="#96501e"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","dIN_na":"ina_dIN_na","dIN_kSlow":"ik_dIN_kSlow","dIN_ca":"ica_dIN_ca"}

		# init cell		
		Cell.__init__(self, cap=1.0*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# parameters for dINs
		self.parameters = {"erev_lk": -52.0, 
		                   "gmax_lk": 1.405*molt_gmax,
		                   "erev_na": 50.0,
		                   "gmax_na": 240.5*molt_gmax, # Down to 210.5*molt_gmax reduces mid-cycle dINs
		                   "erev_k": -81.5,
		                   "gmax_ks": 150.0*molt_gmax,
		                   "T_ca": 300.0,
		                   "Sin_ca": 100e-9,
		                   "Sout_ca": 10e-6,
		                   "perm_ca": 1.425e-10/1e-5} # this gives the same results as Bob's dINs

		# ion channels
		self.soma.cm=1.0*molt_cap
		self.soma.g_pas=self.parameters["gmax_lk"]
		self.soma.e_pas=self.parameters["erev_lk"]
		self.soma.gmax_dIN_na=self.parameters["gmax_na"]
		self.soma.gmax_dIN_kSlow=self.parameters["gmax_ks"]
		self.soma.T_dIN_ca=self.parameters["T_ca"]
		self.soma.perm_dIN_ca=self.parameters["perm_ca"]
		self.soma.Sin_dIN_ca=self.parameters["Sin_ca"]
		self.soma_Sout_dIN_ca=self.parameters["Sout_ca"]
		self.soma.ena=self.parameters["erev_na"]
		self.soma.ek=self.parameters["erev_k"]

class dINcell_test(Cell):
	def __init__(self,RecAll=2,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="dIN"
		self.color="#96501e"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","dIN_na":"ina_dIN_na","dIN_kSlow":"ik_dIN_kSlow","dIN_ca":"ica_dIN_ca"}

		# init cell		
		Cell.__init__(self, cap=1.0*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# parameters for dINs
		self.parameters = {"erev_lk": -52.0, 
		                   "gmax_lk": 1.405*molt_gmax,
		                   "erev_na": 50.0,
		                   "gmax_na": 240.5*molt_gmax, # Down to 210.5*molt_gmax reduces mid-cycle dINs
		                   "erev_k": -80.0,
		                   "gmax_ks": 150.0*molt_gmax,
		                   "T_ca": 300.0,
		                   "Sin_ca": 100e-9,
		                   "Sout_ca": 10e-6,
		                   "perm_ca": 1.425e-10/1e-5} # this gives the same results as Bob's dINs

		# ion channels
		self.soma.cm=1.0*molt_cap
		self.soma.g_pas=self.parameters["gmax_lk"]
		self.soma.e_pas=self.parameters["erev_lk"]
		self.soma.gmax_dIN_na=self.parameters["gmax_na"]
		self.soma.gmax_dIN_kSlow=self.parameters["gmax_ks"]
		self.soma.T_dIN_ca=self.parameters["T_ca"]
		self.soma.perm_dIN_ca=self.parameters["perm_ca"]
		self.soma.Sin_dIN_ca=self.parameters["Sin_ca"]
		self.soma_Sout_dIN_ca=self.parameters["Sout_ca"]
		self.soma.ena=self.parameters["erev_na"]
		self.soma.ek=self.parameters["erev_k"]



class dINcell_sautois(Cell):
	def __init__(self,RecAll=1,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="dIN"
		self.color="#96501e"

		# channel properties and cell init
		self.channels = {"pas":"i_pas","MN_na":"ina_MN_na","MN_kFast":"ik_MN_kFast","MN_kSlow":"ik_MN_kSlow"}

		# init cell		
		Cell.__init__(self, cap=0.4*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# parameters for dINs
		self.parameters = {"erev_lk": -51.0, 
		                   "gmax_lk": 3.6765*molt_gmax,
		                   "erev_na": 50.0,
		                   "gmax_na": 210.00*molt_gmax, # Down to 210.5*molt_gmax reduces mid-cycle dINs
		                   "erev_k": -80,
		                   "gmax_kf": 0.5*molt_gmax,
		                   "gmax_ks": 3.0*molt_gmax} # this gives the same results as Bob's dINs

		# ion channels
		self.soma.cm=0.4*molt_cap
		self.soma.g_pas=self.parameters["gmax_lk"]
		self.soma.e_pas=self.parameters["erev_lk"]
		self.soma.gmax_MN_na=self.parameters["gmax_na"]
		self.soma.gmax_MN_kFast=self.parameters["gmax_kf"]
		self.soma.gmax_MN_kSlow=self.parameters["gmax_ks"]
		self.soma.ena=self.parameters["erev_na"]
		self.soma.ek=self.parameters["erev_k"]
		

		# Na constants
		self.soma.alpha_A_m_MN_na = 13.01
		self.soma.alpha_B_m_MN_na = 0.0
		self.soma.alpha_C_m_MN_na = 4.0
		self.soma.alpha_D_m_MN_na = -1.01
		self.soma.alpha_E_m_MN_na = -12.56
	    	
		self.soma.beta_A_m_MN_na = 5.73 
		self.soma.beta_B_m_MN_na = 0.0
		self.soma.beta_C_m_MN_na = 1.0
		self.soma.beta_D_m_MN_na = 9.01
		self.soma.beta_E_m_MN_na = 9.69

		self.soma.alpha_A_h_MN_na = 0.06 
		self.soma.alpha_B_h_MN_na = 0.0
		self.soma.alpha_C_h_MN_na = 0.0
		self.soma.alpha_D_h_MN_na = 30.88
		self.soma.alpha_E_h_MN_na = 26.0
	    	
		self.soma.beta_A_h_MN_na = 3.06
		self.soma.beta_B_h_MN_na = 0.0
		self.soma.beta_C_h_MN_na = 1.0
		self.soma.beta_D_h_MN_na = -7.09
		self.soma.beta_E_h_MN_na = -10.21 

		# Kf constants
		self.soma.alpha_A_MN_kFast = 3.1
		self.soma.alpha_B_MN_kFast = 0.0
		self.soma.alpha_C_MN_kFast = 1.0
		self.soma.alpha_D_MN_kFast = -31.5
		self.soma.alpha_E_MN_kFast = -9.3
	    	
		self.soma.beta_A_MN_kFast = 0.44
		self.soma.beta_B_MN_kFast = 0.0
		self.soma.beta_C_MN_kFast = 1.0
		self.soma.beta_D_MN_kFast = 4.98
		self.soma.beta_E_MN_kFast = 16.19

		# Ks constants
		self.soma.alpha_A_MN_kSlow = 0.2
		self.soma.alpha_B_MN_kSlow = 0.0
		self.soma.alpha_C_MN_kSlow = 1.0 
		self.soma.alpha_D_MN_kSlow = -6.96
		self.soma.alpha_E_MN_kSlow = -7.74
	    	
		self.soma.beta_A_MN_kSlow = 0.05
		self.soma.beta_B_MN_kSlow = 0.0
		self.soma.beta_C_MN_kSlow = 2.0
		self.soma.beta_D_MN_kSlow = -18.07
		self.soma.beta_E_MN_kSlow = 6.1

class dINcell_hull2015(Cell):
	def __init__(self,RecAll=2,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="dIN"
		self.color="#96501e"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","dIN_na":"ina_dIN_na","dIN_kFast":"ik_dIN_kFast","dIN_kSlow":"ik_dIN_kSlow","dIN_ca":"ica_dIN_ca"}

		# init cell		
		Cell.__init__(self, cap=1.0*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# parameters for dINs
		self.parameters = {"erev_lk": -52.0, 
		                   "gmax_lk": 2.5*molt_gmax,
		                   "erev_na": 50.0,
		                   "gmax_na": 300*molt_gmax, 
		                   "erev_k": -81.5,
		                   "gmax_kf": 25.0*molt_gmax,
		                   "gmax_ks": 20.0*molt_gmax,
		                   "T_ca": 300.0,
		                   "Sin_ca": 100e-9,
		                   "Sout_ca": 10e-6,
		                   "perm_ca": 1.6e-10/1e-5} 

		# ion channels
		self.soma.cm=1.0*molt_cap
		self.soma.g_pas=self.parameters["gmax_lk"]
		self.soma.e_pas=self.parameters["erev_lk"]
		self.soma.gmax_dIN_na=self.parameters["gmax_na"]
		self.soma.gmax_dIN_kFast=self.parameters["gmax_kf"]
		self.soma.gmax_dIN_kSlow=self.parameters["gmax_ks"]
		self.soma.T_dIN_ca=self.parameters["T_ca"]
		self.soma.perm_dIN_ca=self.parameters["perm_ca"]
		self.soma.Sin_dIN_ca=self.parameters["Sin_ca"]
		self.soma_Sout_dIN_ca=self.parameters["Sout_ca"]
		self.soma.ena=self.parameters["erev_na"]
		self.soma.ek=self.parameters["erev_k"]



class MNcell(Cell):
	def __init__(self,RecAll=2,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="mn"
		self.color="#00963c"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","MN_na":"ina_MN_na","MN_kFast":"ik_MN_kFast","MN_kSlow":"ik_MN_kSlow"}

		# init cell		
		Cell.__init__(self, cap=1.0*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# parameters for dINs
		self.parameters = {"L": 100.0,
		                  "erev_lk": -61.0, 
		                  "gmax_lk": 2.4691*molt_gmax, 
		                  "erev_na": 50.0,
		                  "gmax_na": 110*molt_gmax,
		                  "erev_k": -80.0,
		                  "gmax_kf": 8*molt_gmax,
		                  "gmax_ks": 1*molt_gmax}

		# ion channels
		self.soma.g_pas=self.parameters["gmax_lk"]
		self.soma.e_pas=self.parameters["erev_lk"]
		self.soma.gmax_MN_na=self.parameters["gmax_na"]
		self.soma.gmax_MN_kFast=self.parameters["gmax_kf"]
		self.soma.gmax_MN_kSlow=self.parameters["gmax_ks"]
		self.soma.ena=self.parameters["erev_na"]
		self.soma.ek=self.parameters["erev_k"]


class RBcell(Cell):
	def __init__(self,RecAll=2,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="rb"
		self.color="#ffd232"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","MN_na":"ina_MN_na","MN_kFast":"ik_MN_kFast","MN_kSlow":"ik_MN_kSlow"}

		# init cell		
		Cell.__init__(self, cap=.4*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, theta=theta, channels=self.channels)

		# Na constants
		self.soma.cm = 0.4*molt_cap
		self.soma.g_pas = 4.3573*molt_gmax
		self.soma.e_pas = -70.0
		
		self.soma.gmax_MN_na=120.0*molt_gmax 
		self.soma.gmax_MN_kFast=1.5*molt_gmax
		self.soma.gmax_MN_kSlow=8.0*molt_gmax

		# Na constants
		self.soma.alpha_A_m_MN_na = 13.01
		self.soma.alpha_B_m_MN_na = 0.0
		self.soma.alpha_C_m_MN_na = 1.0
		self.soma.alpha_D_m_MN_na = -1.01
		self.soma.alpha_E_m_MN_na = -12.56
	    	
		self.soma.beta_A_m_MN_na = 5.73 
		self.soma.beta_B_m_MN_na = 0.0
		self.soma.beta_C_m_MN_na = 1.0
		self.soma.beta_D_m_MN_na = 6.01
		self.soma.beta_E_m_MN_na = 9.69

		self.soma.alpha_A_h_MN_na = 0.06 
		self.soma.alpha_B_h_MN_na = 0.0
		self.soma.alpha_C_h_MN_na = 0.0
		self.soma.alpha_D_h_MN_na = 29.88
		self.soma.alpha_E_h_MN_na = 26.0
	    	
		self.soma.beta_A_h_MN_na = 2.04
		self.soma.beta_B_h_MN_na = 0.0
		self.soma.beta_C_h_MN_na = 1.0
		self.soma.beta_D_h_MN_na = -8.09
		self.soma.beta_E_h_MN_na = -10.21 

		# Kf constants
		self.soma.alpha_A_MN_kFast = 3.1
		self.soma.alpha_B_MN_kFast = 0.0
		self.soma.alpha_C_MN_kFast = 1.0
		self.soma.alpha_D_MN_kFast = -32.5
		self.soma.alpha_E_MN_kFast = -9.3
	    	
		self.soma.beta_A_MN_kFast = 0.44 
		self.soma.beta_B_MN_kFast = 0.0
		self.soma.beta_C_MN_kFast = 1.0
		self.soma.beta_D_MN_kFast = 3.98
		self.soma.beta_E_MN_kFast = 16.19

		# Ks constants
		self.soma.alpha_A_MN_kSlow = 0.2
		self.soma.alpha_B_MN_kSlow = 0.0
		self.soma.alpha_C_MN_kSlow = 1.0
		self.soma.alpha_D_MN_kSlow = -7.96
		self.soma.alpha_E_MN_kSlow = -7.74
	    	
		self.soma.beta_A_MN_kSlow = 0.05
		self.soma.beta_B_MN_kSlow = 0.0
		self.soma.beta_C_MN_kSlow = 2.0
		self.soma.beta_D_MN_kSlow = -19.07
		self.soma.beta_E_MN_kSlow = 6.1


class tINcell(Cell):
	def __init__(self,RecAll=1,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		# type specific properties
		self.whatami="tIN"
		self.color="#008080"
		
		# channel properties and cell init
		self.channels = {"pas":"i_pas","MN_na":"ina_MN_na","MN_kFast":"ik_MN_kFast","MN_kSlow":"ik_MN_kSlow"}

		# init cell		
		Cell.__init__(self, cap=.4*molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)

		# Na constants
		self.soma.cm = 0.4*molt_cap
		self.soma.g_pas = 2.1786*molt_gmax
		self.soma.e_pas = -55.0 
		
		self.soma.gmax_MN_na = 680.0*molt_gmax
		self.soma.gmax_MN_kFast = 40.0*molt_gmax
		self.soma.gmax_MN_kSlow = 20.0*molt_gmax

		d=15
		# Na constants
		self.soma.alpha_A_m_MN_na = 8.67
		self.soma.alpha_B_m_MN_na = 0.0
		self.soma.alpha_C_m_MN_na = 0.5
		self.soma.alpha_D_m_MN_na = -13.01-d
		self.soma.alpha_E_m_MN_na = -18.56
	    	
		self.soma.beta_A_m_MN_na = 5.73 -2
		self.soma.beta_B_m_MN_na = 0.0
		self.soma.beta_C_m_MN_na = 1.0
		self.soma.beta_D_m_MN_na = -2.99-d
		self.soma.beta_E_m_MN_na = 9.69

		self.soma.alpha_A_h_MN_na = 0.04 
		self.soma.alpha_B_h_MN_na = 0.0
		self.soma.alpha_C_h_MN_na = 0.0
		self.soma.alpha_D_h_MN_na = 15.8-d
		self.soma.alpha_E_h_MN_na = 26.0
	    	
		self.soma.beta_A_h_MN_na = 4.08
		self.soma.beta_B_h_MN_na = 0.0
		self.soma.beta_C_h_MN_na = 0.001
		self.soma.beta_D_h_MN_na = -19.09-d
		self.soma.beta_E_h_MN_na = -10.21 

		# Kf constants
		self.soma.alpha_A_MN_kFast = 3.1
		self.soma.alpha_B_MN_kFast = 0.0
		self.soma.alpha_C_MN_kFast = 1.0
		self.soma.alpha_D_MN_kFast = -35.5-d
		self.soma.alpha_E_MN_kFast = -9.3
	    	
		self.soma.beta_A_MN_kFast = 1.1 
		self.soma.beta_B_MN_kFast = 0.0
		self.soma.beta_C_MN_kFast = 1.0
		self.soma.beta_D_MN_kFast = 0.98-d
		self.soma.beta_E_MN_kFast = 16.19

		# Ks constants
		self.soma.alpha_A_MN_kSlow = 0.2
		self.soma.alpha_B_MN_kSlow = 0.0
		self.soma.alpha_C_MN_kSlow = 1.0 
		self.soma.alpha_D_MN_kSlow = -10.96-d
		self.soma.alpha_E_MN_kSlow = -7.74
	    	
		self.soma.beta_A_MN_kSlow = 0.05
		self.soma.beta_B_MN_kSlow = 0.0
		self.soma.beta_C_MN_kSlow = 1.0
		self.soma.beta_D_MN_kSlow = -22.07-d
		self.soma.beta_E_MN_kSlow = 6.1


class xIN_one_comp(Cell):
	
	def __init__(self,RecAll=2,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):
		# type specific properties
		self.whatami="xin"
		self.color="#008080"
		# channel properties and cell init
		self.channels = {"pas":"i_pas","MN_na":"ina_MN_na","MN_kFast":"ik_MN_kFast","MN_kSlow":"ik_MN_kSlow"}
		# init cell		
		Cell.__init__(self, cap=molt_cap, RecAll=RecAll, varDt=varDt, atol=atol, rtol=rtol, channels=self.channels)
		# parameters 
		self.parameters = {"L": 100.0,
		                  "erev_lk": -61.0, 
		                  "gmax_lk": 2.4691*molt_gmax, 
		                  "erev_na": 50.0,
		                  "gmax_na": 110*molt_gmax,
		                  "erev_k": -80.0,
		                  "gmax_kf": 8*molt_gmax,
		                  "gmax_ks": 1*molt_gmax}
		self.soma.cm=1*molt_cap/2.0
		self.soma.insert("pas")
		self.soma.g_pas=self.parameters["gmax_lk"]
		self.soma.e_pas=self.parameters["erev_lk"]
		self.soma.insert("MN_na")
		self.soma.gmax_MN_na=5*self.parameters["gmax_na"]
		self.soma.insert("MN_kFast")
		self.soma.gmax_MN_kFast=5*self.parameters["gmax_kf"]
		self.soma.insert("MN_kSlow")
		self.soma.gmax_MN_kSlow=5*self.parameters["gmax_ks"]
		self.soma.ena=self.parameters["erev_na"]
		self.soma.ek=self.parameters["erev_k"]

class HHcell:
	def __init__(self,RecAll=1,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):
		# type specific properties
		self.whatami="hh"
		self.color="k"
		self.index="None"
		self.pos=0.0
		self.theta=theta
		self.RecAll=RecAll

		# soma properties
		self.soma=h.Section(cell=self)
		self.soma.nseg=1
		self.soma.L=100.0 
		self.surface_area=1e3
		self.soma.diam=self.surface_area/(np.pi*self.soma.L)
		self.soma.cm=1.0
		
		self.soma.insert("hh")
		self.soma.el_hh=-54.3 # like Roman's model
		#self.soma.gnabar_hh=self.soma.gnabar_hh*molt_gmax
		#self.soma.gkbar_hh=self.soma.gkbar_hh*molt_gmax
		#self.soma.gl_hh=self.soma.gl_hh*molt_gmax

		# set current clamp and voltage clamp 
		self.CC = h.IClamp(self.soma(0.5))
		self.VC = h.VClamp(self.soma(0.5))

		# vector of connections
		self.connlist = []

		# record spike time stamps
		self.threshold=0.0
	        self.record = {}
	        self.nc_spike = h.NetCon(self.soma(0.5)._ref_v,None,self.threshold,0.0,0.0,sec=self.soma)
	        self.record['spk'] = h.Vector()
	        self.nc_spike.record(self.record['spk'])

		if varDt:
			# setup the variable time step integrator
			self.Hines = h.CVode()
			self.Hines.active(1)
			#self.Hines.use_local_dt(1)
			self.Hines.atol(atol)
			self.Hines.rtol(rtol)

		if RecAll>0:
			self.record['t']=h.Vector()
			self.record['vm']=h.Vector()
			if varDt:
				self.Hines.record(self.soma(0.5)._ref_v,self.record['vm'],self.record['t'],sec=self.soma)
			else:
				self.record['t'].record(h._ref_t)
				self.record['vm'].record(self.soma(0.5)._ref_v)

	def connect(self,dest,InputSynName=None,w=0.0,delay=0.0,gmax=0.0):
		connection(self,dest=dest,InputSynName=InputSynName,w=w,delay=delay,gmax=gmax)

	# destroy
	def destroy(self):
		del self.record
		del self.soma


class xINcell:
	def __init__(self,RecAll=2,varDt=False,atol=1e-5,rtol=1e-5,theta=0.0):

		'''
		varDt: True set the use of local time step integrator, False uses predefined step dt
		RecAll: 0 save nothing, 1 saves time and voltage, 2 saves also synaptic currents
        	'''

		self.RecAll = RecAll
		self.side=0 # 0=unassigned, 1=left, 2=right

		# parameters 
		self.parameters = {"L": 100.0,
		                  "erev_lk": -61.0, 
		                  "gmax_lk": 2.4691*molt_gmax, 
		                  "erev_na": 50.0,
		                  "gmax_na": 110*molt_gmax,
		                  "erev_k": -80.0,
		                  "gmax_kf": 8*molt_gmax,
		                  "gmax_ks": 1*molt_gmax}

		# type specific properties
		self.whatami="xin"
		self.color='#FF00FF'
		self.index='None'
		self.pos=0.0
		
		# soma properties
		self.soma=h.Section(cell=self)
		self.soma.nseg=1
		self.soma.L=100.0
		self.soma_surface_area=1e3 
		self.soma.diam=self.soma_surface_area/(np.pi*self.soma.L)
		self.soma.cm=molt_cap/2.0
		self.soma.Ra=100*1e2/(4*np.pi) # 100 is MOhm of total axial resistance [R=Ra*L/(pi*(d/2)^2)]
		
		self.soma.insert("pas")
		self.soma.g_pas=self.parameters["gmax_lk"]
		self.soma.e_pas=self.parameters["erev_lk"]
		self.soma.insert("MN_na")
		self.soma.gmax_MN_na=self.parameters["gmax_na"]
		self.soma.insert("MN_kFast")
		self.soma.gmax_MN_kFast=self.parameters["gmax_kf"]
		self.soma.insert("MN_kSlow")
		self.soma.gmax_MN_kSlow=self.parameters["gmax_ks"]
		self.soma.ena=self.parameters["erev_na"]
		self.soma.ek=self.parameters["erev_k"]

		# axon properties
		self.axon=h.Section(cell=self)
		self.axon.nseg=1
		self.axon.L=100.0
		self.axon_surface_area=1e3 
		self.axon.diam=self.axon_surface_area/(np.pi*self.soma.L)
		self.axon.cm=molt_cap/2.0
		self.axon.Ra=100*1e2/(4*np.pi) # 100 is MOhm of total axial resistance [R=Ra*L/(pi*(d/2)^2)]

		self.axon.insert("pas")
		self.axon.g_pas=self.parameters["gmax_lk"]
		self.axon.e_pas=self.parameters["erev_lk"]
		self.axon.insert("MN_na")
		self.axon.gmax_MN_na=5*self.parameters["gmax_na"]
		self.axon.insert("MN_kFast")
		self.axon.gmax_MN_kFast=5*self.parameters["gmax_kf"]
		self.axon.insert("MN_kSlow")
		self.axon.gmax_MN_kSlow=5*self.parameters["gmax_ks"]
		self.axon.ena=self.parameters["erev_na"]
		self.axon.ek=self.parameters["erev_k"]

		self.soma.connect(self.axon,0,0)

		# set current clamp and voltage clamp 
		self.CC = h.IClamp(self.soma(0.5))
		self.VC=h.VClamp(self.soma(0.5))

		# vector of connections
		self.connlist = []

		# record spike time stamps
		self.threshold=0.0
	        self.record = {}
	        self.nc_spike = h.NetCon(self.axon(0.5)._ref_v,None,self.threshold,0.0,0.0,sec=self.axon)
	        self.record['spk'] = h.Vector()
	        self.nc_spike.record(self.record['spk'])

		if varDt:
			# setup the variable time step integrator
			self.Hines = h.CVode()
			self.Hines.active(1)
			#self.Hines.use_local_dt(1)
			self.Hines.atol(atol)
			self.Hines.rtol(rtol)

		if RecAll>0:
			self.record['t']=h.Vector()
			self.record['t_axon']=h.Vector()
			self.record['vm']=h.Vector()
			self.record['vm_axon']=h.Vector()
			if varDt:
				self.Hines.record(self.soma(0.5)._ref_v,self.record['vm'],self.record['t'],sec=self.soma)
				self.Hines.record(self.axon(0.5)._ref_v,self.record['vm_axon'],self.record['t_axon'],sec=self.axon)
			else:
				self.record['t'].record(h._ref_t)
				self.record['vm'].record(self.soma(0.5)._ref_v)
				self.record['t_axon'].record(h._ref_t)
				self.record['vm_axon'].record(self.axon(0.5)._ref_v)
			if RecAll>1:
				self.record['IAMPA']=h.Vector()
				self.record['IAMPA_STD']=h.Vector()
				self.record['INMDA']=h.Vector()
				self.record['INMDA_STD']=h.Vector()
				self.record['INH']=h.Vector()
				self.record['INH_STD']=h.Vector()
		
	# set the connections previously appensed
	def connect(self,dest,InputSynName,w=0.0,delay=0.0,gmax=0.0,seed=None):
		if InputSynName=='ampa':
			if hasattr(dest,"syn_ampa")==0:
				dest.syn_ampa = h.syn_ampa(dest.soma(0.5))
			self.connlist.append(h.NetCon(self.axon(0.5)._ref_v,dest.syn_ampa,0.0,delay,w,sec=self.axon))
			if self.RecAll==2:
				dest.record['IAMPA'].record(dest.syn_ampa._ref_i,sec=self.soma)
		if InputSynName=='ampa_std':
			if hasattr(dest,"syn_ampa_std")==0:
				dest.syn_ampa_std = h.syn_ampa_std(dest.soma(0.5))
			self.connlist.append(h.NetCon(self.axon(0.5)._ref_v,dest.syn_ampa_std,0.0,delay,w,sec=self.axon))
			if self.RecAll==2:
				dest.record['IAMPA_STD'].record(dest.syn_ampa._ref_i,sec=self.soma)
		if InputSynName=='nmda':
			if hasattr(dest,"syn_nmda")==0:
				dest.syn_nmda = h.syn_nmda(dest.soma(0.5))
			self.connlist.append(h.NetCon(self.axon(0.5)._ref_v,dest.syn_nmda,0.0,delay,w,sec=self.axon))
			if self.RecAll==2:
				dest.record['INMDA'].record(dest.syn_nmda._ref_i,sec=self.soma)
		if InputSynName=='nmda_std':
			if hasattr(dest,"syn_nmda_std")==0:
				dest.syn_nmda_std = h.syn_nmda_std(dest.soma(0.5))
			self.connlist.append(h.NetCon(self.axon(0.5)._ref_v,dest.syn_nmda_std,0.0,delay,w,sec=self.axon))
			if self.RecAll==2:
				dest.record['INMDA_STD'].record(dest.syn_nmda_std._ref_i,sec=self.soma)
		if InputSynName=='ampa_var':
			if hasattr(dest,"syn_ampa_var")==0:
				dest.syn_ampa_var = h.syn_ampa_var(dest.soma(0.5))
				dest.rnd_ampa=h.Random(seed)
				dest.rnd_ampa.uniform(0,1)
				dest.syn_ampa_var.noiseFromRandom(dest.rnd_ampa)
				#dest.syn_ampa_var.seed(seed)
			self.connlist.append(h.NetCon(self.axon(0.5)._ref_v,dest.syn_ampa_var,0.0,delay,w,sec=self.axon))
			if self.RecAll==2:
				dest.record['IAMPA_VAR'].record(dest.syn_ampa_var._ref_i,sec=self.soma)
		if InputSynName=='nmda_var':
			if hasattr(dest,"syn_nmda_var")==0:
				dest.syn_nmda_var = h.syn_nmda_var(dest.soma(0.5))
				dest.rnd_nmda=h.Random(seed)
				dest.rnd_nmda.uniform(0,1)
				dest.syn_nmda_var.noiseFromRandom(dest.rnd_nmda)
				#dest.syn_nmda_var.seed(seed)
			self.connlist.append(h.NetCon(self.axon(0.5)._ref_v,dest.syn_nmda_var,0.0,delay,w,sec=self.axon))
			if self.RecAll==2:
				dest.record['INMDA_VAR'].record(dest.syn_nmda_var._ref_i,sec=self.soma)
		if InputSynName=='inh':
			if hasattr(dest,"syn_inh")==0:
				dest.syn_inh = h.syn_inh(dest.soma(0.5))
			self.connlist.append(h.NetCon(self.axon(0.5)._ref_v,dest.syn_inh,0.0,delay,w,sec=self.axon))
			if self.RecAll==2:
				dest.record['INH'].record(dest.syn_inh._ref_i,sec=self.soma)
		if InputSynName=='inh_std':
			if hasattr(dest,"syn_inh_std")==0:
				dest.syn_inh_std = h.syn_inh_std(dest.soma(0.5))
			self.connlist.append(h.NetCon(self.axon(0.5)._ref_v,dest.syn_inh_std,0.0,delay,w,sec=self.axon))
			if self.RecAll==2:
				dest.record['INH_STD'].record(dest.syn_inh._ref_i,sec=self.soma)
		if InputSynName=='nmda_sat':
			if hasattr(dest,"syn_nmda_sat")==0:
				dest.syn_nmda_sat = h.syn_nmda_sat(dest.soma(0.5))
			dest.connlist.append(h.NetCon(self.axon(0.5)._ref_v,dest.syn_nmda_sat,0.0,delay,w,sec=self.axon))
			if dest.RecAll==2:
				dest.record['i_nmda_sat']=h.Vector() 
				dest.record['i_nmda_sat'].record(dest.syn_nmda_sat._ref_i,sec=self.soma)
		if InputSynName=='ampa_sat':
			if hasattr(dest,"syn_ampa_sat")==0:
				dest.syn_ampa_sat = h.syn_ampa_sat(dest.soma(0.5))
			dest.connlist.append(h.NetCon(self.axon(0.5)._ref_v,dest.syn_ampa_sat,0.0,delay,w,sec=self.axon))
			if dest.RecAll==2:
				dest.record['i_ampa_sat']=h.Vector() 
				dest.record['i_ampa_sat'].record(dest.syn_nmda_sat._ref_i,sec=self.soma)
	# destroy
	def destroy(self):
		del self.record
		del self.connlist
		del self.soma
		del self.axon

















		
