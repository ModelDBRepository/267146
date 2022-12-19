
from nrnTemplate.CellTypes.CellTemplateCommon import dINcell, MNcell, aINcell, xINcell, RBcell, tINcell, MHRcell, HHcell, dINcell_hull2015
import matplotlib.pylab as plt
import numpy as np
import random
import numpy.random as rnd
from mpi.mp_util import UniqueProcessMap
from util import plot_matrix_detailed, probabilistic_model_extension, inj_current, inj_current2, plotLeftRightSpikeTrain, plotLeftRightVoltageExt, plotLeftRightVoltageOffset, TimeVoltageTrace, classification_name, classify_behaviour
import time, datetime, os
import param
from shutil import copyfile
from neuron import h

# create a directory in figures/ named after date and time to store output files
path_tmp = "figures/"
today = datetime.datetime.now()
todaystr = today.isoformat()
os.mkdir(path_tmp+todaystr)
save_path = path_tmp+todaystr+"/"

# save files of specification for parameters
file_main = "main.py"
if ~os.path.isfile(save_path+file_main):
	copyfile(file_main,save_path+file_main)
file_util = "util.py"
if ~os.path.isfile(save_path+file_util):
	copyfile(file_util,save_path+file_util)
file_param = "param.py"
if ~os.path.isfile(save_path+file_param):
	copyfile(file_param,save_path+file_param)

par = param.create_params()

RecAll = par.RecAll
varDt = par.varDt
atol = par.atol
rtol = par.rtol
dt = par.dt
time_end = par.time_end

cell_types = par.cell_types
num_types = par.num_types
halves = par.halves
vect_index = par.vect_index
left_index = par.left_index
right_index = par.right_index
colors = par.colors

fixed_delay = par.fixed_delay
sa_prop = par.sa_prop
var_delay = par.var_delay
std_on = par.std_on

pos = par.pos
w = par.w
alpha = par.alpha
beta1 = par.beta1
beta2 = par.beta2

v = par.v
w = par.w
n_active = par.n_active
pos_active = par.pos_active
delay = par.delay
duration = par.duration
amplitude_mean = par.amplitude_mean
amplitude_std = par.amplitude_std


# Cells Creation 
def CellCreation(RecAll=0,varDt=False):
	cellist=[]
	for i in xrange(len(cell_types)):
		for j in vect_index[i]:
			if cell_types[i] in ("rb","dla","cin","mn"):
				cellist.append(MNcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol))
			elif cell_types[i] == "dlc":
				cellist.append(MNcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol))
			elif cell_types[i] == "ain":
				cellist.append(MNcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol))
			elif cell_types[i] == "ecin":
				cellist.append(MNcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol))
			elif cell_types[i] == "din":
				cellist.append(dINcell_hull2015(RecAll=1,varDt=varDt,atol=atol,rtol=rtol,theta=0.0))
			elif cell_types[i] == "xin":
				cellist.append(xINcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol))
			elif cell_types[i] == "tst":
				cellist.append(MNcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol))
			elif cell_types[i] == "tin":
				cellist.append(tINcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol))
			elif cell_types[i] == "tsp":
				cellist.append(RBcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol))			
			elif cell_types[i] == "mhr":
				cellist.append(MHRcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol))
			elif cell_types[i] == "dinr":
				cellist.append(dINcell(RecAll=1,varDt=varDt,atol=atol,rtol=rtol,theta=0.0))
			else:
				print cell_types[i]
				raise Exception("Cell created is not included in the list of available cells")
			cellist[-1].whatami=cell_types[i]
			cellist[-1].color=colors[i]
			cellist[-1].index=len(cellist)-1
			cellist[-1].type_num=i
			cellist[-1].pos=pos[j] 
			if j in left_index[i]:
				cellist[-1].body_side=1
			else:
				cellist[-1].body_side=2
	return cellist


# create synapse between neural types
def connection(pre,post,dist):
	if pre.whatami=="xin" and post.whatami=="xin":
		key = pre.whatami + " -> " + post.whatami
		if key in w.viewkeys():
			specs = w[key]
			tmp = rnd.choice(v)
			for spec in specs:
				syn_type = spec[0]
				w_mean = spec[1]
				w_std = spec[2]
				if pre.body_side==post.body_side:
					pre.connect(post, syn_type, w=w_mean*tmp, delay=spec[3])
				else:
					pre.connect(post, syn_type, w=w_mean*tmp, delay=2.0)

			if hasattr(post,"syn_ampa_std"):
				post.syn_ampa_std.alpha = alpha
			if hasattr(post,"syn_nmda_std"):
				post.syn_nmda_std.alpha = alpha


	else:
		key = pre.whatami + " -> " + post.whatami
		if key in w.viewkeys():
			specs = w[key]
			for spec in specs:
				if spec != None:
					syn_type = spec[0]
					w_mean = spec[1]
					w_std = spec[2]
					if w_std != 0.0:
						weight = rnd.normal(w_mean,w_std*std_on)
					else:
						if pre.whatami=="xin":
							molt=rnd.choice(v) #rnd.uniform(0,1)
							weight=molt*w_mean
						else:
							weight = w_mean	

					if spec[3]=="distance":
						distance = fixed_delay+dist*var_delay
					else:
						distance = spec[3]
					pre.connect(post, syn_type, w=np.abs(weight), delay=distance)

			if hasattr(post,"syn_nmda_sat") and pre.whatami=="din" and post.whatami=="din":
				post.syn_nmda_sat.alpha = beta1
			if hasattr(post,"syn_nmda_std") and pre.whatami=="din" and post.whatami=="din":
				post.syn_nmda_std.alpha = beta1
			if hasattr(post,"syn_ampa_std") and pre.whatami=="xin" and post.whatami!="xin":
				post.syn_ampa_std.alpha = beta2
			if hasattr(post,"syn_nmda_std") and pre.whatami=="xin" and post.whatami!="xin":
				post.syn_nmda_std.alpha = beta2	


def CreateConfigAdjExtended(cellist,sim_num):
	seed1 = int(sim_num)
	rnd.seed(seed1)
	with open(save_path+"seed", "wt") as f_seed:
		f_seed.write(str(seed1))
	A = probabilistic_model_extension(sim_num)
	#plot_matrix_detailed(A) # to plot the adjacency matrix
	#plt.savefig(save_path+"A"+str(sim_num)+".png")

	seed2 = random.getrandbits(32) # int(sim_num) 
	rnd.seed(seed2)
	with open(save_path+"seed", "wt") as f_seed:
		f_seed.write(str(seed2))
	n=len(cellist)
	for i in xrange(n):
		for j in xrange(n):
			if A[i,j]:
				connection(cellist[i],cellist[j],abs(cellist[i].pos-cellist[j].pos))

	# dIN gj
	gj_strength=0.2e-3

	for i in left_index[4]:
		for j in left_index[4]:
			if abs(cellist[i].pos-cellist[j].pos)<100.0: 
				cellist[i].connect(cellist[j],"gap",gmax=gj_strength)

	for i in right_index[4]:
		for j in right_index[4]:
			if abs(cellist[i].pos-cellist[j].pos)<100.0: 
				cellist[i].connect(cellist[j],"gap",gmax=gj_strength)


# run one swimming simulation
def SwimmingSimulation(tstop,sim_num):
	t_start = time.time()
	print "Running Entire Simulation..."
	
	print "Create Cells ..."
	cellist=CellCreation(RecAll=RecAll,varDt=varDt)
	print "Cells Created"

	print "Create Connectivity ..."
	CreateConfigAdjExtended(cellist,sim_num)
	print "Connectivity Created"

	print "Run Numerical Integration..."

	inj_current(cellist[pos_active:pos_active+4],delay,duration,amp_mean=amplitude_mean,amp_std=amplitude_std) # trunk skin touch
	#inj_current(cellist[2046:2046+2],delay,duration,amp_mean=amplitude_mean,amp_std=amplitude_std) # head skin touch 
	#inj_current(cellist[2216:2216+13],450.0,400.0,amp_mean=0.2,amp_std=0.02) # head pressure (13 tSps activated)
	#mhr_stop_protocol_perrins_2002(cellist[2296:2297], delay=1500.0, dur=30.0, amp_mean=0.3,amp_std=0.0, num_impulses=5) # protocol of injection of a single MHR

	all_times = RunSim(tstop=tstop,dt=dt)

	# ===== PLOTTING ======	

	print "End of the Integration"

	start_plot = 50
	# plot rostral spike trains
	plt.figure(1,figsize=(20,10))
	plotLeftRightSpikeTrain(cellist,[start_plot,tstop])
	plt.subplot(2,1,1)
	plt.xlim([start_plot,tstop])
	plt.subplot(2,1,2)
	plt.xlim([start_plot,tstop])


	# plot selected cells
	plot_idxs=[]
	plot_types=["rb","tsp","tst","mhr","tin","dlc","dla","xin","din","cin","ain","mn"]
	for type_id in plot_types:
		cells=[cell for cell in cellist if cell.whatami==type_id and cell.body_side==1]
		tmp=[len(cell.record["spk"]) for cell in cells]
		if len(tmp)>0:
			if type_id=="xin" or type_id=="din":
				for i in xrange(3):
					idx=rnd.choice([cell.index for cell in cells])
					plot_idxs.append(idx)
			else:
				idx=tmp.index(max(tmp))
				plot_idxs.append(cells[idx].index)
		else:
			plot_idxs.append(rnd.choice([cell.index for cell in cells]))
	for type_id in plot_types:
		cells=[cell for cell in cellist if cell.whatami==type_id and cell.body_side==2]
		tmp=[len(cell.record["spk"]) for cell in cells]
		if len(tmp)>0:
			if type_id=="xin" or type_id=="din":
				for i in xrange(3):
					idx=rnd.choice([cell.index for cell in cells])
					plot_idxs.append(idx)
			else:
				idx=tmp.index(max(tmp))
				plot_idxs.append(cells[idx].index)


	# plot rostral dINs voltage
	plt.figure(4,figsize=(15,8))
	plotLeftRightVoltageOffset([cellist[i] for i in plot_idxs],[start_plot,tstop],offset=50.0)
	plt.savefig(save_path+"tst_dlcr_tin_dinr_volt"+str(sim_num)+".png")

	# saving voltages 
	tmp = []
	tmp.append([cellist[idx].color for idx in plot_idxs])

	(t,v) = TimeVoltageTrace(cellist[idx])
	tmp.append(t)
	for idx in plot_idxs:
		(t,v) = TimeVoltageTrace(cellist[idx])
		tmp.append(v)
	np.save(save_path+"voltages_one_for_each"+str(sim_num)+".npy",tmp)

	np.save(save_path+"pos_one_for_each"+str(sim_num)+".npy",[cellist[idx].pos for idx in plot_idxs])

	tmp = []
	for idx in vect_index[4]:
		if pos[idx]<1000:
			(t,v) = TimeVoltageTrace(cellist[idx])
			tmp.append(v)
	np.save(save_path+"voltages_hdIN"+str(sim_num)+".npy",tmp)

	
	# algorithm for behavioural classification
	(out,tL,tR,tstar)=classify_behaviour([np.array(cell.record["spk"]) for cell in cellist],tmp,sim_num,t)

	plt.figure(1)
	plt.title(classification_name(out)+', tL='+str(tL)+', tR='+str(tR)+', tstar='+str(tstar))
	plt.savefig(save_path+"spk_train"+str(sim_num)+".png")


	# plot rostral dINs voltage
	plt.figure(5,figsize=(15,8))
	plotLeftRightVoltageExt([cellist[i] for i in [index for index in vect_index[4] if pos[index]<1000]],[start_plot,tstop]) 
	plt.subplot(2,1,1)
	plt.ylim([-60,0])
	plt.subplot(2,1,2)
	plt.ylim([-60,0])
	plt.savefig(save_path+"din_volt"+str(sim_num)+".png")

	plt.figure(5,figsize=(20,10))
	plotLeftRightSpikeTrain(cellist,[start_plot,tstop])
	plt.subplot(2,1,1)
	plt.xlim([tstop-300,tstop])
	plt.subplot(2,1,2)
	plt.xlim([tstop-300,tstop])
	#plt.savefig(save_path+"last_spk_train"+str(sim_num)+".png")


	# save spike trains
	np.save(save_path+"spikes"+str(sim_num)+".npy",[np.array(cell.record["spk"]) for cell in cellist])

	t_end = time.time()
	print "Simulation Took {0}s.".format(t_end-t_start)


	spk_mns=[]
	pos_mns=[]
	side_mns=[]
	for cell in cellist:
		if cell.whatami=="mn":
			spk_mns.append(np.array(cell.record["spk"]))
			pos_mns.append(cell.pos)
			side_mns.append(cell.body_side)
	np.save(save_path+'spk_mns'+str(sim_num),spk_mns)
	np.save(save_path+'pos_mns'+str(sim_num),pos_mns)
	np.save(save_path+'side_mns'+str(sim_num),side_mns)
	
	t_end = time.time()
	print "Simulation Took {0}s.".format(t_end-t_start)

	# return reaction times
	if classification_name(out)==1 or classification_name(out)==2:
		react_time=tstar
	else:
		react_time=None

	destroy(cellist)

def destroy(cellist):
	for cell in cellist:
		cell.destroy()

def SimForAll(sim_num):
	seed = int(random.getrandbits(32))
	rnd.seed(seed)
	with open("seed", "wt") as f_seed:
		f_seed.write(str(seed))
	output=SwimmingSimulation(time_end,sim_num)
	return output

def RunManySim():
	num_process = 5 # number of cores to use
	processes = UniqueProcessMap(num_process)
	I = range(1,101)
	out=processes.map(SimForAll,I)
	np.save(save_path+"out.npy",out)
	return out

# run simulation using Euler/CVode method
def RunSim(v_init=-80.0,tstop=0.0,dt=0.01):
	all_times=[]
	t_start = time.time()
	h.dt = dt
	h.t = 0.0
	h.finitialize(v_init)
	while h.t<tstop:
		all_times.append(h.t)
		h.fadvance()
	return all_times


def mhr_stop_protocol_perrins_2002(cellist, delay=0.0, dur=0.0, amp_mean=0.0,amp_std=0.0, num_impulses=0):
	for cell in cellist:
		cell.impulses = []
		time_interval = 20
		for i in xrange(num_impulses):
			cell.impulses.append(h.IClamp(cell.soma(0.5)))
			cell.impulses[-1].delay = delay + (dur+time_interval)*i
			cell.impulses[-1].dur = dur
			cell.impulses[-1].amp = rnd.normal(amp_mean,amp_std)



























