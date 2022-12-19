
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm
import numpy.random as rnd
import random
from scipy.interpolate import interp1d
import param 
par = param.create_params()

cell_types = par.cell_types
num_types = par.num_types
halves = par.halves
vect_index = par.vect_index
left_index = par.left_index
right_index = par.right_index
tot_num_x_type = par.tot_num_x_type
colors = par.colors

pos = par.pos

pos_old = np.loadtxt("data/pos.txt")

def tmfr(cellist,tw,dtF,I):
	st=[]
	for cell in cellist:
		st.append(np.array(cell.record["spk"]))
	st=np.concatenate(st)
	t = np.arange(I[0],I[1],dtF)

	tmfr=[]
	for ti in t:
		tmfr.append(1.0*sum((st>ti)&(st<(ti+tw)))/len(cellist))
	return (t,tmfr)

def plot_spikes_data(st,I=''):
	st=np.array([train.tolist() for train in st])
	plt.subplot(2,1,2)
	for i in xrange(len(left_index)):
		if i!=3 and i!=4:
			plot_spk_train(pos[left_index[i]],st[left_index[i]],I,colors=colors[i])
			plt.ylim([500,2000])
			plt.ylabel("Left")
	plt.subplot(2,1,1)
	for i in xrange(len(right_index)):
		if i!=3 and i!=4:
			plot_spk_train(pos[right_index[i]],st[right_index[i]],I,colors=colors[i])
			plt.ylim([500,2000])
			plt.ylabel("Right")

def SpikeTrain(cellist):
	return [np.array(cell.record['spk']) for cell in cellist]

def TimeVoltageTrace(cell):
	return (np.array(cell.record['t']),np.array(cell.record['vm']))

def plot_spk_train(pos,st,I='',colors=''):
	for (i,spike_train) in enumerate(st):
		plt.plot(spike_train,pos[i]*np.ones_like(spike_train),marker='.',markersize=5,markerfacecolor=colors,markeredgecolor=colors,linestyle='None')
	if len(I):
		plt.xlim(I)

def inj_current(cells,delay,dur,amp_mean=0,amp_std=0):
	for cell in cells:
		cell.CC.delay = delay
		cell.CC.dur = dur
		cell.CC.amp = rnd.normal(amp_mean,amp_std)


def inj_current2(cells,delay,dur,amp_mean=0,amp_std=0):
	for cell in cells:
		cell.CC2.delay = delay
		cell.CC2.dur = dur
		cell.CC2.amp = rnd.normal(amp_mean,amp_std)

def plotLeftRightSpikeTrainType(cellist,cell_type,I):
	st=np.array([train.tolist() for train in SpikeTrain(cellist)])
	plt.subplot(2,1,2)
	for i in xrange(len(left_index)):
		plot_spk_train([cell.pos for cell in [cellist[x] for x in left_index[i]]],st[left_index[i]],I,colors=colors[i])
		plt.ylim([500,3620])
		plt.ylabel("Left")
	plt.subplot(2,1,1)
	for i in xrange(len(right_index)):
		plot_spk_train([cell.pos for cell in [cellist[x] for x in right_index[i]]],st[right_index[i]],I,colors=colors[i])
		plt.ylim([500,3620])
		plt.ylabel("Right")


def plotLeftRightSpikeTrain(cellist,I):
	st=np.array([train.tolist() for train in SpikeTrain(cellist)])
	plt.subplot(2,1,2)
	for i in xrange(len(left_index)):
		if i!=25: # no mn plot for alan
			plot_spk_train([cell.pos for cell in [cellist[x] for x in left_index[i]]],st[left_index[i]],I,colors=colors[i])
			plt.ylim([-131,3620])
			plt.ylabel("Left")
	plt.subplot(2,1,1)
	for i in xrange(len(right_index)):
		if i!=25: # no mn plot for alan
			plot_spk_train([cell.pos for cell in [cellist[x] for x in right_index[i]]],st[right_index[i]],I,colors=colors[i])
			plt.ylim([-131,3620])
			plt.ylabel("Right")

def plotLeftRightVoltageOffset(cellist,I,offset=0.0):
	n = len(cellist)
	for i in xrange(n):
		if i<n/2:
			plt.subplot(2,1,2)
			(t,v) = TimeVoltageTrace(cellist[i])
			plt.plot(t,v+offset*i,color=cellist[i].color,linewidth=1.0)
		else:
			plt.subplot(2,1,1)
			(t,v) = TimeVoltageTrace(cellist[i])
			plt.plot(t,v+offset*(i-n/2),color=cellist[i].color,linewidth=1.0)
	plt.subplot(2,1,1)
	plt.ylim([-80,40+offset*n/2])
	plt.xlim(I)
	plt.ylabel("Right")
	plt.subplot(2,1,2)
	plt.xlim(I)
	plt.ylabel("Left")
	plt.ylim([-80,40+offset*n/2])

def plotLeftRightVoltage(cellist,I):
	n = len(cellist)
	for i in xrange(n):
		if cellist[i].body_side==1:
			plt.subplot(4,1,4)
			(t,v) = TimeVoltageTrace(cellist[i])
			plt.plot(t,v,color=cellist[i].color,linewidth=1.0)
		elif cellist[i].body_side==2:
			plt.subplot(4,1,1)
			(t,v) = TimeVoltageTrace(cellist[i])
			plt.plot(t,v,color=cellist[i].color,linewidth=1.0)
		else:
			print "Error no body side declaration"
	plt.subplot(4,1,4)
	plt.xlim(I)
	plt.ylabel("Left")
	plt.ylim([-80,40])
	plt.subplot(4,1,1)
	plt.xlim(I)
	plt.ylabel("Right")
	plt.ylim([-80,40])

def plotLeftRightVoltageExt(cellist,I):
	n = len(cellist)
	for i in xrange(n):
		if i<n/2:
			plt.subplot(2,1,2)
			(t,v) = TimeVoltageTrace(cellist[i])
			plt.plot(t,v,color=cellist[i].color,linewidth=1.0)
		else:
			plt.subplot(2,1,1)
			(t,v) = TimeVoltageTrace(cellist[i])
			plt.plot(t,v,color=cellist[i].color,linewidth=1.0)
	plt.subplot(2,1,2)
	plt.ylabel("Left")
	plt.ylim([-80,40])
	plt.xlim([I[0],I[1]])
	plt.subplot(2,1,1)
	plt.xlim([I[0],I[1]])
	plt.ylim([-80,40])
	plt.ylabel("Right")



def probabilistic_model_extension(num):

	P = np.load('data/P1000.npy')
	rnd.seed(num)

	n = num_types[-1]
	Q = np.zeros((sum(tot_num_x_type),sum(tot_num_x_type)))

	# dla->xIN
	for idx1 in left_index[6]:
		for idx2 in left_index[7]:
			Q[idx1,idx2] = 0.03 

	for idx1 in right_index[6]:
		for idx2 in right_index[7]:
			Q[idx1,idx2] = 0.03 

	# dlc->xIN
	for idx1 in left_index[1]:
		for idx2 in right_index[7]:
			Q[idx1,idx2] = 0.014 

	for idx1 in right_index[1]:
		for idx2 in left_index[7]:
			Q[idx1,idx2] = 0.014  


	# xIN->CPG
	prob = 0.05 # decreasing this param reduces sync->rest transitions
	dist_tresh = 1000 
	rostral_cpg_left = np.concatenate([[int(idx) for idx in left_index[numb] if pos[idx]<dist_tresh] for numb in [4,2,3,5] if len([int(idx) for idx in left_index[numb] if pos[idx]<dist_tresh])>0])
	rostral_cpg_right = np.concatenate([[int(idx) for idx in right_index[numb] if pos[idx]<dist_tresh] for numb in [4,2,3,5] if len([int(idx) for idx in right_index[numb] if pos[idx]<dist_tresh])>0])
	
	for i in left_index[7]:
		for j in rostral_cpg_left:
			Q[i,j] = prob
		for j in rostral_cpg_right:
			if rnd.rand()<0.:
				Q[i,j] = prob

	for i in right_index[7]:
		for j in rostral_cpg_right:
			Q[i,j] = prob
		for j in rostral_cpg_left:
			if rnd.rand()<0.:
				Q[i,j] = prob

	# xIN->xIN
	cell_num = 5 # mn 
	Psmall = P[np.ix_(par.left_index_old[cell_num],par.left_index_old[cell_num])] + P[np.ix_(par.left_index_old[cell_num],par.left_index_old[cell_num])].transpose()
	molt=2.5
	n_xin = len(left_index[7])
	for i in xrange(n_xin):
		for j in range(n_xin):
			Q[i+left_index[7][0],j+left_index[7][0]] = molt*Psmall[i,j]
			if rnd.rand()<0.33:
				Q[i+left_index[7][0],j+right_index[7][0]] = molt*Psmall[i,j]

	n_xin=len(right_index[7])
	for i in xrange(n_xin):
		for j in range(n_xin):
			Q[i+right_index[7][0],j+right_index[7][0]] = molt*Psmall[i,j]
			if rnd.rand()<0.33:
				Q[i+right_index[7][0],j+left_index[7][0]] = molt*Psmall[i,j]

	# tSt->tIN
	for i in left_index[8]:
		for j in left_index[9]:
			Q[i,j] = GeneralizeData(P[rnd.choice(par.left_index_old[0]),par.left_index_old[6]])

	for i in right_index[8]:
		for j in right_index[9]:
			Q[i,j] = GeneralizeData(P[rnd.choice(par.right_index_old[0]),par.right_index_old[6]])

	# tSt->rdlc
	for i in left_index[8]:
		for j in [idx for idx in left_index[1] if pos[idx]<700]:
			Q[i,j] = GeneralizeData(P[rnd.choice(par.left_index_old[0]),par.left_index_old[1]])

	for i in right_index[8]:
		for j in [idx for idx in right_index[1] if pos[idx]<700]:
			Q[i,j] = GeneralizeData(P[rnd.choice(par.right_index_old[0]),par.right_index_old[1]])


	# tIN->xIN
	for idx1 in left_index[9]:
		for idx2 in left_index[7]:
			Q[idx1,idx2] = 0.03 

	for idx1 in right_index[9]:
		for idx2 in right_index[7]:
			Q[idx1,idx2] = 0.03 

	
	# tIN->dIN
	for post_num in range(4,5):
		pre = par.left_index_old[6] # dla
		post = par.left_index_old[post_num]
		f = distance_dependent_prob(pre,post)
		for i in left_index[9]:
			for j in left_index[post_num]:
				x = pos[j]-pos[i]
				Q[i,j] = f(-x)


		pre = par.right_index_old[6] # dla
		post = par.right_index_old[post_num]
		f = distance_dependent_prob(pre,post)
		for i in right_index[9]:
			for j in right_index[post_num]:
				x = pos[j]-pos[i]
				Q[i,j] = f(-x)
	

	# tSp->MHR
	for i in left_index[10]:
		for j in left_index[11]:
			Q[i,j] = GeneralizeData(P[rnd.choice(par.left_index_old[0]),par.left_index_old[6]])

	for i in right_index[10]:
		for j in right_index[11]:
			Q[i,j] = GeneralizeData(P[rnd.choice(par.right_index_old[0]),par.right_index_old[6]])

	# MHR-> ipsilaterally to CPG
	for post_num in range(2,6):
		# left MHRs
		pre = par.left_index_old[1][-2:-1]
		post = par.right_index_old[post_num]
		f = distance_dependent_prob(pre,post)
		# contralateral
		for i in left_index[11]:
			for j in right_index[post_num]:
				x = pos[j]-pos[i]
				Q[i,j] = f(-x)
		# 20% ipsilateral
		for i in rnd.random_integers(left_index[11][0],left_index[11][-1],len(left_index[11])*2/10):
			for j in left_index[post_num]:
				x = pos[j]-pos[i]
				Q[i,j] = f(-x)

		# right MHRs
		pre = par.right_index_old[1][-2:-1]
		post = par.left_index_old[post_num]
		f = distance_dependent_prob(pre,post)
		# contralateral
		for i in right_index[11]:
			for j in left_index[post_num]:
				x = pos[j]-pos[i]
				Q[i,j] = f(-x)
		# 20% ipsilateral
		for i in rnd.random_integers(right_index[11][0],right_index[11][-1],len(right_index[11])*2/10):
			for j in right_index[post_num]:
				x = pos[j]-pos[i]
				Q[i,j] = f(-x)




	Aconnectome = np.load("data/A_connectome"+str(num)+".npy")
	A = np.zeros((n,n))
	for i in xrange(n):
		for j in xrange(n):
			if rnd.rand()<Q[i,j]:
				A[i,j] = 1

	# anatomical connectome
	A[np.ix_(xrange(sum(tot_num_x_type[0:8])),xrange(sum(tot_num_x_type[0:8])))]=Aconnectome;

	return A

def plot_matrix_detailed(Q):

	n=len(Q)
	num_x_type=[num_types[i+1]-num_types[i] for i in xrange(len(num_types)-1)]
	num_x_type.insert(0,0)
	halves=np.divide(num_x_type,2)
	types=cell_types
	#order_best=[0,8,10,6,1,9,11,7,2,3,4,13,5,12]
	order_best=[0,8,10,6,1,9,11,7,2,3,4,5]

	grid_color1="grey"
	grid_color2="grey"
	lw=0.5
	tr=0.5
	alpha=1

	Q=Q.transpose()
	Qnew=np.zeros(Q.shape)

	v_ind=[]
	for i in xrange(len(num_x_type)+1):
		if i!=0:
			v_ind.append(np.sum(num_x_type[0:i]))

	new_idx=[]
	for i in xrange(len(order_best)):
		if v_ind[order_best[i]+1]-v_ind[order_best[i]]>0:
			new_idx.append(range(v_ind[order_best[i]],v_ind[order_best[i]+1]))

	# === test with different colors ===
	idx_inh = np.concatenate([range(v_ind[2],v_ind[3]), range(v_ind[3],v_ind[4]), range(v_ind[11],v_ind[12])]) # aIN, cIN, MHR indexes
	Q[:,idx_inh] = - Q[:,idx_inh]

	new_idx=np.concatenate((new_idx)).tolist()
	Qnew=Q[new_idx,:]
	Qnew=Qnew[:,new_idx]

	colors_new = [colors[i] for i in order_best]
	num_x_type = [num_x_type[i+1] for i in order_best]
	num_x_type.insert(0,0)
	halves = np.divide(num_x_type,2)
	types = [cell_types[i] for i in order_best]

	fig, ax = plt.subplots(figsize=(14,14))
	ax.matshow(Qnew.transpose(),norm=mpl.colors.Normalize(vmin=-1.5, vmax=1.5), cmap='seismic')
	#plt.matshow(Qnew.transpose(),cmap='Greys_r')

	plt.plot([-100,n-0.5],[n-0.5,n-0.5],color=grid_color1,alpha=alpha)
	plt.plot([-100,n-0.5],[-100,-100],color=grid_color1,alpha=alpha)
	plt.plot([n-0.5,n-0.5],[-100,n-0.5],color=grid_color1,alpha=alpha)
	plt.plot([-100,-100],[-100,n-0.5],color=grid_color1,alpha=alpha)
	#plt.text(-60,-40,'P', fontsize=20)

	for i in xrange(len(num_x_type)):
		line_half=np.sum(num_x_type[0:i])+halves[i]-tr
		line=np.sum(num_x_type[0:i])-tr
		plt.plot([0,n-tr],[line_half,line_half],'--',color=grid_color2,linewidth=lw,alpha=alpha)
		plt.plot([-100,n-tr],[line,line],'-',color=grid_color1,linewidth=lw,alpha=alpha)
		plt.plot([line_half,line_half],[0,n-tr],'--',color=grid_color2,linewidth=lw,alpha=alpha)
		plt.plot([line,line],[-100,n-tr],'-',color=grid_color1,linewidth=lw,alpha=alpha)
	
	for i in xrange(len(types)):
		plt.text(np.sum(num_x_type[0:i+1])+halves[i+1]-tr-25,-40,types[i],fontsize=12,fontweight='bold',color=colors_new[i],rotation=65)
		plt.text(-80,np.sum(num_x_type[0:i+1])+halves[i+1]-tr,types[i],fontsize=12,fontweight='bold',color=colors_new[i],rotation=20)

	plt.xlim([-105,num_types[-1]+10])
	plt.ylim([-105,num_types[-1]+10])
	plt.axis("off")

def GeneralizeData(input_data):
	data=np.sort(input_data)
	cdf = sm.distributions.ECDF(data)
	p=np.unique(cdf(data)).tolist()
	p.insert(0,0)
	x=np.unique(data).tolist()
	x.insert(0,x[0])
	w=rnd.rand()
	r=next(tmp[0] for tmp in enumerate(p[1:]) if w<tmp[1])
	if p[r+1]-p[r] != 0:
		return x[r]+(x[r+1]-x[r])*(w-p[r])/(p[r+1]-p[r])
	else:
		return x[r]


def probability_visualization(P):
	plt.figure()
	plt.imshow(P,cmap='Greys_r')
	plt.show()


def period(spk,tstop=1500):
	T=[]
	bound1=50
	bound2=70
	for i in vect_index[5]:
		if len(spk[i])>20:
			tmp=spk[i][-1]-spk[i][-2]
			if tmp>bound1 and tmp<bound2:
				T.append(tmp)

	return np.median(T)


def indexes(A):
	idx=np.where(sum(A[np.ix_(vect_index[4],vect_index[4])])>=thresh1)
	conn=sum(A[np.ix_(750+idx[0],vect_index[4])])
	din_mean=np.mean(conn)
	din_std=np.std(conn)

	idx=np.where(sum(A[np.ix_(vect_index[4],vect_index[3])])>=thresh2)
	conn=sum(A[np.ix_(366+idx[0],vect_index[4])])
	cin_mean=np.mean(conn)
	cin_std=np.std(conn)

	y_hat=din_mean-cin_mean
	y_mdl=79.13-1.33*y_hat

	return (din_mean,din_std,cin_mean,cin_std)


def distance_dependent_prob(pre,post):
	P = np.load("data/P1000.npy")
	half_pos = 3000
	histo = range(-half_pos,half_pos+1,40)
	all_prob=[]
	for i in pre:
		y_histo = np.zeros(len(histo))
		count = np.zeros(len(histo))
		for j in post:
			if P[i,j]>0:
				h = pos_old[j]-pos_old[i]
				find_pos = [r for r in xrange(len(histo)-1) if histo[r]<h and histo[r+1]>=h]
				y_histo[find_pos] = y_histo[find_pos]+P[i,j]
				count[find_pos] = count[find_pos]+1
	
		tmp=[]
		for k in xrange(len(y_histo)):
			if count[k]!=0:	
				tmp.append(y_histo[k]/count[k])
			else:
				tmp.append(0)
	
		all_prob.append(tmp)

	prob_mean = np.mean(all_prob,axis=0)
	f = interp1d(np.multiply(histo,1),prob_mean,kind='linear',fill_value="extrapolate")
	return f


def first_firing_times(spk,tstop):
	n=len(spk)
	spk_left=spk[0:n/2]
	spk_right=spk[n/2+1:n]

	start_left=[]
	for train in spk_left:
		if len(train)>=3:
			start_left.append(train[0])

	start_right=[]
	for train in spk_right:
		if len(train)>=3:
			start_right.append(train[0])

	return (start_left,start_right)



def average_dIN_volt(v):
	ave_v=[]
	for i in xrange(v.shape[1]):
		ave_v.append(np.mean(v[:,i]))
	return ave_v


def classify_behaviour(spk,hdin_volt,m,t): # only works with fixed time step integration
	din_spks_l=[spk[j] for j in left_index[4] if pos[j]<1000.0]
	din_spks_r=[spk[j] for j in right_index[4] if pos[j]<1000.0]

	swim_l_start=np.median([x[0] for x in din_spks_l if len(x)>3])
	swim_r_start=np.median([x[0] for x in din_spks_r if len(x)>3])

	if np.all([~np.isnan(swim_l_start), ~np.isnan(swim_r_start)]):
		if swim_l_start<=swim_r_start:
			tstar=np.mean([x[0] for x in din_spks_l if len(x)>0])
		else:
			tstar=np.mean([x[0] for x in din_spks_r if len(x)>0])
	else: 
		tstar=0

	n=len(hdin_volt)
	hdin_v_l=hdin_volt[0:n/2]
	hdin_v_r=hdin_volt[n/2:-1]

	adin_l=average_dIN_volt(np.array(hdin_v_l))
	adin_r=average_dIN_volt(np.array(hdin_v_r))


	thresh=-27
	idxL=[k for k in xrange(len(adin_l)-1) if adin_l[k+1]>=thresh and adin_l[k]<thresh]
	if len(idxL)>0:
		tL=t[idxL[0]]
	else:
		tL=None
	idxR=[k for k in xrange(len(adin_r)-1) if adin_r[k+1]>=thresh and adin_r[k]<thresh]
	if len(idxR)>0:
		tR=t[idxR[0]]
	else:
		tR=None


	plt.figure()
	plt.plot(t,adin_l,'b',label='left',linewidth=2.0)
	plt.plot(t,adin_r,'r',label='right',linewidth=2.0)
	plt.legend(loc='upper left',fontsize=18)
	plt.ylabel('<v$_{hdIN}$>')
	if tstar is not 0:
		plt.fill_between(t, -60, -10, where=t>tstar,facecolor='green', alpha=0.25)
	plt.ylim([-80,-10])
	if len(idxL)>0:
		idx_l=idxL[0]
		plt.plot(t[idx_l], adin_l[idx_l],'bo',markersize=8) 
	idxR=[k for k in xrange(len(adin_r)-1) if adin_r[k+1]>=thresh and adin_r[k]<thresh]
	if len(idxR)>0:
		idx_r=idxR[0]
		plt.plot(t[idx_r], adin_r[idx_r],'ro',markersize=8) 

	if tL!=None and tR!=None:
		if abs(tL-tR)>3.0:
			if tL<tR:
				out=1 # swim left
			else:
				out=2 # swim right
		else:
			out=3 # sync
	else:
		if tL!=None or tR!=None:
			out = 4 # one sided
		else:
			out = 5 # no swim

	return (out,tL,tR,tstar)


def classification_name(idx):
	if idx==0:
		return "Undetected"
	if idx==1:
		return "Swim left"
	if idx==2:
		return "Swim right"
	if idx==3:
		return "Sync"
	if idx==4:
		return "One sided"
	if idx==5:
		return "No swim"




























