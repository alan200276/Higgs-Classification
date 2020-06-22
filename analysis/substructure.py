# Importing
import numpy as np
import math
import random
import os
import re
import matplotlib
import time
import pickle 

from itertools import groupby
from io import StringIO
# from mpl_toolkits.axes_grid1 import make_axes_locatable

import pyjet
from pyjet import cluster, DTYPE_PTEPM
from pyjet.testdata import get_event

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.pyplot import cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import PatchCollection

import scipy
import scipy.optimize as opt
from scipy.interpolate import griddata
from scipy import interpolate
from scipy.optimize import least_squares


# This file contains some substructure variables

def find_new_var_3(reclustered,clustered,pt_cut = 1):
	#The point is that many jets only have two subjets so the above is not really that efficient... hm
	new_var = []
	for sqn in range(len(reclustered)):
		net_R = 0
		jet = clustered[sqn][0]
		con = jet.constituents()
		sub_1 = reclustered[sqn][0]
		try:
			sub_2 = reclustered[sqn][1]
		except:
			#print('ah')
			sub_2 = sub_1
		for i in range(len(con)):
			if con[i].pt < pt_cut:
				continue
			#net_R = net_R + (con[i].pt/jet.pt)**(-1)*(R(sub_1,con[i])*R(sub_2,con[i]))
			net_R = net_R + (R(sub_1,con[i])*R(sub_2,con[i]))
		new_var.append(-net_R)
	return new_var

def find_new_var_3_norm(reclustered,clustered,pt_cut=1):
	new_var = []
	for sqn in range(len(reclustered)):
		net_R = 0
		jet = clustered[sqn][0]
		con = jet.constituents()
		sub_1 = reclustered[sqn][0]
		max_rad_1 = 0
		max_rad_2 = 0
		for i in range(0,len(con)):
			if R(con[i],sub_1) > max_rad_1:
				max_rad_1 = R(con[i],sub_1)
		try:
			sub_2 = reclustered[sqn][1]
		except:
			print('ah')
			sub_2 = sub_1
		for i in range(0,len(con)):
				if R(con[i],sub_2) > max_rad_2:
					max_rad_2 = R(con[i],sub_2)

		I = np.zeros((2,2))
		basis_1 = np.array([jet.px,jet.py,jet.pz])/np.linalg.norm(np.array([jet.px,jet.py,jet.pz]))
		basis_2 = np.array([-jet.py,jet.px,0])/np.linalg.norm(np.array([-jet.py,jet.px,0]))
		basis_3 = np.cross(basis_1,basis_2)
		transform = np.array([basis_1,basis_2,basis_3])
		
		norm_factor = 0
		
		for i in range(0,len(con)):
			if con[i].pt < pt_cut:
				continue
			new_par = transform.dot(np.array([con[i].px,con[i].py,con[i].pz]))
			rel_pt = (new_par[1]**2+new_par[2]**2)**(1/2)
			norm_factor = norm_factor + rel_pt
			net_R = net_R + (rel_pt)*(R(sub_1,con[i])*R(sub_2,con[i]))/(max_rad_1*max_rad_2)
		net_R = net_R/norm_factor
		new_var.append(net_R)
	return new_var
		
def find_new_var_N_2(reclustered,clustered):
	new_var = []
	for i in clustered:
		jcon = i[0].constituents()
		#Takes jcon, jet constituents
		p_total = np.sum([con.pt for con in jcon])
		v_1e2 = 0
		for i in range(len(jcon)):
			for j in range(i+1,len(jcon)):
				v_1e2 = v_1e2+ jcon[i].pt*jcon[j].pt*R(jcon[i],jcon[j])/(p_total**2)
		#v_1e2 = np.sum([[(con1.pt/p_total)*(con2.pt/p_total)*R(con1,con2) for con1 in jcon] for con2 in jcon])/2
		v_2e3 = 0
		for i in range(len(jcon)):
			for j in range(i+1,len(jcon)):
				for k in range(j+1,len(jcon)):
					v_2e3 = v_2e3 + jcon[i].pt*jcon[j].pt*jcon[j].pt*min(R(jcon[i],jcon[j]),R(jcon[j],jcon[k]),
																			R(jcon[i],jcon[k]))/(p_total**3)
		new_var.append(v_2e3/(v_1e2**2))
	return(new_var)
		
def find_new_var_pf(reclustered,clustered):
	#https://arxiv.org/pdf/0807.0234.pdf
	#Planar Flow
	new_var = []
	for sqn in range(len(reclustered)):
		
		I = np.zeros((2,2))
		jet = clustered[sqn][0]
		
		basis_1 = np.array([jet.px,jet.py,jet.pz])/np.linalg.norm(np.array([jet.px,jet.py,jet.pz]))
		basis_2 = np.array([-jet.py,jet.px,0])/np.linalg.norm(np.array([-jet.py,jet.px,0]))
		basis_3 = np.cross(basis_1,basis_2)
		
		transform = np.array([basis_1,basis_2,basis_3])
		printed = 0
		for par in clustered[sqn][0].constituents():
			#OK so we need to transfer into coordinates defined by the jet itself
			#new_par = np.dot(transform,np.array(par.px,par.py,par.pz))
			new_par = transform.dot(np.array([par.px,par.py,par.pz]))
			
			I = I + [[new_par[1]**2/par.e,new_par[1]*new_par[2]/par.e],[new_par[1]*new_par[2]/par.e,
																		new_par[2]**2/par.e]]
		I = I/clustered[sqn][0].mass
			#if printed == 0:
				#print(par.e,new_par)
				#print(new_par)
				#print([par.px,par.py,par.pz])
				#print(jet.px,jet.py,jet.pz)
			#	printed = 0
		#print(I)
		#print(transform)
		#new_var.append(4*(I[0,0]*I[1,1]-I[0,1]*I[1,0])/((I[0,0]+I[1,1])**2))
		new_var.append(4*np.linalg.det(I)/(np.trace(I)**2))
	return new_var

def find_new_var_radius(reclustered,clustered):
	#Finds the radius by looping over particles and finding the diameter and dividing by 2
	new_var = []
	for sqn in range(len(reclustered)):
		max_rad = 0
		con = clustered[sqn][0].constituents()
		for i in range(0,len(con)):
			for j in range(i,len(con)):
				if R(con[i],con[j]) > max_rad:
					max_rad = R(con[i],con[j])
		new_var.append(max_rad/2)
		#if new_var[-1] < 0.2:
		#	print(con)
	return new_var

def find_new_var_int(reclustered,clustered):
	#First, we want to start working in the transverse plane of the jet, with [jet, b-direction, other direction]
	#as our axis. So:
	new_var = []
	alpha = 1
	for sqn in range(len(reclustered)):
		good = 0
		jet = clustered[sqn][0]
		b_jet_1 = reclustered[sqn][0]
		
		basis_1 = np.array([jet.px,jet.py,jet.pz])/np.linalg.norm(np.array([jet.px,jet.py,jet.pz]))
		basis_2 = np.array([b_jet_1.px,b_jet_1.py,b_jet_1.pz])
		basis_2 = basis_2 - np.dot(basis_2,basis_1)*basis_1
		basis_2 = basis_2/np.linalg.norm(basis_2)
		basis_3 = np.cross(basis_1,basis_2)
		
		transform = np.array([basis_1,basis_2,basis_3])
		b_jet_1_newcoord = transform.dot(np.array([b_jet_1.px,b_jet_1.py,b_jet_1.pz]))
		b_jet_1_newcoord = b_jet_1_newcoord/b_jet_1_newcoord[0]
		#print(b_jet_1_newcoord)
		
		try:
			b_jet_2 = reclustered[sqn][1]
			b_jet_2_newcoord = transform.dot(np.array([b_jet_2.px,b_jet_2.py,b_jet_2.pz]))
			b_jet_2_newcoord = b_jet_2_newcoord/b_jet_2_newcoord[0]
		except:
			b_jet_2_newcoord = - b_jet_1_newcoord
		
		num = 0
		shift = (b_jet_1_newcoord+b_jet_2_newcoord)/2
		for par in clustered[sqn][0].constituents():
			new_par = transform.dot(np.array([par.px,par.py,par.pz]))
			new_par = new_par/new_par[0]
			new_par = new_par + shift
			#good = good + (1-new_par[1]**2/b_jet_1_newcoord[1]**2)
			#good = good + (-new_par[1]**3+b_jet_1_newcoord[1]**3)*(new_par[1]**3+b_jet_1_newcoord[1]**3)/b_jet_1_newcoord[1]**6
			#good = good + np.e**(-(new_par[1]/b_jet_1_newcoord[1])**2-(new_par[2]/b_jet_1_newcoord[1])**2)
			#good = good + 1/((new_par[1]/b_jet_1_newcoord[1])**2+1)
			#good = good + (abs(new_par[1]) < alpha*b_jet_1_newcoord[1])*(abs(new_par[2]) < alpha*b_jet_1_newcoord[1])
			#good = good + (new_par[1]**2+new_par[2]**2)/b_jet_1_newcoord[1]**2
			good = good + (1-new_par[1]**4/b_jet_1_newcoord[1]**4)*np.e**(-0.6*(new_par[1]**2+new_par[2]**2)/b_jet_1_newcoord[1]**2)
			num = num + 1
		#new_var.append(good)
		new_var.append(good/num)
	return new_var

def find_new_var_beta_3(reclustered,clustered,pt_cut = 1):
	new_var = []
	for sqn in range(len(reclustered)):
		tau_1 = np.zeros(3)
		tau_2 = np.zeros(2)
		d_0 = 0
		jet = clustered[sqn][0]
		sub_1 = reclustered[sqn][0]
		try:
			sub_2 = reclustered[sqn][1]
		except:
			sub_2 = sub_1
		for k in jet:
			if k.pt < pt_cut:
				continue
			d_0 = d_0 + k.pt*0.8
			tau_1 = tau_1 + k.pt*np.array([R(k,jet)**0.5,R(k,jet)**1,R(k,jet)**2])
			d_2 = min(R(sub_1,k),R(sub_2,k))
			tau_2 = tau_2 + k.pt*np.array([d_2**1,d_2**2])
		tau_1 = tau_1/d_0
		tau_2 = tau_2/d_0
		c = 0; d = 0.5; e = -1; a = 0; b = 0;
		#new_var.append(tau_1[0]**a*tau_1[1]**b*tau_1[2]**c*tau_2[0]**d*tau_2[1]**e)
		new_var.append(tau_1[0]**2*tau_2[0]**0.5*tau_2[1]**-1)
	return(new_var)

def find_new_var_beta_3_kT(reclustered,clustered,pt_cut = 1):
	new_var = []
	for sqn in range(len(reclustered)):
		tau_1 = np.zeros(3)
		tau_2 = np.zeros(2)
		d_0 = 0
		jet = clustered[sqn][0]
		
		#Now, we need to recluster jet into two kT subjets
		sequence_Cluster = cluster(jet, R=0.2,p=1)
		jets_Cluster = sequence_Cluster.inclusive_jets()
		
		sub_1 = jets_Cluster[0]
		try:
			sub_2 = jets_Cluster[1]
		except:
			sub_2 = sub_1
		for k in jet:
			if k.pt < pt_cut:
				continue
			d_0 = d_0 + k.pt*0.8
			tau_1 = tau_1 + k.pt*np.array([R(k,jet)**0.5,R(k,jet)**1,R(k,jet)**2])
			d_2 = min(R(sub_1,k),R(sub_2,k))
			tau_2 = tau_2 + k.pt*np.array([d_2**1,d_2**2])
		tau_1 = tau_1/d_0
		tau_2 = tau_2/d_0
		c = 0; d = 0.5; e = -1; a = 0; b = 0;
		#new_var.append(tau_1[0]**a*tau_1[1]**b*tau_1[2]**c*tau_2[0]**d*tau_2[1]**e)
		new_var.append(tau_1[0]**2*tau_2[0]**0.5*tau_2[1]**-1)
	return(new_var)
	
def find_new_var_beta_rb(reclustered,clustered,pt_cut=1):
	new_var = []
	for sqn in range(len(reclustered)):
		tau_1 = np.zeros(3)
		tau_2 = np.zeros(2)
		d_0 = 0
		jet = clustered[sqn][0]
		sub_1 = reclustered[sqn][0]
		try:
			sub_2 = reclustered[sqn][1]
		except:
			sub_2 = sub_1
		for k in jet:
			if k.pt < pt_cut:
				continue
			d_0 = d_0 + k.pt*0.8
			tau_1 = tau_1 + k.pt*np.array([R(k,jet)**0.5,R(k,jet)**1,R(k,jet)**2])
			d_2 = min(R(sub_1,k),R(sub_2,k))
			tau_2 = tau_2 + k.pt*np.array([d_2**1,d_2**2])
		tau_1 = tau_1/d_0
		tau_2 = tau_2/d_0
		c = 0; d = 0.5; e = -1; a = 0; b = 0;
		#new_var.append(tau_1[0]**a*tau_1[1]**b*tau_1[2]**c*tau_2[0]**d*tau_2[1]**e)
		#new_var.append(tau_1[0]**2*tau_2[0]**0.5*tau_2[1]**-1)
		
		net_R = 0
		jet = clustered[sqn][0]
		con = jet.constituents()
		sub_1 = reclustered[sqn][0]
		try:
			sub_2 = reclustered[sqn][1]
		except:
			#print('ah')
			sub_2 = sub_1
		for i in range(len(con)):
			if con[i].pt < pt_cut:
				continue
			#net_R = net_R + (con[i].pt/jet.pt)**(-1)*(R(sub_1,con[i])*R(sub_2,con[i]))
			net_R = net_R + (R(sub_1,con[i])*R(sub_2,con[i]))
		new_var.append(np.log(tau_1[0]**2*tau_2[0]**0.5*tau_2[1]**-1*20-1*net_R))
	return(new_var)
		
def find_new_var_kmeans(reclustered,clustered):
	#We use k_means with no pT dependence to cluster, and see how many clusteres
	new_var = []
	for sqn in range(len(clustered)):
		data = np.array([[x.eta,x.phi] for x in clustered[sqn][0].constituents()])
		mean = np.average(data,axis=0)
		badness_1 = np.sum([(x[0]-mean[0])**2+(x[1]-mean[1])**2 for x in data])
		try:
			mean = np.array([[reclustered[sqn][0].eta,reclustered[sqn][0].phi],
					[reclustered[sqn][1].eta,reclustered[sqn][1].phi]])
		except:
			mean = np.array([[reclustered[sqn][0].eta,reclustered[sqn][0].phi],
					[reclustered[sqn][0].eta+0.4,reclustered[sqn][0].phi]])
		badness = 0
		old_badness = -1
		for i in range(100):
			cluster = [[],[]]
			badness_c1 = 0; badness_c2 = 0;
			for k in data:
				if np.linalg.norm(k-mean[0],2) < np.linalg.norm(k-mean[1],2):
					cluster[0].append(k)
					badness_c1 = badness_c1 + np.linalg.norm(k - mean[0],2)
				else:
					cluster[1].append(k)
					badness_c2 = badness_c2 + np.linalg.norm(k - mean[1],2)
			if len(cluster[0]) != 0:
				badness_c1 = badness_c1/len(cluster[0])
			if len(cluster[1]) != 0:
				badness_c2 = badness_c2/len(cluster[1])
			badness = np.average([badness_c1,badness_c2])
			mean[0] = np.average(cluster[0])
			mean[1] = np.average(cluster[1])
			if old_badness == badness:
				break
			old_badness = badness
		new_var.append(badness_1/badness)
	return(new_var)

def find_IRC_kmeans(reclustered,clustered):
	class angle:
		def __init__(self,eta,phi):
			self.eta = eta
			self.phi = phi
		def __repr__(self):
			return "(" + str(self.eta) + "," + str(self.phi)+")"
	def average(sqn):
		#For a sequence of particles, finds the pT weighted average in eta,phi space
		pt_sum = np.sum([x.pt for x in sqn])
		avg_eta = np.sum([x.eta*(x.pt/pt_sum) for x in sqn])
		#print(str(avg_eta) + " avg_eta")
		#avg_phi = np.min(map(lambda phi : np.sum([dphi(x.phi, phi)**2*x.pt for x in sqn]),
		#					 np.arange(0,2*math.pi,0.01)))
		#Average phi is really annoying actually. Convention is -pi < phi < pi. To do this, 
		#let us split the points into phi > 0, phi < 0. Take the averages. If the two averages are
		#closer (say closer than pi) than take the average of the two with no shenanigans. If not, wrap.
		#print([x.phi < 0 for x in sqn])
		#print(np.array([x.phi < 0 for x in sqn]).shape)
		#print(sqn)
		try:
			under_0 = np.array(sqn)[np.array([x.phi < 0 for x in sqn])];
		except:
			under_0 = []
		#print([x.phi for x in under_0])
		try:
			over_0 = np.array(sqn)[np.array([x.phi >= 0 for x in sqn])];
		except:
			over_0 = []
		under_pt = np.sum([x.pt for x in under_0]); over_pt = np.sum([x.pt for x in over_0]);
		if under_pt == 0:
			return angle(avg_eta,np.sum([x.phi*x.pt for x in over_0])/over_pt)
		elif over_pt == 0:
			return angle(avg_eta,np.sum([x.phi*x.pt for x in under_0])/under_pt)
		under_avg = np.sum([x.phi*x.pt for x in under_0])/under_pt; 
		over_avg = np.sum([x.phi*x.pt for x in over_0])/over_pt;
		#print("avgs " + str(under_avg) + " " + str(over_avg))
		if abs(over_avg - under_avg) < 3.1415:
			avg_phi = (over_avg*over_pt+under_avg*under_pt)/(over_pt + under_pt)
		else:
			avg_phi = (over_avg*over_pt + (under_avg+math.pi*2)*under_pt)/(over_pt + under_pt)
			if avg_phi > math.pi:
				avg_phi = avg_phi - math.pi*2
		
		return angle(avg_eta,avg_phi)
	new_var = []
	for sqn in range(len(clustered)):
		#print("hmm")
		parts = clustered[sqn][0].constituents()
		mean = average(parts)
		badness_1 = np.sum([R(x,mean)**2*x.pt for x in parts])/np.sum([x.pt for x in parts])
		try:
			mean = [reclustered[sqn][0],reclustered[sqn][1]]
		except:
			mean = [reclustered[sqn][0],angle(reclustered[sqn][0].eta+0.1,reclustered[sqn][0].phi+0.1)]
		badness = 0
		old_badness = -1
		#print("hah")
		for i in range(100):
			#print(mean)
			#print("hoh")
			cluster = [[],[]]
			badness_c1 = 0; badness_c2 = 0;
			for x in parts:
				if R(x,mean[0]) < R(x,mean[1]):
					cluster[0].append(x)
				else:
					cluster[1].append(x)
			pT_1 = np.sum([x.pt for x in cluster[0]])
			pT_2 = np.sum([x.pt for x in cluster[1]])
			badness_c1 = np.sum([R(x,mean[0])**2*(x.pt/pT_1) for x in cluster[0]])
			badness_c2 = np.sum([R(x,mean[1])**2*(x.pt/pT_1) for x in cluster[1]])
			badness = np.average([badness_c1,badness_c2])
			#print("huh")
			mean[0] = average(cluster[0])
			mean[1] = average(cluster[1])
			#print("hih")
			if old_badness == badness:
				break
			old_badness = badness
		new_var.append(badness_1/badness)
	return(new_var)

def find_jet_pull(reclustered,clustered):
	
	"""
		Returns variables related to the jet pull vector
		Reclustered : Contains the 0.2 R anti-kT b-tagged subjets
		Clustered   : Contains the 0.8 R anti-kT 'higgs' jet (which includes the 0.2R subjets)
	"""
	
	def y(p):
		""" Returns the rapidity"""
		return ((1/2)*math.log((p.e+p.pz)/(p.e-p.pz)))
	
	new_var = []
	for sqn in range(len(reclustered)):
		con_1 = reclustered[sqn][0].constituents()
		J1 = reclustered[sqn][0]
		try:
			con_2 = reclustered[sqn][1].constituents()
			J2 = reclustered[sqn][1]
		except:
			new_var.append([0,0,0])
			continue
		
		r_phi_1  = np.array([dphi(x.phi,J1.phi) for x in con_1])  
		r_y_1	= np.array([y(x)-y(J1) for x in con_1])
		r_norm_1 = (r_phi_1**2+r_y_1**2)**(1/2)
		pt_1	 = np.array([x.pt for x in con_1])
		
		"""This is the pull vector of the second b-jet acting on the first b-jet (sorted by pT)"""
		pull_1  = np.array([np.sum([pt_1*r_y_1*r_norm_1]),
							np.sum([pt_1*r_phi_1*r_norm_1])])/J1.pt
		pull_1_norm = (pull_1[0]**2+pull_1[1]**2)**(1/2)
		
		r_phi_2  = np.array([dphi(x.phi,J2.phi) for x in con_2])  
		r_y_2	= np.array([y(x)-y(J2) for x in con_2])
		r_norm_2 = (r_phi_2**2+r_y_2**2)**(1/2)
		pt_2	 = np.array([x.pt for x in con_2])
		
		"""Similarly the pull vector of the first on the second"""
		pull_2  = np.array([np.sum([pt_2*r_y_2*r_norm_2]),
							np.sum([pt_2*r_phi_2*r_norm_2])])/J2.pt
		pull_2_norm = (pull_2[0]**2+pull_2[1]**2)**(1/2)
		
		"""The vector from one b-subjet to another"""
		J1_J2 = [y(J2)-y(J1),dphi(J2.phi,J1.phi)]
		J2_J1 = [y(J1)-y(J2),dphi(J1.phi,J2.phi)]
		J_norm = (J1_J2[0]**2+J1_J2[1]**2)**(1/2)
		
		"""The angles the pull vectors make with the vector from one subjet to another"""
		A_1 = (pull_1[0]*J1_J2[0]+pull_1[1]*J1_J2[1])/(pull_1_norm*J_norm)
		A_2 = (pull_2[0]*J2_J1[0]+pull_2[1]*J2_J1[1])/(pull_2_norm*J_norm)
		
		"""The angle between the pull vectors (just a test)"""
		pA =  (pull_1[0]*pull_2[0]+pull_1[1]*pull_2[1])/(pull_1_norm*pull_2_norm)
		new_var.append([math.acos(A_1),math.acos(A_2),math.acos(pA)])
		
	return new_var

def find_connexion(reclustered,clustered):
	#The idea is that given a metric between particles (e.g. the anti_kT metric),
	#we want to find the minimum epsilon such that if each particle is connected
	#to all other particles closer than epsilon away, then the resultant graph
	#is fully connected
	new_var = []
	class node():
		def add_connection(self,node):
			self.connected.append(node)
		def __init__(self,pt,eta,phi,e,pz,index):
			self.connected = []
			self.pt = pt
			self.eta = eta
			self.phi = phi
			self.e = e
			self.pz = pz
			self.to_look = 0
			self.index = index
		def distance(self,node):
			#return R_y(self,node)**2
			return max(self.pt**(-2),node.pt**(-2))*R_y(self,node)**2 #The R^2 factor don't matter
	for x in range(len(clustered)):
		parts = clustered[x][0].constituents()
		node_list = []
		for y in range(len(parts)):
			node_list.append(node(parts[y].pt,parts[y].eta,parts[y].phi,parts[y].e,parts[y].pz,y))
		alive_nodes = [node_list[0]]
		seen_nodes = [node_list[0]]
		epsilon = 0
		ghost_edges = []
		while True:
			#time.sleep(1)
			#print(ghost_edges)
			while len(alive_nodes) > 0:
				while alive_nodes[0].to_look < len(node_list):
					i_temp = alive_nodes[0].index
					i_look = alive_nodes[0].to_look
					if i_look == i_temp:
						alive_nodes[0].to_look = i_look + 1
						continue
					d_temp = alive_nodes[0].distance(node_list[i_look])
					if d_temp < epsilon:
						if node_list[i_temp] not in seen_nodes:
							alive_nodes.append(node_list[i_temp])
							seen_nodes.append(node_list[i_temp])
					else:
						where_in = 0
						for k in range(len(ghost_edges)):
							if d_temp > ghost_edges[k][2]:
								where_in = where_in + 1
							else:
								break
						ghost_edges.insert(where_in,[i_temp,i_look,d_temp])
					alive_nodes[0].to_look = i_look + 1
				alive_nodes.pop(0)
			if len(seen_nodes) == len(node_list):
				break
			while True:
				new_edge = ghost_edges.pop(0)
				epsilon = new_edge[2]
				if node_list[new_edge[1]] not in seen_nodes:
					alive_nodes.append(node_list[i_temp])
					seen_nodes.append(node_list[i_temp])
					break
			if len(seen_nodes) == len(node_list):
				break
		new_var.append(epsilon)
	return(new_var) 

	# Correct phi range
def fix_phi(phi):
	while phi > math.pi:
		phi = phi - 2*math.pi
	while phi < -math.pi:
		phi = phi + 2*math.pi
	return phi

# Returns the difference in phi between phi, and phi_center
# as a float between (-PI, PI)
def dphi(phi,phi_c):
	
	dphi_temp = phi - phi_c
	while dphi_temp > np.pi:
		dphi_temp = dphi_temp - 2*np.pi
	while dphi_temp < -np.pi:
		dphi_temp = dphi_temp + 2*np.pi
	return (dphi_temp)

# Rapidity
def y(p):
	return ((1/2)*math.log((p.e+p.pz)/(p.e-p.pz)))

def R(con1,con2):
	return (((con1.eta-con2.eta)**2+dphi(con1.phi,con2.phi)**2)**(1/2))

def R_y(con1,con2):
	return (((y(con1)-y(con2))**2+dphi(con1.phi,con2.phi)**2)**(1/2))

def N_2(jcon):
	#Takes jcon, jet constituents
	p_x_total = np.sum([con.px for con in jcon])
	p_y_total = np.sum([con.py for con in jcon])
	p_total = (p_x_total**2+p_y_total**2)**(1/2)
	
	v_1e2 = 0
	for i in range(len(jcon)):
		for j in range(i+1,len(jcon)):
			v_1e2 = v_1e2+ jcon[i].pt*jcon[j].pt*R(jcon[i],jcon[j])/(p_total**2)
	v_2e3 = 0
	for i in range(len(jcon)):
		for j in range(i+1,len(jcon)):
			for k in range(j+1,len(jcon)):
				v_2e3 = v_2e3 + jcon[i].pt*jcon[j].pt*jcon[j].pt*min(R(jcon[i],jcon[j])*R(jcon[i],jcon[k]),
																	 R(jcon[j],jcon[k])*R(jcon[i],jcon[j]),
														R(jcon[i],jcon[k])*R(jcon[j],jcon[k]))/(p_total**3)
	return v_2e3/(v_1e2**2)