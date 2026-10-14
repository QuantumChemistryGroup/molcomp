#!/usr/bin/python
import sys, string, os, glob, re
import subprocess
from math import *
from multiprocessing import Pool
import multiprocessing
import random
from func_otmol_sort import compute_otmol_rmsd

file_path = os.path.abspath(os.path.dirname(__file__))
encoding = 'utf-8'  

def GETE (j):
	E = 'NaN'
	f = open(j, 'r')
	s = f.readline()
	s = f.readline()
	p = re.compile(r'[+-]?\d+\.\d+')
	d = [float(i) for i in p.findall(s)]  # Convert strings to float
	if len(d) !=0:
		E = float(d[0])
	f.close()
	if E == 'NaN':
		print ('No energy in ', j, '. Stop.')
		#exit(0)
	return E


def parallel_chunk_otmol(args):
	j, option1, value = args
	out = compute_otmol_rmsd(j[0], j[1])
	m = [j[0],j[1]]
	# compare with criterion
	# rmsd-based
	# get energy
	E1 = GETE(j[0])
	E2 = GETE(j[1])
	dE = round(abs(E2-E1)*627.51,2)
	out1 = float(out)
	#print (j)
	#print (out1)
	#print (float(value))
	if out1 < value:
		m.append(0)
	else:
		m.append(1)
	# append rmsd, avdev, lardev	
	m.append(round(float(out1),3))
	# append dE
	m.append(dE)
	#print (m)
	return m


if __name__ == '__main__':
	#encoding='utf-8'
	# definitions
	option1 = ''
	option2 = ''
	value = ''
	maxdist = ''
	rmsd = ''
	mode = sys.argv[1].lower()
	topo = ''
	if sys.argv[-1].lower=='readtopo':
		topo=' --readtopo=true'
	top = glob.glob('topo.txt')
	if (topo == ' --readtopo=true') and len(top) == 0:
		print ("Your topo.txt file is missing. Stop.\n")
		exit(0)
	if mode=='otmol':
		option1=' '
		value=0.125
	else: 
		print ("Your mode is not clear. Stop.\n")
		exit(0)
	filelist = glob.glob('*.xyz')
	#filelist = glob.glob('*B97-3c.out')
	filelist.sort()
	outputs=[]
	chunks = []
	for i in range(len(filelist)):
		for j in range(i+1,len(filelist),1):
			pair = [filelist[i],filelist[j]]
			chunks.append([pair, option1, value])
	cpu_count = multiprocessing.cpu_count()
#	for i in chunks:
#		print (i)
#		parallel_chunk_otmol(i)
	with Pool(cpu_count) as pool:
		if ('otmol' in mode):
			result = pool.map(parallel_chunk_otmol,chunks)

	f=open(mode+'.txt','w')
	for i in result:
		for j in range(len(i)):
			f.write(str(i[j])+'   ')
		f.write('\n')
	f.close()
