#!/usr/bin/python
import sys, string, os, glob, re
import subprocess
from math import *
from multiprocessing import Pool
import multiprocessing
import random
from func_RMSD_RDKit import align_xyz_files

mode = ''
option1 = ''
option2 = ''
value = ''
maxdist = ''
rmsd = ''
file_path = os.path.abspath(os.path.dirname(__file__))

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


def parallel_chunk_rdkit(j):
	out, aligned_mol = align_xyz_files(j[0], j[1], print_output=False)
	m = [j[0],j[1]]
	# compare with criterion
	# rmsd-based
	# get energy
	E1 = GETE(j[0])
	E2 = GETE(j[1])
	dE = round(abs(E2-E1)*627.51,2)
	out1 = float(out)
	if out1 < value:
		m.append(0)
	else:
		m.append(1)
	# append rmsd, avdev, lardev	
	m.append(round(float(out1),3))
	# append dE
	m.append(dE)
	return m

if __name__ == '__main__':
	multiprocessing.set_start_method('fork')
	encoding='utf-8'
	# definitions
	mode = sys.argv[1].lower()
	if mode=='rdkit-125':
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
			chunks.append([filelist[i],filelist[j]])
	cpu_count = multiprocessing.cpu_count()
	with Pool(cpu_count) as pool:
		if ('irmsd' not in mode) and ('rdkit' not in mode):
			result = pool.map(parallel_chunk,chunks)
		elif ('irmsd' in mode):
			result = pool.map(parallel_chunk_irmsd,chunks)
		elif ('rdkit' in mode):
			result = pool.map(parallel_chunk_rdkit,chunks)

	f=open(mode+'.txt','w')
	for i in result:
		for j in range(len(i)):
			f.write(str(i[j])+'   ')
		f.write('\n')
	f.close()
