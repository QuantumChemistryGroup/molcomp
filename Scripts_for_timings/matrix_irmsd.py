#!/usr/bin/python
import sys, string, os, glob, re
import subprocess
from math import *
from multiprocessing import Pool
import multiprocessing
import random
from irmsd import read_structures, get_irmsd_molecule
#import irmsd
#print(dir(irmsd)) 

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
		pass
		#print ('No energy in ', j, '. Stop.')
		#exit(0)
	return E


def parallel_chunk_irmsd(j):
	mol_a_list = read_structures(j[0])
	mol_b_list = read_structures(j[1])

	# 2. Extract the molecule objects from the lists
	mol_a = mol_a_list[0]
	mol_b = mol_b_list[0]

	# 3. Calculate the permutation-invariant RMSD
	irmsd_value = get_irmsd_molecule(mol_a, mol_b)
	out1 = irmsd_value[0]
	m = [j[0],j[1]]
	# compare with criterion
	# rmsd-based
	# get energy
	E1 = GETE(j[0])
	E2 = GETE(j[1])
	if ((E1 != 'NaN') and (E2 != 'NaN')):
		dE = round(abs(E2-E1)*627.51,2)
	else:
		dE = 'NaN'
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
	topo = ''
	if sys.argv[-1].lower=='readtopo':
		topo=' --readtopo=true'
	top = glob.glob('topo.txt')
	if (topo == ' --readtopo=true') and len(top) == 0:
		print ("Your topo.txt file is missing. Stop.\n")
		exit(0)
	if mode=='irmsd-def':
		option1=' compare '
		value=0.125
	elif mode=='irmsd-500':
		option1=' compare '
		value=0.500
	elif mode=='irmsd-cre':
		option1=' compare --quaternion '
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
		if ('irmsd' in mode):
			result = pool.map(parallel_chunk_irmsd,chunks)

	f=open(mode+'.txt','w')
	for i in result:
		for j in range(len(i)):
			f.write(str(i[j])+'   ')
		f.write('\n')
	f.close()
