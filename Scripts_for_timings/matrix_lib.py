#!/usr/bin/python
import sys, string, os, glob, re
import subprocess
from math import *
from multiprocessing import Pool
import multiprocessing
import random
from ctypes import CDLL, c_char_p, c_int, POINTER, c_double, cast, create_string_buffer

mode = ''
option1 = ''
option2 = ''
value = ''
maxdist = ''
rmsd = ''
file_path = os.path.abspath(os.path.dirname(__file__))

lib = CDLL(file_path + "/libmolcomp.so")
lib.main.argtypes = [c_int, POINTER(c_char_p)]
lib.main.restype = c_int


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
		exit(0)
	return E

def parallel_chunk(task_args):
	j, config = task_args
	
	mode = config['mode']
	value = config['value']
	rmsd = config['rmsd']
	maxdist = config['maxdist']
	option1_bytes = config['option1_bytes']
	encoding = config['encoding']

	res_array = (c_double * 3)()

	c_args_buffers = [
		create_string_buffer((file_path + "/molcomp").encode(encoding)),
		create_string_buffer(j[0].encode(encoding)),
		create_string_buffer(j[1].encode(encoding))
	]
	for opt in option1_bytes:
		c_args_buffers.append(create_string_buffer(opt))
		
	argc = len(c_args_buffers)
	
	argv_type = c_char_p * (argc + 1)
	argv = argv_type()
	for idx, buf in enumerate(c_args_buffers):
		argv[idx] = cast(buf, c_char_p)
		
	argv[argc] = cast(res_array, c_char_p)

	status = lib.main(argc, argv)

	if status != 0:
		return [j[0], j[1], 1, 0.0, 0.0, 0.0, 'NaN']

	out = [str(res_array[0]), str(res_array[1]), str(res_array[2])]
	#print (out)
	m = [j[0],j[1]]
	# compare with criterion
	# rmsd-based
	# get energy
	E1 = GETE(j[0])
	E2 = GETE(j[1])
	dE = round(abs(E2-E1)*627.51,2)
	if (mode=='rmsd-graph' or mode == 'rmsd-hungarian' or mode == 'rmsd-kabsch'):
		out1=''
		out1 = float(out[0])
		if out1 < value:
			m.append(0)
		else:
			m.append(1)
	elif (mode == 'angle'):
		out1 = float(out[0])
		out2 = float(out[1])
		out3 = float(out[2])
		if (out1 < rmsd and out2 < maxdist and out3 < value):
			m.append(0)
		else:
			m.append(1)
	else:
		print ("Your mode (dist?) is not implemented yet. Stop.")
		exit(0)
	# append rmsd, avdev, lardev	
	for i in out:
		m.append(round(float(i),3))
	# append dE
	m.append(dE)
	return m


if __name__ == '__main__':
	multiprocessing.set_start_method('fork', force=True)
	encoding='utf-8'
	# definitions
	mode = sys.argv[1].lower()
	topo = ''
	
	if sys.argv[-1].lower() == 'readtopo':
		topo=' --readtopo=true'
	top = glob.glob('topo.txt')
	if (topo == ' --readtopo=true') and len(top) == 0:
		print ("Your topo.txt file is missing. Stop.\n")
		exit(0)
	if mode=='rmsd-kabsch':
		option1=' --connect=false --alg=kabsch --keepenantio=false '
		value=0.125
	elif mode=='rmsd-hungarian':
		option1=' --connect=false --alg=std --keepenantio=false '
		value=0.125
	elif mode=='rmsd-graph':
		option1=' --connect=true --alg=std --keepenantio=false --compare=angle'+topo
		value=0.125
	elif mode=='dist':
		option1=' --connect=true --alg=std --compare=dist --keepenantio=false'+topo
		value=0.125
		distsm=0.75
		rat=4.0
		distlr = 1.0
		exit(0)
	elif mode=='angle':
		option1=' --connect=true --alg=std --keepenantio=false --compare=angle'+topo
		value=30.0
		rmsd = 0.5
		maxdist = 1.5
	else: 
		print ("Your mode is not clear. Stop.\n")
		exit(0)
	filelist = glob.glob('*.xyz')
	filelist.sort()
	option1_bytes = [opt.encode(encoding) for opt in option1.split() if opt]

	config_payload = {
		'mode': mode,
		'value': value,
		'rmsd': rmsd,
		'maxdist': maxdist,
		'option1_bytes': option1_bytes,
		'encoding': encoding
	}

	outputs=[]
	chunks = []
	for i in range(len(filelist)):
		for j in range(i+1,len(filelist),1):
			chunks.append(([filelist[i],filelist[j]], config_payload))
			
	cpu_count = multiprocessing.cpu_count()
	with Pool(cpu_count) as pool:
		result = pool.map(parallel_chunk, chunks)
	f=open(mode+'.txt','w')
	for i in result:
		for j in range(len(i)):
			f.write(str(i[j])+'   ')
		f.write('\n')
	f.close()
