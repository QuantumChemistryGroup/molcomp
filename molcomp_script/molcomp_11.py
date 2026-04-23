#!/usr/bin/python
import sys, string, os, glob, re
import subprocess
from math import *
from multiprocessing import Pool
import multiprocessing
import random

mode = ''
option1 = ''
option2 = ''
value = ''
topo= ''
maxdist = ''
rmsd = ''

def POPOPT (mode):
	global option1
	global option2
	global value
	global topo
#	global distsm
#	global distlr
#	global rat
	global rmsd
	global maxdist
	topo = ''
	if sys.argv[-1].lower()=='readtopo':
		topo=' --readtopo=true'
	top = glob.glob('topo.txt')
	if (topo == ' --readtopo=true') and len(top) == 0:
		print ("Your topo.txt file is missing. Stop.\n")
		exit(0)
	if mode=='rmsd-kabsch':
		option2=' --connect=false --alg=kabsch --keepenantio=false --compare=dist'
		value=0.125
	elif mode=='rmsd-hungarian':
		option2=' --connect=false --alg=std --keepenantio=false --compare=dist'
		value=0.125
	elif mode=='rmsd-graph':
		option2=' --connect=true --alg=std --keepenantio=false --compare=angle'+topo
		value=0.125
	elif mode=='angle':
		option2=' --connect=true --alg=std --keepenantio=false --compare=angle'+topo
		value=30.0
		rmsd = 0.5
		maxdist = 1.5
	elif mode=='dist':
		option2=' --connect=true --alg=std --dist=true --keepenantio=false'+topo
		value=0.125
		distsm=0.75
		rat=4.0
		distlr=1.0
		print ("Dist is not implemented yet. Exit.")
		exit(0)


############################
# CHUNK-based options
############################

file_path = os.path.abspath(os.path.dirname(__file__))
# conformers in chunk
n=24
#iterations
iter1=3000
# permutations
perm=8
# also calculate BO Nuclear Repulsion energy
NR = True
# or calculate RMSD with respect to the first structure
CROSSRMSD = True

MIX = False # xyz are from different methods & software

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


def IsFloat(s):
	try: 
		float(s)
		return True
	except ValueError:
		return False

# function to populate options and value

elements = {
       'H':'1',
       'He':'2',
       'Li':'3',
       'Be':'4',
       'B':'5',
       'C':'6',
       'N':'7',
       'O':'8',
       'F':'9',
       'Ne':'10',
       'Na':'11',
       'Mg':'12',
       'Al':'13',
       'Si':'14',
       'P':'15',
       'S':'16',
       'Cl':'17',
       'Ar':'18',
       'K':'19',
       'Ca':'20',
       'Sc':'21',
       'Ti':'22',
       'V':'23',
       'Cr':'24',
       'Mn':'25',
       'Fe':'26',
       'Co':'27',
       'Ni':'28',
       'Cu':'29',
       'Zn':'30',
       'Ga':'31',
       'Ge':'32',
       'As':'33',
       'Se':'34',
       'Br':'35',
       'Kr':'36',
       'Rb':'37',
       'Sr':'38',
       'Y':'39',
       'Zr':'40',
       'Nb':'41',
       'Mo':'42',
       'Tc':'43',
       'Ru':'44',
       'Rh':'45',
       'Pd':'46',
       'Ag':'47',
       'Cd':'48',
       'In':'49',
       'Sn':'50',
       'Sb':'51',
       'Te':'52',
       'I':'53',
       'Xe':'54',
       'Cs':'55',
       'Ba':'56',
       'La':'57',
       'Ce':'58',
       'Pr':'59',
       'Nd':'60',
       'Pm':'61',
       'Sm':'62',
       'Eu':'63',
       'Gd':'64',
       'Tb':'65',
       'Dy':'66',
       'Ho':'67',
       'Er':'68',
       'Tm':'69',
       'Yb':'70',
       'Lu':'71',
       'Hf':'72',
       'Ta':'73',
       'W':'74',
       'Re':'75',
       'Os':'76',
       'Ir':'77',
       'Pt':'78',
       'Au':'79',
       'Hg':'80',
       'Tl':'81',
       'Pb':'82',
       'Bi':'83',
       'Po':'84',
       'At':'85',
       'Rn':'86',
       'Fr':'87',
       'Ra':'88',
       'Ac':'89',
       'Th':'90',
       'Pa':'91',
       'U':'92',
       'Np':'93',
       'Pu':'94',
       'Am':'95',
       'Cm':'96',
       'Bk':'97',
       'Cf':'98',
       'Es':'99',
       'Fm':'100',
       'Md':'101',
       'No':'102',
       'Lr':'103',
       'Rf':'104',
       'Db':'105',
       'Sg':'106',
       'Bh':'107',
       'Hs':'108',
       'Mt':'109',
       'Ds':'110',
       'Rg':'111',
       'Cn':'112',
       'Uut':'113',
       'Fl':'114',
       'Uup':'115',
       'Lv':'116',
       'Uuh':'117',
       'Uuh':'118',
}

def NucRep (Conf_E):
	coord = []
	f = open(Conf_E[0], 'r')
	for line in f:
		d = str.split(line)
		if (len(d) == 4 and IsFloat(d[1]) == True and IsFloat(d[2]) == True and IsFloat(d[3]) == True):
			coord.append([d[0],d[1],d[2],d[3]])
	f.close()
	E = 0
	for i in range(len(coord)-1):
		for j in range(i+1, len(coord), 1):
			dx = float(coord[i][1])-float(coord[j][1])
			dy = float(coord[i][2])-float(coord[j][2])
			dz = float(coord[i][3])-float(coord[j][3])
			d = sqrt(dx*dx + dy*dy + dz*dz)*1.8897261259077824
			q1 = float(elements[coord[i][0]])
			q2 = float(elements[coord[j][0]])
			E = E + (q1*q2)/d
			#E = E + (332.0538*(q1*q2)/d)/627.5095
	#print (Conf_E[0],E)
	Conf_E.append(E)
	return Conf_E

def RMSD(j):
	Conf_E = j
	Conf_E0 = refstr
	#print (Conf_E, Conf_E0)
	option3 = ' --connect=false --alg=std --keepenantio=false'
	E = ''
	out = subprocess.check_output(file_path+"/molcomp "+Conf_E[0]+' '+Conf_E0+option3, shell=True)
	out = out.decode(encoding)
	out = re.findall(r"\d*\.?\d+", out)
	E = float(out[0])
	Conf_E.append(E)						
	return Conf_E

def DUP (dE, out):
	M = ''
	if (mode=='rmsd-graph' or mode == 'rmsd-hungarian' or mode == 'rmsd-kabsch'):
		out1=''
		out1 = float(out[0])
		if out1 < value:
			M = 0
		else:
			M = 1
	elif (mode == 'angle'):
		out1 = float(out[0])
		out2 = float(out[1])
		out3 = float(out[2])
		if (out1 < rmsd and out2 < maxdist and out3 < value):
			M = 0
		else:
			M = 1
	else:
		print ("Your mode (dist?) is not implemented yet. Stop.")
		exit(0)
	return M
	
def parallel_chunk(j):
	duppair = []
	dup = set()
	un = 0
	pairs = []
	for i in range(len(j)-1):
		if (j[i][0] not in dup):
			for m in range(i+1,len(j),1):
				if j[m][0] not in dup:
					dE = abs(j[m][1] - j[i][1])
					#if dE < Ethld:
					tmp_pair = ''
					if j[m][1] > j[i][1]:
						tmp_pair = [j[i][0],j[m][0]]
					elif j[m][1] < j[i][1]:
						tmp_pair = [j[m][0],j[i][0]]
					else:
						tmp_pair=[j[i][0],j[m][0]]
						tmp_pair.sort(key=natural_keys)
						#print ("Energies are identical: ",tmp_pair)
					#tmp_pair = [j[i][0],j[m][0]]
					#tmp_pair.sort()
					tmp_pair = tuple(tmp_pair)
					if ( tmp_pair not in all_pairs ):
						un=un+1
						out = subprocess.check_output(file_path+"/molcomp "+tmp_pair[0]+' '+tmp_pair[1]+option2, shell=True)
						out = out.decode(encoding)
						out = re.findall(r"\d*\.?\d+", out)
						# compare with criterion
						out1 = ''
						out1 = DUP (dE*627.51, out)
						if out1 == 0:
							pair = []
							if tmp_pair[1] == j[m][0]:
								dup.add(j[m][0])
								for l in range(len(j[i])):
									pair.append(j[i][l])
								for l in range(len(j[m])):
									pair.append(j[m][l])
								for l in range(len(out)):
									pair.append(round(float(out[l]),3))
								pair.append(round(dE*627.51,2))
							elif tmp_pair[1] == j[i][0]:
								dup.add(j[i][0])
								for l in range(len(j[m])):
									pair.append(j[m][l])
								for l in range(len(j[i])):
									pair.append(j[i][l])
								for l in range(len(out)):
									pair.append(round(float(out[l]),3))
								pair.append(round(dE*627.51,2))
								break
							else:
								exit(0)
							duppair.append(pair)
						else:
							pairs.append(tmp_pair)
	#print ("un: ", un)
	return [duppair,pairs]



	

def POPE (filelist):
	conf_E = []
	for i in filelist:
		f = open(i, 'r')
		E = 'NaN'
		s = f.readline()
		s = f.readline()
		d = str.split(s)
		if len(d) == 1 and IsFloat(d[0]) == True:
			E = float(d[0])
		else:
			E = 'NaN'
		if E !='NaN':
			conf_E.append([i,E])
		else:
			print ("Warning: E= ",E, " in ", i)
		f.close()
	#conf_E.sort(key=lambda x: x[1])
	return conf_E
all_pairs = set()
refstr = ''
if __name__ == '__main__':
	encoding='utf-8'
	# definitions
	SingleFile = False
	if '.xyz' in sys.argv[1].lower():
		SingleFile=True
		mode = sys.argv[2].lower()
	else:
		mode = sys.argv[1].lower()
	if (SingleFile == False):
		# check if there are file/folders
		if os.path.isdir('unique-'+mode):
			print ('Folder ', 'unique-'+mode, 'exists. Remove it first.')
			exit(0)
		#if os.path.isfile(mode+'.txt'):
		#	print ('File ', mode+'.txt', 'exists. Remove it first.')
		#	exit(0)
		if os.path.isfile('dup-'+mode+'.txt'):
			print ('File ', 'dup-'+mode+'.txt', 'exists. Remove it first.')
			exit(0)
	else:
		if os.path.isdir('tmp'):
			print ('Folder ', 'tmp', 'exists. Remove it first.')
			exit(0)
		# check if there are file/folders
		if os.path.isfile(sys.argv[1][0:-4]+'-unique-'+mode+'.xyz'):
			print ('File ', sys.argv[1][0:-4]+'-unique-'+mode+'.xyz', 'exists. Remove it first.')
			exit(0)
		if os.path.isfile(sys.argv[1][0:-4]+'-'+mode+'.txt'):
			print ('File ', sys.argv[1][0:-4]+'-'+mode+'.txt', 'exists. Remove it first.')
			exit(0)
		if os.path.isfile(sys.argv[1][0:-4]+'-dup-'+mode+'.txt'):
			print ('File ', sys.argv[1][0:-4]+'-dup-'+mode+'.txt', 'exists. Remove it first.')
			exit(0)
	path = os.getcwd()
	POPOPT(mode)
	if SingleFile == True:
		os.mkdir('tmp')
		os.system('cp '+sys.argv[1]+' tmp')
		if topo==' --readtopo=true':
			os.system('cp '+'topo.txt'+' tmp')
		os.chdir(path+'/tmp')
		os.system('python3 '+file_path+'/split_xyz_E.py '+sys.argv[1])
		os.system('rm '+sys.argv[1])	
	filelist = glob.glob('*.xyz')
	conf_Etmp = POPE(filelist[0:1])
	#print (conf_Etmp)
	if len(conf_Etmp) == 0:
		MIX=True
	# return sorted array
	if MIX == False:
		conf_E = POPE(filelist)
		refstr = conf_E[0][0]
	else:
		if (NR == False):
			print ('If MIX == True, then NR should be True. Stop.')
			exit(0)
		conf_E = []
		for i in filelist:
			conf_E.append([i])
		refstr = conf_E[0][0]
	cpu_count = multiprocessing.cpu_count()
	if NR == True:
		# calculate Nuclear repulsion
		with Pool(cpu_count) as pool:
			conf_E = pool.map(NucRep,conf_E)
	if CROSSRMSD == True:
		# calculate Nuclear repulsion
		with Pool(cpu_count) as pool:
			conf_E = pool.map(RMSD,conf_E)
		if MIX == True:
			refstr = conf_E[-1][0]
			with Pool(cpu_count) as pool:
				conf_E = pool.map(RMSD,conf_E)
	conf_E.sort(key=lambda x: x[1])
	# if t==2, then only E
	# if t==3, then E + (NR || RMSD)
	# if t==4, then E + (NR && RMSD)
	t = len(conf_E[0])
	dup_all = set()
	pm=0
	N = n	
	print ('Parallel starts\n')
	for it in range (iter1):
		dup_all_it = []
		dup_all_it = set()
		#print (conf_E)
		print ("Iteration: ", it, "permuation: ", pm, " (",n,")")
		chunks = []
		a = len(dup_all)
		for i in range(0,len(conf_E), n):
			chunks.append(conf_E[i:i + n])
			#print (conf_E[i:i + n])
		with Pool(cpu_count) as pool:
			res = pool.map(parallel_chunk,chunks)
			# res[0] - duplicates
#			for i in range(len(res)):
#				for j in range(len(res[i])):
#					print (res[i][j])
			for mm in res:
				for m1 in range(len(mm)):
					if len(mm[m1]) != 0:
						if m1==0:
							for m2 in mm[m1]:
								m2 = tuple(m2)
								dup_all.add(m2)
								dup_all_it.add(m2)
						else:
							for m2 in mm[m1]:
								all_pairs.add(m2)
#			for mm in res[1]:
#				if len(mm)!=0:
#					for m1 in mm:
#						all_pairs.append(m1)
			print (len(all_pairs))
		print ("LEN: ", len(dup_all))
		if len(dup_all) == a and pm >= perm:
			break 
		elif len(dup_all) == a and pm < perm:
			#n=int(n+n*(0.1)*pm)
			#random.shuffle(conf_E)
			# flip a coin
			coin=1
			if t !=2:
				coin = random.randint(1, t-1)
			else:
				coin = 1
			step = int((128-N)/perm)+coin
			n=n+step
			if n % 2 !=0:
				conf_E.sort(key=lambda x: x[coin], reverse=False)
			else: 
				conf_E.sort(key=lambda x: x[coin], reverse=True)
			pm=pm+1
		else:
			for j in dup_all_it:
				tmp = j[t]
				#print (tmp)
				for m in conf_E:
					if tmp == m[0]: 
						conf_E.remove(m)
			pm=0
			n=N
	if (pm !=perm):
		print ("Warning: increse number of iter="+str(iter1)+'\n')
		
	# sort again based on energies
	conf_E.sort(key=lambda x: x[1])
	print ('UINQUE = ',len(conf_E))
	dup_all = list(dup_all)
	dup_all.sort(key=lambda x: x[t+1])
	if SingleFile == False:
		f=open('dup-'+mode+'.txt', 'w')
		f1=open('DUPXYZ-'+mode+'.txt', 'w')
		f.write('DUP = '+str(len(dup_all))+'\n')
		f1.write('DUP = '+str(len(dup_all))+' (2nd is the DUP)\n')
		# use t!!!
		for i in dup_all:
			f1.write(i[0]+' '+i[t]+'\n')
			for j in range(len(i)):
				if t==2:
					if j!=1 and j!=3:
						f.write(str(i[j])+' ')
				elif t == 3:
					if j!=1 and j!=2 and j!=4 and j!=5:					
						f.write(str(i[j])+' ')
				elif t == 4:
					if j!=1 and j!=2 and j!=3 and j!=5 and j!=6 and j!=7:					
						f.write(str(i[j])+' ')
			f.write('\n')
		f1.close()
		f.close()
		# copying unique structures
		os.mkdir('unique-'+mode)
		# conf_E after program contains only unique structures
		for i in conf_E:
#			Copy = True
#			for j in dup_all:
#				if i == j[2]:
#					Copy = False
#			if Copy == True:
#				os.system('cp '+i+' '+'unique-'+mode)
				os.system('cp '+i[0]+' '+'unique-'+mode)
	else:
		f=open(sys.argv[1][0:-4]+'-dup-'+mode+'.txt', 'w')
		f.write('DUP = '+str(len(dup_all))+'\n')
		f1=open(sys.argv[1][0:-4]+'-DUPXYZ-'+mode+'.txt', 'w')
		f1.write('DUP = '+str(len(dup_all))+' (2nd is the DUP)\n')
		for i in dup_all:
			f1.write(i[0]+' '+i[t]+'\n')
			for j in range(len(i)):
				if t==2:
					f.write(str(i[j])+' ')
				elif t == 3:
					if j!=2 and j!=5:					
						f.write(str(i[j])+' ')
				elif t == 4:
					if j!=2 and j!=3 and j!=6 and j!=7:					
						f.write(str(i[j])+' ')
			f.write('\n')
		f.close()
		f1.close()
		os.system('mv '+sys.argv[1][0:-4]+'-dup-'+mode+'.txt ..')
		os.system('mv '+sys.argv[1][0:-4]+'-DUPXYZ-'+mode+'.txt ..')
		f1=open(sys.argv[1][0:-4]+'-unique-'+mode+'.xyz', 'w')
		for i in conf_E:
#			Copy = True
#			for j in dup_all:
#				if i[0] == j[2]:
#					Copy = False
#			if Copy == True:
				f0 = open(i[0], 'r')
				for line in f0:
					f1.write(line)
				f0.close()
		f1.close()
		os.system('mv '+sys.argv[1][0:-4]+'-unique-'+mode+'.xyz ..')
		os.chdir(path)
		os.system('rm -rf tmp')
