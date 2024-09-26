# GOAL: reconstruct the hopcount tool we use in Tufts Galaxy (galaxy.med.tufts.edu) that is not present on usegalaxy.org.  

import glob 

sams = glob.glob('/Users/tailob/Desktop/jax/repositories/tnseq-meta-integration/multi_fasta_SAM/*/*.sam')

for s in sams:
	datadict = {}
	for line in open(s,'r').readlines():
		if not line.startswith('@') and not line.startswith('/'):
			info = line.split('\t')
			flag = info[1]; ref = info[2]; site = int(info[3]); seq = info[9]
			if flag == '0': 
				if (ref,site) in datadict: 
					datadict[(ref,site)][0]+=1
				else: datadict[(ref,site)]=[1,0]
			elif flag == '16':
				if (ref,site+15) in datadict: 
					datadict[(ref,site+15)][1]+=1
				else: 
					datadict[(ref,site+15)]=[0,1]
					
					
	output = open(s[:-3]+'tabular','w')
	output.write('Reference\tPosition\tLocus\tGene\tPlusCount\tMinusCount\tTotalCount\tProduct\tProteinID\tNote\tSequence\n')
	datakeys = sorted(datadict.keys(),key=lambda x: x[1])
	for k in datakeys: 
		output.write(k[0]+'\t'+str(k[1])+'\t\t\t')
		output.write('%i\t%i\t%i\t\t\t\t\n' %(datadict[k][0],datadict[k][1],datadict[k][0]+datadict[k][1]))