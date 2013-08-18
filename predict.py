#!/usr/bin/python
import re,sys,os,subprocess
class FormatAll:
	def __init__(self,fname):
		#self.pwd=os.path.dirname(os.path.abspath(__file__))+"/results/"
		self.frname=fname+".fasta"
		self.fwname=fname+".arff"
		self.seq={}
		self.d={}
		for x in open(self.frname):
			x=x.rstrip()
			m=re.search(">(.*)",x)
			if m:
				p=m.group(1)
				self.d[p]={}
				self.seq[p]=""
				continue
			self.seq[p]+=x		
	def FormatProsite(self):				
		l=[]
		for x in open("motifs.list"):
			x=x.rstrip()
			l.append(x)
		for p in self.d.keys():
			for x in l:
				self.d[p][x]=0
		fw=open("temp_"+self.frname,"w")
		#subprocess.call(["./ps_scan.pl",self.frname],stdout=fw)
		os.system("./ps_scan.pl %s >%s"%(self.frname,"temp_"+self.frname))#careful
		lr=[]
		for x in open("temp_"+self.frname):#careful
			x=x.rstrip()
			m=re.search(">(.*?)\s:\s(.*?)\s",x)
			if m:
				p=m.group(1)
				if not p in lr:
					lr.append(p)
				s=m.group(2)
				continue
			self.d[p][s]+=1
		#print len(lr)
		#os.system("rm temp")
	def hodc(self,seq,n):
		l=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
		self.dc={}
		for x in l:
			for y in l:
				self.dc[x+"x"*(n-1)+y]=0.
		i=0
		while i<=len(seq)-n-1:
			s=seq[i]+"x"*(n-1)+seq[i+n]
			if s in self.dc.keys():
				self.dc[s]+=1
			i+=1
		for x in self.dc.keys():
			self.dc[x]=self.dc[x]/(len(seq)-n)
	def FormatCksaap(self,n):
		for x in self.d:
			for i in range(1,n):
				self.hodc(self.seq[x],i)
				self.d[x].update(self.dc)
	def FormatArff(self):
		f=open(self.fwname,"w")
		l=['IxxxxE', 'IxxxD', 'KxxD', 'IxxxxxK', 'PxxxxS', 'KxxxD', 'IxD', 'IxxK', 'PxxS', 'LxxxxK', 'QI', 'IxxxxxD', 'VxK', 'DxxxxV', 'LxxxxxK', 'KxxxA', 'AxxxK', 'SxP', 'PS00007', 'DxxxxxK', 'QxxxxP', 'AxM', 'KxxxxI', 'AN', 'ExxxxV', 'HxxxD', 'PxxxS', 'PxP', 'TxD', 'GxxxI', 'VxxxK', 'HxxxxS', 'KxxxxV', 'QxxP', 'YxxxxxQ', 'PS00016', 'EI', 'AH', 'HxY', 'MxxH', 'LxxE', 'ExxK', 'DxL', 'AxxC', 'SP', 'IxW', 'AxxxE', 'QL', 'GxxA', 'DxxxxT', 'DK', 'DxxxxxL', 'NxxxxxM', 'VxE', 'PP', 'PS51233', 'PS51211', 'NxxxA', 'DxxxI', 'DxxxxK', 'ID', 'IxV', 'ExR', 'PxxH', 'KxL', 'PxxxR', 'HxxxxR', 'MxxxxN', 'LxxW', 'MH', 'AxP', 'KxV', 'NxxA', 'LxxxM', 'DxxxxxI', 'RxI', 'HxF', 'SxW', 'KxxV', 'CxD', 'FxP', 'CxxxxxE', 'DxxxxW', 'DxxI', 'PS00342', 'PS00995', 'RxxV', 'AK', 'VxxA', 'AI', 'ExxxxxI', 'KxxxxY', 'AxxE', 'KxY', 'KxxxxE', 'CxxxxE', 'AxxxxxM', 'GxxxxD', 'WxxxxxE', 'PxxxxT', 'VxxxxD', 'IxxxxA', 'ND', 'VxM', 'LxxxxxW', 'TxxE', 'ExxxxW', 'IxxxM', 'EA', 'HxL', 'KxT', 'QxxxV', 'QxxxR', 'AxxxH', 'SxxxxP', 'VxxxxxK', 'DA', 'GxxI', 'CxxxxY', 'AxxxxP', 'IxxG', 'PS00109', 'PS50853', 'GV', 'PxxxP', 'DxxxxL', 'PxS', 'LxxxxxD', 'LC', 'AxxxxxD', 'AxxxxxK', 'PT', 'WxH', 'KxxxxxN', 'YxxxxR', 'SxA', 'DL', 'AxxxV', 'DxxG', 'PS50010', 'AxxK', 'ExxxxxY', 'ExxxxxF', 'YxxxxxF', 'KxxF', 'AP', 'HxxxxxP', 'LxK', 'DxxxM', 'GxI', 'SxxP', 'EL', 'PS00299', 'PS50212', 'PS00720', 'PS50009', 'CS', 'TxxxxxC', 'AxxD', 'GxN', 'IxxxK', 'PxxxxA', 'PS00005', 'KxxxxA', 'HxxxP', 'SxxxQ', 'RxxxxP', 'RxT', 'PxxxxP', 'QxxxxxR', 'YxxD', 'WP', 'SxxxH', 'ExG', 'PS50835', 'KxxxxxI', 'TxxK', 'YxxxxQ', 'VN', 'FxE', 'QxxxxV', 'KxxxxC', 'MxxxxxH', 'AxxxxI', 'KxxT', 'KA', 'PS50053', 'GxxxxxI', 'IxE', 'WxxxxxT', 'SxQ', 'FK', 'PxH', 'ExxxxN', 'KxxxxxM', 'PxA', 'GxxxxxV', 'KxxxE', 'DxxxN', 'QxxxA']
		f.write("@relation "+self.fwname+"\n\n")
		for x in open("arff"):#care arff
			f.write(x)
		f.write("\n")
		f.write("@data\n")
		#p=self.d.keys()[0]
		for y in self.d.keys():
			fm=y+","
			for x in l:
				fm+=str(self.d[y][x])+","
			fm+="yes"
			f.write(fm+"\n")
		f.close()
	def EasyFormat(self):
		self.FormatProsite()
		self.FormatCksaap(7)
		self.FormatArff()
class RandomForestPredict:
	def __init__(self,fname):
		d={}
		for i in range(1,11):
			p=os.popen("java weka.classifiers.trees.RandomForest -T %s -l %s -p 1"%(fname+".arff","pe_"+str(i)+"_cksaap_prosite_200.model"))
			p=p.readlines()
			for x in p:
				x=x.rstrip()
				m=re.search("\+",x)
				if m:
					m2=re.search("\+\s+(.*?)\s+?\((.*?)\)",x)
					if m2:
						#print m2.group(1),m2.group(2)
						if not m2.group(2) in d.keys():
							d[m2.group(2)]=1-float(m2.group(1))
						else:
							d[m2.group(2)]+=1-float(m2.group(1))
				else:
					m2=re.search("1:yes      1:yes\s+(.*?)\s\((.*?)\)",x)
					if m2:
						#print m2.group(1),m2.group(2)
						if not m2.group(2) in d.keys():
							d[m2.group(2)]=float(m2.group(1))
						else:
							d[m2.group(2)]+=float(m2.group(1))
		l=[]
		for x in d:
			l.append((x,d[x]/10))
		#print l
		l=sorted(l,key=lambda l:l[1])
		#print l
		self.r=l
if __name__=="__main__":
	if len(sys.argv)==3:
		fname=sys.argv[1]
		fw=sys.argv[2]
		p=FormatAll(fname)
		p.EasyFormat()
		p=RandomForestPredict(fname)
		r=p.r
		f=open(fw,"w")
		for x in r:
			f.write(x[0]+"\t"+str(x[1])+"\n")
		f.close()
	else:
		print "\npredict.py <FastaFile> <ResultsFile>\n"
