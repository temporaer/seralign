import re, string, os, glob

graphs = {}

atomRe = re.compile('atm\((d\d+),(d\d+\_\d+),(\S+),(\S+),(\S+)\)')
bondRe = re.compile('bond\((d\d+),(d\d+\_\d+),(d\d+\_\d+),(\d+)\)')

fh = open('atom_bond.pl')
data = fh.readlines()
fh.close()

def setDoubleHashVal (hash,key1,key2,val,listSortingOn=0,useDict=False):
    try:
        hash[key1]
    except:
        if useDict and dictAvailable:
            hash[key1] = Dict()
        else:
            hash[key1] = {}
    hash[key1][key2] = val

def setTripleHashVal (hash,key1,key2,key3,val,listSortingOn=0,useDict=False):
    try:
        hash[key1]
    except:
        if useDict and dictAvailable:
            hash[key1] = Dict()
        else:
            hash[key1] = {}
        pass
    try:
        hash[key1][key2]
    except:
        if useDict and dictAvailable:
            hash[key1][key2] = Dict()
        else:
            hash[key1][key2] = {}
        pass
    hash[key1][key2][key3] = val

def setBond(graphs,example,atomsName1,atomsName2,bondType):
    try:
        graphs[examp]['atoms'][atomsName1]
        graphs[examp]['atoms'][atomsName2]
    except:
        raise "Atoms not known : %s %s" % (atomsName1,atomsName2)
    setDoubleHashVal(graphs[examp]['bonds'],atomsName1,atomsName2,bondType)
    setDoubleHashVal(graphs[examp]['bonds'],atomsName2,atomsName1,bondType)

def getBond(graphs,example,atomsName1,atomsName2):
    try:
        return graphs[example]['bonds'][atomsName1][atomsName2]
    except:
        pass
    return '0'
    

def toCSV(graphs,example,sorting=None):
    if not sorting:
        sorting = graphs[example]['atoms'].keys()
    retString  = example + '\n'
    retString += string.join(map(lambda x: '%s' % x,['None','None','None','None']+sorting),',')+'\n'
    #retString = 'None,' + string.join(map(lambda x: '%s_%s' % (graphs[example]['atoms'][x][2],x),sorting),',')+'\n'
    for atom in sorting:
        retVec = [atom, graphs[example]['atoms'][atom][2], graphs[example]['atoms'][atom][3], graphs[example]['atoms'][atom][4]]
        retVec += map(lambda x: getBond(graphs,example,atom,x),sorting)
        retString += string.join(retVec,',')+'\n'
    return retString

for line in data:
    atomMa = atomRe.search(line)
    if atomMa:
        (examp,atomsName,atomType,atomSomething,atomWeight) = atomMa.groups()
        try:
            graphs[examp]
        except:
            graphs[examp] =  {'atoms': {},'bonds': {}}
            pass
        atomsCount = len(graphs[examp]['atoms'].keys())
        graphs[examp]['atoms'][atomsName] = (atomsCount,atomsName,atomType,atomSomething,atomWeight)
    else:
        bondMa = bondRe.search(line)
        if bondMa:
            (examp,atomsName1,atomsName2,bondType) = bondMa.groups()
            setBond(graphs,examp,atomsName1,atomsName2,bondType)

files = glob.glob("../188/*.py")
for file in files:
	s = string.join(open(file,"r").readlines(),"\n")
	base = os.path.splitext(os.path.basename(file))[0]
	posfile = open(base+"-pos.csv", "w")
	negfile = open(base+"-neg.csv", "w")
	exec(s)
	for k in graphs.keys():
		try:
			p = pos[k]
			if p:
				posfile.write(toCSV(graphs,k))
				posfile.write("\n")
			else:
				negfile.write(toCSV(graphs,k))
				negfile.write("\n")
		except:
			continue
