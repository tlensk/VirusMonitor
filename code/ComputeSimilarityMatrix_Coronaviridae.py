"""
Updated on 03/18/20
@author: Tatiana Lenskaia
"""



def CountGC(seq):
    seq = seq.lower();
    n_seq = len(seq)
    
    n_a = seq.count("a");
    n_c = seq.count("c");
    n_g = seq.count("g");
    n_t = seq.count("t");
    
    n_bases = n_a + n_c + n_g + n_t
    
    n_other = n_seq - n_bases;
    
    if n_bases != 0:
        gc = round(1.0*(n_c+n_g)/n_bases*100,2)
    else:
        gc = "N/A"
    
    return [gc, n_other]



def GetListFromFile(fInName):
    t = []
    fIn = open(fInName, "r")
    for line in fIn:
        line = line.strip()
        if line != "":
            if line not in t:
                t.append(line)
    #print(len(t))
    fIn.close()
    return t


def GetDictFromFile(fInName, sep, header):
    fIn = open(fInName,"r")
    lines = fIn.readlines()
    if header == 1:
        header_line = lines[0]
        lines = lines[1:]
    
    d = {}
    t = []
    
    n = len(lines[0].split(sep))
    
    for line in lines:
        line = line.strip()
        if line != "":
            t_line = line.split(sep)
            if len(t_line) != n:
                print("Check format!", line, t_line)
                if t_line[0] not in d:
                    t.append(t_line[0])
                    d[t_line[0]] = t_line[1:]
                
            else:
                if t_line[0] not in d:
                    t.append(t_line[0])
                    d[t_line[0]] = t_line[1:]
                else:
                    print("First colum has non-unique values!")
    fIn.close()
    return [t,d, header_line]



def GetText(finName):
    '''Extracts text from a single fasta file'''
    fin = open(finName, 'r')
    text = ''
    for line in fin:
        line = line.strip()
        if (line != "") and (line[0] != '>'):
            text = text + line
    fin.close()
    return text



def Rev_cmp_upper(st):
    '''Computes the reverse complement of a string'''
    st = st.upper()
    cmp_st = st.translate(str.maketrans("ACGT","TGCA"))
    rev_cmp_st = cmp_st[::-1]
    return rev_cmp_st;

def CreateDictLocD_upper(text, m, gtp = "l"):
    if (m <= 0) or (m > len(text)):
        print("(n = "+str(m)+") is not a valid window size for this genome!!!");
        return {}
	
    d_g = dict()
    nn = len(text);
    gtype = gtp.lower();
    gtype = gtype[0];
	
    if gtype == "c":
        text = text + text[0:(m-1)];
        lastpos = nn;
    elif gtype == "l":
        lastpos = nn-m+1;
    else:
        print("Is this genome linear or circular?");
        return d_g;
		
    for ii in range (lastpos):
        bl = text[ii:(ii+m)]; 
        bl = bl.upper();
        if bl in d_g :
            dd = d_g[bl];
            dd[ii] = ""
            d_g[bl] = dd;
        else:
            dd = {}
            dd[ii] = ""
            d_g[bl] = dd;
    return d_g;	


def FindIntersection(d_g11, d_g22):
	
	t_ints = list()
	nn1 = len(d_g11);
	nn2 = len(d_g22);
	
	if nn1 <= nn2:
		d_first = d_g11;
		d_last = d_g22;
	else:
		d_first = d_g22;
		d_last = d_g11;
	
	for s in d_first:
		if s in d_last:
			t_ints.append(s);
	return t_ints;




def PrintMatrix(k, pref, row_index, col_index, matrix, sep = ","):
    fOut = open(str(k)+pref+".csv","w");
    fOut.write(str(k)+sep+sep.join(str(x) for x in col_index)+"\n")
    for i in range(len(matrix)):
        s = str(row_index[i])
        for el in matrix[i]:
            s = s+sep+str(el)
        fOut.write(s+"\n")
        
    fOut.close()
    return

def PrintSubMatrix(k, pref, row_names, row_index, col_names, col_index, matrix, sep = ","):
    if row_index != [] and col_index != []:
        fOut = open(str(k)+"_"+pref+".csv","w");
        #fOut.write(str(k)+sep+sep.join(str(x) for x in col_index)+"\n")
        s = str(k)
        for jj in range(len(col_index)):
            s = s+sep+col_names[col_index[jj]]
        fOut.write(s+"\n")

        
        for ii in range(len(row_index)):
            i = row_index[ii]
            s = str(row_names[row_index[ii]])
            
            for jj in range(len(col_index)):
                j = col_index[jj]
                
                if i == j:
                    sign = "-"
                else:
                    sign = ""
                    
                s = s+sep+sign+str(matrix[i][j])
            
            fOut.write(s+"\n")
            
        fOut.close()
    return


#============================================




fInName = "3168_desc.txt"
res = GetDictFromFile(fInName,"\t",1)
d_host= res[1]
print(len(d_host))

t_list = res[0][0:2]

n_list = len(t_list)
pref = "3168_"

# Window length
m = 40

matrix = [[0 for i in range(n_list)] for j in range(n_list)] 

for i in range(n_list):
    fInNameA = t_list[i].split(".")[0]+".fasta"
    textA = GetText(fInNameA)
    dictA =  CreateDictLocD_upper(textA, m, "l")

    
    j = i-1
    #print(i)
    while(j < n_list-1):
        j = j+1
        fInNameB = t_list[j].split(".")[0]+".fasta"
        textB = GetText(fInNameB)
        dictB =  CreateDictLocD_upper(textB, m, "l")
        t_ints = FindIntersection(dictA,dictB)
        n_ints = len(t_ints)
        if n_ints != 0:
            matrix[i][j] = n_ints
            matrix[j][i] = n_ints

PrintMatrix(m,"_corona(3168)", t_list, t_list, matrix, sep = ",")



