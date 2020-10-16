"""
Updated on 03/18/20
@author: Tatiana Lenskaia
"""

import time as tm
import sys
import os
import copy
"""
if len(sys.argv) == 2:
    path = sys.argv[1]+"\\"
    cur_dir = sys.argv[1] + " directory"
else:
    path = ""
    cur_dir = "the current directory"



if not os.path.exists(path+'Output_files'):
    os.makedirs(path+'Output_files')


def CheckInputFiles(path):
    if os.path.exists(path+"list.txt"):
        fList = open(path+"list.txt")
        
        flag = 0
        
        for line in fList:
            line = line.strip()
            if line != "":
                if (os.path.exists(path+"Input_files\\FASTA\\"+line+".fasta")==False) or (os.path.exists(path+"Input_files\\CRISPRCasdb\\"+line+".csv") == False):
                    flag = 1
                    break
        fList.close()
        
        if flag == 1:
            print("Not all input files are supplied")
            print("Please check that .fasta and .csv files are provided for each replicon accession number from list.txt")
    else:
        
        
        print("The list of replicons is missing in", cur_dir)
        print("Please provide a list of replicon accession numbers (list.txt)")
    return
CheckInputFiles(path)
"""
#=============================================

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
    return [t,d]


#Updated: May 4, 2019    
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


def SearchOrganism_plm_cas(fRname, name,kind, count, ddd, d_rep_sps, d_fname):
    '''Analyzes self-targeting events in a given organism'''
    ct_chr = 0
    ct_plm = 0
    ct_other = 0
       
    d_sps = {}

    for seq in ddd:
        
        if "chromo" in ddd[seq].lower():
            ct_chr = ct_chr + 1
        elif "plasmid" in ddd[seq].lower():
            ct_plm = ct_plm + 1
        else:
            ct_other = ct_other+1
        
        
        if seq in d_rep_sps:
 
            for sp in d_rep_sps[seq]:
                nn = len(sp)
                if nn not in d_sps:
                    dd = {}
                    dd[sp] = [seq]
                    d_sps[nn] = dd
                else:
                    dd = d_sps[nn]
                    if sp in dd:
                        dd[sp].append(seq)
                    else:
                        dd[sp] = [seq]
                    d_sps[nn] = dd
                  
    # d_sps contains a dictionary of spacers for a genome
    #print(d_sps,"\n")
    
    d_res = {}
    d_gcas = {}
    d_gcrs = {}
    
    sep2 = ","
    
    for rep in d_gens[name] :
        fIn1 = open(path+"Input_files\\CRISPRCasdb\\"+rep+".csv","r")
        lines1 = fIn1.readlines()
        d_crs = {}
        d_cas = {}
        for ln in lines1:
            ln = ln.strip()
            if (ln != ""):
                t_ln = ln.split(sep2)
        
                r_id = t_ln[0]
                r_class = t_ln[1]
                r_name = t_ln[2]
                
        
                if r_class.lower() == "locuscrispr":
                    if r_id not in d_crs:
                        d_crs[r_id] = t_ln[2:]
                        
                if r_class.lower() == "clustercas":
                    if r_id not in d_cas:
                        d_cas[r_id] = t_ln[2:]
                            
                        
        fIn1.close()
    
        d_gcas[rep] = d_cas
        d_gcrs[rep] = d_crs

    
    
    for nn in d_sps:
        
        for rep in d_gens[name] :
            
            d_crs = d_gcrs[rep]
            d_cas = d_gcas[rep]
            
            fInName = path+"Input_files\\FASTA\\"+rep+".fasta"
            
            text = GetText(fInName)
            
            if len(text) != int(d_fname[rep][1]):
                print("Mismatch!!!")
            else:
                d_b = CreateDictLocD_upper(text,nn,d_fname[rep][2])

           
                for sp in d_sps[nn]:
                    sp_rc = Rev_cmp_upper(sp)
                    sp_len = len(sp)
                                  
                    # sp spacer
                    # nn length of a spaccer                  
                    
                    c_dir = 0
                    c_rev = 0
                    
                    c_in = 0
                    c_out = 0
                    
                    
                    
                    if sp in d_b:
                        c_dir = len(d_b[sp])
                        
                        t = sorted(d_b[sp].keys())
                        
                        
                        for pos in t:
                            fl = 0
                            
                            for it in d_crs:
                                ttt = d_crs[it]
                                if (pos >= int(ttt[1])) and (pos+sp_len <= int(ttt[1])+int(ttt[2])):
                                    fl = 1
                                    break
                
                            if fl == 1:
                                d_b[sp][pos] = "F"+"|"+ttt[0]+"_x"
                                c_in = c_in+1
                            else:
                                d_b[sp][pos] = "F"+"|"+"-"
                                c_out = c_out+1
                           
                    
                    if sp_rc in d_b:
                        c_rev = len(d_b[sp_rc])
                        
                        t = sorted(d_b[sp_rc].keys())


                    
                        for pos in t:
                            fl = 0
                            
                            for it in d_crs:
                                ttt = d_crs[it]
                                
                                if (pos >= int(ttt[1])) and (pos+sp_len <= int(ttt[1])+int(ttt[2])):
                                    fl = 1
                                    break
                
                            if fl == 1:
                                d_b[sp_rc][pos] = "R"+"|"+ttt[0]+"_x"
                                c_in = c_in+1
                            else:
                                d_b[sp_rc][pos] = "R"+"|"+"-"
                                c_out = c_out+1

                    
                    if sp in d_res:
                        if rep not in d_res[sp]:
                            d_res[sp][rep] = [c_dir,c_rev, c_in, c_out]
                        else:
                            print("Duplicate seq", rep)
                    else:
                        d_res[sp] = {}
                        if rep not in d_res[sp]:
                            d_res[sp][rep] = [c_dir,c_rev, c_in, c_out]
                        else:
                            print("Duplicate seq", rep)
                d_b.clear()
                
                
    
    num_sp = len(d_res)
    num_stsp = 0
    
    
    for sp in d_res:

        sp_total = 0
        sp_in = 0
        sp_out = 0
        sp_dir = 0
        sp_rev = 0
        
        sp_len = len(sp)
        
        for rep in d_res[sp]:
            total = d_res[sp][rep][0]+d_res[sp][rep][1]
            sp_dir = sp_dir + d_res[sp][rep][0]
            sp_rev = sp_rev + d_res[sp][rep][1]
            sp_total = sp_total+total
            sp_in = sp_in + d_res[sp][rep][2]
            sp_out = sp_out +d_res[sp][rep][3]

        
        if sp_out > 0:
            num_stsp = num_stsp + 1
        
        
        
        

        cas_ind = ""
        st_ind = ""
        plm_ind = ""
        
        if sp_len in d_sps:
            d_len_sps = d_sps[sp_len]
    

            
            
            t_cas = []
            fl_plm = 0
            for it in d_len_sps[sp]:
                
                d_cas = d_gcas[it]
                t_cas.append(d_cas)
                if "plasmid" in d_gens[name][it].lower():
                    fl_plm = 1
            
            if fl_plm == 0:
                plm_ind = "No"
            else:
                plm_ind = "Yes"
                
                
             
            flag = 0
            for jj in range(len(t_cas)):
                if len(t_cas[jj]) > 0:
                    flag = 1
                    
            if flag ==1:
                cas_ind = "Yes"
            else:
                cas_ind = "No"
        
        if sp_out > 0:
            st_ind = "Yes"
        else:
            st_ind = "No"
            
        
        fOutSP = open(fRname,"a")
        print(kind,name,sp,len(sp), sp_total, sp_in, sp_out, st_ind,plm_ind,cas_ind, file = fOutSP, sep = "\t")
        fOutSP.close()
        

    return [ct_chr, ct_plm,ct_other, num_sp, num_stsp]
    
'''
def PrintMatrix(k, pref, row_index, col_index, matrix, sep = ","):
    fOut = open(str(k)+"_"+pref+".csv","w");
    fOut.write(str(k)+sep+sep.join(str(x) for x in col_index)+"\n")
    
    n_cols = len(col_index)
    n_rows = len(row_index)
    
    
    for i in range(n_rows):
        s = row_index[i]
        print(i,matrix[i])
    
        for el in matrix[i]:
            s = s+sep+str(el)
        
        fOut.write(s+"\n")
        
    fOut.close()
    return
'''

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
'''
t = ['1','2','3']
n = len(t)

matrix = [[0 for i in range(n)] for j in range(n)]
print(matrix)

n_r = len(matrix[0])
n_c = len(matrix)

#print(n_r,n_c)
#PrintMatrix(n, "test", t, t, matrix, sep = ",")

fOut = open(str(m)+"_corona.csv","w")
for row in matrix:
    print(row)
    
    st = ",".join(row)
    print(st)


'''




# Summary for replicons in genomes
fInName = "3168_summary.txt"
res = GetDictFromFile(fInName,"\t",0)
d_fname = res[1]
#print(d_fname)
print(len(d_fname))


"""
fInName1 = "3168_allinfo.txt"
res1 = GetDictFromFile(fInName1,"\t",1)
d_host = res1[1]
#print(len(d_host))
#print(d_host)
"""


fInName = "3168_complete genomes_NCBIVirus.txt"
res = GetDictFromFile(fInName,"\t",1)
d_host= res[1]
print(len(d_host))


#Full list
#t_list =  sorted(d_fname.keys())

# Reduced list from a file
t_list = GetListFromFile("test_sample.txt")
#t_list = GetListFromFile("block_ids.txt")



n_list = len(t_list)
pref = "test_"

# Window length
m = 40

matrix = [[0 for i in range(n_list)] for j in range(n_list)] 

for i in range(n_list):
    fInNameA = t_list[i].split(".")[0]+".fasta"
    textA = GetText(fInNameA)
    #print(fInNameA, len(textA))

    dictA =  CreateDictLocD_upper(textA, m, d_fname[t_list[i]][3])
    j = i-1
    print(i)
    while(j < n_list-1):
        j = j+1
        fInNameB = t_list[j].split(".")[0]+".fasta"
        textB = GetText(fInNameB)
        dictB =  CreateDictLocD_upper(textB, m, d_fname[t_list[j]][3])
        t_ints = FindIntersection(dictA,dictB)
        n_ints = len(t_ints)
        if n_ints != 0:
            #print(i,j,len(t_ints))
            matrix[i][j] = n_ints
            matrix[j][i] = n_ints
        #print(fInNameB, len(textB))
        #print(i,j)
    #print("\n")

#print(matrix)

PrintMatrix(m,"_corona(test)", t_list, t_list, matrix, sep = ",")


"""
t_left = copy.copy(t_list)
t_clusters = []

while(t_left != []):
    
    print("t_left",t_left)
    ind = t_list.index(t_left[0])
    #t = [t_left[0]]
    t = [ind]
    t_cl = [ind]
    print(t)
    t_left.remove(t_left[0])
    while(t != []):
        for i in range(n_list):
            if (matrix[t[0]][i] != 0) and (t[0]!= i):
                if i not in t_cl:
                #if t_list[i] not in t:
                    #t.append(t_list[i])
                    t.append(i)
                    t_cl.append(i)
                    
                    if t_list[i] in t_left:
                        t_left.remove(t_list[i])
        t.remove(t[0])
        print("t t_cl",t,t_cl)
        
    #print(t)
    if t_cl != []:
        t_clusters.append(t_cl)


print(len(t_clusters))
nn = 0
for t in t_clusters:
    nn=nn+1
    print(t)

    for it in t:
        print(t_list[it],d_fname[t_list[it]][0],d_host[t_list[it].split(".")[0]])
    PrintSubMatrix(m, "_cluster(bl)"+str(nn), t_list, t, t_list, t, matrix, sep = ",")
    
    print("\n")
"""
    
    


