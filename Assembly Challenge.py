#Foroutan(2577847)-Khojasteh(2567833),-Assignment5-Exercise1

import sys
import collections
from Bio import SeqIO, Seq, SeqRecord
from Bio.SeqRecord import SeqRecord


#reverse complement Function
def rev_com (list_seq):
    return Seq.reverse_complement(list_seq)

#possible next K_mers
def after(km):
    for x in 'ACGT':
        yield km[1:]+x

#possible before K_mers
def before(km):
    for x in 'ACGT':
        yield x + km[:-1]

def dif(x,y,size):
    d=0
    for i in range(0,size):
        if x[i]!=y[i]:
            d+=1
            if d>1:
                break
    return d

     
#Getting reads from FASTA file
def build_kmers(file_in,size):
    
    print(file_in,size)
    
    reads=SeqIO.parse(open(file_in, mode='r'), 'fasta')
            
#******************************************************************************
    #task a) 
    
    #creating all k_mers of reads
    dic = collections.defaultdict(int) #dictionary of k_mers with the number of occurance (their coverage)
    incoming_nei_dic=collections.defaultdict(list) #the value of dic is the list of incoming k_mer neighbours
    outgoing_nei_dic=collections.defaultdict(list) #the value of dic is the list of outgoing k_mer neighbours
    
    k_mers=[]
    for read in reads:
        #print("read ", read.id,": ",read.seq)
        list_read=str(read.seq) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ .seq
        #print(list_read)
        
        for k in range(0,len(list_read)-size+1):
            k_mers.append(list_read[k:k+size])
        
        #reverse_complement
        #and add its K_mers
        list_read = rev_com(list_read)
        #print("reverse_complement",'\n',list_read)
        for k in range(0,len(list_read)-size+1):
            k_mers.append(list_read[k:k+size])

    #Assign empty list   
    for km in k_mers:
        incoming_nei_dic[km]=[] 
        outgoing_nei_dic[km]=[] 
    
    for km in k_mers:
        dic[km] += 1  #compute coverage of each k_mer

    #omit k_mers with low coverage
    
    limit=1
    dt = [km for km in dic if dic[km] <= limit]
    dt2=[]
    if len(dic)<8000:
        for km in dt:                
            for x in dic:
                if dif(km,x,size)==1:
                    dt2.append(km)
                
    
        for km in dt2:
            del dic[km]
    else:
        for km in dt:
            del dic[km]
            
  #neighbours        
    for km in dic:    
        #find neighbours
        for x in after(km): #possible next neighbour k_mer
            if x in dic and (x not in outgoing_nei_dic[km]):    # if they exist in dic
                outgoing_nei_dic[km].append(x)

        for x in before(km): #possible next neighbour k_mer
            if x in dic and (x not in incoming_nei_dic[km]):    # if they exist in dic
                incoming_nei_dic[km].append(x)
                
    #print("in",incoming_nei_dic)
    #print("out",outgoing_nei_dic)
    
        
   # print(dic)    
    return k_mers,dic,incoming_nei_dic,outgoing_nei_dic
 
    
            
def input_fun(file_in,file_out,size,query_K_mer):
  #creating all k_mers
    dic = collections.defaultdict(int)
    incoming_nei_dic = collections.defaultdict(list) 
    outgoing_nei_dic = collections.defaultdict(list)
    k_mers=[]
    k_mers,dic,incoming_nei_dic,outgoing_nei_dic = build_kmers(file_in,size)
    

  #print neighbours
    if incoming_nei_dic[query_K_mer] == []:
        print("No incoming neighbour")
    else:
        print("incoming - ", ','.join(incoming_nei_dic[query_K_mer]))

    if outgoing_nei_dic[query_K_mer]==[]:
        print("No outgoing neighbour")
    else:
        print("outgoing - ", ','.join(outgoing_nei_dic[query_K_mer]))
        
    '''for n in dic:
        print(n,"....",incoming_nei_dic[n])'''

  
    
  #Contigs
    new_vs=[]
    starts=[] #nodes which have no incoming
    contigs=[]
    for km in dic:
        if len(incoming_nei_dic[km])==0 or len(incoming_nei_dic[km])>1: #indegree(node)==0 or >1
            starts.append(km)
    #print(starts)
               
    for n in starts:
        temp = n   #start of the contig
        if len(outgoing_nei_dic[n])==1:
            next_n=outgoing_nei_dic[n][0]
            
               
        while len(outgoing_nei_dic[n])==1 and len(incoming_nei_dic[next_n])==1:
            temp = temp + next_n[-1] #merge
            n=next_n
            if len(outgoing_nei_dic[n])!=1:
                break
            next_n=outgoing_nei_dic[n][0]
            
                
        contigs.append(temp)
        new_vs.append(temp)
        
        '''if len(outgoing_nei_dic[n])>1:
            for nei in outgoing_nei_dic[n]:
                starts.append(nei)'''

    #print(new_vs)
    
#################################### output
    
    max_val=new_vs[0]
    lengths=[len(new_vs[0])]
    for i in range(1,len(new_vs)):
        lengths.append(len(new_vs[i]))
        if len(max_val) < len(new_vs[i]):
            max_val = new_vs[i]
    #print(lengths)
               
    num=1
    print(">", num,"\n", max_val)
    num+=1
    print(">", num,"\n", Seq.reverse_complement(max_val))
    num+=1

        
    print("\nThe length of the obtained sequence:",len(max_val))

   #writing in a file
    f= open(file_out,"w+")

    #writing assembled seq in the out put file
    num=1
    f.write(">%d\n" %num)
    f.write("%s\n"%max_val)
    num+=1
    f.write(">%d\n" %num)
    f.write("%s\n"%Seq.reverse_complement(max_val))

    '''#writing all contigs in the out put file
    for co in contigs:
        f.write("%s\n"%co)'''
    
    f.close()

    

if __name__ == '__main__':
    input_fun(sys.argv[1], sys.argv[2], int(sys.argv[3]), str(sys.argv[4]))
