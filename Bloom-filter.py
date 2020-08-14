#Foroutan(2577847)-Khojasteh(2567833),-Assignment4-Exercise3

import sys


def h1(s, n):
    hash = 5381
    for c in s:
        hash = (( hash << 5) + hash ) + ord (c)
    return hash % n

def h2(s, n):
    hash = 0
    for c in s:
        hash += ord (c)
    return hash % n


#adds the element (of type string) to the set
def add(element,bloom):
    bloom[h1(element,len(bloom))]=1
    bloom[h2(element,len(bloom))]=1
    #print(bloom)
    
    

#check whether the element is in the set or not
def contains(element,bloom):
    value_h1=h1(element,len(bloom))
    value_h2=h2(element,len(bloom))

    if(bloom[value_h1]==1 and bloom[value_h2]==1):
        print(value_h1,',',value_h2,',','T')
    else:
        print(value_h1,',',value_h2,',','F')

    
#Bloom filter
def bloom_filter(n,s1,s2):
    l1=s1.split() #list of add strings
    l2=s2.split() #list of test strings
    bloom=[0]*n
    #add the each string of iput file to the set
    for word in l1:
        add(word,bloom)

    for word in l2:
        contains(word,bloom)
        

    

if __name__ == '__main__':
    bloom_filter(int(sys.argv[1]), open(sys.argv[2]).read().strip(), open(sys.argv[3]).read().strip())
