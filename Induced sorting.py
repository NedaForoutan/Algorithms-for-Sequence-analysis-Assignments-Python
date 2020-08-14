#Foroutan(2577847)-Khojasteh(2567833).-Assignment3-Exercise3

#reading the text input file
name_text=input("Enter the name of text file:") #The name of the text is got from input

file=open(name_text)
if file.mode=='r':
    text=file.read() # text is read

file.close()

#getting the sorted LMS positions
sorted_LMS=input("Enter the sorted LMS positions:") #the sorted LMS positions are gotten as an input
LMS=sorted_LMS.split(',')
LMS = [int(i) for i in LMS] #List of int of Sorted_LMS

#########################################

s=text    #+'$'  text file is included $
n=len(s)

#compute the type array for s sequence
type_array=[]
type_array.append('S') #type[n-1]='S'

for p in range((n-2),-1,-1):
    if s[p]>s[p+1] :
        type_array.insert(0,'L')
    elif s[p]<s[p+1] :
        type_array.insert(0,'S')
    else : # s[p]==s[p+1]
        type_array.insert(0,type_array[0])


#sorted LMS positions are gotten from input and splitted as LMS

#implementing buckets
buckets=['$']
for i in range (0,n):
    if s[i] not in buckets:
        buckets.append(s[i])
buckets.sort() #['$','a','c','g','t']
begin=[0]   #beginning of each buckets
end=[0]     #end of each buckets
for i in range(1,len(buckets)):
    num=s.count(buckets[i])
    begin.append(1+end[i-1])  #[0,1,7,13,15]
    end.append(num+end[i-1])  #[0,6,12,14,21]

                                
#initializing pos with unknown at each position
pos = ['?']*n

#writing the positions of the LMS-positions at the end of their respective buckets.
for i in range((len(LMS)-1),-1,-1):
    index=buckets.index(s[LMS[i]])
    write=False
    while(write==False):
        for position in range(end[index],begin[index]-1,-1):
            if pos[position]=='?':
                pos[position]=LMS[i]
                write=True
                break

pos[0]=n-1 #pos[$]=len(s)-1 #For example 21

#Sorting the L-positions

for i in range(0,n):    #from the left to the right
    if pos[i]=='?' or type_array[pos[i]-1]=='S':   #skip
        continue
    elif type_array[pos[i]-1]=='L':
        index=buckets.index(s[pos[i]-1])
        write=False
        while(write==False):
            for position in range(begin[index],end[index]+1):
                if pos[position]=='?':
                    pos[position]=pos[i]-1
                    write=True
                    break


#omitting S-positions from pos in order to sort them
for i in range(0,n):
    if pos[i]=='?':
        continue
    if type_array[pos[i]]=='S':
        pos[i]='?'


#Sorting the S-positions

for i in range(n-1,-1,-1): #from right to the left
    if pos[i]=='?' or type_array[pos[i]-1]=='L':   #skip
        continue
    elif type_array[pos[i]-1]=='S':
        index=buckets.index(s[pos[i]-1])
        write=False
        while(write==False):
            for position in range(end[index],begin[index]-1,-1):
                if pos[position]=='?':
                    pos[position]=pos[i]-1
                    write=True
                    break



#print suffix array in a comma-separated list
print("suffix array:")
print(pos[1],end="")
for i in range(2,n):
    print(',',pos[i],end="")
        
























