## Implementation of Hashtable
# Author: Yang Yu
# Date: 2018/10/31
# Description: Hash table/Hash function/collision handling
# https://www.youtube.com/watch?v=KyUTuwz_b7Q

from linked_list import linked_list
import numpy as np

## Functions
class node:
    
    name = [] # Date within the node
    score = -1 # Score
    nextptr = -1 # Pointer
    
    def __init__(self,name=[],score=-1):
        self.name = name
        self.score = score
        
    def copy(self):
        node_new = node(self.name,self.score)
        node_new.nextptr = self.nextptr
    
    def data_print(self):
        print(name,':',score)
        
class linked_list_HT(linked_list):
    
    def __init__(self,node0):
        linked_list.__init__(self,node0)
        
    def find_node(self,s):
        L = linked_list.list_length(self)
        nodetmp = self.node0
        pos = -1
        for ii in range(L):
            if s == nodetmp.name:
               return nodetmp
            else:
                nodetmp = nodetmp.nextptr
        return -1
        
    def del_node(self,s):
        L = linked_list.list_length(self)
        nodetmp = self.node0
        nodeprev = -1
        pos = -1
        for ii in range(L):
            if s == nodetmp.name:
                pos = ii
                if ii == 0:
                    self.node0 = nodetmp.nextptr
                    nodetmp.nextptr = -1
                elif ii == L-1:
                    nodeprev.nextptr = -1
                    nodetmp.nextptr = -1
                else:
                    nodeprev.nextptr = nodetmp.nextptr
                    nodetmp.nextptr = -1
                self.nodeNum = self.nodeNum-1
                break
            else:
                nodeprev = nodetmp
                nodetmp = nodetmp.nextptr
        return pos
        
class Hash_table:
    
    hash_L = -1
    hash_llptrs = [] # Record the ptrs of the first element of each node
    
    def __init__(self,hash_L):
        self.hash_L = hash_L
        self.hash_llptrs = []
        for ii in range(hash_L):
            self.hash_llptrs.append(-1) # Set to -1 for empty

    def hashfunction(self,s):
        # Convert s characters to ASCII value
        s_val = np.array([ord(c) for c in s])
        # Calculate Hash val
        hash_indx = np.mod(np.sum(s_val),self.hash_L)
        
        return hash_indx
    
    def hash_table_add(self,s,score):
        hash_indx = self.hashfunction(s)
        node_new = node(s,score) # Create the new node for insertion
        
        if self.hash_llptrs[hash_indx] == -1: # Empty node
            ll_new = linked_list_HT(node_new)
            self.hash_llptrs[hash_indx] = ll_new
        else:
            self.hash_llptrs[hash_indx].append(node_new)
    
    def hash_table_add_all(self,s_ens,score_ens):
        L = len(s_ens)
        for ii in range(L):
            self.hash_table_add(s_ens[ii],score_ens[ii])
    
    def find_node(self,s):
        hash_indx = self.hashfunction(s) # Find the index
        # Find the s in the linked list
        ll = self.hash_llptrs[hash_indx]
        llnode = ll.find_node(s)
        
        return llnode
        
    def del_node(self,s):
        hash_indx = self.hashfunction(s) # Find the index
        # Find the s in the linked list
        ll = self.hash_llptrs[hash_indx]
        pos = ll.del_node(s)
        
        return pos
        
## Routines
# s_ens = ['Jan','Tim','Mia','Sam','Leo','Ted','Bea','Lou','Ada','Max','Zoe']
s_ens = ['Mia','Tim','Bea','Zoe','Sue','Len','Moe','Lou','Rae','Max','Tod'] # Collision case
scores = np.random.randint(0,100,len(s_ens))
hash_L = 11 # For collision test
        
HT = Hash_table(hash_L)

# Test of hash value
s_ens_L = len(s_ens)
Hash_val = -1*np.ones(s_ens_L)
for ii in range(s_ens_L):
    Hash_val[ii] = HT.hashfunction(s_ens[ii])
print(Hash_val)

# Add one element
s0 = s_ens[0]
score0 = scores[0]
HT.hash_table_add(s0,score0)
s1 = s_ens[4]
score1 = scores[1]
HT.hash_table_add(s1,score1)

# Add all elements
HT = Hash_table(hash_L)
HT.hash_table_add_all(s_ens,scores)

# Find given name's score
score = HT.find_node(s1).score
print(score)

scores_r = np.zeros(s_ens_L)
for ii in range(s_ens_L):
    scores_r[ii] = HT.find_node(s_ens[ii]).score
print(scores)
print(scores_r)

# Del a node
stodel = 'Sue'
score_prev = HT.find_node(stodel).score
pos = HT.del_node(stodel)
stodel_node = HT.find_node(stodel) # Shall be -1
if stodel_node == -1:
    score_del = -1
print(score_prev,score_del)