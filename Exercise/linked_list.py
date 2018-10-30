## Python linked list
# Author: Yang Yu
# Date: 2018/10/30
# Description: linked list

class node:
    
    data = [] # Date within the node
    nextptr = -1 # Pointer
    
    def __init__(self,data=[]):
        self.data = data
        
    def copy(self):
        node_new = node(self.data)
        node_new.nextptr = self.nextptr
    
    def data_print(self):
        print(data)
    
class linked_list:
    
    node0 = node() # First node
    nodeNum = 0 # Linked list length
    
    def __init__(self,node0):
        self.node0 = node0
        # Check initilized list lengh
        self.nodeNum = self.list_length()
    
    def copy(self):
        ll_new = linked_list(self.node0)
    
    def list_length(self):
        cnt = 1
        nodetmp = self.node0
        while nodetmp.nextptr != -1: # Not the final node
            cnt=cnt+1
            nodetmp = nodetmp.nextptr
        
        return cnt
    
    def gotonode(self,idx): # Got to idx node
        if idx+1 > self.nodeNum: # Longer than the list
            print("Error: node index off-boundary")
            return -1
        else:
            nodetmp = self.node0
            for ii in range(idx):
                nodetmp = nodetmp.nextptr
            
            return nodetmp
            
    def append(self,nodeE):
        nodeZ = self.gotonode(self.nodeNum-1) # Go to the final node 
        nodeZ.nextptr = nodeE # Point to the nodeE
        self.nodeNum = self.nodeNum+1
    
    def insert(self,nodeI,insertidx):
        nodeI_1 = self.gotonode(insertidx) # Go to insert idx
        nodeI1 = self.gotonode(insertidx+1) # Go to insert idx+1
        
        nodeI_1.nextptr = nodeI
        nodeI.nextptr = nodeI1 
        self.nodeNum = self.nodeNum+1
        
    def delete(self,deleteidx=-1):
        
        if deleteidx == -1:
            deleteidx = self.nodeNum-1
        
        if deleteidx == 0: # First
            nodeD = self.gotonode(deleteidx)
            nodeD1 = self.gotonode(deleteidx+1)
            
            nodeD.nextptr = -1     
            self.node0 = nodeD1
            
        elif deleteidx == self.nodeNum-1: # Last
            nodeD_1 = self.gotonode(deleteidx-1)
            
            nodeD_1.nextptr = -1
            
        else:
            nodeD_1 = self.gotonode(deleteidx-1)
            nodeD = self.gotonode(deleteidx)
            nodeD1 = self.gotonode(deleteidx+1)
            
            nodeD_1.nextptr = nodeD1
            nodeD.nextptr = -1
            
        self.nodeNum = self.nodeNum-1
    
    def print1(self,idx=0):
        
        if idx < 0:
            idx = 0
        elif idx > self.nodeNum:
            idx = self.nodeNum-1
        
        nodeP = self.gotonode(idx)
        print(nodeP.data)
    
    def printA(self):
        nodetmp = self.node0
        print(nodetmp.data)
        while nodetmp.nextptr != -1: # Not the final node
            nodetmp = nodetmp.nextptr
            print(nodetmp.data)

## Routines
# Construct linked list
ll_length = 10

node0 = node(0)
node_prev = node0
for ii in range(ll_length-1):
    node_new = node(ii+1)
    node_prev.nextptr = node_new
    node_prev = node_new

ll = linked_list(node0)
ll.print1(5)
ll.printA()

node_new = node(100)
ll.append(node_new)

node_insert = node(1000)
newidx = 7
ll.insert(node_insert,newidx)
ll.printA()

ll.delete(8)
ll.printA()

ll.delete(-1)
ll.delete(0)
ll.printA()