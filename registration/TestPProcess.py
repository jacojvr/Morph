# -*- coding: utf-8 -*-
import time
import numpy as np
import pprocess

t=time.time()
MM = np.array([[1],[2],[3],[4],[5]])


def calculate(t):
  i,j=t
  #time.sleep(3)
  return np.array([sum(MM),i*N+j])
     
for inc in range(1,3):
    global N
    N=2*inc

    sequence = []
    for i in range(0,N):
      for j in range(0,N):
	sequence.append([i,j])
    #sequence = np.array(sequence)

    results = pprocess.pmap(calculate,sequence,limit=1)
    results = np.array([results])
    for i in range(0,N):
      for result in results[i*N:i*N+N]:
	print result,
    print
  
print "Time taken:", time.time()-t