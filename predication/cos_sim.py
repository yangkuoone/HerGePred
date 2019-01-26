# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 22:11:26 2017

@author: yk
"""

def cos(vector1, vector2):
    dot_product = 0.0
    normA = 0.0
    normB = 0.0
    for a, b in zip(vector1, vector2):
        dot_product += a*b  
        normA += a**2  
        normB += b**2  
    if normA == 0.0 or normB == 0.0:
        return None  
    else:  
        return dot_product / ((normA*normB)**0.5)  
def calSim(inFile, outFile):
    
    fr = open(inFile, 'r')
    fw = open(outFile, 'w')
    fw.truncate()
    
    count = 0
    dic = {}
    for line in fr:
        arr = line.strip().split(' ')
        arr1 = [int(i) for i in arr]
        dic[count] = arr1
        count += 1
    
    dicNum = len(dic)
    for i in range(0, dicNum):
        for j in range(i+1, dicNum):
            v1 = dic[i]
            v2 = dic[j]
            sim = cos(v1, v2)
            fw.write(str(i)+'\t'+str(j)+'\t'+str(sim)+'\n')
    fr.close()
    fw.flush()
    fw.close()
