#-*- coding:utf-8 -*-
# 根据位置判断节点编号与竖向位移
from odbAccess import *
import numpy as np
import os

file='HM.odb'
o=openOdb(file)
nodes=o.rootAssembly.instances['SOIL-1'].nodes
#确定坐标
loc=[20,70,110,170]
Node={}   #记录节点号
f=open('Object2.txt','w')
for n in nodes:
    coor=n.coordinates
    if abs(coor[2]-24.0)>=0.01:   ##纵坐标轴线
        continue
    for i in range(len(loc)):
        if abs(coor[1]-loc[i])<=0.01:
            Node[str(loc[i])]=n.label
vall=o.steps['kw'].frames[-1].fieldOutputs['U'].values
vall_label=[]
for v in vall:
    vall_label.append(v.nodeLabel)
for l in Node.keys():
    u=vall[vall_label.index(Node[l])].data[2]    #竖向位移
    a='Node{:d}:loc={:d},U={:f}\n'.format(Node[l],int(l),u)
    f.write(a)
o.close()
f.close()