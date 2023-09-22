#-*- coding:utf-8 -*-
#2022.11.29隧道受压反演算例
import numpy as np
import time
import os
from abaqus import *
from abaqusConstants import *
from odbAccess import *
from Function.funcs import *
import shutil
from Function.struct import Object
import copy

##运行目标模型
P=Object()
file1='HM.inp'
f=open(file1,'r')
A=f.readlines()
for (i,a) in enumerate(A):
    if '*Elset, elset=soil' in a:
        k0=i
    if '*Nset, nset=r1' in a:
        k1=i
soils=[]
for i in range(k0+1,k1):
    a=A[i][:-1].split(',')
    for am in a:
        soils.append(int(am))
f.close()
pwd=os.getcwd()
Inv_num=100
Length=200.0
ELASTIC=np.zeros((Inv_num,1))
#划分目标强度分布

a=np.log(5)/2500
for i in range(ELASTIC.shape[0]):
    #ela=5+0.2*i
    ela=25*np.exp(-a*(i-49.5)**2)
    ELASTIC[i][0]=1e9*ela
    
'''
d=int(Inv_num/4)
D=[0,d,d*2,d*3,d*4]    ##分为四段
for i in range(len(D)-1):
    E=5e9*(i+1)
    for j in range(D[i],D[i+1]):
        ELASTIC[j][0]=E
 '''       
        
        
        
        
f=open('Object.txt','w')
for i in range(ELASTIC.shape[0]):
    f.write('{:8e}\n'.format(ELASTIC[i][0]))
f.close()
Ele_ori=ELASTIC.copy()
#
P.add_elset('soil',np.array(soils))
P=Position('HM.odb',Length,Inv_num,0,ELASTIC,P)
#计算
os.mkdir('Origin')
file2='Origin/HM.inp'
modify(file1,file2,P)
os.chdir('Origin')
mdb.JobFromInputFile(name='HM',inputFileName='HM.inp',parallelizationMethodExplicit=DOMAIN,
                     numDomains=8,multiprocessingMode=DEFAULT, numCpus=8)
mdb.jobs['HM'].submit()
os.chdir(pwd)
time.sleep(10)
finished('Origin')
time.sleep(2)
##确定目标值
Monitor={}
Label=[]
L=[]
file='Origin/HM.odb'

o=openOdb('Origin/HM.odb')
nodes=o.rootAssembly.instances['SOIL-1'].nodes
loc=[12,60,100,140,188]
Node={}   #记录节点号
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
    Monitor[str(Node[l])]=[int(l),u]
    Label.append(str(Node[l]))
    L.append(l)
o.close()
Label=np.array(Label)
Label=Label[np.argsort(L)]

## 运行初始正模型
if os.path.exists('Krun'):
    shutil.rmtree('Krun')
os.mkdir('Krun')
ELASTIC=10e9*np.ones((Inv_num,1))
file2='Krun/HM.inp'
Chang_Material(P,ELASTIC,Inv_num)
modify(file1,file2,P)
os.chdir('Krun')
mdb.JobFromInputFile(name='HM',inputFileName='HM.inp',parallelizationMethodExplicit=DOMAIN,
                     numDomains=8,multiprocessingMode=DEFAULT, numCpus=8)
mdb.jobs['HM'].submit()
os.chdir(pwd)
time.sleep(10)
finished('Krun')
U0=read_out('Krun',Monitor)
Fun0=[]
F_tar=[]
for k in Label:
    Fun0.append(U0[k])
    F_tar.append(Monitor[k][1])
Fun0=np.array(Fun0).reshape(-1,1)
F_tar=np.array(F_tar).reshape(-1,1)
Error={}
a='Iter0:\n'
for k in U0.keys():
    Error[k]=(U0[k]-Monitor[k][1])/Monitor[k][1]
    a=a+'point{:s}={:6f}  '.format(k,Error[k])
f=open('Error.txt','w')
f.write(a+'\n')
f.close()
Q=get_Q(Length,Inv_num,Monitor,0.0750,50.0,0.30)
[U,S,V]=np.linalg.svd(Q)
if os.path.exists('plot'):
    shutil.rmtree('plot') 
os.mkdir('plot')
Plot(ELASTIC,Length,0)
Rate=[]
r=sum(abs(Ele_ori-ELASTIC)/Ele_ori)/ELASTIC.shape[0]
Rate.append(r)


##反演循环
#在每一个主成分，反演参数加入相应的分量
Princ_num=5 #主成分分析的个数
Iter_num=8 #循环次数
lamda=2e8
Dm=3.0e9
for Iter in range(Iter_num):
    if os.path.exists('running'):
        shutil.rmtree('running')
    os.mkdir('running')
    for i in range(Princ_num):
        path='running/run{:d}'.format(i+1)
        os.mkdir(path)
        u=U[:,i].reshape(-1,1)
        Elastic=ELASTIC.copy()
        k=lamda/np.max(abs(u))
        Elastic=Elastic+k*u
        Q=copy.deepcopy(P)
        Chang_Material(Q,Elastic,Inv_num)
        file1='HM.inp'
        file2=path+'/HM.inp'
        modify(file1,file2,Q)
        os.chdir(path)
        mdb.JobFromInputFile(name='HM',inputFileName='HM.inp',parallelizationMethodExplicit=DOMAIN,
                     numDomains=8,multiprocessingMode=DEFAULT, numCpus=8)
        mdb.jobs['HM'].submit() 
        os.chdir(pwd)
        del mdb.jobs['HM']
        del Q
    time.sleep(10)
    ZETA=np.zeros((Label.shape[0],Inv_num))   ##储存变量ZETA
    for i in range(Princ_num):
        path='running/run{:d}'.format(i+1)
        finished(path)
        u=read_out(path,Monitor)
        Fun1=[]
        for k in Label:
            Fun1.append(u[k])
        Fun1=np.array(Fun1).reshape(-1,1)
        zeta=(Fun1-Fun0)/lamda
        ZETA[:,i]=zeta.reshape(-1) 
    #参数更新
    HQ=np.zeros((Label.shape[0],Inv_num))
    HQHT=np.zeros((Label.shape[0],Label.shape[0]))
    for i in range(Princ_num):
        HQ=HQ+np.dot(ZETA[:,i].reshape(-1,1),U[:,i].reshape(1,-1))
        HQHT=HQHT+np.dot(ZETA[:,i].reshape(-1,1),ZETA[:,i].reshape(1,-1))
    cst=0.1*np.max(abs(np.diag(HQHT)))
    HQHT=HQHT+cst*np.diag(np.diag(HQHT))
    HQHT=np.linalg.inv(HQHT)
    Epcy=np.dot(HQHT,HQ).T
    Delta=np.dot(Epcy,Fun0-F_tar)
    Dm=Dm*0.95
    k=Dm/np.max(abs(Delta))
    ELASTIC=ELASTIC-k*Delta
    Chang_Material(P,ELASTIC,Inv_num)
    r=sum(abs(Ele_ori-ELASTIC)/Ele_ori)/ELASTIC.shape[0]
    Rate.append(r)
    #更新后的计算
    if os.path.exists('Krun'):
        shutil.rmtree('Krun')
    os.mkdir('Krun')
    file1='HM.inp'
    file2='Krun/HM.inp'
    modify(file1,file2,P)
    os.chdir('Krun')
    mdb.JobFromInputFile(name='HM',inputFileName='HM.inp',parallelizationMethodExplicit=DOMAIN,
                         numDomains=8,multiprocessingMode=DEFAULT, numCpus=8)
    mdb.jobs['HM'].submit()
    os.chdir(pwd)
    time.sleep(10)
    finished('Krun')
    U0=read_out('Krun',Monitor)
    Fun0=[]
    for k in Label:
        Fun0.append(U0[k])
    Fun0=np.array(Fun0).reshape(-1,1)
    a='Iter{:d}:\n'.format(Iter+1)
    for k in U0.keys():
        Error[k]=(U0[k]-Monitor[k][1])/Monitor[k][1]
        a=a+'point{:s}={:6f}  '.format(k,Error[k])
    f=open('Error.txt','a')
    f.write(a+'\n')
    f.close()
    Plot(ELASTIC,Length,Iter+1)
f=open('Rate.txt','w')
for i in range(len(Rate)):
    f.write('Rate={:6e}\n'.format(Rate[i]))
f.close()            
 