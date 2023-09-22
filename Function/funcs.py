#-*- coding:utf-8 -*-
#2022.11.30    反演需要的函数
import numpy as np
import time
import os
import shutil
from .struct import Object
from odbAccess import *
import matlab.engine

def Position(odbfile,length,num,Init_p,ELASTIC,P):   #对隧道全长进行定位分区
    ELASTIC=ELASTIC.reshape(-1)
    o=openOdb(odbfile)
    nodes=o.rootAssembly.instances['SOIL-1'].nodes
    node_label=[]
    for n in nodes:
        node_label.append(n.label)
    node_arg=np.argsort(node_label)
    elements=o.rootAssembly.instances['SOIL-1'].elements
    Center=np.zeros((len(elements),3))
    d=length/num
    for E in elements:
        N=E.connectivity
        Y=[]
        for n in N:
            y=nodes[node_arg[n-1]].coordinates[1]
            Y.append(y)
        y=np.mean(Y)
        l=int((y-Init_p)//d)+1
        c=np.array([y,l,E.label]).reshape(1,-1)
        if Center.shape==(0,):
            Center=c
        else:
            Center=np.concatenate((Center,c),axis=0)
    for i in range(1,num+1):
        eles=Center[:,-1][Center[:,l]==i]
        name='rock{:d}'.format(i)
        P.add_elset(name, eles)
        elastic=ELASTIC[i-1]
        stress=[5e5,1e6]
        strain=[0,5e-4]
        P.add_material(name, elastic, stress, strain)
        P.add_section(name, name, name)
    o.close()
    return P
   
def modify(file1,file2,P):
    f1=open(file1,'r')
    f2=open(file2,'w')
    A=f1.readlines()
    for (i,a) in enumerate(A):
        if '*Nset, nset=soil' in a:
            k0=i
        if '*End Part'  in a:
            k1=i
        if '** MATERIALS' in a:
            k2=i
        if '** PREDEFINED FIELDS'in a:
            k3=i 
    for i in range(k0):
        f2.write(A[i])
    for k in P.elsets.keys():
        a='*Elset, elset={:s}\n'.format(k)
        f2.write(a)
        elements=P.elsets[k]
        l=np.arange(0,elements.shape[0],16)
        l=list(l)
        l.append(elements.shape[0])
        for i in range(len(l)-1):
            a=','.join('%6d'%m for m in elements[l[i]:l[i+1]])+'\n'
            f2.write(a) 
    a='*Nset, nset=wole, generate\n      1,  133665,       1\n'
    a=a+'*Elset, elset=wole, generate\n      1,  124800,       1\n'
    f2.write(a)
    for k in P.sections.keys():
        a='** Section: {:s}\n'.format(k)
        a=a+'*Solid Section, elset={:s}, material={:s}\n,\n'.format(P.sections[k].el_name,P.sections[k].m_name)
        f2.write(a)
    for i in range(k1,k2+2):
        f2.write(A[i])
    for k in P.materials.keys():
        a='*Material, name={:s}\n*Density\n2500.,\n*Elastic\n'.format(k)
        a=a+'{:8e}, 0.25\n*Mohr Coulomb\n45., 0.1\n*Mohr Coulomb Hardening\n'.format(P.materials[k].elastic)
        for i in range(len(P.materials[k].stress)):
            a=a+'{:8e},{:8e}\n'.format(P.materials[k].stress[i],P.materials[k].strain[i])
        f2.write(a)
    for i in range(k3,len(A)):
        f2.write(A[i])
    f1.close()
    f2.close()
    
def finished(path):
    filename=path+'/HM.lck'
    while 1:
        if os.path.exists(filename):
            time.sleep(1)
        else:
            break
    time.sleep(2)

def read_out(path,Monitor):
    file=path+'\HM.odb'
    o=openOdb(file)
    vall=o.steps['kw'].frames[-1].fieldOutputs['U'].values
    vall_label=[]
    U={}
    for v in vall:
        vall_label.append(v.nodeLabel)
    for k in Monitor.keys():
        Id=int(k)
        u=vall[vall_label.index(Id)].data[2]    #读取竖向位移
        U[k]=u
    o.close()
    return U

def get_Q(Length,num,Monitor,d1,d2,lamda):  #d1,d2表示距离和模拟位移影响的系数,lamda表示位置和位移的影响系数
    Q1=np.zeros((num,num))
    Q2=np.zeros((num,num))
    sqr=Monitor.keys()
    for i in range(num):
        for j in range(num):
            Q1[i,j]=np.exp(-abs(i-j)*d1)
    U=np.zeros((num,))
    Distance=Length/num   #每一块的纵向长度
    for i in range(num):
        D=[]   #记录各个距离
        for j in range(len(sqr)):
            p=(i+0.5)*Distance
            d=abs(p-Monitor[sqr[j]][0])
            if d==0:
                D.append(0)
                continue
            else:
                D.append(np.exp(-d*0.5))
        if 0 in D:
            k=D.index(0)
            U[i]=Monitor[sqr[k]][1]
        else:
            u=0
            for k in range(len(sqr)):
                u=u+D[k]*Monitor[sqr[k]][1]/sum(D)
            U[i]=u
    for i in range(num):
        for j in range(num):
            Q2[i,j]=np.exp(-abs(U[i]-U[j])*d2)
    Q=lamda*Q1+(1-lamda)*Q2
    return Q

def Chang_Material(P,Elastic,num):
    Elastic=Elastic.reshape(-1)
    for i in range(num):
        name='rock{:d}'.format(i+1)
        elastic=Elastic[i]
        P.change_material(name,elastic)  
        
def Plot(ELASTIC,Length,Iter):
    eng=matlab.engine.start_matlab()
    Elastic=ELASTIC.reshape(-1)
    f=open('Elastic.dat','w')
    for i in range(Elastic.shape[0]):
        a='{:8e}\n'.format(Elastic[i])
        f.write(a)
    f.close()
    filename='Elastic_{:d}.dat'.format(Iter)
    eng.Plot_material(Length,Iter)
    shutil.copyfile('Elastic.dat', filename)