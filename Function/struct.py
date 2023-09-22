#-*- coding:utf-8 -*-
##2022.11.30  创建结构体
class Material():
    def __init__(self,elastic,stress,strain):   #弹性模量，应力，应变 
        self.elastic=elastic
        self.stress=stress
        self.strain=strain
class Section():
    def __init__(self,el_name,m_name):
        self.m_name=m_name
        self.el_name=el_name
class Object():
    def __init__(self):
        self.elsets={}
        self.materials={}
        self.sections={}
    def add_elset(self,name,elements):
        self.elsets[name]=elements
    def add_material(self,name,elastic,stress,strain):
        self.materials[name]=Material(elastic, stress, strain)
    def change_material(self,name,elastic):    #改变材料的弹性模量
        self.materials[name].elastic=elastic
    def add_section(self,name,el_name,m_name,):
        self.sections[name]=Section(el_name,m_name)
    
    