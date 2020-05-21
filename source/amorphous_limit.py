#!/usr/bin/python
# -*- coding: UTF-8 -*-
import os,yaml
import numpy as np
from scipy import integrate
import scipy.constants as C
import csv


pi = C.pi
k_B = C.k
h = C.h
hbar = h/(2*pi)
N_A = C.N_A

class Cahill_Model(object):
    def __init__(self,name):
        self.__sample_name = name
        self.__number_atoms_set = 0
        self.__number_atoms_in_cell_set = 0
        self.__volume_cell_set = 0
        self.__sound_modes_set = 0
        self.__density_set = 0
        self.__relative_atomic_mass_set = 0
        self.__average_atomic_volume_set = 0
                
        
    @property
    def sample_name(self):
        return self.__sample_name
    

    # 读写纵波声速
    @property
    def sound_velocity_l(self):
        if self.__sound_modes_set == 1:
            return self.__sound_modes[0]
        else:
            return self.__sound_velocity_l
    @sound_velocity_l.setter
    def sound_velocity_l(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("logitudinal sound velocity must be a number")
        self.__sound_velocity_l = value
    
    #读写横波声速
    @property
    def sound_velocity_t(self):
        if self.__sound_modes_set == 1:
            return (self.__sound_modes[1]+self.__sound_modes[2])/2
        else:
            return self.__sound_velocity_t
    @sound_velocity_t.setter
    def sound_velocity_t(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("transverse sound velocity must be a float or int")
        self.__sound_velocity_t = value
    
    #返回声学模式
    @property
    def sound_modes(self):
        if self.__sound_modes_set == 0:
            value=np.array([self.sound_velocity_l,
                self.sound_velocity_t,
                self.sound_velocity_t])
            return value
        else:
            return self.__sound_modes
    @sound_modes.setter
    def sound_modes(self,value):
        self.__sound_modes_set = 1
        self.__sound_modes = np.array(value)
    

    #读写物质的单胞体积
    @property
    def volume_cell(self):
        return self.__volume_cell
    @volume_cell.setter
    def volume_cell(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("volume_cell must be a float or int")
        self.__volume_cell_set = 1
        self.__volume_cell = value
    
    #读写样品的密度
    @property
    def density(self):
        return self.__density
    @density.setter
    def density(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("density must be a float or int")  
        self.__density_set = 1
        self.__density = value
    
    # 读写物质中的相对原子质量。
    @property
    def relative_atomic_mass(self):
        return self.__relative_atomic_mass
    @relative_atomic_mass.setter
    def relative_atomic_mass(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("relative atomic mass must be a number")
        self.__relative_atomic_mass_set = 1
        self.__relative_atomic_mass = value    
    
    #读写物质的温度，temperature,T均表示温度
    @property
    def temperature(self):
        return self.__temperature
    @temperature.setter
    def temperature(self,value):
#        if not (isinstance(value,float) or isinstance(value,int)):
#            raise ValueError("temperature must be a float or int")
        self.__temperature = value
    @property
    def T(self):
        return self.__temperature
    
    #读写化学式中的原子数
    @property
    def number_atoms(self):
        return self.__number_atoms
    @number_atoms.setter
    def number_atoms(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("number atoms must be a float or int")
        self.__number_atoms_set = 1
        self.__number_atoms = value  
    
    #读写单胞中原子数目
    @property
    def number_atoms_in_cell(self):
        return self.__number_atoms_in_cell
    @number_atoms_in_cell.setter
    def number_atoms_in_cell(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("number atoms in cell must be an int")
        self.__number_atoms_in_cell_set = 1
        self.__number_atoms_in_cell = value     

    #读写平均原子体积
    @property
    def average_atomic_volume(self):
        if self.__average_atomic_volume_set == 0:
            if self.__density_set == 1 and self.__relative_atomic_mass_set == 1 and self.__number_atoms_set == 1:
                return self.relative_atomic_mass /(self.number_atoms*self.density*N_A)*1E24
            elif self.__volume_cell_set == 1 and self.__number_atoms_in_cell_set == 1:
                return self.volume_cell/self.number_atoms_in_cell
        elif self.__average_atomic_volume_set == 1:
            return self.__average_atomic_volume
    @average_atomic_volume.setter
    def average_atomic_volume(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("number atoms must be a float or int")
        self.__average_atomic_volume_set = 1
        self.__average_atomic_volume = value
    
    
    #读取原子数密度，相当于平均原子体积的倒数。
    @property
    def number_denisty_of_atom(self):
        if self.__density_set == 1 and self.__relative_atomic_mass_set == 1 and self.__number_atoms_set == 1:
            value = self.number_atoms*self.density*N_A/self.relative_atomic_mass*1E6
        elif self.__volume_cell_set == 1 and self.__number_atoms_in_cell_set == 1:
            value = self.number_atoms_in_cell/self.volume_cell*1E30
        elif self.__average_atomic_volume_set == 1:
            value = 1/self.average_atomic_volume*1E30
        return value
    

        
    #读写截至频率
    @property
    def cutoff_frequency(self):
        value = self.sound_modes*(h/k_B) \
            *np.power(3*self.number_denisty_of_atom/(4*pi),1/3)
        return value
    
    @property    
    def frequency_integral(self):
        a = np.array(self.cutoff_frequency/self.temperature)
        value= [integrate.quad(lambda x:np.power(x,3)*np.exp(-x)/np.power(1-np.exp(-x),2),0,a[i])[0] for i in range(3)]
        return value
    
    @property
    def const_amorphous_lattice_thermal(self):
        value = 1/2*np.power(pi/6,1/3)*k_B  \
            *np.power(self.number_denisty_of_atom,2/3)  \
            *(2*self.sound_velocity_t+self.sound_velocity_l)
        return value
    
    @property
    def amorphous_lattice_thermal(self):

        klmin = np.power(pi/6,1/3)*k_B   \
                *np.power(self.number_denisty_of_atom,2/3) \
                *np.sum(self.sound_modes*np.power(self.temperature/self.cutoff_frequency,2)*self.frequency_integral)
        return klmin
        
def read():
    CurrentPath=os.getcwd()
    YamlFile=os.path.join(CurrentPath,"input.yaml")

    with open(YamlFile,"r") as f:
        value = yaml.load(f,Loader=yaml.FullLoader)
    return value
def calculate():
    parameter = read()
    s = Cahill_Model(parameter["Sample_Name"])
    if "Logitudinal_Sound_Velocity"  and "Transverse_Sound_Velocity" in parameter.keys():
        s.sound_velocity_l = float(parameter["Logitudinal_Sound_Velocity"])
        s.sound_velocity_t = float(parameter["Transverse_Sound_Velocity"])
    elif "Sound_Modes" in parameter.keys():
        s.sound_modes = np.array((parameter["Sound_Modes"]).split(","),float)
    
    if "Sample_Density" and "Relative_Atomic_Mass" and "Number_Atoms" in parameter.keys():
        s.density = float(parameter["Sample_Density"])
        s.relative_atomic_mass = float(parameter["Relative_Atomic_Mass"])
        s.number_atoms = float(parameter["Number_Atoms"])
    elif "Volume_Cell" and "Number_Atoms_in_Cell" in parameter.keys():
        s.volume_cell = float(parameter["Volume_Cell"])
        s.number_atoms_in_cell = float(parameter["Number_Atoms_in_Cell"])
    elif "Average_Atomic_Volume" in parameter.keys():
        s.average_atomic_volume = float(parameter["Average_Atomic_Volume"])
    
    start_temperature = float(parameter["Temperature"]["Start_Temperature"])
    end_temperature = float(parameter["Temperature"]["End_Temperature"])
    interval_temperature = float(parameter["Temperature"]["Interval_Temperature"])
    
    print("-"*36+"基本信息"+"-"*36)
    print(" "*4+"本程序由13skeleton编写,如有任何问题，请直接联系邮箱。(J.Pei@foxmail.com)")
    print("""    参考文献：
    1. Cahill, D. G., Watson, S. K. & Pohl, R. O. Lower limit to the thermal con
       ductivity of disordered crystals. Phys. Rev. B 46, 6131–6140 (1992).
    2. Huang, B.-L. & Kaviany, M. Ab initio and molecular dynamics predictions f
       or electron and phonon transport in bismuth telluride. Phys. Rev. B 77, 
       125209(2008).
    3. He, Y. et al. High Thermoelectric Performance in Non-Toxic Earth-Abundant 
       Copper Sulfide. Adv. Mater. 26, 3974–3978 (2014).
    4. Dong, J. et al. Reducing Lattice Thermal Conductivity of MnTe by Se Alloying 
       toward High Thermoelectric Performance. ACS Appl. Mater. Interfaces 11, 28221–
       28227 (2019).

      """)
    
    print("-"*36+"输入参数"+"-"*36)
    print("样品名称",s.sample_name)
    if "Logitudinal_Sound_Velocity"  and "Transverse_Sound_Velocity" in parameter.keys():
        print("纵波声速",s.sound_velocity_l)
        print("横波声速",s.sound_velocity_t)
    elif "Sound_Modes" in parameter.keys():
        print("声速模式",s.sound_modes)
        
    if "Sample_Density" and "Relative_Atomic_Mass" and "Number_Atoms" in parameter.keys():
        print("样品密度",s.density)
        print("相对原子质量",s.relative_atomic_mass)
        print("化学式中原子数目",s.number_atoms)
    elif "Volume_Cell" and "Number_Atoms_in_Cell" in parameter.keys():
        print("单胞体积",s.volume_cell)
        print("单胞中原子数目",s.number_atoms_in_cell)
    elif "Average_Atomic_Volume" in parameter.keys():
        print("平均原子体积",s.average_atomic_volume)
        
    print(" ")
    print(" ")
    print(" ")
    
    results_temperature = []
    for i in np.arange(start_temperature,end_temperature,interval_temperature):
        s.temperature = i
        list_for_temperature = (format(s.temperature,".1f"),format(s.amorphous_lattice_thermal,".6f"))
        results_temperature.append(list_for_temperature)
    
    with open("out.csv","w",encoding="utf-8",newline="") as csvfile:
        myinput = csv.writer(csvfile)
        myinput.writerow(["#输入参数"])
        myinput.writerow(["样品名称",s.sample_name])
        if "Logitudinal_Sound_Velocity"  and "Transverse_Sound_Velocity" in parameter.keys():
            myinput.writerow(["纵波声速",s.sound_velocity_l])
            myinput.writerow(["横波声速",s.sound_velocity_t])
        elif "Sound_Modes" in parameter.keys():
            myinput.writerow(["声速模式",s.sound_modes[0],s.sound_modes[1],s.sound_modes[2]])

        if "Sample_Density" and "Relative_Atomic_Mass" and "Number_Atoms" in parameter.keys():
            myinput.writerow(["样品密度",s.density])
            myinput.writerow(["相对原子质量",s.relative_atomic_mass])
            myinput.writerow(["化学式中原子数目",s.number_atoms])
        elif "Volume_Cell" and "Number_Atoms_in_Cell" in parameter.keys():
            myinput.writerow(["单胞体积",s.volume_cell])
            myinput.writerow(["单胞中原子数目",s.number_atoms_in_cell])
        elif "Average_Atomic_Volume" in parameter.keys():
            myinput.writerow(["平均原子体积",s.average_atomic_volume])        
        myinput.writerow([" "," "])
        
    with open("out.csv","a",encoding="utf-8",newline="") as csvfile:
        myoutput = csv.writer(csvfile)
        myoutput.writerow(["#输出结果"])
        if "Sample_Density" and "Relative_Atomic_Mass" and "Number_Atoms" in parameter.keys():
            myoutput.writerow(["平均原子体积(A^3)",s.average_atomic_volume])
        elif "Volume_Cell" and "Number_Atoms_in_Cell" in parameter.keys():
            myoutput.writerow(["平均原子体积(A^3)",s.average_atomic_volume])
            
        myoutput.writerow(["原子数密度(A^-3)",s.number_denisty_of_atom/1E30])
        myoutput.writerow(["非晶极限恒定值版(W/m/K)",s.const_amorphous_lattice_thermal])
        myoutput.writerow(["温度","非晶极限(W/m/K)"])
        myoutput.writerows(results_temperature)
    
    print("-"*36+"#输出结果"+"-"*36)
    if "Sample_Density" and "Relative_Atomic_Mass" and "Number_Atoms" in parameter.keys():
        print("平均原子体积(A^3)",s.average_atomic_volume)
    elif "Volume_Cell" and "Number_Atoms_in_Cell" in parameter.keys():
        print("平均原子体积(A^3)",s.average_atomic_volume)
        
    print("原子数密度(A^-3)",format(s.number_denisty_of_atom/1E30,".3f"))
    
    print("非晶极限恒定值版(W/m/K)",format(s.const_amorphous_lattice_thermal,".3f"))
    print("温度","非晶极限(W/m/K)")
    for i in range(3):
        print(results_temperature[i][0],results_temperature[i][1])
    
    print("...")
    print("...")
    print("...")
    
    print("计算完成，输出结果请查看out.csv文件")
    
    
    
if __name__ == "__main__":
    calculate()
    a=input("按任意键退出")

