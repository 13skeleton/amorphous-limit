## Cahill 的非晶极限
import scipy.constants as C
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
import  argparse
## 定义常量


pi = C.pi
kB = C.k
h = C.h
A = 1E-10

prompt = ">>>"


## 截止频率
def cutoff_frequency(v,n):
    theta = np.array(v)*(h/kB)*np.power(3*n/(4*pi),1/3)
    return theta

    
def frequency_integral(cutoff_frequency,temperature):
    integrat_datas = []
    for i in range(3):
        a = np.array(cutoff_frequency/temperature)[i]
        data,err = integrate.quad(lambda x:x**(3)*np.exp(-x)/np.power((1-np.exp(-x)),2),0,a)
        integrat_datas.append(data)
    return integrat_datas
    
def amorphous_lattic_thermal(v,n,temperature):
    sound_modes = np.array(v)
    theta = cutoff_frequency(v,n)
    integrat_datas = frequency_integral(theta,temperature)
    
    klmin = np.power(pi/6,1/3)*kB*np.power(n,2/3)*np.sum(sound_modes*np.power(temperature/theta,2)*integrat_datas)
    return klmin
    

def plot_lattic_thermal(x,y):
    plt.figure(figsize=(8,6))
    plt.plot(x,y,label="amorphous limitation")
    plt.xlabel("$T \; (K)$")
    plt.ylabel("$\kappa_{L} \; (Wm^{-1}K^{-1})$")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    ## 解析
    parser = argparse.ArgumentParser()
    parser.description="Cahill非晶极限理论模型"

    parser.add_argument("-p","--plot",help="绘制非晶极限图",action="store_true")
    parser.add_argument("-o","--out",help="模拟结果保存到dat文件",action="store_true")
    args = parser.parse_args()
    
    
    sound_modes =[]
    for i in range(3):
        velocity = float(input("请输入第%d种模式的声速。(单位:m/s)."%(i+1)))
        sound_modes.append(velocity)
    volume_unit_cell = float(input("请输入晶胞体积。(单位:A^3)."))
    number_atoms_unit_cell = int(input("请输入晶胞内原子数目。(单位： 个)"))
    
    number_density_of_atoms = float(number_atoms_unit_cell/volume_unit_cell/np.power(A,3))
    temperature_min= float(input("请输入初始温度。(单位:K)"))
    temperature_max = float(input("请输入终止温度。(单位:K)"))
    temperatures = np.arange(temperature_min,temperature_max+0.1,0.1)
    lattic_thermal_amorphous_list = []
    for i in temperatures:
        lattic_thermal_amorphous = amorphous_lattic_thermal(sound_modes,number_density_of_atoms,i)
        lattic_thermal_amorphous_list.append(lattic_thermal_amorphous)
    
    ## 输入结果文件
    print("*"*30+"输入结果"+"*"*30)
    print("声速分别为：",sound_modes[0],sound_modes[1],sound_modes[2])
    print("样品原子密度为：",number_density_of_atoms)
    print("初始温度为：",temperature_min)
    print("终止温度为：",temperature_max)
    
    ## 输入结果文件保存
    with open("Cahill非晶极限-模拟参数.in","w") as f:
        f.write("声速分别为：%8.3f %8.3f %8.3f"%(sound_modes[0],sound_modes[1],sound_modes[2])+"\n")
        f.write("晶胞体积为: %8.3f"%(volume_unit_cell)+"\n")
        f.write("晶胞中原子个数为: %8.3f"%(number_atoms_unit_cell)+"\n")
        f.write("初始温度为: %8.3f"%(temperature_min)+"\n")
        f.write("终止温度为: %8.3f"%(temperature_max)+"\n")

    ### 绘图
    if args.plot:
        plot_lattic_thermal(temperatures,lattic_thermal_amorphous_list)
    if args.out:
        ## 输出结果到文件
        print("*"*30+"结果文件保存中"+"*"*30)
        results = np.array([temperatures,lattic_thermal_amorphous_list]).T
        np.savetxt("Cahill非晶极限-模拟结果.dat",results,fmt="%8.3f",delimiter=" ",header="T(K)  k_L(W/m/K)",comments="#")
        print("保存成功")
        
    a=input("按任意键退出")



