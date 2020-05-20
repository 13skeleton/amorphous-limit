[TOC]

# Cahill模型求解非晶极限

## Cahill模型原始版(与温度相关)


$$
\kappa_{l,min}=\left(\frac{\pi}{6} \right)^{1/3}k_{B}\cdot n^{2/3}\sum_{i} v_{i}\left(\frac{T}{\theta_{i}}\right)^{2}\int_{0}^{\theta_{i}/T}\frac{x^{3}e^{x}}{(e^x-1)^{2}}dx
$$

$$
\theta_{i}=v_{i}(h/k_{B})(3 n/4 \pi)^{1/3}
$$

$$
n=\frac{1}{V_{0}}
$$
或
$$
n=\frac{\rho \cdot N_{A}}{m_{0}}
$$
其中,$k_{B}$表示玻尔兹曼常数; $n$表示原子数密度; $v_{i}$表示声频支声速的三种模式，一种纵波声速，两种横波声速; $\theta_{i}$表示三种声速模式对应的截止频率; $T$表示温度; $x$表示振动频率; $V_{0}$表示单胞体积; $m_{0}$表示相对分子质量;  $\rho$表示块体样品的密度;  $N_{A}$表示阿伏加德罗常数。

### 参考文献:

1. Cahill, D. G., Watson, S. K. & Pohl, R. O. Lower limit to the thermal conductivity of disordered crystals. *Phys. Rev. B* **46**, 6131–6140 (1992).
2. May, A. F. & Snyder, G. J. Introduction to Modeling Thermoelectric Transport at High Temperatures. in *Materials, Preparation, and Characterization in Thermoelectrics* (ed. Rowe, D. M.) 207–224 (CRC Press, 2017). doi:[10.1201/b11891-11](https://doi.org/10.1201/b11891-11).
3. Dong, J. *et al.* Reducing Lattice Thermal Conductivity of MnTe by Se Alloying toward High Thermoelectric Performance. *ACS Appl. Mater. Interfaces* **11**, 28221–28227 (2019).

> 董金峰论文中的公式有错误，特别注意。

## Cahill模型简化版(与温度无关)

$$
\kappa_{l,min} = \frac{1}{2}\left(\frac{\pi}{6}\right)^{1/3}k_{B}V_{0}^{-2/3}(2v_{t}+v_{l})
$$

其中，$v_{t}$表示横波声速，$v_{l}$表示纵波声速

### 参考文献:

1. He, Y. *et al.* High Thermoelectric Performance in Non-Toxic Earth-Abundant Copper Sulfide. *Adv. Mater.* **26**, 3974–3978 (2014).
2. Tan, G. *et al.* Extraordinary role of Hg in enhancing the thermoelectric performance of p-type SnTe. *Energy Environ. Sci.* **8**, 267–277 (2015).

