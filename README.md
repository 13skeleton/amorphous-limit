# amorphous-limit

## 程序简介

本程序采用python3编写，具体用到的第三方库有numpy, scipy, pyyaml。如有问题，及时联系J.Pei(J.Pei@foxmail.com)。

本程序源代码托管在github上面，如需要查看最新版本程序，请移步至: [https://github.com/13skeleton/amorphous-limit](https://github.com/13skeleton/amorphous-limit)

## 程序下载方法

**Windows:**

1. 打开 [https://github.com/13skeleton/amorphous-limit](https://github.com/13skeleton/amorphous-limit) 链接，点击"clone or download"按钮，将zip文件下载至本地

2. 解压缩.zip文件，打开bin文件夹，双击运行".exe"文件。

**linux:**

> 自己解决吧，不会的，联系我。

## 程序使用方法

1. 在程序运行目录下准备一个“input.yaml”文件。

具体格式如下:

```yaml
#基本设置
Sample_Name: Bi0.3Sb1.7Te3-cross-plane #样品名称，无实际意义，仅为区分
#声频支声速总共有三种模式，分为1个纵波声速和2个横波声速。常规声速测试测得的横波声速为平均值纵波声速和横波声速 (默认)
Logitudinal_Sound_Velocity: 2372 #纵波声速 unit: m/s
Transverse_Sound_Velocity: 1483 #横波声速 unit: m/s

#声频支声速总共有三种模式，分为1个纵波声速和2个横波声速。常规声速测试测得的横波声速为平均值，如可以知道两种横波各自的值的话可以设置如下结果。(可选)
#Sound_Modes: 7160.305,4250.755,4250.755   #分别代表纵波声速，横波声速1,横波声速2.


#求取原子数密度(平均原子体积的倒数)

#对于求原子数密度可以通过密度和相对原子质量求得，这种方案是推荐的解决思路。(默认)
Sample_Density: 4.85   #样品的密度 unit: g/cm^3
Relative_Atomic_Mass: 652.469 #化学式中的相对原子质量 g/mol
Number_Atoms: 5  #化学式中的原子数目

#如果确实不知道密度和相对原子质量，如果知道单胞体积，也可求得原子数密度。(可选)
#Volume_Cell: 169.23  #单胞体积 unit: A^3
#Number_Atoms_in_Cell: 5  #单胞中原子数目

#如果知道平均原子体积，后续的样品密度，化学式中的相对原子质量，单胞体积等都可设置为注释。
#Average_Atomic_Volume: 33.846  #单位: A^3


Temperature:
    Start_Temperature: 300  #起始温度 单位: K
    End_Temperature: 800    #终止温度 单位: K
    Interval_Temperature: 1 # 温度间隔 单位: K

```

> 该输入文件遵循yaml的书写规范。可自行调整。

2. 运行程序

```bash
python amorphous_limit.py
```

