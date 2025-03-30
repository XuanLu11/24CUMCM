'''这是第4题中用于计算调头曲线相关参数的程序,
这些参数在program4_2.py,program5.py中作为已知量被使用'''
#调用库和方法
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy.optimize import root_scalar

#阿基米德螺旋线的参数
B=1.7/2/np.pi #曲线r=B*theta的参数计算
R=9/2 #调头空间的圆的半径

def polarToRec(r:float,theta:float)->tuple:
    '''将极坐标转换成直角坐标的函数(返回值x,y)'''
    return r*np.cos(theta),r*np.sin(theta)
def lengthPolar(r1:float,r2:float,theta1:float,theta2:float)->float:
    '''计算极坐标中两点之间的欧式距离的方法(返回值d)'''
    a=r1*r1+r2*r2-2*r1*r2*math.cos(theta1-theta2)
    return math.sqrt(a)
def getCurveTangent(theta:float)->float:
    '''计算曲线上某点的切线的斜率(返回值k)'''
    #该曲线的切线过程的化简表达式
    t=np.cos(theta)-theta*np.sin(theta)
    if t==0:
        return np.nan
    return (theta*np.cos(theta)+np.sin(theta))/t

theta1=R/B
x0,y0=polarToRec(R,theta1)
k1=-1/getCurveTangent(theta1)
d=abs(k1*x0-y0)/math.sqrt(k1*k1+1)
l1=2*math.sqrt(R*R-d*d)
r1=4*R*R/3/l1
r2=r1/2
#此处计算圆的圆心位置的公式来自于对圆与螺旋线的切点和切线方向的观察,从而确定了它们的位置
x1=x0+1/math.sqrt(1+k1*k1)*r1
y1=y0+k1/math.sqrt(1+k1*k1)*r1
x2=-x0-1/math.sqrt(1+k1*k1)*r2
y2=-y0-k1/math.sqrt(1+k1*k1)*r2

print('X0='+'%.6f'%x0)
print('Y0='+'%.6f'%y0)
print('X1='+'%.6f'%x1)
print('Y1='+'%.6f'%y1)
print('R1='+'%.6f'%r1)
print('X2='+'%.6f'%x2)
print('Y2='+'%.6f'%y2)
print('R2='+'%.6f'%r2)
#两圆弧的切点
x3=(x1+2*x2)/3
y3=(y1+2*y2)/3
#计算两圆弧对应的圆心角
vector1=np.array([x1-x0,y1-y0])
vector2=np.array([x1-x3,y1-y3])
alpha1=np.arccos(np.dot(vector1,vector2)/(np.linalg.norm(vector1)\
    *np.linalg.norm(vector1)))
vector1=np.array([x2+x0,y2+y0])
vector2=np.array([x2-x3,y2-y3])
alpha2=np.arccos(np.dot(vector1,vector2)/(np.linalg.norm(vector1)\
    *np.linalg.norm(vector1)))

print('X3='+'%.6f'%x3)
print('Y3='+'%.6f'%y3)
print('ALPHA1='+'%.6f'%alpha1)
print('ALPHA2='+'%.6f'%alpha2)

#计算两圆弧的弧长
s1=alpha1*r1
s2=alpha2*r2

print('S1='+'%.6f'%s1)
print('S2='+'%.6f'%s2)

#用参数方程表示两圆弧时的参数边界取值
beta1=np.arctan((y1-y0)/(x1-x0))+np.pi
beta2=np.arctan((y3-y1)/(x3-x1))
beta3=np.arctan((y2-y3)/(x2-x3))-np.pi
beta4=np.arctan((y2+y0)/(x2+x0))

print('BETA1='+'%.6f'%beta1)
print('BETA2='+'%.6f'%beta2)
print('BETA3='+'%.6f'%beta3)
print('BETA4='+'%.6f'%beta4)
'''
计算结果为:
X0=-2.711856
Y0=-3.591078
X1=-0.760009
Y1=-1.305726
R1=3.005418
X2=1.735932
Y2=2.448402
R2=1.502709
X3=0.903952
Y3=1.197026
ALPHA1=3.021487
ALPHA2=3.021487
S1=9.080830
S2=4.540415
BETA1=4.005538
BETA2=0.984051
BETA3=-2.157542
BETA4=0.863945
'''