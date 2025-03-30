'''用于计算第1题的程序'''
#调用库和方法
import numpy as np
import math
import pandas as pd
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['KaiTi']
plt.rcParams['axes.unicode_minus'] = False
#曲线的参数
B=.55/2/np.pi #曲线r=B*theta的参数计算
#求根的数值解方法
ROOT_METHOD='brentq'
#龙头速度
HEAD_VECOCITY=1
#龙的孔的间距
LENGTH_LONG=(341-2*27.5)/100
LENGTH_SHORT=(220-2*27.5)/100

def polarToRec(r:float,theta:float)->float:
    '''将极坐标转换成直角坐标的函数(返回值x,y)'''
    return r*np.cos(theta),r*np.sin(theta)
def lengthPolar(r1:float,r2:float,theta1:float,theta2:float)->float:
    '''计算极坐标中两点之间的欧式距离的方法(返回值d)'''
    a=r1*r1+r2*r2-2*r1*r2*math.cos(theta1-theta2)
    return math.sqrt(a)
def lengthToCenter(theta:float,b:float)->float:
    '''求出阿基米德螺旋线某点到其中心的弧长(返回值s)'''
    return b/2*(theta*math.sqrt(1+theta**2)+math.log(theta+math.sqrt(1+theta**2)))
def headDthetaByTime(theta:float,b:float,length:float)->float:
    '''求某一点沿曲线往中心移动了距离length后所处的位置(返回值theta)'''
    dtheta=0.1 #小步长长度
    f=lambda x:lengthToCenter(theta,b)-lengthToCenter(x,b)-length 
    #需要求根的表达式
    theta1=theta
    while(theta1>=0):#小步长寻找零点存在的区间
        theta1-=dtheta
        s=lengthToCenter(theta,b)-lengthToCenter(theta1,b)
        if s>=length:
            root=root_scalar(f,bracket=[theta1,theta],method=ROOT_METHOD) 
            #调用区间内求根的数值解方法
            return root.root
    return np.nan
def curveCalculateTheta(theta:float,length:float,b:float)->float:
    '''求曲线上离某点欧氏距离恰好为length的点(此点离中心更远)的位置(返回值theta)'''
    dtheta=0.1 #小步长长度
    r=b*theta
    f=lambda x:lengthPolar(r,b*x,theta,x)-length 
    #需要求根的表达式
    theta1=theta
    while(True):#小步长寻找零点存在的区间
        theta1+=dtheta
        r1=b*theta1
        if lengthPolar(r,r1,theta,theta1)>=length:
            root=root_scalar(f,bracket=[theta,theta1],method=ROOT_METHOD) 
            #调用区间内求根的数值解方法
            return root.root
def getLineTangent(theta1:float,theta2:float)->float:
    '''计算曲线上两点的割线的斜率(返回值k)'''
    #该曲线的割线过程的化简表达式
    t=theta1*np.cos(theta1)-theta2*np.cos(theta2)
    if t==0:
        return np.nan
    return (theta1*np.sin(theta1)-theta2*np.sin(theta2))/t
def getCurveTangent(theta:float)->float:
    '''计算曲线上某点的切线的斜率(返回值k)'''
    #该曲线的切线过程的化简表达式
    t=np.cos(theta)-theta*np.sin(theta)
    if t==0:
        return np.nan
    return (theta*np.cos(theta)+np.sin(theta))/t
def getNextVelocity(theta1:float,theta2:float,v1:float)->float:
    '''根据上一个点的速率计算紧邻着的下一个点的速率'''
    k0=getLineTangent(theta1,theta2) #割线斜率
    k1=getCurveTangent(theta1) #切线1斜率
    k2=getCurveTangent(theta2) #切线2斜率
    #cos=1/sqrt(1+tan^2)
    t1=1+k0*k1
    t2=1+k0*k2
    if t1==0:
        c1=0
    else:
        t1=(k0-k1)/t1
        c1=1/np.sqrt(1+t1**2)
    if t2==0:
        c2=0
    else:
        t2=(k0-k2)/t2
        c2=1/np.sqrt(1+t2**2)
    if c2==0:
        return np.nan
    #上述为计算两个夹角的余弦值的过程
    return v1*c1/c2 #计算下一个点的速度
theta=np.zeros((224,301))
for i in range(301): #计算不同位置在不同时间下的位置(theta)
    print('theta turn:'+str(i))
    if i>0:
        theta[0,i]=headDthetaByTime(theta[0,i-1],B,HEAD_VECOCITY)
    else:
        theta[0,i]=16*2*np.pi #龙头的起始点的位置
    theta[1,i]=curveCalculateTheta(theta[0,i],LENGTH_LONG,B)
    for j in range(2,224):
        theta[j,i]=curveCalculateTheta(theta[j-1,i],LENGTH_SHORT,B)
v=np.zeros((224,301))
for i in range(301): #计算不同位置在不同时间下的速度(velocity)
    print('velocity turn:'+str(i))
    v[0,i]=HEAD_VECOCITY #龙头的速度始终是确定的
    for j in range(1,224):
        v[j,i]=getNextVelocity(theta[j-1,i],theta[j,i],v[j-1,i])
r=B*theta #计算不同位置在不同时间下的位置(r)
x=np.zeros((224,301))
y=np.zeros((224,301))
for i in range(224): #计算不同位置在不同时间下的位置(x,y)
    for j in range(301):
        x[i,j],y[i,j]=polarToRec(r[i,j],theta[i,j])
result=np.zeros((2*224,301))
for i in range(224): #排列成可以输出的格式
    result[2*i]=x[i]
    result[2*i+1]=y[i]

result=np.round(result,6) #保留6位小数
v=np.round(v,6) #保留6位小数
data1=pd.DataFrame(result)
data1.to_excel('data1.xlsx')#保存了result1.xlsx中的sheet1的数值
data2=pd.DataFrame(v)
data2.to_excel('data2.xlsx')#保存了result1.xlsx中的sheet2的数值

x0=[i for i in range(301)]
plt.plot(x0,v[0],':',color='#FF0000',label='龙头')
plt.plot(x0,v[1],'-',color='#00FF00',label='第1节龙身')
plt.plot(x0,v[51],'-',color='#0000FF',label='第51节龙身')
plt.plot(x0,v[101],'-',color='#FF8800',label='第101节龙身')
plt.plot(x0,v[151],'-',color='#FF00FF',label='第151节龙身')
plt.plot(x0,v[201],'-',color='#00FFFF',label='第201节龙身')
plt.plot(x0,v[223],':',color='#000000',label='龙尾(后)')
plt.xlabel('时间(s)')
plt.ylabel('速度(m/s)')
plt.title('速度-时间')
plt.legend()
plt.grid(True)
plt.show()
x0=[i for i in range(224)]
plt.plot(x0,v[:,0],'-',color='#FF0000',label='第0s')
plt.plot(x0,v[:,75],'-',color='#00FF00',label='第75s')
plt.plot(x0,v[:,150],'-',color='#0000FF',label='第150s')
plt.plot(x0,v[:,225],'-',color='#FF8800',label='第225s')
plt.plot(x0,v[:,300],'-',color='#FF00FF',label='第300s')
plt.xlabel('位置(第X个)')
plt.ylabel('速度(m/s)')
plt.title('速度-位置')
plt.legend()
plt.grid(True)
plt.show()
