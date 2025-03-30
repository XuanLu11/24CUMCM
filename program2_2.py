'''计算碰撞时的位置和速度'''
#调用库和方法
import numpy as np
import math
import pandas as pd
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['KaiTi']
plt.rcParams['axes.unicode_minus'] = False
#阿基米德螺旋线的参数
B=.55/2/np.pi #曲线r=B*theta的参数计算
#求根的数值解方法
ROOT_METHOD='brentq'
#龙头速度
HEAD_VECOCITY=1
#龙的孔的间距
LENGTH_LONG=(341-2*27.5)/100
LENGTH_SHORT=(220-2*27.5)/100
#时间设置
T=414.56 #由test2.1.py多次测量得到的时间值

def polarToRec(r:float,theta:float)->tuple:
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
    dtheta=0.1
    f=lambda x:lengthToCenter(theta,b)-lengthToCenter(x,b)-length
    theta1=theta
    while(theta1>=0):#小步长寻找零点存在的区间
        theta1-=dtheta
        s=lengthToCenter(theta,b)-lengthToCenter(theta1,b) 
        #需要求根的表达式
        if s>=length:
            root=root_scalar(f,bracket=[theta1,theta],method=ROOT_METHOD) 
            #调用区间内求根的数值解方法
            return root.root
    return np.nan
def curveCalculateTheta(theta:float,length:float,b:float)->float:
    '''求曲线上离某点欧氏距离恰好为length的点(此点离中心更远)的位置(返回值theta)'''
    dtheta=0.1
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

theta=np.zeros(224) #计算此时的位置(theta)
theta[0]=headDthetaByTime(16*2*np.pi,B,HEAD_VECOCITY*T)
theta[1]=curveCalculateTheta(theta[0],LENGTH_LONG,B)
for i in range(2,224):
    theta[i]=curveCalculateTheta(theta[i-1],LENGTH_SHORT,B)
r=B*theta
x=np.zeros(224)
y=np.zeros(224)
v=np.zeros(224)
for i in range(224): #计算此时的位置(x,y)
    x[i],y[i]=polarToRec(r[i],theta[i])
v[0]=HEAD_VECOCITY
for j in range(1,224):
    v[j]=getNextVelocity(theta[j-1],theta[j],v[j-1])
result=np.array([x,y,v]).T
result=np.round(result,6)
data1=pd.DataFrame(result)
data1.to_excel('data3.xlsx') #保存了result2.xlsx中的sheet1的数值

x0=[i for i in range(224)]
plt.plot(x0,v,'-',color='#0000FF',label='碰撞时的速度')
plt.scatter([0,1,51,101,151,201,223],\
    [v[0],v[1],v[51],v[101],v[151],v[201],v[223]],color='#FF0000',label='特殊点')
#特殊点包含:龙头,第1节龙身,第51节龙身,第101节龙身,第151节龙身,第201节龙身,龙尾(后)
plt.xlabel('位置(第X个)')
plt.ylabel('速度(m/s)')
plt.title('速度-位置')
plt.legend()
plt.grid(True)
plt.show()

alpha=np.linspace(0,36*np.pi,5000)
a=B*alpha
x0=a*np.cos(alpha)
y0=a*np.sin(alpha)
plt.plot(x0,y0,color='#0000FF',label='盘入螺旋线')
plt.scatter([x[0],x[1],x[51],x[101],x[151],x[201],x[223]],\
    [y[0],y[1],y[51],y[101],y[151],y[201],y[223]],color='#FF0000',s=20,label='特殊点')
#特殊点包含:龙头,第1节龙身,第51节龙身,第101节龙身,第151节龙身,第201节龙身,龙尾(后)
plt.scatter(x[2:51],y[2:51],color='#00FF00',s=10,label='非特殊点')
plt.scatter(x[52:101],y[52:101],color='#00FF00',s=10)
plt.scatter(x[102:151],y[102:151],color='#00FF00',s=10)
plt.scatter(x[152:201],y[152:201],color='#00FF00',s=10)
plt.scatter(x[202:223],y[202:223],color='#00FF00',s=10)
plt.grid(True)
plt.legend()
plt.axis('equal')
plt.xlabel('X坐标')
plt.ylabel('Y坐标')
plt.title('碰撞发生时所处的位置')
plt.show()