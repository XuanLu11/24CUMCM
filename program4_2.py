'''这是用于计算第4题的程序,使用了program4_1.py得到的结果'''
#调用库和方法
import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['KaiTi']
plt.rcParams['axes.unicode_minus'] = False
#曲线的参数
B=1.7/2/np.pi
R=9/2
GAMMA1=R/B
#求根的数值解方法
ROOT_METHOD='brentq'
#龙头速度
HEAD_VECOCITY=1
#龙的孔的间距
LENGTH_LONG=(341-2*27.5)/100
LENGTH_SHORT=(220-2*27.5)/100
#两段圆弧和切点的参数(这些参数由test4.1.py得到)
#盘入螺旋线与大圆弧的切点
X0=-2.711856
Y0=-3.591078
#大圆弧的圆心和半径
X1=-0.760009
Y1=-1.305726
R1=3.005418
#小圆弧的圆心和半径
X2=1.735932
Y2=2.448402
R2=1.502709
#两圆弧的切点
X3=0.903952
Y3=1.197026
#两圆弧的圆心角
ALPHA1=3.021487
ALPHA2=3.021487
#两圆弧的弧长
S1=9.080830
S2=4.540415
#大圆弧的参数方程的参数的取值范围
BETA1=4.005538
BETA2=0.984051
#小圆弧的参数方程的参数的取值范围
BETA3=-2.157542
BETA4=0.863945

def getLargeCircleXY(theta:float)->tuple:
    '''由大圆弧上的某点对应的参数得到它的坐标'''
    x=X1+R1*np.cos(theta)
    y=Y1+R1*np.sin(theta)
    return x,y
def getSmallCircleXY(theta:float)->tuple:
    '''由小圆弧上的某点对应的参数得到它的坐标'''
    x=X2+R2*np.cos(theta)
    y=Y2+R2*np.sin(theta)
    return x,y
def EuclideanRec(x1:float,x2:float,y1:float,y2:float)->float:
    '''计算直角坐标系中两点之间的欧式距离的方法(返回值d)'''
    return math.sqrt((x1-x2)**2+(y1-y2)**2)
def polarToRec(r:float,theta:float)->tuple:
    '''将极坐标转换成直角坐标的函数(返回值x,y)'''
    return r*np.cos(theta),r*np.sin(theta)
def lengthPolar(r1:float,r2:float,theta1:float,theta2:float)->float:
    '''计算极坐标中两点之间的欧式距离的方法(返回值d)'''
    a=r1*r1+r2*r2-2*r1*r2*math.cos(theta1-theta2)
    return math.sqrt(a)
def helixLengthToCenter(theta:float,b:float)->float:
    '''求出阿基米德螺旋线某点到其中心的弧长(返回值s)'''
    return b/2*(theta*math.sqrt(1+theta**2)+math.log(theta+math.sqrt(1+theta**2)))
def helixHeadDthetaByTime(theta:float,b:float,length:float)->float:
    '''求某一点沿曲线往中心移动了距离length后所处的位置(返回值theta)'''
    dtheta=0.1 #小步长长度
    f=lambda x:helixLengthToCenter(theta,b)-helixLengthToCenter(x,b)-length 
    #需要求根的表达式
    theta1=theta
    while(theta1>=0):#小步长寻找零点存在的区间
        theta1-=dtheta
        s=helixLengthToCenter(theta,b)-helixLengthToCenter(theta1,b)
        if s>=length:
            root=root_scalar(f,bracket=[theta1,theta],method=ROOT_METHOD) 
            #调用区间内求根的数值解方法
            return root.root
    return np.nan
def curveCalculateTheta(theta:float,length:float,b:float)->float:
    '''求曲线上离某点欧氏距离恰好为length的点(此点离中心更远)的位置(返回值theta)'''
    dtheta=0.1 #小步长长度
    r=b*theta
    f=lambda x:lengthPolar(r,b*x,theta,x)-length #需要求根的表达式
    theta1=theta
    while(True):#小步长寻找零点存在的区间
        theta1+=dtheta
        r1=b*theta1
        if lengthPolar(r,r1,theta,theta1)>=length:
            root=root_scalar(f,bracket=[theta,theta1],method=ROOT_METHOD) 
            #调用区间内求根的数值解方法
            return root.root
def helixAwayDthetaByTime(theta:float,b:float,length:float)->float:
    '''求某一点沿螺旋线远离中心移动了距离length后所处的位置(返回值theta)'''
    dtheta=0.1 #小步长长度
    f1=lambda x:helixLengthToCenter(x,b)-helixLengthToCenter(theta,b)-length 
    #需要求根的表达式
    theta1=theta
    while(True):#小步长寻找零点存在的区间
        theta1+=dtheta
        s=helixLengthToCenter(theta1,b)-helixLengthToCenter(theta,b)
        if s>=length:
            root=root_scalar(f1,bracket=[theta,theta1],method=ROOT_METHOD) 
            #调用区间内求根的数值解方法
            return root.root
def getNextPointPlace(type:int,length:float,nowParameter:float,b:float)->tuple:
    '''得到下一个点所在的曲线类型和在这个曲线内的参数值'''
    if type==0: #前驱在盘入螺旋线上,后继只能在盘入螺旋线上
        typeReturn=0
        f1=lambda x:lengthPolar(b*nowParameter,b*x,nowParameter,x)-length
        theta1=nowParameter
        dtheta=0.1
        while(True):#小步长寻找零点存在的区间
            r1=b*theta1
            if lengthPolar(b*nowParameter,r1,nowParameter,theta1)>=length:
                root=root_scalar(f1,bracket=[nowParameter,theta1],method=ROOT_METHOD) 
                #调用区间内求根的数值解方法
                nextParameter=root.root
                break
            theta1+=dtheta
    elif type==1: #前驱在大圆弧上,后继可能在大圆弧上,也可能在盘入螺旋线上
        dParameter=2*np.arcsin(length/2/R1)
        if dParameter+nowParameter>BETA1:
            #后继在盘入螺旋线上
            x1,y1=getLargeCircleXY(nowParameter)
            typeReturn=0
            theta1=GAMMA1
            dtheta=0.1
            while(True):
                x2,y2=polarToRec(b*theta1,theta1)
                if EuclideanRec(x1,x2,y1,y2)>length:
                    def f2(x):
                        x2,y2=polarToRec(b*x,x)
                        return EuclideanRec(x1,x2,y1,y2)-length
                    root=root_scalar(f2,bracket=[GAMMA1,theta1],method=ROOT_METHOD) 
                    #调用区间内求根的数值解方法
                    nextParameter=root.root
                    break
                theta1+=dtheta
        else:
            #后继在大圆弧上
            typeReturn=1
            nextParameter=dParameter+nowParameter
    elif type==2: #前驱在小圆弧上,后继可能在小圆弧上,也可能在大圆弧上
        dParameter=2*np.arcsin(length/2/R2)
        if nowParameter-dParameter<BETA3:
            #后继在大圆弧上
            typeReturn=1
            x1,y1=getSmallCircleXY(nowParameter)
            def f4(x):
                x2,y2=getLargeCircleXY(x)
                return EuclideanRec(x1,x2,y1,y2)-length
            root=root_scalar(f4,bracket=[BETA2,BETA1],method=ROOT_METHOD) 
            #调用区间内求根的数值解方法
            nextParameter=root.root
        else:
            #后继在小圆弧上
            typeReturn=2
            nextParameter=nowParameter-dParameter
    else: #前驱在盘出螺旋线上,后继可能在盘出螺旋线上,也可能在小圆弧上
        para1=curveCalculateTheta(GAMMA1,length,b)
        if nowParameter<para1:
            #后继在小圆弧上
            typeReturn=2
            x1,y1=polarToRec(b*nowParameter,nowParameter)
            x1,y1=-x1,-y1
            def f6(x):
                x2,y2=getSmallCircleXY(x)
                return EuclideanRec(x1,x2,y1,y2)-length
            root=root_scalar(f6,bracket=[BETA3,BETA4],method=ROOT_METHOD) 
            #调用区间内求根的数值解方法
            nextParameter=root.root
        else:
            #后继在盘出螺旋线上
            typeReturn=3
            f7=lambda x:lengthPolar(b*nowParameter,b*x,nowParameter,x)-length
            theta1=nowParameter
            dtheta=0.1
            while(True):#小步长寻找零点存在的区间
                theta1-=dtheta
                theta1=GAMMA1 if theta1<GAMMA1 else theta1
                r1=b*theta1
                if lengthPolar(b*nowParameter,r1,nowParameter,theta1)>=length:
                    root=root_scalar(f7,bracket=[theta1,nowParameter],method=ROOT_METHOD) 
                    #调用区间内求根的数值解方法    
                    nextParameter=root.root
                    break
    return typeReturn,nextParameter
def getHeadPlace(length:float)->float:
    '''得到龙头所在的曲线位置'''
    typeReturn=1
    parameterReturn=0
    #length为龙头沿着进入调头空间的方向到盘入曲线与调头空间的交点的相对距离
    if length<=0:
        typeReturn=0
        parameterReturn=helixAwayDthetaByTime(GAMMA1,B,-length)
    elif length<=S1:
        typeReturn=1
        parameterReturn=BETA1-length/R1
    elif length<=S1+S2:
        typeReturn=2
        parameterReturn=BETA3+(length-S1)/R2
    else:
        typeReturn=3
        parameterReturn=helixAwayDthetaByTime(GAMMA1,B,length-S1-S2)
    return typeReturn,parameterReturn
def getTangent(type:int,parameter:float)->float:
    '''计算某个点所在曲线的切线斜率值'''
    if type==0:
        t=np.cos(parameter)-parameter*np.sin(parameter)
        if t==0:
            return np.nan
        return (parameter*np.cos(parameter)+np.sin(parameter))/t
    elif type==1:
        return -1/np.tan(parameter)
    elif type==2:
        return -1/np.tan(parameter)
    else:
        t=np.cos(parameter)-parameter*np.sin(parameter)
        if t==0:
            return np.nan
        return (parameter*np.cos(parameter)+np.sin(parameter))/t
def getNextVelocity(k1:float,k2:float,k:float,v:float)->float:
    '''得到下一个点的速度'''
    #t表示正切值,c表示余弦值
    t1=(k-k1)/(1+k*k1)
    t2=(k-k2)/(1+k*k2)
    c1=1/math.sqrt(1+t1*t1)
    c2=1/math.sqrt(1+t2*t2)
    return v*c1/c2

type=np.zeros((224,201))
#0:在盘入螺旋线上,1:在大圆弧上,2:在小圆弧上,3:在盘出螺旋线上
parameter=np.zeros((224,201))
#点在曲线上的位置所对应的参数
k=np.zeros((224,201))
v=np.zeros((224,201))
x=np.zeros((224,201))
y=np.zeros((224,201))
for i in range(201):
    t=i-100
    print('t:'+str(t))
    type[0,i],parameter[0,i]=getHeadPlace(t*HEAD_VECOCITY)
    type[1,i],parameter[1,i]=\
        getNextPointPlace(type[0,i],LENGTH_LONG,parameter[0,i],B)
    for j in range(2,224):
        type[j,i],parameter[j,i]=\
            getNextPointPlace(type[j-1,i],LENGTH_SHORT,parameter[j-1,i],B)
    for j in range(224):
        k[j,i]=getTangent(type[j,i],parameter[j,i])
for i in range(201):
    for j in range(224):
        p=parameter[j,i]
        if type[j,i]==0:
            x[j,i],y[j,i]=polarToRec(B*p,p)
        elif type[j,i]==1:
            x[j,i],y[j,i]=getLargeCircleXY(p)
        elif type[j,i]==2:
            x[j,i],y[j,i]=getSmallCircleXY(p)
        else:
            x[j,i],y[j,i]=polarToRec(B*p,p)
            x[j,i],y[j,i]=-x[j,i],-y[j,i]
for i in range(201): #计算不同位置在不同时间下的速度(velocity)
    print('velocity turn:'+str(i))
    v[0,i]=HEAD_VECOCITY
    for j in range(1,224):
        if (type[j,i]==1 and type[j-1,i]==1) or (type[j,i]==2 and type[j-1,i]==2): 
            #在同一段圆弧上的点的速率一样
            v[j,i]=v[j-1,i]
        else:
            v[j,i]=getNextVelocity(k[j-1,i],k[j,i],\
                (y[j,i]-y[j-1,i])/(x[j,i]-x[j-1,i]),v[j-1,i])
result=np.zeros((2*224,201))
for i in range(224):
    result[2*i]=x[i]
    result[2*i+1]=y[i]
result=np.round(result,6)
v=np.round(v,6)
data1=pd.DataFrame(result)
data1.to_excel('data4.xlsx') #保存了result4.xlsx中的sheet1的数值
data2=pd.DataFrame(v)
data2.to_excel('data5.xlsx') #保存了result4.xlsx中的sheet2的数值

x0=[i-100 for i in range(201)]
plt.plot(x0,v[0],'-',color='#FF0000',label='龙头')
plt.plot(x0,v[1],'-',color='#00FF00',label='第1节龙身')
plt.plot(x0,v[101],'-',color='#0000FF',label='第101节龙身')
plt.plot(x0,v[223],'-',color='#FF0088',label='龙尾(后)')
plt.xlabel('时间(s)')
plt.ylabel('速度(m/s)')
plt.title('速度-时间(步长为1s)')
plt.legend()
plt.grid(True)
plt.show()

x0=[i-100 for i in range(101)]
plt.plot(x0,v[0,:101],'-',color='#FF0000',label='龙头')
plt.plot(x0,v[1,:101],'-',color='#00FF00',label='第1节龙身')
plt.plot(x0,v[101,:101],'-',color='#0000FF',label='第101节龙身')
plt.plot(x0,v[223,:101],'-',color='#FF0088',label='龙尾(后)')
plt.xlabel('时间(s)')
plt.ylabel('速度(m/s)')
plt.title('速度-时间(步长为1s)')
plt.legend()
plt.grid(True)
plt.show()

x0=[i for i in range(101)]
plt.plot(x0,v[0,100:],'-',color='#FF0000',label='龙头')
plt.plot(x0,v[1,100:],'-',color='#00FF00',label='第1节龙身')
plt.plot(x0,v[101,100:],'-',color='#0000FF',label='第101节龙身')
plt.plot(x0,v[223,100:],'-',color='#FF0088',label='龙尾(后)')
plt.xlabel('时间(s)')
plt.ylabel('速度(m/s)')
plt.title('速度-时间(步长为1s)')
plt.legend()
plt.grid(True)
plt.show()

x0=[i for i in range(224)]
plt.plot(x0,v[:,0],'-',color='#FF0000',label='t=-100s')
plt.plot(x0,v[:,34],'-',color='#00FF00',label='t=-66s')
plt.plot(x0,v[:,67],'-',color='#0000FF',label='t=-33s')
plt.plot(x0,v[:,100],'-',color='#FF8800',label='t=0s')
plt.xlabel('位置(第X个)')
plt.ylabel('速度(m/s)')
plt.title('速度-位置')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(x0,v[:,100],'-',color='#FF0000',label='t=0s')
plt.plot(x0,v[:,111],'-',color='#FF0088',label='t=11s')
plt.plot(x0,v[:,114],'-',color='#00FF00',label='t=14s')
plt.plot(x0,v[:,150],'-',color='#0000FF',label='t=50s')
plt.plot(x0,v[:,200],'-',color='#FF8800',label='t=100s')
plt.xlabel('位置(第X个)')
plt.ylabel('速度(m/s)')
plt.title('速度-位置')
plt.legend()
plt.grid(True)
plt.show()