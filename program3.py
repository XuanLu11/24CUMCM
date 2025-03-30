'''这是在第3题中用于计算最小螺距的程序'''
#调用库和方法
import numpy as np
import math
from scipy.optimize import root_scalar

#求根的数值解方法
ROOT_METHOD='brentq'
#龙的孔的间距
LENGTH_LONG=(341-2*27.5)/100
LENGTH_SHORT=(220-2*27.5)/100
#调头空间的圆的半径
R=9/2

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
    f=lambda x:lengthPolar(r,b*x,theta,x)-length #需要求根的表达式
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
def projection(point:np.array,start:np.array,end:np.array)->np.array:
    '''得到某点在某两点确定的直线上的投影点'''
    v=end-start #线段的方向向量
    w=point-start #点到线段某一连线的向量
    mu=np.dot(w,v)/np.dot(v, v)
    return start+mu*v
def lineOverlap(points1:np.array,points2:np.array)->bool:
    '''判断两组点所对应的线段是否有重叠'''
    for i in range(4):
        t=points2[i]
        for j in range(4):
            for k in range(j+1,4):
                if np.dot(t-points1[j],t-points1[k])<=0: 
                    #向量的数量积为非正数意味着线段有重叠
                    return True
    return False
def rectOverlap(X1:np.array,X2:np.array,Y1:np.array,Y2:np.array)->bool:
    '''判断两个矩形是否有重叠'''
    proj1=np.zeros([4,2])
    proj2=np.zeros([4,2])
    for i in range(2):#只需要考虑不平行的两条轴即可
        x1=X1[i]
        x2=X1[i+1]
        y1=Y1[i]
        y2=Y1[i+1]
        for j in range(4):
            proj1[j]=projection(np.array([X1[j],Y1[j]]),\
                np.array([x1,y1]),np.array([x2,y2]))
            proj2[j]=projection(np.array([X2[j],Y2[j]]),\
                np.array([x1,y1]),np.array([x2,y2]))
        if lineOverlap(proj1,proj2)==False:
            return False
        x1=X2[i]
        x2=X2[i+1]
        y1=Y2[i]
        y2=Y2[i+1]
        for j in range(4):
            proj1[j]=projection(np.array([X1[j],Y1[j]]),\
                np.array([x1,y1]),np.array([x2,y2]))
            proj2[j]=projection(np.array([X2[j],Y2[j]]),\
                np.array([x1,y1]),np.array([x2,y2]))
        if lineOverlap(proj1,proj2)==False:
            return False
    return True
def getVertexList(x:np.array,y:np.array)->tuple:
    '''获取矩形的顶点列表,这里的函数与program2_1.py中的函数
    存在略微的差别,此函数只能解决单个时刻的坐标计算'''
    rectX=np.zeros((223,4)) #矩阵的坐标排序符合顺时针方向
    rectY=np.zeros((223,4))
    #把手与板凳边缘的距离
    l1=.275
    l2=.15
    for i in range(223):
        if x[i+1]>x[i]:#为了方便讨论,在此进行的一次判断
            x1=x[i]
            x2=x[i+1]
            y1=y[i]
            y2=y[i+1]
        else:
            x2=x[i]
            x1=x[i+1]
            y2=y[i]
            y1=y[i+1]
        k=(y2-y1)/(x2-x1) #计算斜率
        kx=1/math.sqrt(1+k**2)
        ky=k/math.sqrt(1+k**2)
        #计算得到矩形各点的坐标的过程
        rectX[i,0]=x1-kx*l1
        rectY[i,0]=y1-ky*l1
        rectX[i,1]=rectX[i,0]-ky*l2
        rectY[i,1]=rectY[i,0]+kx*l2
        rectX[i,0]=rectX[i,0]+ky*l2
        rectY[i,0]=rectY[i,0]-kx*l2
        rectX[i,2]=x2+kx*l1
        rectY[i,2]=y2+ky*l1
        rectX[i,3]=rectX[i,2]+ky*l2
        rectY[i,3]=rectY[i,2]-kx*l2
        rectX[i,2]=rectX[i,2]-ky*l2
        rectY[i,2]=rectY[i,2]+kx*l2
    return rectX,rectY

'''
第一批参数为:
d_start=.55
d_end=.3
dd=.01
->
最小螺距(m)=0.460
第二批参数为:
d_start=.46
d_end=.45
dd=.001
->
最小螺距(m)=0.451
'''
d_start=.55 #由第2题的结果可知此时符合题目要求
d_end=.3 #取板凳的宽度作为估计值下界
dd=.01 #螺距的变化步长
d_best=d_start
d=d_start
theta_length=np.pi #龙头的初始位置(调头空间+一段冗余的长度)
N=40
dtheta=theta_length/N #龙头位置的变化步长
flag=False
while ((d>=d_end) and (flag==False)):
    print('d='+'%.3f'%(d))
    b=d/2/np.pi
    theta=np.zeros(224)
    theta[0]=R/b+theta_length
    while ((theta[0]>=R/b) and (flag==False)): 
        #模拟了龙在进入调头空间前的一段距离
        #以防止漏掉之前发生碰撞而在进入调头空间这一瞬间未发生碰撞的情形
        theta[1]=curveCalculateTheta(theta[0],LENGTH_LONG,b)
        for i in range(2,224):
            theta[i]=curveCalculateTheta(theta[i-1],LENGTH_SHORT,b)
        r=b*theta
        x=np.zeros(224)
        y=np.zeros(224)
        for i in range(224):
            x[i],y[i]=polarToRec(r[i],theta[i])
        rectX,rectY=getVertexList(x,y)
        for i in range(2): #只需要考虑龙头和第1条龙身
            for k in range(i+2,223):
                if rectOverlap(rectX[i],rectX[k],rectY[i],rectY[k])==True:
                    flag=True
                    print('Hit!'+'%.3f'%(d))
                    d_best=d+dd
                    break
            if flag==True:
                break
        theta[0]-=dtheta
    d-=dd #计算这一轮的螺距
print('最小螺距(m)='+'%.3f'%(d_best))
#由程序运行结果可知最小螺距(m)为0.451