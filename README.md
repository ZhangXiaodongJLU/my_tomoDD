# Our tomoDD..
## 简介
-------
本程序大部分将使用**Fortran90**、**C**以及**Matlab**编写。
=======
## tomoDD@H,J Zhang 文件说明
------
### aprod.f(已更新为F90)
#### 版本 hypoDD,Zhang修改
#### 类型 subroutine aprod(mode, m, n, x, y, leniw, lenrw, iw, rw)
#### 文件说明
本子程序用于计算矩阵a与列向量x、y的乘积，据称用于**LSQR**算法。
#### 变量说明
	变量名  变量类型                         变量说明
	1.mode  INTEGER                   选择计算模式
                                1      计算y=y+a*x，不改变x
                                2      计算x=x+a(转置)*y，不改变y
	2.m     INTEGER                   a矩阵行数  
	3.n     INTEGER                   a矩阵列数
	4.x     REAL,DIMENSION(n)         x列向量
	5.y     REAL,DIMENSION(m)         y列向量
	6.leniw INTEGER                   iw向量长度，为矩阵a非零元素个数*2+1
	7.lenrw INTEGER                   rw向量长度，为矩阵a非零元素个数
	8.iw    INTEGER,DIMENSION(leniw)  
                               iw[1]           矩阵a非零元素个数
                           iw[2:(iw[1]+1)]     各非零元素行指标
                         iw[iw[1]+2:2*iw[1]+1] 各非零元素列指标
	9.rw    REAL,DIMENSION(lenrw)          各非零元素(行优先)
------
### atoangle.c

### atoangle.c

### chtof.c

### cluster_tomoDD.f

### compat.h

### covar.f

### datetime_.c

### datum.f(已更新为F90)
#### 版本 Zhang版本
#### 类型 subroutine datum(itf, iyr, imo, idy, ihr, imn)
#### 文件说明
本子程序用于输入相对事件参考年的分钟数转化为绝对事件(精确到分)，事件参考年由输入参数决定。
#### 变量说明
	变量名    变量类型                         变量说明
	1.itf     INTEGER                   输入参量，事件相对分钟数
	2.iyr     INTEGER                   输入参量，事件参考年  
	3.imo     INTEGER                   输出参量，事件绝对年
	4.idy     INTEGER                   输出参量，事件绝对月
	5.ihr     INTEGER                   输出参量，事件绝对时
	6.imn     INTEGER                   输出参量，事件绝对分
#### 算例
	call datum(190946,2008,0,0,0,0)
	输出：190946 0 5 12 14 26
------
### delaz.f

### delaz2.f

### direct1.f

### dist.f(已更新为F90，未测试)
#### 版本 Zhang版本
#### 类型 subroutine dist(xlat, xlon, xkm, ykm)
#### 文件说明
本子程序用于将经纬度转化为本地笛卡尔坐标。
#### 变量说明
	变量名      变量类型                         变量说明
	1.xlat      DOUBLE PRECISION             输入参量，纬度
	2.xlon      DOUBLE PRECISION             输入参量，经度

	3.xlat       REAL                        输出参量，x坐标
	4.xlon       REAL                        输出参量，y坐标
#### 算例
缺失

------
### dtres_tomoDD.f

### exist.f(已更新为F90)
#### 版本 Zhang版本
#### 类型 subroutine exist(fn)
#### 文件说明
本子程序用于判断输入文件是否存在。若存在，无反馈继续执行主程序；若不存在，主程序报错退出。
#### 变量说明
	变量名      变量类型                         变量说明
	1.fn         CHARACTER(80)                   输入字符串(文件名)
#### 算例
	fn="123.dat"
	call exist(fn)
	输出：FILE DOES NOT EXIST / CHECK IDATA,IPHASE: 123.dat
------

### f77types.h

### freeunit.f(已更新为F90)
#### 版本 Zhang版本
#### 类型 subroutine freeunit(iunit)
#### 文件说明
本子程序用于寻找10-999第一个自由通道号。
#### 变量说明
	变量名      变量类型                         变量说明
	1.iunit     INTEGER                   输出参量，自由通道号
------
### geocoord.inc(已更新为F90版本，未测试)
#### 文件说明
本文件为变量定义文件。
#### 调用关系
	dist.f      include ‘geocoord.inc’
	redist.f    include ‘geocoord.inc’
	sdc2.f      include 'geocoord.inc'
------

### geometry.h

### getdata_tomoDD.f

### getinp_tomoDD.f

### hypot_.c

### ifindi.f

### indexxi.f

### juliam.f(已更新为F90)
#### 版本 Zhang版本
#### 类型 integer function juliam(iyr, imo, idy, ihr, imn)
#### 文件说明
本子函数用于输入事件绝对事件(精确到分)转化为事件相对分钟数。
#### 变量说明
	变量名  变量类型                         变量说明
	1.iyr     INTEGER                   输入参量，事件参考年  
	2.imo     INTEGER                   输入参量，事件绝对年
	3.idy     INTEGER                   输入参量，事件绝对月
	4.ihr     INTEGER                   输入参量，事件绝对时
	5.imn     INTEGER                   输入参量，事件绝对分

	6.juliam  INTEGER                   返回值，事件相对分钟数
#### 算例
	itf=juliam(0,5,12,14,26)
	print *,itf
	输出：190946
------

### IsfitH_tomoDD_Isqr.f

### IsfitHV_tomoDD_Isqr.f

### IsfitV_tomoDD_Isqr.f

### Isqr.f

### Makefile

### Makefile.syn

### matmult1.f

### matmult2.f

### matmult3.f

### mdian1.f

### normlz.f

### partials_tomoDD.f

### ran.f(已更新为F90)
#### 版本 Zhang版本
#### 类型 subroutine ran(jlo, jhi, j)
#### 文件说明
本子程序用生成伪随机数。
#### 变量说明
	变量名      变量类型                         变量说明
	1.jlo        REAL                       输入参量，随机数下限
	2.jhi        REAL                       输入参量，随机数上限

	3.j          REAL                       输出参量，伪随机数
#### 算例
	real j
	call ran(100,200,j)
	输出：121.132172
------

### ray_common.inc

### Ray3VD.f

### redist.f(部分更新为F90，未测试)
#### 版本 Zhang版本
#### 类型 subroutine redist(xkm, ykm, xlat, xlon)
#### 文件说明
本子程序用于将本地笛卡尔坐标转化为经纬度。
#### 变量说明
	变量名      变量类型                         变量说明
	1.xkm        REAL                       输入参量，x坐标
	2.ykm        REAL                       输入参量，y坐标

	3.xlat       DOUBLE PRECISION           输出参量，纬度
	4.xlon       DOUBLE PRECISION           输出参量，经度
#### 算例
缺失

------
### refract.f

### resstat_tomoDD.f

### rpad_.c

### scopy.f

### sdc2.f(已更新为F90，未测试)
#### 版本 Zhang版本
#### 类型 subroutine sdc2(x, y, xlat, xlon, i)
#### 文件说明
本子程序用于控制本地笛卡尔坐标和经纬度的转化。
#### 变量说明
	变量名      变量类型                         变量说明
	1.x          REAL                             x坐标
	2.y          REAL                             y坐标
	3.xlat       DOUBLE PRECISION                 纬度
	4.xlon       DOUBLE PRECISION                 经度
	1.i          INTEGER                   选择模式
                                +1      call redist(x,y,xlat,xlon)
                                -1      call dist(xlat,xlon,x,y)

#### 算例
缺失

------
### setorg.f

### skip_tomoDD.f

### snrm2.f

### sort.f

### sorti.f

### sscal.f

### sscanf3_.c

### sscanf4_.c

### svd.f

### syn_time.f

### tiddid.f

### tomoDD_syn.f

### tomoDD.f

### tomoDD.inc

### trialsrc_tomoDD.f

### trimlen.f(已更新为F90)
#### 版本 Zhang版本
#### 类型 integer function trimlen(t)
#### 文件说明
本子函数用于计算字符串不包括尾随空格的长度。
#### 变量说明
	变量名      变量类型                         变量说明
	1.t         CHARACTER(*)                   输入字符串

	2.trimlen   INTEGER                  字符串不包括尾随空格的长度
#### 算例
	t="abcdefg   "
	length=juliam(t)
	print *,length
	输出：7
------

### ttime.f

### vmodel.f

### weighting_tomoDD.f
