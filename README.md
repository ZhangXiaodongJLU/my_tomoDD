# Our tomoDD
## 简介
  本程序大部分将使用**Fortran90**编写 **c**编写。
## tomoDD@H,J Zhang 文件说明
------
### aprod.f
#### 版本 hypoDD,Zhang修改
#### 类型 subroutine aprod(mode, m, n, x, y, leniw, lenrw, iw, rw)
#### 程序说明
本子函数用于计算矩阵a与列向量x、y的乘积，据称用于**LSQR**算法。
#### 变量说明
	变量名  变量类型                         变量说明
	1.mode  INTEGER                   选择计算模式
                                1      计算y=y+a*x，不改变x
                                2      计算x=x+a（转置）*y，不改变y
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

### datum.f

### delaz.f

### delaz2.f

### direct1.f

### dist.f

### dtres_tomoDD.f

### exist.f

### f77types.h

### freeunit.f

### geocoord.inc

### geometry.h

### getdata_tomoDD.f

### getinp_tomoDD.f

### hypot_.c

### ifindi.f

### indexxi.f

### juliam.f
	integer function juliam(iyr, imo, idy, ihr, imn)
	
	implicit none
	
	Parameters:
	integer iyr, imo, idy, ihr, imn	!(input)
					!iyr < 4000 for 32-bit int
	Local variables:
	integer	kl
	integer kmo(12)
	integer ky, km, kd
	integer ky0
	integer ky1
	integer	ky4
	integer	l
	integer	leap
	
	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/juliam.f,v1.4 2001/02/19 01:31:06 julian Exp julian $"/
	save rcsid
	
	data kmo/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
	data leap/1/
	
	ky=iyr
	km=imo
	kd=idy
	if(km.le.0) km=1
	juliam=365*ky
	kd=kmo(km)+kd
	ky4=ky/4
	ky1=ky/100
	ky0=ky/1000
	kl=leap*(ky4-ky1+ky0)
	l=0
	if(ky4*4.eq.ky.and.(ky1*100.ne.ky.or.ky0*1000.eq.ky)) l=leap
	if(l.ne.0.and.km.lt.3)kl=kl-leap
	juliam=juliam+kd+kl
	juliam=juliam*24+ihr
	juliam=juliam*60+imn
	return
	end ! of integer function juliam	
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

### ran.f

### ray_common.inc

### Ray3VD.f

### redist.f

### refract.f

### resstat_tomoDD.f

### rpad_.c

### scopy.f

### sdc2.f

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

### trimlen.f

### ttime.f

### vmodel.f

### weighting_tomoDD.f
