# Our tomoDD
## 简介
  本程序大部分将使用**Fortran90**、**C**编写。
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

