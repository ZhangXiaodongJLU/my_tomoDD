# Our tomoDD
## 简介
  本程序大部分将使用**Fortran90**编写 **c**编写。
## tomoDD@H,J Zhang 文件说明
------
### aprod.f
#### 版本 hypoDD,Zhang修改
#### 类型 subroutine aprod(mode, m, n, x, y, leniw, lenrw, iw, rw)
#### 程序说明
本子函数用于计算矩阵a与列向量x、y的乘积，主要使用**LSQR**和**SVD**算法。
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
	subroutine getinp_tomoDD (MAXEU,log,fn_inp,
      & fn_cc, fn_ct, fn_sta, fn_eve, fn_vel, fn_abs,
      & fn_loc, fn_reloc, fn_res, fn_stares, fn_srcpar,
      & idata, iphase,
      & minobs_cc, minobs_ct,
      & amaxres_cross, amaxres_net, amaxdcc, amaxdct,
      & noisef_dt, maxdist,
      & awt_ccp, awt_ccs, awt_ctp, awt_cts, awt_ctd, adamp,
      & istart, maxiter, isolv, niter, aiter, ajoint, threshold,
      & iclust, ncusp, icusp,
      & lat_Orig, lon_Orig, dep_Orig, iorig,
      & rota, stepl, CC_format,
      & weight1, weight2, weight3, air_dep)
      	
	implicit none
	
	include'ray_common,inc'
	
	Parameters:
	integer		maxev		! Array dimension
	integer		log		! Log-file indentifier
	character	fn_inp*80	! File of control info.
	character	fn_cc*80	! File of cross-corr. times
	character	fn_ct*80	! File of catalog times
	character	fn_sta*80	! Station file
	character	fn_eve*80	! Event file
	character	fn_loc*80	! Output file of original locs
	character	fn_reloc*80	! Output file of final locs
	character	fn_res*80	! Output residual file
	character	fn_stares*80	! Output station file
	character	fn_srcpar*80	! Output source-parameter file
	character 	fn_abs*80	! Absolute travel time data file
	character	fn_vel*80	! The velocity from each inversion
	integer		idata		! 0:Synthetics
					! 1:Cross-correlation
					! 2:catalog
					! 3:Both
	integer		iphase		! 1: P; 2: S; 3: Both
	integer		minobs_cc	! Min. obs./pair for ccor. data
	integer		minbos_ct	! Min. obs./pair for cat. data
	real		amaxres_cross(20)! [1..nither] Ccor. res. thresh.
	real		amaxres_net(20)	! [1..nither] Cat. res. thresh.
	real		amaxdcc(20)	! [1..nither] Ccor. link-dist. limit
	real		amaxdct(20)	! [1..nither] Cat. link-dist. limit
	real		noisef_dt	! Synthetic noise
	real		maxdist		! Max. cluster-station distance
	real		threshold(20)	! Determine the threshold
	real 		awt_ctd(20)	! relative weighting between abs and diff cat
	real		awt_ccp(20)	! [1..niter] Wts. for ccor. P
	real		awt_ccs(20)	! [1..niter] Wts. for ccor. S
	real		awt_ctp(20)	! [1..niter] Wts. for cat. P
	real		awt_ctp(20)	! [1..niter] Wts. for cat. S
	real		adamp(20)	! [1..niter] Damping (lsqr only)
	integer		adjoint(20)	! Joint inversion or not
	integer		istart		! 1:From single source
					! 2:From network sources
	integer		maxiter		
	integer 	isolv		! 1:SVD; 2: LSQR
	integer		niter		! No. of iteration sets
	integer		aiter(0:20)	! [1..niter] Iterations/set
	integer		icluster	! Cluster to relocate (0:all)
	integer		ncusp		! No. of event keys in icusp[]
	integer		icusp(maxev)	! [1..niter] Events to relocate
	real 		stepl
	real		lat_Orig
	real		lon_Orig
	real		dep_Orig			
	real		rota
	integer		iorig
	integer		CC_format	! CC data format:1 for hypoDD,2 for simple fmt
	real		weight1
	real		weight2
	real		weight3
	real		air_dep
	
	
 	c		Local varibles
	integer		fu_inp
	integer		i
	integer		ii
	integer		l
	character	line*80
	integer		trimlen
	
	c--- newest format: 083000 with iteration step dependent weighting
	c-- open input file:
	call freeunit(fu_inp)
	open (fu_inp,status='unknown',file=fn_inp,err=998)
	ncusp= 0
	niter= 0  ! number of iteration blocks
	l=1
	ii=1
	
	
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
