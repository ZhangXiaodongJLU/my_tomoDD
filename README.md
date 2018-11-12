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

### atoangle_.c

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

### f77types.h 头文件
此函数用于F77编译器读取C语言数据类型时的一些作用函数
#### 函数说明
	  函数名		函数作用
	  F77INT	  默认整数数据类型
	  F77LOG	  默认逻辑型数据类型
	  F77SLA	  长字符串参数类型
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
	
	c-- Loop to read each parameter lines, skipping comments
	210 read (fu_inp,'(a)',end=220)line
	if (line(1:1).eq.'*' .or.line(2:2).eq.'*') goto 210
	if (l.eq.1) read (line,'(a)',err=999) fn_cc
	if (l.eq.2) read (line,'(a)',err=999) fn_ct
	if (l.eq.3) read (line,'(a)',err=999) fn_eve
	if (l.eq.4) read (line,'(a)',err=999) fn_sta
	if (l.eq.5) read (line,'(a)',err=999) fn_loc
	if (l.eq.6) read (line,'(a)',err=999) fn_reloc
	if (l.eq.7) read (line,'(a)',err=999) fn_stares
	if (l.eq.8) read (line,'(a)',err=999) fn_res
	if (l.eq.9) read (line,'(a)',err=999) fn_srcpar
	if (l.eq.10) read (line,'(a)',err=999) fn_vel
	if (l.eq.11) read (line,'(a)',err=999) fn_abs
	if (l.eq.12) read (line,*,err=999) idata, iphase, maxdist
	if (l.eq.13) read (line,*,err=999) minobs_cc,minobs_ct,CC_format
	if (l.eq.14) then
	    read (line,*,err=999) istart, isolv, niter,
      &	weight1, weight2, weight3, air_dep
        end if
	if (l.eq.15) then ! parameters for pesudo-bending ray tracing
		 read (line,*,err=999) i3d, delt1, nidp, iskip, scale1,
      &	       scale2, iusep, iuses, iuseq
        end if
	if (l.eq.16) then
		 read (line,*,err=999) invdel, ifixl, xfac, tlim,
      &	       nitpb(1), nitpb(2), stepl
        end if
	if (l.eq.17) then 
	     read (line,*,err=999) lat_Orig, lon_Orig, dep_Orig, iorig, rota
        end if

       

	c--Read iteration instructions
      	if (l.ge.18 .and. l.le.17+niter) then
         i=l-17
         read (line,*,err=999) aiter(i),
      & awt_ccp(i), awt_ccs(i), amaxres_cross(i), amaxdcc(i),
      & awt_ctp(i), awt_cts(i), amaxres_net(i), amaxdct(i), awt_ctd(i),
      & adamp(i),ajoint(i),threshold(i)
        endif

	c--Read specific clusters/events to relocate
      	if (l.eq.18+niter) read (line,*,err=999) iclust
      	if (l.ge.19+niter) then
         read (line,*,err=999,end=230) (icusp(i),i=ii,ii+7)
	230      ii= i
      	endif
      	l= l+1
        goto 210
	220   close (fu_inp)
        ncusp= ii-1

	c- rearrange aiter:
        do i=2,niter
          aiter(i)= aiter(i-1)+aiter(i)
        enddo

	c- check files
        call exist (fn_eve)
        call exist (fn_sta)
        call exist (fn_abs)

      	if ((idata.eq.1 .or.idata.eq.3).and.trimlen(fn_cc).gt.1)
      &  call exist(fn_cc)
        if ((idata.eq.2 .or.idata.eq.3).and.trimlen(fn_ct).gt.1)
      &  call exist (fn_ct)

        maxiter= aiter(niter)
	c synthetic noise:
        noisef_dt= 0.002

	c write log output: of newest format
	600   if (trimlen(fn_loc).lt.2) fn_loc= 'tomoDD.loc'
        if (trimlen(fn_reloc).lt.2) fn_reloc= 'tomoDD.reloc'
        write (6,'("INPUT FILES:",/,
       &"cross dtime data: ",a,/,"catalog dtime data: ",a,/,
       &"events: ",a,/,"stations: ",a,/,"OUTPUT FILES:",/,
       &"initial locations: ",a,/,"relocated events: ",a,/,
       &"event pair residuals: ",a,/,"station residuals: ",a,/,
       &"source parameters: ",a)')
       &fn_cc(1:trimlen(fn_cc)),
       &fn_ct(1:trimlen(fn_ct)),
       &fn_eve(1:trimlen(fn_eve)),
       &fn_sta(1:trimlen(fn_sta)),fn_loc(1:trimlen(fn_loc)),
       &fn_reloc(1:trimlen(fn_reloc)),fn_res(1:trimlen(fn_res)),
       &fn_stares(1:trimlen(fn_stares)),fn_srcpar(1:trimlen(fn_srcpar))

        write (log,'("Input parameters: (from ",a,")",/,
       &"  cross dtime file: ",a,/,"  catalog dtime file: ",a,/,
       &"  station file: ",a,/,"  event file: ",a,/,
       &"  initial locations: ",a,/,"  relocated events: ",a)')
       &fn_inp(1:trimlen(fn_inp)),
       &fn_cc(1:trimlen(fn_cc)),
       &fn_ct(1:trimlen(fn_ct)),
       &fn_sta(1:trimlen(fn_sta)),
       &fn_eve(1:trimlen(fn_eve)),fn_loc(1:trimlen(fn_loc)),
       &fn_reloc(1:trimlen(fn_reloc))

        write (log,'(
       &"  event pair file: ",a,/,"  station residual file: ",a,/,
       &"  source parameter file: ",a,/,
       &"  IDATA= ",i2,2X,"IPHASE= ",i2,2x,"MAXDIST= ",f5.0,/,
       &"  MINOBS_CC= ",i3,2x,"MINOBS_CT= ",i3,/,"  ISTART= ",i1,2x,
       &"ISOLV= ",i1,2x)')
       &fn_res(1:trimlen(fn_res)),
       &fn_stares(1:trimlen(fn_stares)),fn_srcpar(1:trimlen(fn_srcpar)),
       &idata,iphase,maxdist,minobs_cc,minobs_ct,istart,isolv

        aiter(0)=0
        write (log, '("  ITER ",i2,"-",i2,
       &": DAMP= "f5.1,/,"    WT_CCP= ",f7.4,2X,"WT_CCS= ",f7.4,2x,
       &"MAXR_CC= ",f7.4,2X,"MAXD_CC= ",f7.2,2X,/,
       &"    WT_CTP= ",f7.4,2x,"WT_CTS= ",f7.4,2x,"MAXR_CT= ",f7.4,2x,
       &"MAXD_CT= ",f7.2, 2x, "JOINT= ", i3)')
       &(aiter(i-1)+1,aiter(i),adamp(i),awt_ccp(i),awt_ccs(i),
       & amaxres_cross(i),
       & amaxdcc(i),awt_ctp(i),awt_cts(i), amaxres_net(i), amaxdct(i),
       & ajoint(i),
       & i=1,niter)

	c--- write the pseudo-bending parameters
      	write(log,*)'i3d=',i3d
      	write(log,*)'delt1=',delt1
      	write(log,*)'ndip=',ndip
      	write(log,*)'iskip=',iskip
      	write(log,*)'scale1=',scale1
      	write(log,*)'scale2=',scale2
      	write(log,*)'iusep=',iusep
      	write(log,*)'iuses=',iuses
     	write(log,*)'iuseq=',iuseq
	c--- smoothing constraint
      	write(log,*)'Smoothing applied.....'
      	write(log,*)'Weight1=',weight1
      	write(log,*)'Weight2=',weight2
      	write(log,*)'Weight3=',weight3
	c--- CC data format
      	if(CC_format.eq.1) then
	  write(log,*)'hypoDD CC format is used'
        elseif(CC_format.eq.2) then
	  write(log,*)'Another CC format is used! It has the following format:'
	  write(log,*)'EveID1 EveID2 Station Diff_time CC_coef CC_phase'
        else
	  write(log,*)'Neither of the format is chosen! Stop the program!'
	  write(log,*)'CC_format must be 1 or 2!'
	  stop
        endif	

	c--Repeat number of clusters, events to relocate
      	if (iclust.eq.0) then
          write (*,*) 'Relocate all clusters'
          write (log,*) 'Relocate all clusters'
        else
          write (*,*) 'Relocate cluster number ',iclust
          write (log,*) 'Relocate cluster number ',iclust
        end if

      	if (ncusp.eq.0) then
          write (*,*) 'Relocate all events'
          write (log,*) 'Relocate all events'
        else
          write (*,*) 'Relocate ',ncusp,' events'
          write (log,*) 'Relocate ',ncusp,' events'
        end if
        return

	c--Input error handling
	998   write(*,*)'>>> ERROR OPENING CONTROL PARAMETER FILE'
              goto 1000

	999   write (*,*)'>>> ERROR READING CONTROL PARAMETERS IN LINE ',l
              write (*,*) line
	1000  stop 'Program run aborted.'
      	      end  ! of subroutine getinp
	
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
	subroutine skip_tomoDD(log, kiter, minwght,
      & ndt, nev, nsrc, nsta,
      &	ev_cusp, ev_date, ev_time, ev_mag,
      &	ev_lat, ev_lon, ev_dep, ev_x, ev_y, ev_z,
      &	ev_herr, ev_zerr, ev_res, ev_type,
      &	src_cusp, src_lat, src_lon, src_dep,
      &	src_lat0, src_lon0, src_type,
      &	src_x, src_y, src_z, src_t, src_x0, src_y0, src_z0, src_t0,
      &	sta_lab, sta_lat, sta_lon, sta_dist, sta_az,
      &	sta_rmsc, sta_rmsn, sta_np, sta_ns, sta_nnp, sta_nns,
      &	dt_sta, dt_c1, dt_c2, dt_idx, dt_dt, dt_qual, dt_cal,
      &	dt_ista, dt_ic1, dt_ic2,
      &	dt_res, dt_wt, dt_offs,
      &	tmp_ttp, tmp_tts, tmp_xp, tmp_yp, tmp_zp, 
      & tmp_xs, tmp_ys, tmp_zs,
      & tmp_vp, tmp_vs,
      & tmp_vp_index, tmp_vs_index, nct, ncc)

	implicit none

	include'tomoDD.inc'
	include'ray_common.inc'

      c Parameters:
	integer		log
	integer		kiter
	real		minwght
	integer		ndt
	integer		nev
	integer		nsrc
	integer		nsta
	integer		ev_cusp(MAXEVE)	! [1..MAXEVE]
	integer		ev_date(MAXEVE)	! [1..MAXEVE]
	integer		ev_time(MAXEVE)	! [1..MAXEVE]
	real		ev_mag(MAXEVE)	! [1..MAXEVE]
	real		ev_lat(MAXEVE)	! [1..MAXEVE]
	real		ev_lon(MAXEVE)	! [1..MAXEVE]
	real		ev_dep(MAXEVE)	! [1..MAXEVE]
	real		ev_x(MAXEVE)	! [1..MAXEVE]
	real		ev_y(MAXEVE)	! [1..MAXEVE]
	real		ev_z(MAXEVE)	! [1..MAXEVE]
	real		ev_herr(MAXEVE)	! [1..MAXEVE]
	real		ev_zerr(MAXEVE)	! [1..MAXEVE]
	real		ev_res(MAXEVE)	! [1..MAXEVE]
	integer         ev_type(MAXEVE) ! [1..MAXEVE]
	integer		src_cusp(MAXEVE)! [1..MAXEVE]
	doubleprecision	src_lat(MAXEVE)	! [1..MAXEVE]
	doubleprecision	src_lon(MAXEVE)	! [1..MAXEVE]
	real		src_dep(MAXEVE)	! [1..MAXEVE]
	real		src_lat0(MAXEVE)! [1..MAXEVE]
	real		src_lon0(MAXEVE)! [1..MAXEVE]
	real		src_x(MAXEVE)	! [1..MAXEVE]
	real		src_y(MAXEVE)	! [1..MAXEVE]
	real		src_z(MAXEVE)	! [1..MAXEVE]
	real		src_t(MAXEVE)	! [1..MAXEVE]
	real		src_x0(MAXEVE)	! [1..MAXEVE]
	real		src_y0(MAXEVE)	! [1..MAXEVE]
	real		src_z0(MAXEVE)	! [1..MAXEVE]
	real		src_t0(MAXEVE)	! [1..MAXEVE]
	integer         src_type(MAXEVE)! [1..MAXEVE]
	character	sta_lab(MAXSTA)*7! [1..MAXSTA]
	real		sta_lat(MAXSTA)	! [1..MAXSTA]
	real		sta_lon(MAXSTA)	! [1..MAXSTA]
	real		sta_dist(MAXSTA)! [1..MAXSTA]
	real		sta_az(MAXSTA)	! [1..MAXSTA]
	real		sta_rmsc(MAXSTA)! [1..MAXSTA]
	real		sta_rmsn(MAXSTA)! [1..MAXSTA]
	integer		sta_np(MAXSTA)	! [1..MAXSTA]
	integer		sta_ns(MAXSTA)	! [1..MAXSTA]
	integer		sta_nnp(MAXSTA)	! [1..MAXSTA]
	integer		sta_nns(MAXSTA)	! [1..MAXSTA]
	character	dt_sta(MAXDATA)*7! [1..MAXDATA]
	integer		dt_c1(MAXDATA)	! [1..MAXDATA]
	integer		dt_c2(MAXDATA)	! [1..MAXDATA]
	integer		dt_idx(MAXDATA)	! [1..MAXDATA]
	real		dt_dt(MAXDATA)	! [1..MAXDATA]
	real		dt_qual(MAXDATA)! [1..MAXDATA]
	real		dt_cal(MAXDATA)	! [1..MAXDATA]
	integer		dt_ista(MAXDATA)! [1..MAXDATA]
	integer		dt_ic1(MAXDATA)	! [1..MAXDATA]
	integer		dt_ic2(MAXDATA)	! [1..MAXDATA]
	real		dt_res(MAXDATA)	! [1..MAXDATA]
	real		dt_wt(MAXDATA)	! [1..MAXDATA]
	real		dt_offs(MAXDATA)! [1..MAXDATA]
	real		tmp_ttp(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_tts(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_xp(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_yp(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_zp(MAXSTA,MAXEVE)! [1..MAXSTA,1..MAXEVE]
	real		tmp_xs(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_ys(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_zs(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	integer         tmp_vp_index(MAXSTA,MAXEVE,8*MAXNZ+1)
	real            tmp_vp(MAXSTA,MAXEVE,8*MAXNZ)
        integer         tmp_vs_index(MAXSTA,MAXEVE,8*MAXNZ+1)
	real            tmp_vs(MAXSTA,MAXEVE,8*MAXNZ)
	integer		nct
	integer		ncc

      c	Local variables:
	integer		i, j, k, m
	integer		icusp(MAXEVE)	! [1..nev] Event keys
	integer		ifindi
	integer		iicusp(MAXEVE)	! [1..nev] Index table into ev_cusp[]
	integer		nccold
	integer		nctold
	integer		ndtold
	integer		sta_itmp(MAXSTA)


      	write(log,'("skipping data...")')

     c  Skip data with large resiudals
      	if (kiter.eq.1) then
          ndtold = ndt
          nccold = 0
          nctold = 0
      	endif
      	ncc = 0
      	nct = 0
      	j = 1
      	do i=1,ndt
           if (kiter.eq.1) then
              if (dt_idx(i).le.2) then
                 nccold = nccold+1
              else
                 nctold = nctold+1
              endif
           endif
           if (dt_wt(i).ge.minwght) then
              dt_sta(j) = dt_sta(i)
              dt_c1(j) = dt_c1(i)
              dt_c2(j) = dt_c2(i)
              dt_idx(j) = dt_idx(i)
              dt_qual(j) = dt_qual(i)
              dt_dt(j) = dt_dt(i)
              dt_cal(j) = dt_cal(i)
              dt_res(j) = dt_res(i)
              dt_wt(j) = dt_wt(i)
              dt_offs(j) = dt_offs(i)
              if (dt_idx(i).le.2) then
                  ncc = ncc+1
              else
                  nct = nct+1
              endif
              j = j+1
           endif
      	enddo
      	ndt = j-1
       	write(log,'("# obs = ",i9," (",f5.1,"%)")')
       &ndt, (ndt*100.0/ndtold)
        if (nccold.gt.0.and.nctold.gt.0) then
           write(log,'("# obs cc = ",i9," (",f5.1,"%)")')
       &   ncc, (ncc*100.0/nccold)
           write(log,'("# obs ct = ",i9," (",f5.1,"%)")')
       &   nct, (nct*100.0/nctold)
        endif

      c	Skip events
      	do i=1,ndt
           dt_ic1(i) = dt_c1(i)   !dt_ic1 is just a workspace array here!
           dt_ic2(i) = dt_c2(i)   !dt_ic2 is just a workspace array here!
        enddo
      	call sorti(ndt, dt_ic1)
      	call sorti(ndt, dt_ic2)
      	k = 1
      	do i=1,nev
           if (ifindi(ndt, dt_ic1, ev_cusp(i)).gt.0 .or.
       &       ifindi(ndt, dt_ic2, ev_cusp(i)).gt.0) then
              ev_date(k) = ev_date(i)
              ev_time(k) = ev_time(i)
              ev_cusp(k) = ev_cusp(i)
              ev_lat(k) = ev_lat(i)
              ev_lon(k) = ev_lon(i)
              ev_dep(k) = ev_dep(i)
              ev_mag(k) = ev_mag(i)
              ev_herr(k) = ev_herr(i)
              ev_zerr(k) = ev_zerr(i)
              ev_res(k) = ev_res(i)
              ev_x(k) = ev_x(i)
              ev_y(k) = ev_y(i)
              ev_z(k) = ev_z(i)
	      ev_type(k)=ev_type(i)
              k = k+1
           endif
      	enddo
      	nev = k-1
      	write(log,'("# events = ",i9)') nev

      c	Skip sources
      c Uses sorted dt_ic[12] arrays from above
      	if (nsrc.ne.1) then
           k = 1
           do i=1,nsrc
              if (ifindi(ndt, dt_ic1, src_cusp(i)).gt.0 .or.
        &         ifindi(ndt, dt_ic2, src_cusp(i)).gt.0) then
                 src_cusp(k) = src_cusp(i)
                 src_lat(k) = src_lat(i)
               	 src_lon(k) = src_lon(i)
                 src_lat0(k) = src_lat0(i)
                 src_lon0(k) = src_lon0(i)
                 src_dep(k) = src_dep(i)
                 src_x(k) = src_x(i)
                 src_y(k) = src_y(i)
                 src_z(k) = src_z(i)
                 src_t(k) = src_t(i)
	         src_type(k)=src_type(i)
                 src_x0(k) = src_x0(i)
                 src_y0(k) = src_y0(i)
                 src_z0(k) = src_z0(i)
                 src_t0(k) = src_t0(i)
                 do j=1,nsta
                    tmp_ttp(j,k) = tmp_ttp(j,i)
                    tmp_tts(j,k) = tmp_tts(j,i)
                    tmp_xp(j,k) = tmp_xp(j,i)
                    tmp_yp(j,k) = tmp_yp(j,i)
                    tmp_zp(j,k) = tmp_zp(j,i)
                    tmp_xs(j,k) = tmp_xs(j,i)
                    tmp_ys(j,k) = tmp_ys(j,i)
                    tmp_zs(j,k) = tmp_zs(j,i)
                    do m=1, tmp_vp_index(j,i,1)
		       tmp_vp(j,k,m)=tmp_vp(j,i,m)
		       tmp_vp_index(j,k,m+1)=tmp_vp_index(j,i,m+1)
		    enddo
		    tmp_vp_index(j,k,1)=tmp_vp_index(j,i,1)
		    if (iuses.eq.2) then
	            do m=1, tmp_vs_index(j,i,1)
		       tmp_vs(j,k,m)=tmp_vs(j,i,m)
		       tmp_vs_index(j,k,m+1)=tmp_vs_index(j,i,m+1)
		    enddo
		    tmp_vs_index(j,k,1)=tmp_vs_index(j,i,1)	
		    endif
                  enddo
                  k = k+1
              endif
           enddo
           nsrc = k-1
      	endif

     c	Clean stations
        do i=1,nsta
            sta_itmp(i) = 0
        enddo
        do j=1,ndt
            do i=1,nsta
               if (dt_sta(j).eq.sta_lab(i)) then
                  sta_itmp(i) = 1
                  goto 200	! break
               endif
            enddo
     200    continue
      	 enddo
         k = 1
      	 do i=1,nsta
            if (sta_itmp(i).eq.1) then
               sta_lab(k) = sta_lab(i)
               sta_lat(k) = sta_lat(i)
               sta_lon(k) = sta_lon(i)
               sta_dist(k) = sta_dist(i)
               sta_az(k) = sta_az(i)
               sta_np(k) = sta_np(i)
               sta_ns(k) = sta_ns(i)
               sta_nnp(k) = sta_nnp(i)
               sta_nns(k) = sta_nns(i)
               sta_rmsc(k) = sta_rmsc(i)
               sta_rmsn(k) = sta_rmsn(i)
               do j=1,nsrc
                  tmp_ttp(k,j) = tmp_ttp(i,j)
               	  tmp_tts(k,j) = tmp_tts(i,j)
                  tmp_xp(k,j) = tmp_xp(i,j)
                  tmp_yp(k,j) = tmp_yp(i,j)
                  tmp_zp(k,j) = tmp_zp(i,j)
                  tmp_xs(k,j) = tmp_xs(i,j)
                  tmp_ys(k,j) = tmp_ys(i,j)
                  tmp_zs(k,j) = tmp_zs(i,j)
                  do m=1, tmp_vp_index(i,j,1)
		     tmp_vp(k,j,m)=tmp_vp(i,j,m)
		     tmp_vp_index(k,j,m+1)=tmp_vp_index(i,j,m+1)
	          enddo
	          tmp_vp_index(k,j,1)=tmp_vp_index(i,j,1)
	          if (iuses.eq.2) then
	          do m=1, tmp_vs_index(i,j,1)
		     tmp_vs(k,j,m)=tmp_vs(i,j,m)
		     tmp_vs_index(k,j,m+1)=tmp_vs_index(i,j,m+1)
	          enddo
	          tmp_vs_index(k,j,1)=tmp_vs_index(i,j,1)	
	          endif
               enddo
               k = k+1
            endif
         enddo
         nsta = k-1
         write(log,'("# stations = ",i9)') nsta

    c	 Index station labels and cuspids
      	 call indexxi(nev, ev_cusp, iicusp)
         do i=1,nev
            icusp(i) = ev_cusp(iicusp(i)) !icusp is just a workspace array here!
         enddo
         do i=1,ndt
            do j=1,nsta
               if (dt_sta(i).eq.sta_lab(j)) then
                  dt_ista(i) = j
                  dt_ic1(i) = iicusp(ifindi(nev, icusp, dt_c1(i)))
                  dt_ic2(i) = iicusp(ifindi(nev, icusp, dt_c2(i)))
                  goto 300	! continue 2
               endif
            enddo
            write(*,'("FATAL ERROR (indexing).")') 
            stop   
    300     continue
         enddo

         end !of subroutine skip_tomoDD




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
	subroutine weighting_tomoDD(log, ndt, mbad, amcusp, idata,
      & kiter, ineg,
      &	maxres_cross, maxres_net, maxdcc, maxdct, minwght,
      &	wt_ccp, wt_ccs, wt_ctp, wt_cts, wtdd,
      &	dt_c1, dt_c2, dt_idx, dt_qual, dt_res, dt_offs,
      &	dt_wt)

        implicit none

        include "tomoDD.inc"

     c	Parameters:
        integer	log
        integer	ndt
        integer	mbad
        integer	amcusp(1000)
        integer	idata
        integer	kiter
        integer	ineg
        real	maxres_cross
        real	maxres_net
        real	maxdcc
        real	maxdct
        real	minwght
        real    wtdd		! relative weight of abs time over diff data
        real	wt_ccp
        real	wt_ccs
        real	wt_ctp
        real	wt_cts
        integer	dt_c1(MAXDATA)	! (1..MAXDATA)
        integer	dt_c2(MAXDATA)	! (1..MAXDATA)
        integer	dt_idx(MAXDATA)	! (1..MAXDATA)
        real	dt_qual(MAXDATA) ! (1..MAXDATA)
        real	dt_res(MAXDATA)	! (1..MAXDATA)
        real	dt_offs(MAXDATA) ! (1..MAXDATA)
        real	dt_wt(MAXDATA)	! (1..MAXDATA)
	
     c	Local variables:
        character	dattim*25
        real	dt_tmp(MAXDATA)
        real    dt_cat(MAXDATA)
        integer	i, j, k
        real	mad_cc
        real	mad_ct
        real    mad_cat
        real	maxres_cc
        real	maxres_ct
        real    maxres_cat
        real	med_cc
        real	med_ct
        real    med_cat
	integer ncat
	integer	ncc
	integer	nct
	integer nncat
	integer	nncc
	integer nnct
	
        call datetime(dattim)
	
     c synthetics:
	if(idata.eq.0) then
	   do i=1,ndt
	      dt_wt(i)=1 	
	   enddo
	endif

     c--- get a priori data weights:
     c intial (a priori) weights:
     c s=linspace(0.0,100.0,101); ss= (s/100).^2; plot(s/100,ss); % coherency
     c s=linspace(0.0,2.0,101); ss= (1./(2.^s)); plot(s,ss); % pick qual

     c all the quality transf is done in getdata. old format listed qualities,
     c new format list weights directly.
	ineg= 0  		!flag, =1 if neg weights exist
	do i=1,ndt
	   if(dt_idx(i).eq.1)
      & dt_wt(i)= wt_ccp * dt_qual(i) ! compat. with new format
	   if(dt_idx(i).eq.2)
      & dt_wt(i)= wt_ccs * dt_qual(i)	! compat. with new format
	   if(dt_idx(i).eq.3)
      & dt_wt(i)= wt_ctp * dt_qual(i)    ! compatib with new format 17/01/00
          if(dt_idx(i).eq.4)
      & dt_wt(i)= wt_cts * dt_qual(i)    ! compatib with new format 17/01/00

     c Applying relative weighting to absolute catalog data (02/03)
     c When setting wtdd=0, absolute data will be removed totally.
	  if(dt_c1(i).eq.dt_c2(i)) ! Absolute data
      & dt_wt(i)= wtdd*dt_wt(i)

        do j=1,mbad
             if(dt_c1(i).eq.amcusp(j).or.dt_c2(i).eq.amcusp(j)) then
		dt_wt(i)= 0.0
		ineg= 1
             endif
          enddo
	enddo

     c--- re-weighting: :
	if(((idata.eq.1.or.idata.eq.3).and.
     &    (maxres_cross.ne.-9.or.maxdcc.ne.-9)).or.
     &   ((idata.eq.2.or.idata.eq.3).and.
     &    (maxres_net.ne.-9.or.maxdct.ne.-9))) then
	   write(log,'("re-weighting ... ", a)') dattim
	  
     c---    get median and MAD of residuals
	   if(idata.eq.3) then
	      if(maxres_cross.ge.1) then
     c cross data:
		 k= 1
		 do i=1,ndt
		    if(dt_idx(i).le.2) then
		       dt_tmp(k)= dt_res(i)
		       k= k+1
		    endif
		 enddo
		 call mdian1(dt_tmp,k-1,med_cc)
     c 071200...
		 do i=1,k-1
		    dt_tmp(i)= abs(dt_tmp(i)-med_cc)
		 enddo
		 call mdian1(dt_tmp,k-1,mad_cc)
		 mad_cc= mad_cc/0.67449	!MAD for gaussian noise
	      endif
	      if(maxres_net.ge.1) then
     c- catalog data:
		 k= 1
		 j= 1
		 do i=1,ndt
     cz- residual reweighting set up for the absolute time 
		    if(dt_idx(i).ge.3) then 
		       if(dt_c1(i).ne.dt_c2(i)) then !diff. catalog time
			  dt_tmp(k)= dt_res(i)
			  k= k+1
		       else
			  dt_cat(j)= dt_res(i) ! abs. catalog time
			  j= j+1
		       endif
		    endif
		 enddo
      c       for diff. catalog time
		 if(k.gt.1) then 
		    call mdian1(dt_tmp,k-1,med_ct)
		    do i=1,k-1
		       dt_tmp(i)= abs(dt_tmp(i)-med_ct)
		    enddo
		    call mdian1(dt_tmp,k-1,mad_ct)
		    mad_ct= mad_ct/0.67449 !MAD for gaussian noise
		 endif   
      c       for absolute catalog time
		 if(j.gt.1) then
		    call mdian1(dt_cat,j-1,med_cat)
		    do i=1,j-1
		       dt_cat(i)= abs(dt_cat(i)-med_cat)
		    enddo
		    call mdian1(dt_cat,j-1,mad_cat)
		    mad_cat= mad_cat/0.67449 !MAD for gaussian noise	 
		 endif      
	      endif
	   elseif((idata.eq.1.and.maxres_cross.ge.1).or.
         &	      (idata.eq.2.and.maxres_net.ge.1)) then
	      k=1
	      j=1
	      do i=1,ndt
		 if(dt_c1(i).ne.dt_c2(i)) then ! for difference time
		    dt_tmp(k)= dt_res(i)
		    k=k+1
		 else
		    dt_cat(j)= dt_res(i) !for Abs data
		    j=j+1
		 endif
	      enddo
       c       for difference time
	      if(k.gt.1) then 
		 call mdian1(dt_tmp,k-1,med_cc)
		 do i=1,k-1
		    dt_tmp(i)= abs(dt_tmp(i)-med_cc)
		 enddo
		 call mdian1(dt_tmp,k-1,mad_cc)
		 mad_cc= mad_cc/0.67449 !MAD for gaussian noise
	      endif   
        c       for absolute catalog time
	      if(j.gt.1) then
		 call mdian1(dt_cat,j-1,med_cat)
		 do i=1,j-1
		    dt_cat(i)= abs(dt_cat(i)-med_cat)
		 enddo
		 call mdian1(dt_cat,j-1,mad_cat)
		 mad_cat= mad_cat/0.67449 !MAD for gaussian noise	 
	      endif  
	      
	      if(idata.eq.2) mad_ct= mad_cc
	   endif
	endif
     c--- define residual cutoff value:
	maxres_cc= maxres_cross ! absolute cutoff value for CC data
	maxres_ct= maxres_net	! absolute cutoff value for diff. CT
	maxres_cat=maxres_net	! absolute cutoff value for Abs. CT
	if(maxres_cross.ge.1) maxres_cc= mad_cc*maxres_cross
	if(maxres_net.ge.1) then 
	   maxres_ct= mad_ct*maxres_net 
	   maxres_cat= mad_cat*maxres_net
	endif

     c--- apply residual/offset dependent weights to a priori weights
	nncc= 0
	nnct= 0
	nncat=0
	ncc= 0
	nct= 0
	ncat=0
	do i=1,ndt
	   if(dt_idx(i).le.2) then
     c---    cross data:
	      ncc= ncc+1

     c bi ^5 offset weighting for cross data:
     c    exp needs to be uneven so weights become negative for offsets larger
     c    than 2 km. 2km is hardwired, >>>not anymore since 03/23/00
     c s=linspace(0,2.2,30);ss=(1-(s/2).^5).^5;plot(s,ss);axis([0 2.0 -0.0 1]);
	      if(maxdcc.ne.-9)
       & dt_wt(i)= dt_wt(i) * (1 - (dt_offs(i)/(maxdcc*1000))**5)**5


     c bi-cube residual weighting:
     c     needs to be cube so that res > cutoff become negative.
     c s=linspace(-0.2,0.2,101); ss= (1- (abs(s)/0.1).^3).^3;
     c plot(abs(s),ss);  axis([0 0.11 -0.1 1]);

	      if(maxres_cross.gt.0.and.dt_wt(i).gt.0.000001)
     & dt_wt(i)= dt_wt(i) * (1- (abs(dt_res(i))/maxres_cc)**3)**3
	      if(dt_wt(i).lt.minwght) nncc= nncc+1

	   else	
     c--- catalog data:
	      if(dt_c1(i).ne.dt_c2(i)) then !for diff. data
		 nct= nct+1
	      else		! for abs. data
		 ncat=ncat+1
	      endif

     c bi ^3 offset weighting for catalog data:
     c    exp needs to be uneven so weights become negative for offsets larger
     c    than 10 km. 10km is hardwired. not anymore since 03/23/00
     c s=linspace(0,11,100);ss=(1-(s/10).^3).^3;plot(s,ss);axis([0 11 -0.1 1]);

     cz Note that the inter-event distance for absolute catalod data is zero
     ca so there is no distance weighting for absolute data

	if(maxdct.ne.-9)
     &  dt_wt(i)= dt_wt(i) * (1 - (dt_offs(i)/(maxdct*1000))**3)**3
		  
     c bi-cube residual weighting:
     c     needs to be cube so that res > cutoff become negative.
     c s=linspace(-0.2,0.2,101); ss= (1- (abs(s)/0.1).^3).^3;
     c plot(abs(s),ss);  axis([0 0.11 -0.1 1]);
     		if(dt_wt(i).gt.0.000001 .and. maxres_net.gt.0) then
			if(dt_c1(i) .eq. dt_c2(i)) then 
				dt_wt(i)= dt_wt(i) * (1- (abs(dt_res(i))/maxres_cat)**3)**3
					else
						dt_wt(i)= dt_wt(i) * (1- (abs(dt_res(i))/maxres_ct)**3)**3
			endif

		endif
	      
		if(dt_wt(i).lt.minwght) then
			if(dt_c1(i).ne.dt_c2(i)) then
		    		nnct= nnct+1
		 	else
		    		nncat=nncat+1
		 	endif
		 endif
	   endif
	enddo

 
    c--- check if neg residuals exist
	ineg= 0
	do j=1,ndt
	   if(dt_wt(j).lt.minwght) then
	      ineg= 1
	      goto  100
	   endif
	enddo
    100	continue

	if(idata.eq.1.or.idata.eq.3) then
	   write(log,'(" cc res/dist cutoff:",
     & f7.3,"s/",f6.2,"km (",f5.1,"%)")')
     & maxres_cc,maxdcc,(nncc*100.0/ncc)
	endif
	if(idata.eq.2.or.idata.eq.3) then
	   if(nct.gt.0) write(log,'(" diff ct res/dist cutoff [s]:",
     & f7.3,"s/",f6.2,"km (",f5.1,"%)")')
     & maxres_ct,maxdct,(nnct*100.0/nct)
	   if(ncat.gt.0) write(log,'(" abs ct res/dist cutoff [s]:",
     & f7.3,"s/",f6.2,"km (",f5.1,"%)")')
     & maxres_cat,maxdct,(nncat*100.0/ncat)
	endif			! re-weighting
	
	if(ineg.gt.0) kiter= kiter+1
	end			! of subroutine weighting

