*==============================================================================
* Copyright (c) Enrico Forestieri <[name].[surname]@santannapisa.it
*
* Calcola
*            N-1
*            __
*            \        ±j*2*pi*k*n/N
*    X(k) =  /  x(n)*e                    k=0,1,2,...,N-1
*            ¯¯
*            n=0
*
* dove gli x(n) sono complessi, mediante FFT.
*
* L'array x può essere dichiarato dal chiamante sia double complex x(0:nn-1)
* che double precision x(0:2*nn-1). In quest'ultimo caso, la parte reale e
* immaginaria dell'i-esimo dato complesso, i=0,1,...,nn-1, va immagazzinata
* in x(2*i) e x(2*i+1), rispettivamente.
* I valori in x(0:2*nn-1) vengono sostituiti dalla loro trasformata di Fourier
* discreta, se jsign=-1, oppure da N volte la loro inversa, se jsign=1.
*
      subroutine cfft(x,n,jsign)
      implicit real*8(a-h,o-z)
      dimension x(*)
	
      jsn = sign(2,jsign)
      call cft(x(1),x(2),n,n,n,jsn)

      return
      end

*------------------------------------------------------------------------------
* Copyright (c) Enrico Forestieri <[name].[surname]@santannapisa.it
*
* Calcola
*            N-1
*            __
*            \        ±j*2*pi*k*n/N
*    X(k) =  /  x(n)*e                    k=0,1,2,...,N-1
*            ¯¯
*            n=0
*
* dove gli x(n) sono reali, mediante FFT.
*
* Questo corrisponde a calcolare la trasformata di Fourier discreta della
* sequenza x(n), se jsign=-1, o N volte la trasformata inversa se jsign=1.
* ATTENZIONE: N deve essere pari ma non necessariamente una potenza di due
* (non viene fatto alcun controllo a questo riguardo).
*
* - Dal momento che X(N-k)=conjg(X(k)), vengono ritornati solo i valori
*   per k=0,1,2,...,N/2.
* - Gli X(k) sono in generale complessi tranne X(0) e X(N/2) che sono
*   sempre reali.
* - Se il chiamante aveva dichiarato l'array come double precision x(0:N-1),
*   in x(0) viene ritornato X(0), in x(1) viene ritornato X(N/2), mentre
*   in x(2*k) e x(2*k+1), k=1,2,...,N/2-1, vengono ritornati la parte reale
*   e immaginaria di X(k), rispettivamente.
*
      subroutine rfft(x,n,jsign)
      implicit real*8(a-h,o-z)
      parameter(duepi=6.28318530717958647692d0)
      dimension x(*)

      jsn = sign(1,jsign)
      theta = duepi/dble(jsn*n)
      c1 = 0.5d0
      jsn = 2*jsn
      call cft(x(1),x(2),n/2,n/2,n/2,jsn)
      wpr = -2.d0*sin(0.5d0*theta)**2
      wpi = sin(theta)
      wr = 1.d0+wpr
      wi = wpi
      n2p3 = n+3
      do i=2,n/4+1
        i1 = 2*i-1
	i2 = i1+1
	i3 = n2p3-i2
	i4 = i3+1
	h1r = c1*(x(i1)+x(i3))
	h1i = c1*(x(i2)-x(i4))
	h2r = c1*(x(i2)+x(i4))
	h2i = c1*(x(i3)-x(i1))
	x(i1) = h1r+wr*h2r-wi*h2i
	x(i2) = h1i+wr*h2i+wi*h2r
	x(i3) = h1r-wr*h2r+wi*h2i
	x(i4) = -h1i+wr*h2i+wi*h2r
	wt = wr
	wr = wr*wpr-wi*wpi+wr
	wi = wi*wpr+wt*wpi+wi
      enddo
      h1r = x(1)
      x(1) = h1r+x(2)
      x(2) = h1r-x(2)
      return
      end

*------------------------------------------------------------------------------
      subroutine cft(a,b,ntot,n,nspan,isn)
      implicit real*8(a-h,o-z)
c     Copyright (c) CERN https://cernlib.web.cern.ch/cernlib/
c     licensed under the GPL 
c
c     Multivariate complex fourier transform, computed in place
c     using mixed-radix Fast Fourier Transform algorithm.
c     By R. C. Singleton, Stanford Research Institute, Oct. 1968.
c     Arrays A and B originally hold the real and imaginary
c     components of the data, and return the real and
c     imaginary components of the resulting fourier coefficients.
c     Multivariate data is indexed according to the fortran
c     array element successor function, without limit
c     on the number of implied multiple subscripts.
c     The subroutine is called once for each variate.
c     The calls for a multivariate transform may be in any order.
c     NTOT is the total number of complex data values.
c     N is the dimension of the current variable.
c     NSPAN/N is the spacing of consucutive data values
c     while indexing the current variable.
c     The sign of ISN determines the sign of the complex
c     exponential, and the magnitude of isn is normally one.
c
c     For a single-variate transform,
c     NTOT = N = NSPAN = (number of complex data values), f.g.
c     call cft(a,b,n,n,n,1)
c
c     A tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)
c     is computed by
c     call cft(a,b,n1*n2*n3,n1,n1,1)
c     call cft(a,b,n1*n2*n3,n2,n1*n2,1)
c     call cft(a,b,n1*n2*n3,n3,n1*n2*n3,1)
c
c     The data may alternatively be stored in a single complex
c     array a, then the magnitude of isn changed to two to
c     give the correct indexing increment and the second parameter
c     used to pass the initial address for the sequence of
c     imaginary values, e.g.
c
c        real s(2)
c        equivalence (a,s)
c        ....
c        ....
c        call cft(a,s(2),ntot,n,nspan,2)
c
c     Arrays AT(MAXF), CK(MAXF), BT(MAXF), SK(MAXF), and NP(MAXP)
c     are used for temporary storage. If the available storage
c     is insufficient, the program is terminated by a stop.
c     MAXFF must be .ge. the maximum prime factor of n.
c     MAXPP must be .gt. the number of prime factors of n.
c     In addition, if the square-free portion k cf n has two or
c     more prime factors, then maxp must be .ge. k-1.
c     Array storage in nfac for a maximum of 11 factors of n.
c     If n has more than one square-free factor, the product of the
c     square-free factors must be .le. 210
c
      parameter(MAXFF=127,MAXPP=209)
      dimension a(1),b(1)
      dimension nfac(11),np(MAXPP)
c     array storage for maximum prime factor of MAXFF
      dimension at(MAXFF),ck(MAXFF),bt(MAXFF),sk(MAXFF)
      equivalence (i,ii)
c     The following two constants should agree with the array dimensions
      maxf=MAXFF
      maxp=MAXPP
      if(n .lt. 2) return
      inc=isn
c     The following constants are rad = 2.*pi , s72 = sin(0.4*pi) ,
c     c72 = cos(0.4d0*pi) and s120 = sqrt(0.75d0)
      rad = 6.2831853071796d0
      s72 = 0.95105651629515d0
      c72 = 0.30901699437495d0
      s120 = 0.86602540378444d0
      if(isn .ge. 0) go to 10
      s72=-s72
      s120=-s120
      rad=-rad
      inc=-inc
   10 nt=inc*ntot
      ks=inc*nspan
      kspan=ks
      nn=nt-inc
      jc=ks/n
      radf=rad*dble(jc)*0.5d0
      i=0
      jf=0
c     Determine the factors of n
      m=0
      k=n
      go to 20
   15 m=m+1
      nfac(m)=4
      k=k/16
   20 if(k-(k/16)*16 .eq. 0) go to 15
      j=3
      jj=9
      go to 30
   25 m=m+1
      nfac(m)=j
      k=k/jj
   30 if(mod(k,jj) .eq. 0) go to 25
      j=j+2
      jj=j**2
      if(jj .le. k) go to 30
      if(k .gt. 4) go to 40
      kt=m
      nfac(m+1)=k
      if(k .ne. 1) m=m+1
      go to 80
   40 if(k-(k/4)*4 .ne. 0) go to 50
      m=m+1
      nfac(m)=2
      k=k/4
   50 kt=m
      j=2
   60 if(mod(k,j) .ne. 0) go to 70
      m=m+1
      nfac(m)=j
      k=k/j
   70 j=((j+1)/2)*2+1
      if(j .le. k) go to 60
   80 if(kt .eq. 0) go to 100
      j=kt
   90 m=m+1
      nfac(m)=nfac(j)
      j=j-1
      if(j .ne. 0) go to 90
c     Compute fourier transform
  100 sd=radf/dble(kspan)
      cd=2.d0*sin(sd)**2
      sd=sin(sd+sd)
      kk=1
      i=i+1
      if(nfac(i) .ne. 2) go to 400
c     Transform for factor of 2 (including rotation factor)
      kspan=kspan/2
      k1=kspan+2
  210 k2=kk+kspan
      ak=a(k2)
      bk=b(k2)
      a(k2)=a(kk)-ak
      b(k2)=b(kk)-bk
      a(kk)=a(kk)+ak
      b(kk)=b(kk)+bk
      kk=k2+kspan
      if(kk .le. nn) go to 210
      kk=kk-nn
      if(kk .le. jc) go to 210
      if(kk .gt. kspan) go to 800
  220 c1=1.d0-cd
      s1=sd
  230 k2=kk+kspan
      ak=a(kk)-a(k2)
      bk=b(kk)-b(k2)
      a(kk)=a(kk)+a(k2)
      b(kk)=b(kk)+b(k2)
      a(k2)=c1*ak-s1*bk
      b(k2)=s1*ak+c1*bk
      kk=k2+kspan
      if(kk .lt. nt) go to 230
      k2=kk-nt
      c1=-c1
      kk=k1-k2
      if(kk .gt. k2) go to 230
      ak=c1-(cd*c1+sd*s1)
      s1=(sd*c1-cd*s1)+s1
c     The following three statements compensate for truncation
c     error. If rounded arithmetic is used, they may be deleted.
c     c1=0.5d0/(ak**2+s1**2)+0.5d0
c     s1=c1*s1
c     c1=c1*ak
c     Next statement should be deleted if non-rounded arithmetic is used
      c1=ak
      kk=kk+jc
      if(kk .lt. k2) go to 230
      k1=k1+inc+inc
      kk=(k1-kspan)/2+jc
      if(kk .le. jc+jc) go to 220
      go to 100
c     Transform for factor of 3 (optional code)
  320 k1=kk+kspan
      k2=k1+kspan
      ak=a(kk)
      bk=b(kk)
      aj=a(k1)+a(k2)
      bj=b(k1)+b(k2)
      a(kk)=ak+aj
      b(kk)=bk+bj
      ak=-0.5d0*aj+ak
      bk=-0.5d0*bj+bk
      aj=(a(k1)-a(k2))*s120
      bj=(b(k1)-b(k2))*s120
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      kk=k2+kspan
      if(kk .lt. nn) go to 320
      kk=kk-nn
      if(kk .le. kspan) go to 320
      go to 700
c     Transform for factor of 4
  400 if(nfac(i) .ne. 4) go to 600
      kspnn=kspan
      kspan=kspan/4
  410 c1=1.d0
      s1=0
  420 k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      akp=a(kk)+a(k2)
      akm=a(kk)-a(k2)
      ajp=a(k1)+a(k3)
      ajm=a(k1)-a(k3)
      a(kk)=akp+ajp
      ajp=akp-ajp
      bkp=b(kk)+b(k2)
      bkm=b(kk)-b(k2)
      bjp=b(k1)+b(k3)
      bjm=b(k1)-b(k3)
      b(kk)=bkp+bjp
      bjp=bkp-bjp
      if(isn .lt. 0) go to 450
      akp=akm-bjm
      akm=akm+bjm
      bkp=bkm+ajm
      bkm=bkm-ajm
      if(s1 .eq. 0.d0) go to 460
  430 a(k1)=akp*c1-bkp*s1
      b(k1)=akp*s1+bkp*c1
      a(k2)=ajp*c2-bjp*s2
      b(k2)=ajp*s2+bjp*c2
      a(k3)=akm*c3-bkm*s3
      b(k3)=akm*s3+bkm*c3
      kk=k3+kspan
      if(kk .le. nt) go to 420
  440 c2=c1-(cd*c1+sd*s1)
      s1=(sd*c1-cd*s1)+s1
c     The following three statements compensate for truncation
c     error. If rounded arithmetic is used, they may be deleted.
c     c1=0.5d0/(c2**2+s1**2)+0.5d0
c     s1=c1*s1
c     c1=c1*c2
c     Next statement should be deleted if non-rounded arithmetic is used
      c1=c2
      c2=c1**2-s1**2
      s2=2.d0*c1*s1
      c3=c2*c1-s2*s1
      s3=c2*s1+s2*c1
      kk=kk-nt+jc
      if(kk .le. kspan) go to 420
      kk=kk-kspan+inc
      if(kk .le. jc) go to 410
      if(kspan .eq. jc) go to 800
      go to 100
  450 akp=akm+bjm
      akm=akm-bjm
      bkp=bkm-ajm
      bkm=bkm+ajm
      if(s1 .ne. 0.d0) go to 430
  460 a(k1)=akp
      b(k1)=bkp
      a(k2)=ajp
      b(k2)=bjp
      a(k3)=akm
      b(k3)=bkm
      kk=k3+kspan
      if(kk .le. nt) go to 420
      go to 440
c     Transform for factor of 5 (optional code)
  510 c2=c72**2-s72**2
      s2=2.d0*c72*s72
  520 k1=kk+kspan
      k2=k1+kspan
      k3=k2+kspan
      k4=k3+kspan
      akp=a(k1)+a(k4)
      akm=a(k1)-a(k4)
      bkp=b(k1)+b(k4)
      bkm=b(k1)-b(k4)
      ajp=a(k2)+a(k3)
      ajm=a(k2)-a(k3)
      bjp=b(k2)+b(k3)
      bjm=b(k2)-b(k3)
      aa=a(kk)
      bb=b(kk)
      a(kk)=aa+akp+ajp
      b(kk)=bb+bkp+bjp
      ak=akp*c72+ajp*c2+aa
      bk=bkp*c72+bjp*c2+bb
      aj=akm*s72+ajm*s2
      bj=bkm*s72+bjm*s2
      a(k1)=ak-bj
      a(k4)=ak+bj
      b(k1)=bk+aj
      b(k4)=bk-aj
      ak=akp*c2+ajp*c72+aa
      bk=bkp*c2+bjp*c72+bb
      aj=akm*s2-ajm*s72
      bj=bkm*s2-bjm*s72
      a(k2)=ak-bj
      a(k3)=ak+bj
      b(k2)=bk+aj
      b(k3)=bk-aj
      kk=k4+kspan
      if(kk .lt. nn) go to 520
      kk=kk-nn
      if(kk .le. kspan) go to 520
      go to 700
c     Transform for odd factors
  600 k=nfac(i)
      kspnn=kspan
      kspan=kspan/k
      if(k .eq. 3) go to 320
      if(k .eq. 5) go to 510
      if(k .eq. jf) go to 640
      jf=k
      s1=rad/dble(k)
      c1=cos(s1)
      s1=sin(s1)
      if(jf .gt. maxf) go to 998
      ck(jf)=1.d0
      sk(jf)=0.d0
      j=1
  630 ck(j)=ck(k)*c1+sk(k)*s1
      sk(j)=ck(k)*s1-sk(k)*c1
      k=k-1
      ck(k)=ck(j)
      sk(k)=-sk(j)
      j=j+1
      if(j .lt. k) go to 630
  640 k1=kk
      k2=kk+kspnn
      aa=a(kk)
      bb=b(kk)
      ak=aa
      bk=bb
      j=1
      k1=k1+kspan
  650 k2=k2-kspan
      j=j+1
      at(j)=a(k1)+a(k2)
      ak=at(j)+ak
      bt(j)=b(k1)+b(k2)
      bk=bt(j)+bk
      j=j+1
      at(j)=a(k1)-a(k2)
      bt(j)=b(k1)-b(k2)
      k1=k1+kspan
      if(k1 .lt. k2) go to 650
      a(kk)=ak
      b(kk)=bk
      k1=kk
      k2=kk+kspnn
      j=1
  660 k1=k1+kspan
      k2=k2-kspan
      jj=j
      ak=aa
      bk=bb
      aj=0.d0
      bj=0.d0
      k=1
  670 k=k+1
      ak=at(k)*ck(jj)+ak
      bk=bt(k)*ck(jj)+bk
      k=k+1
      aj=at(k)*sk(jj)+aj
      bj=bt(k)*sk(jj)+bj
      jj=jj+j
      if(jj .gt. jf) jj=jj-jf
      if(k .lt. jf) go to 670
      k=jf-j
      a(k1)=ak-bj
      b(k1)=bk+aj
      a(k2)=ak+bj
      b(k2)=bk-aj
      j=j+1
      if(j .lt. k) go to 660
      kk=kk+kspnn
      if(kk .le. nn) go to 640
      kk=kk-nn
      if(kk .le. kspan) go to 640
c     Multiply by rotation factor (except for factors of 2 and 4)
  700 if(i .eq. m) go to 800
      kk=jc+1
  710 c2=1.d0-cd
      s1=sd
  720 c1=c2
      s2=s1
      kk=kk+kspan
  730 ak=a(kk)
      a(kk)=c2*ak-s2*b(kk)
      b(kk)=s2*ak+c2*b(kk)
      kk=kk+kspnn
      if(kk .le. nt) go to 730
      ak=s1*s2
      s2=s1*c2+c1*s2
      c2=c1*c2-ak
      kk=kk-nt+kspan
      if(kk .le. kspnn) go to 730
      c2=c1-(cd*c1+sd*s1)
      s1=s1+(sd*c1-cd*s1)
c     The following three statements compensate for truncation
c     error. If rounded arithmetic is used, they may
c     be deleted.
c     c1=0.5d0/(c2**2+s1**2)+0.5d0
c     s1=c1*s1
c     c2=c1*c2
      kk=kk-kspnn+jc
      if(kk .le. kspan) go to 720
      kk=kk-kspan+jc+inc
      if(kk .le. jc+jc) go to 710
      go to 100
c     Permute the results to normal order---done in two stages
c     permutation for square factors of n
  800 np(1)=ks
      if(kt .eq. 0) go to 890
      k=kt+kt+1
      if(m .lt. k) k=k-1
      j=1
      np(k+1)=jc
  810 np(j+1)=np(j)/nfac(j)
      np(k)=np(k+1)*nfac(j)
      j=j+1
      k=k-1
      if(j .lt. k) go to 810
      k3=np(k+1)
      kspan=np(2)
      kk=jc+1
      k2=kspan+1
      j=1
      if(n .ne. ntot) go to 850
c     permutation for single-variate transform (optional code)
  820 ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=kspan+k2
      if(k2 .lt. ks) go to 820
  830 k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if(k2 .gt. np(j)) go to 830
      j=1
  840 if(kk .lt. k2) go to 820
      kk=kk+inc
      k2=kspan+k2
      if(k2 .lt. ks) go to 840
      if(kk .lt. ks) go to 830
      jc=k3
      go to 890
c     Permutation for multivariate transform
  850 k=kk+jc
  860 ak=a(kk)
      a(kk)=a(k2)
      a(k2)=ak
      bk=b(kk)
      b(kk)=b(k2)
      b(k2)=bk
      kk=kk+inc
      k2=k2+inc
      if(kk .lt. k) go to 860
      kk=kk+ks-jc
      k2=k2+ks-jc
      if(kk .lt. nt) go to 850
      k2=k2-nt+kspan
      kk=kk-nt+jc
      if(k2 .lt. ks) go to 850
  870 k2=k2-np(j)
      j=j+1
      k2=np(j+1)+k2
      if(k2 .gt. np(j)) go to 870
      j=1
  880 if(kk .lt. k2) go to 850
      kk=kk+jc
      k2=kspan+k2
      if(k2 .lt. ks) go to 880
      if(kk .lt. ks) go to 870
      jc=k3
  890 if(2*kt+1 .ge. m) return
      kspnn=np(kt+1)
c     Permutation for square-free factors of n
      j=m-kt
      nfac(j+1)=1
  900 nfac(j)=nfac(j)*nfac(j+1)
      j=j-1
      if(j .ne. kt) go to 900
      kt=kt+1
      nn=nfac(kt)-1
      if(nn .gt. maxp) go to 998
      jj=0
      j=0
      go to 906
  902 jj=jj-k2
      k2=kk
      k=k+1
      kk=nfac(k)
  904 jj=kk+jj
      if(jj .ge. k2) go to 902
      np(j)=jj
  906 k2=nfac(kt)
      k=kt+1
      kk=nfac(k)
      j=j+1
      if(j .le. nn) go to 904
c     determine the permutation cycles of length greater than 1
      j=0
      go to 914
  910 k=kk
      kk=np(k)
      np(k)=-kk
      if(kk .ne. j) go to 910
      k3=kk
  914 j=j+1
      kk=np(j)
      if(kk .lt. 0) go to 914
      if(kk .ne. j) go to 910
      np(j)=-j
      if(j .ne. nn) go to 914
      maxf=inc*maxf
c     Reorder a and b, following the permutation cycles
      go to 950
  924 j=j-1
      if(np(j) .lt. 0) go to 924
      jj=jc
  926 kspan=jj
      if(jj .gt. maxf) kspan=maxf
      jj=jj-kspan
      k=np(j)
      kk=jc*k+ii+jj
      k1=kk+kspan
      k2=0
  928 k2=k2+1
      at(k2)=a(k1)
      bt(k2)=b(k1)
      k1=k1-inc
      if(k1 .ne. kk) go to 928
  932 k1=kk+kspan
      k2=k1-jc*(k+np(k))
      k=-np(k)
  936 a(k1)=a(k2)
      b(k1)=b(k2)
      k1=k1-inc
      k2=k2-inc
      if(k1 .ne. kk) go to 936
      kk=k2
      if(k .ne. j) go to 932
      k1=kk+kspan
      k2=0
  940 k2=k2+1
      a(k1)=at(k2)
      b(k1)=bt(k2)
      k1=k1-inc
      if(k1 .ne. kk) go to 940
      if(jj .ne. 0) go to 926
      if(j .ne. 1) go to 924
  950 j=k3+1
      nt=nt-kspnn
      ii=nt-inc+1
      if(nt .ge. 0) go to 924
      return
c     Error finish, insufficient array storage
  998 isn=0
      print 999
      stop
  999 format('Array bounds exceeded within subroutine cft')
      end

*==============================================================================
