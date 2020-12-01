      PROGRAM HBV_intracellular_model

        IMPLICIT NONE

        DOUBLE PRECISION :: alph,aa,mu,lon,loff,ss,b1,b2,b3,c1,c2,c3,c4,c5,d,e1,e2,e3,e4,e5,e6,e7,nc,ns,nLd,nLm
        DOUBLE PRECISION :: dm,dl,ds,dx,dp,deltp,deltc,deltl,deltm,delts,deltx,deltpc
        DOUBLE PRECISION :: g,rate(12),kf1,knuc,kf2,s1,muR,t,t1,t2,t3,lambda,k1 !rate=(Con1,Con2,Con3,Don1,Don2,Don3,Coff1,Coff2,Coff3,Doff1,Doff2,Doff3)
        DOUBLE PRECISION :: Ibind,Iunbind,psi,kbind,kunbind,lbind,lunbind,Cbind,Cunbind,kcnuc,kcf2,NCbind,NCunbind,kdis
        DOUBLE PRECISION :: mbind,munbind,Rbind,Runbind,kcleave,eta,ebind,eunbind,STtime,RTtime
        DOUBLE PRECISION :: atot,tout,tmax,time,r1,tau,ftime
        DOUBLE PRECISION :: random,dtime,a(588)

        INTEGER :: i,i2,i3,ii,j,k,kk,iseed,iout
        INTEGER :: rcNCe,rc,ccc,siccc,pol,cccpol1,cccpol2,cccpol3,cccpol4,cccpol5,Cm,LSm,Sm,Xm,PCm
        INTEGER :: rib,Crib1,Crib2,LSrib,Srib1,Srib2,Xrib,PCrib,P,L,M,S,X,PC,rcNC,R(29) !R=(R000,R001,...R222,C,D) R021=C on PS1, D on PS2, and free PS3
        INTEGER :: V(29,54) !Array of possible changes of Rijk's and C and D
        INTEGER :: nuc,pgelong(117),emelong(119),sphere(47) !pgelong=(R+4C,R+5C,...,R+120C), emelong=(2C,3C,...,120C), sphere=(2S,3S,...,48S)
        INTEGER :: rcV,RNAV,fil,sph,emV,emVc,san1,san2,san3,san4,san5,san6,san7
        INTEGER :: IFN,cccIFN,GA,PGA,NA,pgNCNA,CAM,CC,CCAM(119),rcCAM,rcCAMe,DD,rcDD,iRNA,CiRNA,eps,Ceps

        INTEGER, PARAMETER :: mxg = 700000000

        INTEGER :: nrun
        rate(:)=0.0d0
        a(:)=0.0d0
        R(:)=0
        pgelong(:)=0
        emelong(:)=0
        sphere(:)=0
        CCAM(:)=0
        V(:,:)=0
        iseed = 555090 

        open (unit=108,file='seed.list',status='unknown')

        do nrun=1,1

        read(108,*)iseed
!        iseed = 384763
        
        time = 0.0d0
        ftime = 0.0d0
        tout = 1.0d-6
        tmax = 600.00d0

        iout = 1
        dtime = 1.0d-6
        
        !!!!!Initial Condition!!!!!
        
        rcNCe=0
        rc = 0
        ccc = 1
        siccc = 0
        pol = 150 !150 !The number of RNA Polymerases, chosen randomly
        cccpol1 = 0
        cccpol2 = 0
        cccpol3 = 0
        cccpol4 = 0
        cccpol5 = 0
        Cm = 0
        LSm = 0
        Sm = 0
        Xm = 0
        PCm = 0
        rib = 10000 !The number of Ribosomes, chosen randomly
        Crib1 = 0
        Crib2 = 0
        LSrib = 0
        Srib1 = 0
        Srib2 = 0
        Xrib = 0
        PCrib = 0
        P = 0
        L = 0
        M = 0
        S = 0
        X = 0
        PC = 0
        rcNC = 0
        rcV = 0
        RNAV = 0
        fil = 0
        sph = 0
        emV = 0
        emVc = 0
        
        san1 = 0
        san2 = 0
        san3 = 0
        san4 = 0
        san5 = 0
        san6 = 0
        san7 = 0
        
        !!drugs
        
        IFN = 0 ! interferon alpha
        cccIFN = 0 ! 
        GA = 0 ! Geldanamycin
        PGA = 0
        NA = 0 ! Nucleos(t)ide analogues
        pgNCNA = 0
        CAM = 0 ! Capsid assembly inhibitors
        CC = 0
        rcCAM = 0 
        rcCAMe = 0
        DD = 0 ! cccDNA inhibitors
        rcDD = 0
        iRNA = 0 ! siRNAs
        CiRNA = 0 
        eps = 0 
        Ceps = 0
        
        
        !!!!!Parameters Values!!!!!
        
        alph = 0.016
        aa = 0.016 !dlog(2)/24, from a paper
        mu = 0.00058 !0.00058 !dlog(2)/(50 * 24), from a paper
        ss = 0.007 !1/150 !Impact of X on desilencing cccDNA and 1 is a reseanable option
        loff = 2.77 !Silencing rate of cccDNA, we don't know its value
        lon = 2.77 !Desilencing rate of cccDNA, we don't know its value
        b1 =  2.0 !0.37 * 6 !* 10000 !0.01728 !RNA Pol and cccDNA binding rate, we don't know its value
        b2 =  2.0 !0.37 * 6
        b3 =  2.0 !0.37 * 6
        nLd = 15 !As surface proteins are 1000-fold higher than HBV DNA, considered this 1000
        nLm = 25
        c1 = 0.43 * 60.0d0 !Production rate of C mRNA, we don't know its value
        c2 = 0.625 * 60.0d0 !Production rate of LS mRNA, we don't know its value
        c3 = 0.71 * 60.0d0 !Production rate of S mRNA, we don't know its value
        c4 = 2.14 * 60.0d0 !Production rate of X mRNA, we don't know its value
        c5 = 0.43 * 60.0d0 !Production rate of PreC mRNA, we don't know its value
        !!!! We can consider c1=c2=c3=c4=c5 !!!!!
        d =  40.0 !40.0 !40.0d0 !0.37 * 60 * 3 !Rib and Cm binding rate, we don't know its value 
        e1 = 0.375 * 60.0d0
        e2 = 0.82 * 60.0d0
        e3 = 0.375 * 60.0d0
        e4 = 0.53 * 60.0d0
        e5 = 0.66 * 60.0d0
        e6 = 1.95 * 60.0d0
        e7 = 1.415 * 60.0d0
        nc = 50 !40 !we don't know its value
        ns = 40 !20 !we don't know its value 
        dm = 0.1386 !dlog(2)/5 from a paper
        dl = 0.2310 !dlog(2)/3 from a paper
        ds = 0.2310 !dlog(2)/3 from a paper
        dx = 0.6931 !dlog(2)/1 By comparing mRNA length, I guessed it.
        dp = 0.1386 !dlog(2)/5 from a paper
        deltp = 0.6931 !dlog(2)/1 from a paper
        deltc = 0.6931 !dlog(2)/1 my assumption
        deltl = 0.6931 !dlog(2)/1 my assumption
        deltm = 0.6931 !dlog(2)/1 my assumption
        delts = 0.6931 !dlog(2)/1 my assumption
        deltx = 0.6931 !dlog(2)/1 from a paper
        deltpc = 0.6931 !dlog(2)/1 from a paper
        g = 5.0 !Using diffusion-limited kinetic kf=10^10
        rate(1) = 0.4 * 60.0 !0.037 * 60 !Con1
        rate(2) = 0.4 * 60.0 !0.037 * 60 !Con2
        rate(3) = 0.4 * 60.0 !0.037 * 60 !Con3
        rate(4) = 0.4 * 60.0d0 !Don1
        rate(5) = 0.4 * 60.0d0 !Don2
        rate(6) = 0.4 * 60.0d0 !Don3 
        rate(7) = 0.4 * 60.0 * 4.0d0 * 0.75 !6.528 !Coff1
        rate(8) = 0.4 * 60.0 * 4.0d0 * 0.75 !6.528 !Coff2
        rate(9) = 0.4 * 60.0 * 4.0d0 * 0.75 !6.528 !Coff3
        rate(10) = 0.4 * 60.0 * 1.2 * 0.75 !Doff1
        rate(11) = 0.4 * 60.0 * 12.3 * 0.75 !Doff2
        rate(12) = 0.4 * 60.0 * 4.53 * 0.75 !Doff3
        kf1 = 60000.0d0 !C-C binding around pgRNA from Zlotnick paper
        knuc = 600.0d0 !C-C binding in an empty capsid from Zlotnick paper
        kf2 = 60000.0d0 !C-C binding in an empty capsid from Zlotnick paper
        nuc = 3 !nucleation number from Zlotnick paper
        
        s1 = 0.058 !dlog(2)/12 from a paper
        muR = 0.029 !dlog(2)/24 from a paper
        t = 0.029 !0.029 !dlog(2)/24 from a paper
        t1 = 0.000029
        t2 = 2.9
        t3 = 2.9 * 3
        lambda = 30.0 * 1.0d-5 !from a paper
        k1 = 60000.0d0
        
        ! Drug parameters !
        Ibind = 2.0 !0.37 * 6 !0 ! binding rate of IFN
        Iunbind = 2.0 * 0.1 * 0.41 !0.37 * 6 * 1 * 0.41 !0 ! unbinding rate of IFN
        psi = 0.2 ! impact of IFN on cccDNA degradation
        kbind = 10.0 ! binding rate of GA
        kunbind = 10.0 * 1200.0 * 0.75 ! unbinding rate of GA
        ebind = 0.37 * 60 ! binding rate of Epsilon
        eunbind = 0.37 * 60 * 0.1 * 0.75 ! binding rate of Epsilon
        lbind = 5.0 !0.4 * 60.0 ! binding rate of NA
        lunbind = 5.0 * 7.4 * 0.75 !0.4 * 60.0 * 7.4 * 0.75 ! unbinding rate of NA
        Cbind = 24.0 ! binding rate of CAM to C
        Cunbind = 24.0 * 4.7 * 0.75 ! unbinding rate of CAM to C
        kcnuc = 600.0 * 10.0 ! nucleation rate of C:CAM
        kcf2 = 60000.0 * 10.0 ! elongation rate of C:CAM
        NCbind = 24.0 ! binding rate of CAM to NC
        NCunbind = 24.0 * 4.7 * 0.75 ! unbinding rate of CAM to NC
        kdis = 0.029 * 10.0 !0.029 * 10 ! disassembly rate of rcNC:CAM
        mbind = 0 ! binding rate of anti DNA repair factors
        munbind = 0 ! unbinding rate of anti DNA repair factors
        Rbind = 40.0 ! binding rate of iRNA 
        Runbind = 40.0 * 0.1 * 0.75 ! unbinding rate of iRNA 
        kcleave =  10.0 ! cleaving rate of iRNA
        eta = 0.0d0 ! viral and SVP releasing inhibitors 
        
        
        STtime = 0.0d0
        RTtime = 168.0d0

        DO WHILE ( time < tmax )

          IF ( rc == 0 .and. rcNC == 0 .and. ccc == 0 .and. siccc == 0 .and. cccIFN == 0) THEN
          IF (cccpol1 == 0 .and. cccpol2 == 0 .and. cccpol3 == 0 .and. cccpol4 == 0 .and. cccpol5 == 0) THEN
             
             GOTO 1
          
          ENDIF
          ENDIF
          
          !!!!!Transition Rates!!!!!
          
          a(387) = alph * DBLE(rcNCe)
          a(388) = muR * DBLE(rcNCe)
          a(389) = mu * DBLE(rc)
          a(1) = aa * DBLE(rc)
          a(2) = mu * (1.0 + psi * DBLE(IFN)) * DBLE(ccc)
          a(3) = loff * DBLE(ccc) / (1.0 + ss * DBLE(X))
          a(4) = lon * DBLE(siccc) * (1.0 + ss * DBLE(X))          
          a(5) = b1 * DBLE(ccc) * DBLE(pol)
          a(582) = DBLE(nLd) * b1 * DBLE(ccc) * DBLE(pol)
          a(583) = DBLE(nLd) * b1 * DBLE(ccc) * DBLE(pol)
          a(584) = b2 * DBLE(ccc) * DBLE(pol)
          a(585) = 0.0d0 !b3 * DBLE(ccc) * DBLE(pol)
          a(6) = c1 * DBLE(cccpol1)
          a(7) = c2 * DBLE(cccpol2)
          a(8) = c3 * DBLE(cccpol3)
          a(9) = c4 * DBLE(cccpol4)
          a(10) = 0.0d0 !c5 * DBLE(cccpol5)

          a(11) = d * DBLE(rib) * DBLE(Cm)
          a(12) = e1 * DBLE(Crib1)
          a(33) = DBLE(nc) * d * DBLE(rib) * DBLE(Cm)
          a(34) = e2 * DBLE(Crib2)
          a(13) = DBLE(nLm) * d * DBLE(rib) * DBLE(LSm)
          a(14) = e3 * DBLE(LSrib)
          a(15) = DBLE(nLm) * d * DBLE(rib) * DBLE(Sm)
          a(16) = e4 * DBLE(Srib1)
          a(390) = DBLE(ns) * DBLE(nLm) * d * DBLE(rib) * DBLE(Sm)
          a(391) = e5 * DBLE(Srib2)
          a(17) = d * DBLE(rib) * DBLE(Xm)
          a(18) = e6 * DBLE(Xrib)
          a(19) = 0.0d0 !d * DBLE(rib) * DBLE(PCm)
          a(20) = 0.0d0 !e7 * DBLE(PCrib)
          
          a(21) = dm * DBLE(Cm)
          a(22) = dl * DBLE(LSm)
          a(23) = ds * DBLE(Sm)
          a(24) = dx * DBLE(Xm)
          a(25) = 0.0d0 !dp * DBLE(PCm) 
          a(26) = deltp * DBLE(P)
          a(27) = deltc * DBLE(R(28))
          a(28) = deltl * DBLE(L)
          a(29) = deltm * DBLE(M)
          a(30) = delts * DBLE(S)
          a(31) = deltx * DBLE(X)
          a(32) = 0.0d0 !deltpc * DBLE(PC)
          
          a(35) = g * DBLE(P) * DBLE(Cm)
          
          DO i=36,89
             ii = i - 35
             IF ( modulo(18,ii) == 0 ) THEN
                CALL PROPENSITY (0,ii,i,R,rate,a,V)
             ENDIF
             
             IF ( modulo(18,ii) /= 0 .and. ii < 9 ) THEN
                j = modulo(ii,3)
                CALL PROPENSITY (j,ii,i,R,rate,a,V)
             ENDIF
             
             IF ( modulo(18,ii) /= 0 .and. ii > 9 .and. ii <= 26 ) THEN
                j = modulo(ii,9)
                CALL PROPENSITY (j,ii,i,R,rate,a,V)
             ENDIF
             
             IF ( ii >= 30 .and. ii <= 52 .and. ii /= 32 .and. ii /= 35 .and. ii /= 44 ) THEN
                k = ii - 26
                IF ( k < 9 ) THEN
                   j = modulo(k,3)
                   CALL PROPENSITY (k-j,k,i,R,rate,a,V)
                ENDIF
                IF ( k > 9 ) THEN
                   kk = modulo(k,9)
                   IF ( kk <= 3 ) THEN
                      CALL PROPENSITY (k-kk,k,i,R,rate,a,V)
                   ENDIF
                   IF ( kk > 3 .and. kk < 6 ) THEN
                      CALL PROPENSITY (k-3,k,i,R,rate,a,V)
                   ENDIF
                   IF ( kk >= 6 ) THEN
                      CALL PROPENSITY (k-6,k,i,R,rate,a,V)
                   ENDIF
                ENDIF
             ENDIF
             
             IF ( ii == 27 ) THEN
                CALL PROPENSITY (12,13,i,R,rate,a,V)
             ENDIF
             
             IF ( ii == 28 ) THEN
                CALL PROPENSITY (12,14,i,R,rate,a,V)
             ENDIF
             
             IF ( ii == 29 ) THEN
                CALL PROPENSITY (15,16,i,R,rate,a,V)
             ENDIF
             
             IF ( ii == 32 ) THEN
                CALL PROPENSITY (15,17,i,R,rate,a,V)
             ENDIF
             
             IF ( ii == 35 ) THEN
                CALL PROPENSITY (21,22,i,R,rate,a,V)
             ENDIF
             
             IF ( ii == 44 ) THEN
                CALL PROPENSITY (21,23,i,R,rate,a,V)
             ENDIF
             
             IF ( ii == 53 ) THEN
                CALL PROPENSITY (24,25,i,R,rate,a,V)
             ENDIF
             
             IF ( ii == 54 ) THEN
                CALL PROPENSITY (24,26,i,R,rate,a,V)
             ENDIF
             
          ENDDO
          
          a(144) = kf1 * DBLE(R(14)) * DBLE(R(28))
          DO i=145,260
             ii = i - 143
             a(i) = kf1 * DBLE(R(28)) * DBLE(pgelong(ii - 1))
          ENDDO
          
          a(261) = s1 * DBLE(pgelong(117))
          a(392) = muR * DBLE(pgelong(117))
          a(262) = muR * DBLE(rcNC)
          a(263) = t * DEXP(-lambda * DBLE(L)) * DBLE(rcNC)
          IF ( S > 79 .and. M > 19 .and. L > 19 ) THEN
              a(264) = (1 - eta) * t * (1.0d0 - DEXP(-lambda * DBLE(L))) * DBLE(rcNC)
          ELSE
              a(264) = 0
          ENDIF
          
          IF ( S > 79 .and. M > 19 .and. L > 19 ) THEN
              a(393) = (1 - eta) * t1 * (1.0d0 - DEXP(-lambda * DBLE(L))) * DBLE(pgelong(117))
          ELSE
              a(393) = 0
          ENDIF
          
          IF ( S > 79 .and. M > 19 .and. L > 19 ) THEN
              a(265) = (1 - eta) * t3 * (1.0d0 - DEXP(-lambda * DBLE(L)))
          ELSE
              a(265) = 0
          ENDIF
          
          a(266) = knuc * DBLE(R(28)) * (DBLE(R(28)) - 1.0d0) / 2
          DO i = 267,267+nuc-3
             ii = i - 265
             a(i) = knuc * DBLE(R(28)) * DBLE(emelong(ii - 1))
          ENDDO
          DO i = 267+nuc-2,384
             ii = i - 265
             a(i) = kf2 * DBLE(R(28)) * DBLE(emelong(ii - 1))
          ENDDO
          
          a(385) = 0.0d0 !muR * DBLE(emelong(119)) / 10
            
          IF ( S > 79 .and. M > 19 .and. L > 19 ) THEN
              a(386) = (1 - eta) * t * (1.0d0 - DEXP(-lambda * DBLE(L))) * DBLE(emelong(119))
          ELSE
              a(386) = 0
          ENDIF
          
          a(394) = k1 * DEXP(-lambda * DBLE(L)) * DBLE(S) * (DBLE(S) - 1.0d0) / 2
          DO i = 395,440
             ii = i - 393
             a(i) = k1 * DBLE(S) * DBLE(sphere(ii - 1)) !k1 * DEXP(-lambda * DBLE(L)) * DBLE(S) * DBLE(sphere(ii - 1))
          ENDDO
          a(441) = (1 - eta) * t2 * DBLE(sphere(47))
          
          !=== Drugs ===!
          
          ! Interferon-alpha
          
          a(442) = Ibind * DBLE(IFN) * DBLE(ccc)
          a(443) = Iunbind * DBLE(cccIFN)
          
          ! Geldanamycin
          
          a(444) = kbind * DBLE(P) * DBLE(GA)
          a(445) = kunbind * DBLE(PGA)
          
          ! Epsilon binding drug
          
          a(587) = ebind * DBLE(Cm) * DBLE(eps)
          a(588) = eunbind * DBLE(Ceps)
          
          ! NAs
          
          a(446) = lbind * DBLE(NA) * DBLE(pgelong(117))
          a(447) = lunbind * DBLE(pgNCNA)
          IF ( S > 79 .and. M > 19 .and. L > 19 ) THEN
              a(448) = (1 - eta) * t1 * (1.0d0 - DEXP(-lambda * DBLE(L))) * DBLE(pgNCNA)
          ELSE
              a(448) = 0
          ENDIF
          
          ! CAMs
          
          a(449) = Cbind * DBLE(CAM) * DBLE(R(28))
          a(450) = Cunbind * DBLE(CC)
!          a(451) = kd * DBLE(CC) * (DBLE(CC) - 1.0d0) / 2
!          DO i = 452,569
!             ii = i - 450
!             a(i) = kd * DBLE(CC) * DBLE(CCAM(ii - 1))
!          ENDDO
          a(451) = kcnuc * DBLE(CC) * (DBLE(CC) - 1.0d0) / 2
          DO i = 452,452+nuc-3
             ii = i - 450
             a(i) = kcnuc * DBLE(CC) * DBLE(CCAM(ii - 1))
          ENDDO
          DO i = 452+nuc-2,569
             ii = i - 450
             a(i) = kcf2 * DBLE(CC) * DBLE(CCAM(ii - 1))
          ENDDO
          a(570) = 0.0d0 !muR * DBLE(CCAM(119))
          IF ( S > 79 .and. M > 19 .and. L > 19 ) THEN
              a(571) = (1 - eta) * t * (1.0d0 - DEXP(-lambda * DBLE(L))) * DBLE(CCAM(119))
          ELSE
              a(571) = 0
          ENDIF
          
          a(572) = NCbind * DBLE(CAM) * DBLE(rcNCe)
          a(573) = NCunbind * DBLE(rcCAMe) 
          a(574) = NCbind * DBLE(CAM) * DBLE(rcNC)
          a(575) = NCunbind * DBLE(rcCAM)
          a(576) = kdis * DBLE(rcCAMe)
          a(577) = kdis * DBLE(rcCAM)
          
          ! Anti DNA repair factors
          
          a(578) = mbind * DBLE(rc) * DBLE(DD)
          a(579) = munbind * DBLE(rcDD) 
          
          ! RNA interference
          
          a(580) = Rbind * DBLE(Cm) * DBLE(iRNA)
          a(581) = kcleave * DBLE(CiRNA)
          a(586) =  Runbind * DBLE(CiRNA)
          
          
          !=== Total Transition Rate ===!

          atot = 0.0d0

          DO i3=1,588
          atot = atot + a(i3)
          ENDDO
          

          !=== Compute Time Increment ===!

          r1 = random(iseed)

          tau = DLOG(1.0d0/r1)
          tau = tau / atot

          time = time + tau

          IF ( time >= tout ) THEN

            WRITE(9,*)tout,ccc+cccpol1+cccpol2+cccpol3+cccpol4+cccpol5+siccc,rcV,rcNC,sph,fil,emV,emVc,RNAV,L         

            iout = iout + 1
            tout = tout + dtime

            IF ( iout == 10 ) THEN
              iout = 1
              IF ( dtime < 1.0d-3 ) dtime = dtime * 10.0d0
            ENDIF

          ENDIF

          !=== Pick Reaction ===!

          r1 = random(iseed)

          r1 = r1 * atot

          i2 = 0
          atot = 0.0d0

          DO WHILE ( atot < r1 )

            i2 = i2 + 1
            atot = atot + a(i2)

          ENDDO
          
          !if(time.gt.500)write(159,*)i2,time

          IF ( i2 == 1 ) THEN

            !=== cccDNA Production ===!

            ccc = ccc + 1
            rc = rc - 1

          ENDIF

          IF ( i2 == 2 ) THEN

            !=== cccDNA Death ===!

            ccc = ccc - 1

          ENDIF

          IF ( i2 == 3 ) THEN

            !=== Scilencing cccDNA ===!

            ccc = ccc - 1
            siccc = siccc + 1

          ENDIF

          IF ( i2 == 4 ) THEN

            !=== Descilencing cccDNA ===!

            ccc = ccc + 1
            siccc = siccc - 1

          ENDIF

          IF ( i2 == 5 ) THEN
          
            !=== Binding RNA Polymerase to cccDNA ===!

            ccc = ccc - 1
            pol = pol - 1
            cccpol1 = cccpol1 + 1

          ENDIF

          IF ( i2 == 6 ) THEN

            !=== C mRNA Production ===!

            cccpol1 = cccpol1 - 1
            Cm = Cm + 1
            pol = pol + 1
            ccc = ccc + 1
            !san3 = san3 + 1
            !WRITE(3,*)time,san3

          ENDIF
          
          IF ( i2 == 582 ) THEN
          
            !=== Binding RNA Polymerase to cccDNA ===!

            ccc = ccc - 1
            pol = pol - 1
            cccpol2 = cccpol2 + 1

          ENDIF

          IF ( i2 == 7 ) THEN

            !=== LS mRNA Production ===!

            cccpol2 = cccpol2 - 1
            LSm = LSm + 1
            pol = pol + 1
            ccc = ccc + 1

          ENDIF
          
          IF ( i2 == 583 ) THEN
          
            !=== Binding RNA Polymerase to cccDNA ===!

            ccc = ccc - 1
            pol = pol - 1
            cccpol3 = cccpol3 + 1

          ENDIF

          IF ( i2 == 8 ) THEN

            !=== S mRNA Production ===!

            cccpol3 = cccpol3 - 1
            Sm = Sm + 1
            pol = pol + 1
            ccc = ccc + 1

          ENDIF
          
          IF ( i2 == 584 ) THEN
          
            !=== Binding RNA Polymerase to cccDNA ===!

            ccc = ccc - 1
            pol = pol - 1
            cccpol4 = cccpol4 + 1

          ENDIF
          
          IF ( i2 == 9 ) THEN

            !=== X mRNA Production ===!

            cccpol4 = cccpol4 - 1
            Xm = Xm + 1
            pol = pol + 1
            ccc = ccc + 1

          ENDIF
          
          IF ( i2 == 585 ) THEN
          
            !=== Binding RNA Polymerase to cccDNA ===!

            ccc = ccc - 1
            pol = pol - 1
            cccpol5 = cccpol5 + 1

          ENDIF
          
          IF ( i2 == 10 ) THEN

            !=== PreC mRNA Production ===!

            cccpol5 = cccpol5 - 1
            PCm = PCm + 1
            pol = pol + 1
            ccc = ccc + 1

          ENDIF
          
          IF ( i2 == 11 ) THEN

            !=== C mRNA Ribosome Binding to produce P ===!

            rib = rib - 1
            Cm = Cm - 1
            Crib1 = Crib1 + 1
            !san2 = san2 + 1
            !WRITE(2,*)time,san2

          ENDIF
          
          IF ( i2 == 12 ) THEN

            !=== P Production ===!

            Crib1 = Crib1 - 1
            P = P + 1
            rib = rib + 1
            Cm = Cm + 1

          ENDIF
          
          IF ( i2 == 13 ) THEN

            !=== LS mRNA Ribosome Binding ===!

            rib = rib - 1
            LSm = LSm - 1
            LSrib = LSrib + 1

          ENDIF
          
          IF ( i2 == 14 ) THEN

            !=== L Production ===!

            LSrib = LSrib - 1
            L = L + 1
            rib = rib + 1
            LSm = LSm + 1

          ENDIF
          
          IF ( i2 == 15 ) THEN

            !=== S mRNA Ribosome Binding to produce M ===!

            rib = rib - 1
            Sm = Sm - 1
            Srib1 = Srib1 + 1

          ENDIF
          
          IF ( i2 == 16 ) THEN

            !=== M Production ===!

            Srib1 = Srib1 - 1
            M = M + 1
            rib = rib + 1
            Sm = Sm + 1

          ENDIF
          
          IF ( i2 == 17 ) THEN

            !=== X mRNA Ribosome Binding ===!

            rib = rib - 1
            Xm = Xm - 1
            Xrib = Xrib + 1

          ENDIF
          
          IF ( i2 == 18 ) THEN

            !=== X Production ===!

            Xrib = Xrib - 1
            X = X + 1
            rib = rib + 1
            Xm = Xm + 1

          ENDIF
          
          IF ( i2 == 19 ) THEN

            !=== PreC mRNA Ribosome Binding ===!

            rib = rib - 1
            PCm = PCm - 1
            PCrib = PCrib + 1

          ENDIF
          
          IF ( i2 == 20 ) THEN

            !=== PC Production ===!

            PCrib = PCrib - 1
            PC = PC + 1
            rib = rib + 1
            PCm = PCm + 1

          ENDIF
          
          IF ( i2 == 21 ) THEN

            !=== C mRNA Degradation ===!

            Cm = Cm - 1
            !san7 = san7 + 1

          ENDIF
          
          IF ( i2 == 22 ) THEN

            !=== LS mRNA Degradation ===!

            LSm = LSm - 1

          ENDIF
          
          IF ( i2 == 23 ) THEN

            !=== S mRNA Degradation ===!

            Sm = Sm - 1

          ENDIF
          
          IF ( i2 == 24 ) THEN

            !=== X mRNA Degradation ===!

            Xm = Xm - 1

          ENDIF
          
          IF ( i2 == 25 ) THEN

            !=== PreC mRNA Degradation ===!

            PCm = PCm - 1

          ENDIF
          
          IF ( i2 == 26 ) THEN

            !=== P Degradation ===!

            P = P - 1

          ENDIF
          
          IF ( i2 == 27 ) THEN

            !=== C Degradation ===!

            R(28) = R(28) - 1

          ENDIF
          
          IF ( i2 == 28 ) THEN

            !=== L Degradation ===!

            L = L - 1

          ENDIF
          
          IF ( i2 == 29 ) THEN

            !=== M Degradation ===!

            M = M - 1

          ENDIF
          
          IF ( i2 == 30 ) THEN

            !=== S Degradation ===!

            S = S - 1

          ENDIF
          
          IF ( i2 == 31 ) THEN

            !=== X Degradation ===!

            X = X - 1

          ENDIF
          
          IF ( i2 == 32 ) THEN

            !=== PC Degradation ===!

            PC = PC - 1

          ENDIF
          
          IF ( i2 == 33 ) THEN

            !=== Cm RNA ribosome Binding to produce C ===!

            rib = rib - 1
            Cm = Cm - 1
            Crib2 = Crib2 + 1
            !san1 = san1 + 1
            !WRITE(1,*)time,san1

          ENDIF
          
          IF ( i2 == 34 ) THEN

            !=== C production ===!

            Crib2 = Crib2 - 1
            R(28) = R(28) + 1
            rib = rib + 1
            Cm = Cm + 1

          ENDIF
          
          IF ( i2 == 35 ) THEN

            !=== P and C mRNA Binding ===!

            Cm = Cm - 1
            P = P - 1
            R(1) = R(1) + 1
            !san4 = san4 + 1
            !WRITE(4,*)time,san4

          ENDIF
          
          IF (i2 >= 36 .and. i2 <= 89) THEN
          
             !=== RNP Binding to C or D ===!
          
             R = R + V(:,i2 - 35)
             
          ENDIF
          
          IF (i2 >= 90 .and. i2 <= 143) THEN
          
             !=== Dissociation of C or D from RNP ===!
          
             R = R - V(:,i2 - 89)
             
          ENDIF
          
          IF ( i2 == 144 ) THEN

            !=== RNP with 4 C Production ===!

            R(14) = R(14) - 1
            R(28) = R(28) - 1
            pgelong(1) = pgelong(1) + 1

          ENDIF
          
          IF (i2 >= 145 .and. i2 <= 260) THEN
            
            !=== RNP with i2-140 C Production ===!
            
            R(28) = R(28) - 1
            pgelong(i2 - 144) = pgelong(i2 - 144) - 1
            pgelong(i2 - 143) = pgelong(i2 - 143) + 1
            
          ENDIF
          
          IF ( i2 == 261 ) THEN

            !=== rcNC Production ===!

            pgelong(117) = pgelong(117) - 1
            rcNC = rcNC + 1

          ENDIF
          
          IF ( i2 == 262 ) THEN

            !=== rcNC Removal ===!

            rcNC = rcNC - 1

          ENDIF
          
          IF ( i2 == 263 ) THEN

            !=== rc DNA Recycling ===!

            rcNC = rcNC - 1
            rc = rc + 1

          ENDIF
          
          IF ( i2 == 264 ) THEN

            !=== Complete Virion Release ===!

            rcNC = rcNC - 1
            S = S - 80
            M = M - 20
            L = L - 20
            rcV = rcV + 1

          ENDIF
          
          IF ( i2 == 265 ) THEN

            !=== SVP Filaments Production ===!

            S = S - 80
            M = M - 20
            L = L - 20
            fil = fil + 1

          ENDIF
          
          IF ( i2 == 266 ) THEN

            !=== 2 C Binding for Empty Virion Production ===!

            R(28) = R(28) - 2
            emelong(1) = emelong(1) + 1

          ENDIF
          
          IF ( i2 >= 267 .and. i2 <= 384 ) THEN
            
            !=== Empty Virion Production ===!
            
            R(28) = R(28) - 1
            emelong(i2 - 266) = emelong(i2 - 266) - 1
            emelong(i2 - 265) = emelong(i2 - 265) + 1
            
          ENDIF
            
          IF ( i2 == 385 ) THEN

            !=== Empty NC Removal ===!

            emelong(119) = emelong(119) - 1

          ENDIF 
          
          IF ( i2 == 386 ) THEN

            !=== Empty Virion Release ===!

            S = S - 80
            M = M - 20
            L = L - 20
            emelong(119) = emelong(119) - 1
            emV = emV + 1

          ENDIF
          
          IF ( i2 == 387 ) THEN

            !=== rcDNA releasing in nucleus ===!

            rcNCe = rcNCe - 1
            rc = rc + 1

          ENDIF 
          
          IF ( i2 == 388 ) THEN

            !=== rcNCe degradation ===!

            rcNCe = rcNCe - 1

          ENDIF 
          
          IF ( i2 == 389 ) THEN

            !=== rcDNA degradation ===!

            rc = rc - 1

          ENDIF
          
          IF ( i2 == 390 ) THEN

            !=== Sm RNA ribosome binding to produce S  ===!

            rib = rib - 1
            Sm = Sm - 1
            Srib2 = Srib2 + 1

          ENDIF 
          
          IF ( i2 == 391 ) THEN

            !=== S production ===!

            Srib2 = Srib2 - 1
            S = S + 1
            rib = rib + 1
            Sm = Sm + 1

          ENDIF 
          
          IF ( i2 == 392 ) THEN

            !=== pgNC degradation ===!

            pgelong(117) = pgelong(117) - 1

          ENDIF 
          
          IF ( i2 == 393 ) THEN

            !=== RNA virion releasing ===!

            pgelong(117) = pgelong(117) - 1
            S = S - 80
            M = M - 20
            L = L - 20
            RNAV = RNAV + 1

          ENDIF 
          
          IF ( i2 == 394 ) THEN

            !=== Sphere elongation ===!

            S = S - 2
            sphere(1) = sphere(1) + 1

          ENDIF 
          
          IF ( i2 >= 395 .and. i2 <= 440 ) THEN
            
            !=== Sphere Production ===!
            
            S = S - 1
            sphere(i2 - 394) = sphere(i2 - 394) - 1
            sphere(i2 - 393) = sphere(i2 - 393) + 1
            
          ENDIF
          
          IF ( i2 == 441 ) THEN

            !=== Sphere elongation ===!

            sphere(47) = sphere(47) - 1
            sph = sph + 1

          ENDIF
          
          !!!!! Drugs !!!!!
          
          IF ( i2 == 442 ) THEN

            !=== IFN binding ===!

            IFN = IFN - 1
            ccc = ccc - 1
            cccIFN = cccIFN + 1

          ENDIF
          
          IF ( i2 == 443 ) THEN

            !=== IFN unbinding ===!

            IFN = IFN + 1
            ccc = ccc + 1
            cccIFN = cccIFN - 1

          ENDIF
          
          IF ( i2 == 444 ) THEN

            !=== GA binding ===!

            GA = GA - 1
            P = P - 1
            PGA = PGA + 1

          ENDIF
          
          IF ( i2 == 445 ) THEN

            !=== GA unbinding ===!

            GA = GA + 1
            P = P + 1
            PGA = PGA - 1

          ENDIF
          
          IF ( i2 == 587 ) THEN

            !=== Epsilon binding ===!

            eps = eps - 1
            Cm = Cm - 1
            Ceps = Ceps + 1
            !san2 = san2 + 1

          ENDIF
          
          IF ( i2 == 588 ) THEN

            !=== Epsilon unbinding ===!

            eps = eps + 1
            Cm = Cm + 1
            Ceps = Ceps - 1
            !san3 = san3 + 1

          ENDIF
          
          IF ( i2 == 446 ) THEN

            !=== NA binding ===!

            NA = NA - 1
            pgelong(117) = pgelong(117) - 1
            pgNCNA = pgNCNA + 1

          ENDIF
          
          IF ( i2 == 447 ) THEN

            !=== NA unbinding ===!

            NA = NA + 1
            pgelong(117) = pgelong(117) + 1
            pgNCNA = pgNCNA - 1

          ENDIF
          
          IF ( i2 == 448 ) THEN

            !=== RNA virion releasing ===!

            pgNCNA = pgNCNA - 1
            S = S - 80
            M = M - 20
            L = L - 20
            RNAV = RNAV + 1

          ENDIF 
          
          IF ( i2 == 449 ) THEN

            !=== CAM binding to C ===!

            CAM = CAM - 1
            R(28) = R(28) - 1
            CC = CC + 1

          ENDIF
          
          IF ( i2 == 450 ) THEN

            !=== CAM unbinding from C ===!

            CAM = CAM + 1
            R(28) = R(28) + 1
            CC = CC - 1

          ENDIF
          
          IF ( i2 == 451 ) THEN

            !=== CC binding to CC ===!

            CC = CC - 2
            CCAM(1) = CCAM(1) + 1

          ENDIF
          
          IF ( i2 >= 452 .and. i2 <= 569 ) THEN
            
            !=== emNC:CAM Production ===!
            
            CC = CC - 1
            CCAM(i2 - 451) = CCAM(i2 - 451) - 1
            CCAM(i2 - 450) = CCAM(i2 - 450) + 1
            
          ENDIF
          
          IF ( i2 == 570 ) THEN

            !=== emNC:CAM degradation ===!

            CCAM(119) = CCAM(119) - 1

          ENDIF
          
          IF ( i2 == 571 ) THEN

            !=== Empty Virion Release ===!

            S = S - 80
            M = M - 20
            L = L - 20
            CCAM(119) = CCAM(119) - 1
            emVc = emVc + 1

          ENDIF
          
          IF ( i2 == 572 ) THEN

            !=== CAM binding to rcNCe ===!

            CAM = CAM - 1
            rcNCe = rcNCe - 1
            rcCAMe = rcCAMe + 1

          ENDIF
          
          IF ( i2 == 573 ) THEN

            !=== CAM unbinding from rcNCe ===!

            CAM = CAM + 1
            rcNCe = rcNCe + 1
            rcCAMe = rcCAMe - 1

          ENDIF
          
          IF ( i2 == 574 ) THEN

            !=== CAM binding to rcNC ===!

            CAM = CAM - 1
            rcNC = rcNC - 1
            rcCAM = rcCAM + 1

          ENDIF
          
          IF ( i2 == 575 ) THEN

            !=== CAM unbinding from rcNC ===!

            CAM = CAM + 1
            rcNC = rcNC + 1
            rcCAM = rcCAM - 1

          ENDIF
          
          IF ( i2 == 576 ) THEN

            !=== rcNCe:CAM disassembly ===!

            CAM = CAM + 1
            rcCAMe = rcCAMe - 1

          ENDIF
          
          IF ( i2 == 577 ) THEN

            !=== rcNC:CAM disassembly ===!

            CAM = CAM + 1
            rcCAM = rcCAM - 1

          ENDIF
          
          IF ( i2 == 578 ) THEN

            !=== DD binding to rcDNA ===!

            rc = rc - 1
            DD = DD - 1
            rcDD = rcDD + 1

          ENDIF
          
          IF ( i2 == 579 ) THEN

            !=== DD unbinding from rcDNA ===!

            rc = rc + 1
            DD = DD + 1
            rcDD = rcDD - 1

          ENDIF
          
          IF ( i2 == 580 ) THEN

            !=== iRNA binding to pgRNA ===!

            Cm = Cm - 1
            iRNA = iRNA - 1
            CiRNA = CiRNA + 1

          ENDIF
          
          IF ( i2 == 581 ) THEN

            !=== pgRNA cleavage ===!

            iRNA = iRNA + 1
            CiRNA = CiRNA - 1

          ENDIF
          
          IF ( i2 == 586 ) THEN

            !=== iRNA unbinding from pgRNA ===!

            Cm = Cm + 1
            iRNA = iRNA + 1
            CiRNA = CiRNA - 1

          ENDIF

        ENDDO

 1      CONTINUE


        write(*,*)nrun,time,rcV,sph,fil,emV,emVc,RNAV,L!,san1,san2,san3,san4,san5,san6,san7

        ENDDO
      
      END PROGRAM
          
      !!!!! Computing Propensity Function for Binding/Unbinding of C and D to/from RNP !!!!!
      
      SUBROUTINE PROPENSITY (j1,ii1,i1,RR1,rate1,a1,V1)
        
        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: j1,ii1,i1,RR1(29)
        DOUBLE PRECISION, INTENT(IN) :: rate1(12)
        INTEGER, INTENT(INOUT) :: V1(29,54)
        DOUBLE PRECISION, INTENT(INOUT) :: a1(588)
        INTEGER :: kk1,kk2,kk3,kk4,kk5
        
        kk1 = ii1 - j1
        kk2 = i1 + 54
        kk3 = i1 - 35
        kk4 = j1 + 1
        kk5 = ii1 + 1
        
        
        IF ( kk1 == 1 ) THEN  !Binding C to PS1
        
           a1(i1) = rate1(1) * DBLE(RR1(kk4)) * DBLE(RR1(28))
           
           a1(kk2) = rate1(7) * DBLE(RR1(kk5))
           
           V1(kk4,kk3) = -1
           
           V1(kk5,kk3) = 1
           
           V1(28,kk3) = -1
           
        ENDIF
        
        IF ( kk1 == 2 ) THEN  !Binding D to PS1
           
           a1(i1) = rate1(4) * DBLE(RR1(kk4)) * DBLE(RR1(29))
           
           a1(kk2) = rate1(10) * DBLE(RR1(kk5))
           
           V1(kk4,kk3) = -1
           
           V1(kk5,kk3) = 1
           
           V1(29,kk3) = -1
           
        ENDIF
        
        IF ( kk1 == 3 ) THEN  !Binding C to PS2
        
           a1(i1) = rate1(2) * DBLE(RR1(kk4)) * DBLE(RR1(28))
           
           a1(kk2) = rate1(8) * DBLE(RR1(kk5))
           
           V1(kk4,kk3) = -1
           
           V1(kk5,kk3) = 1
           
           V1(28,kk3) = -1
           
        ENDIF
        
        IF ( kk1 == 6 ) THEN  !Binding D to PS2
        
           a1(i1) = rate1(5) * DBLE(RR1(kk4)) * DBLE(RR1(29))
           
           a1(kk2) = rate1(11) * DBLE(RR1(kk5))
           
           V1(kk4,kk3) = -1
           
           V1(kk5,kk3) = 1
           
           V1(29,kk3) = -1
           
        ENDIF
        
        IF ( kk1 == 9 ) THEN  !Binding C to PS3
        
           a1(i1) = rate1(3) * DBLE(RR1(kk4)) * DBLE(RR1(28))
           
           a1(kk2) = rate1(9) * DBLE(RR1(kk5))
           
           V1(kk4,kk3) = -1
           
           V1(kk5,kk3) = 1
           
           V1(28,kk3) = -1
           
        ENDIF
        
        IF ( kk1 == 18 ) THEN  !Binding D to PS3
        
           a1(i1) = rate1(6) * DBLE(RR1(kk4)) * DBLE(RR1(29))
           
           a1(kk2) = rate1(12) * DBLE(RR1(kk5))
           
           V1(kk4,kk3) = -1
           
           V1(kk5,kk3) = 1
           
           V1(29,kk3) = -1
           
        ENDIF
        
       RETURN
           
      END SUBROUTINE
      
      !!!!! Creating a Random Number !!!!!
      
      DOUBLE PRECISION FUNCTION RANDOM (ISEED)

        IMPLICIT NONE

        !=== ARGUMENTS ===!

        INTEGER, INTENT(INOUT) :: iseed

        !=== VARIABLES ===!

        INTEGER :: hi,lo,test

        INTEGER, PARAMETER :: ae = 16807
        INTEGER, PARAMETER :: m = 2147483647
        INTEGER, PARAMETER :: q = 127773
        INTEGER, PARAMETER :: re = 2836


        hi = INT(iseed/q)
        lo = MODULO(iseed,q)

        test = ae * lo - re * hi

        IF ( test > 0 ) THEN

          iseed = test

        ELSE

          iseed = test + m

        ENDIF

        RANDOM = DBLE(iseed) / DBLE(m)

        RETURN

      END FUNCTION RANDOM
