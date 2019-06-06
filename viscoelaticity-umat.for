c    
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,
     1 NPT,LAYER,KSPT)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
C
      DO K1=1,111
      STATEV(K1)=0.0
      END DO     
C
      RETURN
      END
C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC) 
C
      INCLUDE 'ABA_PARAM.INC' 
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),PROPS(NPROPS),COORDS(3),DROT(3,3),
     3 DFGRD0(3,3),DFGRD1(3,3)
C
      DIMENSION EG(6),DEG(6),EVN1(6),EVN2(6),EVN3(6),EVN4(6),EVN5(6),
     &EVN6(6), DEN1(6),DEN2(6),DEN3(6),DEN4(6),DEN5(6),DEN6(6),EN(6),
     &TEN1(6),TEN2(6),TEN3(6),TEN4(6),TEN5(6),TEN6(6), DE1(6),DE2(6),
     &DE3(6),DE4(6),DE5(6),DE6(6),DGSTRES(6)
      DOUBLE PRECISION EV0,DEV0,EV,DEV,KVN1,KVN2,KVN3,KVN4,KVN5,KVN6,
     &DVN1,DVN2,DVN3,DVN4,DVN5,DVN6,DV1,DV2,DV3,DV4,DV5,DV6,DKSTRES,GT,
     &KT,K0,G0,VN,TVN1,TVN2,TVN3,TVN4,TVN5,TVN6
    
C
C     prony series----g0 elastic modulus at infintiy, tao1~tao6 is the  relaxation time 
C  
      G0=PROPS(1)
      ENU=PROPS(2)
      G1=PROPS(3)
      G2=PROPS(4)
      G3=PROPS(5)
      G4=PROPS(6)
      G5=PROPS(7)
      G6=PROPS(8)
      TAO1=PROPS(9)
      TAO2=PROPS(10)
      TAO3=PROPS(11)
      TAO4=PROPS(12)
      TAO5=PROPS(13)
      TAO6=PROPS(14)
      K0=G0*2*(1+ENU)/3/(1-2*ENU)
C      
C
C     calucation of deviator strain and  volumetric strain
      EV0=0.0
      DEV0=0.0
      DO I=1,NDI
          EV0=EV0+STRAN(I)
          DEV0=DEV0+DSTRAN(I)
      END DO     
      DEV=DEV0/3
      EV=EV0/3
      DO I=1,NDI
          EG(I)=STRAN(I)-EV
          DEG(I)=DSTRAN(I)-DEV
      END DO
      DO I=NDI+1,NTENS
          EG(I)=STRAN(I)/2
          DEG(I)=DSTRAN(I)/2
      END DO
       
C     initial values
      IF(TIME(1).EQ.0.0)THEN
          DO I=1,NTENS
              STATEV(I)=EG(I)            
          END DO
          STATEV(NTENS*13+1)=EV
      END IF      
C
C      
C     calucation of  deviator stress
C
      DO I=1,NTENS
          EN(I)=STATEV(I)
      END DO
      DO I=1,NTENS
          EVN1(I)=STATEV(NTENS*1+I)
      END DO
      DO I=1,NTENS
          EVN2(I)=STATEV(NTENS*2+I)
      END DO
      DO I=1,NTENS
          EVN3(I)=STATEV(NTENS*3+I)
      END DO
      DO I=1,NTENS
          EVN4(I)=STATEV(NTENS*4+I)
      END DO
      DO I=1,NTENS
          EVN5(I)=STATEV(NTENS*5+I)
      END DO
      DO I=1,NTENS
          EVN6(I)=STATEV(NTENS*6+I)
      END DO
      DO I=1,NTENS
          TEN1(I)=STATEV(NTENS*7+I)
          TEN2(I)=STATEV(NTENS*8+I)
          TEN3(I)=STATEV(NTENS*9+I)
          TEN4(I)=STATEV(NTENS*10+I)
          TEN5(I)=STATEV(NTENS*11+I)
          TEN6(I)=STATEV(NTENS*12+I)
      END DO
c      
      DO I=1,NTENS
          DEN1(I)=DEG(I)/DTIME*(DTIME-TAO1*
     &            (1-EXP(-DTIME/TAO1)))
          DEN2(I)=DEG(I)/DTIME*(DTIME-TAO2*
     &            (1-EXP(-DTIME/TAO2)))
          DEN3(I)=DEG(I)/DTIME*(DTIME-TAO3*
     &            (1-EXP(-DTIME/TAO3)))
          DEN4(I)=DEG(I)/DTIME*(DTIME-TAO4*
     &            (1-EXP(-DTIME/TAO4)))
          DEN5(I)=DEG(I)/DTIME*(DTIME-TAO5*
     &            (1-EXP(-DTIME/TAO5)))
          DEN6(I)=DEG(I)/DTIME*(DTIME-TAO6*
     &            (1-EXP(-DTIME/TAO6)))
      END DO
C        
      DO I=1,NTENS
          DE1(I)=(1-EXP(-DTIME/TAO1))*EN(I)
     &            -(1-EXP(-DTIME/TAO1))*EVN1(I)+DEN1(I)
          DE2(I)=(1-EXP(-DTIME/TAO2))*EN(I)-(1-EXP(-DTIME/TAO2))*EVN2(I)
     &            +DEN2(I)
          DE3(I)=(1-EXP(-DTIME/TAO3))*EN(I)-
     &            (1-EXP(-DTIME/TAO3))*EVN3(I)+DEN3(I)
          DE4(I)=(1-EXP(-DTIME/TAO4))*EN(I)
     &            -(1-EXP(-DTIME/TAO4))*EVN4(I)+DEN4(I)
          DE5(I)=(1-EXP(-DTIME/TAO5))*EN(I)
     &            -(1-EXP(-DTIME/TAO5))*EVN5(I)+DEN5(I)
          DE6(I)=(1-EXP(-DTIME/TAO6))*EN(I)
     &            -(1-EXP(-DTIME/TAO6))*EVN6(I)+DEN6(I)
      END DO
C         increment of  deviator stress
      DO I=1,NTENS
          DGSTRES(I)=2*G0*(DEG(I)-G1*DE1(I)
     &               -G2*DE2(I)-G3*DE3(I)-G4*DE4(I)-G5*DE5(I)-G6*DE6(I))
      END DO
C      update stress
      DO I=1,NTENS
          STRESS(I)=STRESS(I)+DGSTRES(I)
      END DO
C      
C         calucation of volumetric stress

          VN=STATEV(NTENS*13+1)
          KVN1=STATEV(NTENS*13+2)
          KVN2=STATEV(NTENS*13+3)
          KVN3=STATEV(NTENS*13+4)
          KVN4=STATEV(NTENS*13+5)
          KVN5=STATEV(NTENS*13+6)
          KVN6=STATEV(NTENS*13+7)
          TVN1=STATEV(NTENS*13+8)
          TVN2=STATEV(NTENS*13+9)
          TVN3=STATEV(NTENS*13+10)
          TVN4=STATEV(NTENS*13+11)
          TVN5=STATEV(NTENS*13+12)
          TVN6=STATEV(NTENS*13+13)
          DVN1=(DEV)/DTIME*
     &            (DTIME-TAO1*(1-EXP(-DTIME/TAO1)))
          DVN2=(DEV)/DTIME*
     &            (DTIME-TAO2*(1-EXP(-DTIME/TAO2)))
          DVN3=(DEV)/DTIME*
     &            (DTIME-TAO3*(1-EXP(-DTIME/TAO3)))
          DVN4=(DEV)/DTIME*
     &            (DTIME-TAO4*(1-EXP(-DTIME/TAO4)))
          DVN5=(DEV)/DTIME*
     &            (DTIME-TAO5*(1-EXP(-DTIME/TAO5)))
          DVN6=(DEV)/DTIME*
     &            (DTIME-TAO6*(1-EXP(-DTIME/TAO6)))
C         
          DV1=(1-EXP(-DTIME/TAO1))*VN
     &           -(1-EXP(-DTIME/TAO1))*KVN1+DVN1
          DV2=(1-EXP(-DTIME/TAO2))*VN
     &           -(1-EXP(-DTIME/TAO2))*KVN2+DVN2           
          DV3=(1-EXP(-DTIME/TAO3))*VN
     &           -(1-EXP(-DTIME/TAO3))*KVN3+DVN3
          DV4=(1-EXP(-DTIME/TAO4))*VN
     &           -(1-EXP(-DTIME/TAO4))*KVN4+DVN4
          DV5=(1-EXP(-DTIME/TAO5))*VN
     &           -(1-EXP(-DTIME/TAO5))*KVN5+DVN5
          DV6=(1-EXP(-DTIME/TAO6))*VN
     &           -(1-EXP(-DTIME/TAO6))*KVN6+DVN6
c
          DKSTRES=3*K0*(DEV-
     &    G1*DV1-G2*DV2-G3*DV3-G4*DV4-G5*DV5-G6*DV6)
C      update stress
      DO I=1,NDI
          STRESS(I)=STRESS(I)+DKSTRES
      END DO     
C
C      update  Jacobian matrix
C
      GT=G0*(1-G1*(1-TAO1/DTIME*(1-EXP(-DTIME/TAO1)))-
     &   G2*(1-TAO2/DTIME*(1-EXP(-DTIME/TAO2)))-
     & G3*(1-TAO3/DTIME*(1-EXP(-DTIME/TAO3)))
     & -G4*(1-TAO4/DTIME*(1-EXP(-DTIME/TAO4)))-
     & G5*(1-TAO5/DTIME*(1-EXP(-DTIME/TAO5)))-
     & G6*(1-TAO6/DTIME*(1-EXP(-DTIME/TAO6))))
      KT=K0*(1-G1*(1-TAO1/DTIME*(1-EXP(-DTIME/TAO1)))-
     &   G2*(1-TAO2/DTIME*(1-EXP(-DTIME/TAO2)))-
     & G3*(1-TAO3/DTIME*(1-EXP(-DTIME/TAO3)))
     & -G4*(1-TAO4/DTIME*(1-EXP(-DTIME/TAO4)))-
     & G5*(1-TAO5/DTIME*(1-EXP(-DTIME/TAO5)))-
     & G6*(1-TAO6/DTIME*(1-EXP(-DTIME/TAO6))))
      DO I1=1,NTENS
          DO I2=1,NTENS
              DDSDDE(I1,I2)=0.0
          END DO
      END DO
      DO I1=1,NDI
          DO I2=I,NDI
              DDSDDE(I1,I2)=KT-2*GT/3
          END DO
      END DO
      DO I1=1,NDI
          DDSDDE(I1,I1)=4*GT/3+KT
      END DO
      DO I=NDI+1,NTENS
          DDSDDE(I,I)=GT
      END DO

C
C      Update  state variables of deviator strain
      DO I=1,NTENS
          STATEV(I+3)=EN(I)+DEG(I)
      END DO
      DO I=1,NTENS
          STATEV(NTENS*1+I)=EVN1(I)+DE1(I)
      END DO
      DO I=1,NTENS
          STATEV(NTENS*2+I)=EVN2(I)+DE2(I)
      END DO
      DO I=1,NTENS
          STATEV(NTENS*3+I)=EVN3(I)+DE3(I)
      END DO
      DO I=1,NTENS
          STATEV(NTENS*4+I)=EVN4(I)+DE4(I)
      END DO 
      DO I=1,NTENS
          STATEV(NTENS*5+I)=EVN5(I)+DE5(I)
      END DO
      DO I=1,NTENS
          STATEV(NTENS*6+I)=EVN6(I)+DE6(I)
      END DO
      DO I=1,NTENS
          STATEV(NTENS*7+I)=TEN1(I)+
     &                    DEG(I)/DTIME*(DTIME-TAO1*(1-EXP(-DTIME/TAO1)))
          STATEV(NTENS*8+I)=TEN2(I)+
     &       DEG(I)/DTIME*(DTIME-TAO2*(1-EXP(-DTIME/TAO2)))
          STATEV(NTENS*9+I)=TEN3(I)+
     &       DEG(I)/DTIME*(DTIME-TAO3*(1-EXP(-DTIME/TAO3)))
          STATEV(NTENS*10+I)=TEN4(I)+
     &      DEG(I)/DTIME*(DTIME-TAO4*(1-EXP(-DTIME/TAO4)))
          STATEV(NTENS*11+I)=TEN5(I)+
     &       DEG(I)/DTIME*(DTIME-TAO5*(1-EXP(-DTIME/TAO5)))
          STATEV(NTENS*12+I)=TEN6(I)+
     &        DEG(I)/DTIME*(DTIME-TAO6*(1-EXP(-DTIME/TAO6)))
      END DO
          
C      
C     Update  state variables of volumetric strain
          VN=VN+DEV
          KVN1=KVN1+DV1
          KVN2=KVN2+DV2
          KVN3=KVN3+DV3
          KVN4=KVN4+DV4
          KVN5=KVN5+DV5
          KVN6=KVN6+DV6
          TVN1=TVN1+DEV/DTIME*(DTIME-TAO1*(1-EXP(-DTIME/TAO1)))
          TVN2=TVN2+DEV/DTIME*(DTIME-TAO2*(1-EXP(-DTIME/TAO2)))
          TVN3=TVN3+DEV/DTIME*(DTIME-TAO3*(1-EXP(-DTIME/TAO3)))
          TVN4=TVN4+DEV/DTIME*(DTIME-TAO4*(1-EXP(-DTIME/TAO4)))
          TVN5=TVN5+DEV/DTIME*(DTIME-TAO5*(1-EXP(-DTIME/TAO5)))
          TVN6=TVN1+DEV/DTIME*(DTIME-TAO6*(1-EXP(-DTIME/TAO6)))

          STATEV(NTENS*13+1)=VN
          STATEV(NTENS*13+2)=KVN1
          STATEV(NTENS*13+3)=KVN2
          STATEV(NTENS*13+4)=KVN3
          STATEV(NTENS*13+5)=KVN4
          STATEV(NTENS*13+6)=KVN5
          STATEV(NTENS*13+7)=KVN6
          STATEV(NTENS*13+8)=TVN1
          STATEV(NTENS*13+9)=TVN2
          STATEV(NTENS*13+10)=TVN3
          STATEV(NTENS*13+11)=TVN4
          STATEV(NTENS*13+12)=TVN5
          STATEV(NTENS*13+13)=TVN6
C     
      RETURN
      END

