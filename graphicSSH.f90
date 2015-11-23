!INSTITUTE OF PHYSICS - UNIVERSITY OF BRASILIA
!LUIZ A. RIBEIRO JR.
!OCTOBER 07 2012
!PROGRAM TO PLOT DE SSH.F RESULTS USING GNUPLOT

PROGRAM graphicSSH

IMPLICIT NONE
PARAMETER steps=251,N=200,dt=0.01d0,iwritemode=100,IFIELDON=0,IFIELDOFF=50000

!OUTPUT FILES:
!FORT.10 --> LATICE DISPLACEMENT (BOND ORDER)
!FORT.11 --> ELECTRONIC ENERGY SPECTRUM UP 
!FORT.12 --> ELECTRONIC ENERGY SPECTRUM DOWN 
!FORT.13 --> SYSTEM ENERGIES (ELECTRONIC,LATTICE,KINETIC,TOTAL) 
!FORT.14 --> CHARGE DENSITY
!FORT.15 --> SPIN DENSITY
!FORT.16 --> SITE SQUARE VELOCITY (FOR TEMPERATURE CALCULATIONS) 
!FORT.17 --> STITE VELOCITY
!FORT.18 --> ELECTRIC FIELD
!FORT.19 --> OCCUPATION NUMBER FOR ELECTRONS UP AND DOWN

!CALL temperature   (steps,N,dt,iwritemode)
CALL occupation    (steps,N,dt,iwritemode)
CALL bondlength    (steps,N,dt,iwritemode)
CALL chargedensity (steps,N,dt,iwritemode)
!CALL spindensity   (steps,N,dt,iwritemode)
CALL energylevels  (steps,N,dt,iwritemode)
!CALL electricfield (IFIELDOFF,IFIELDON,dt)
!CALL velocities    (steps,N,dt,iwritemode) 


END PROGRAM 

!***********************************************************************************************************
!***********************************************************************************************************
SUBROUTINE temperature(steps,N,dt,iwritemode)
IMPLICIT NONE
INTEGER      :: I,J,K,L
INTEGER      :: steps,N,iwritemode
REAL(KIND=8) :: v2(steps,N),sum1(N),time(steps),auxsum3(steps),y(steps,N),auxy(steps)
REAL(KIND=8) :: dt,sum2,sum3,factor,fcnstK

!TEMPERATURE CALCULATIONS M/Kb = 4K/Kb = factor
factor=9.75d5

OPEN(10,STATUS='OLD')
DO I=1,steps
  READ(10,*)
  READ(10,*)
  READ(10,*) (y(I,J),J=1,N)
END DO
CLOSE(10,STATUS='KEEP')

OPEN(16,STATUS='OLD')
DO I=1,steps
  READ(16,*)
  READ(16,*)
  READ(16,*) (v2(I,J),J=1,N)
END DO
CLOSE(16,STATUS='KEEP')

DO L=1,N
  sum1(L) = 0.0d0
END DO

sum2 = 0.0d0

DO I=1,steps
  time(I) = FLOAT(I-1)*4.0d0*dt*iwritemode
END DO

OPEN(1,FILE='tempxt.dat',STATUS='UNKNOWN')
DO I=1,steps
  DO J=1,N
    DO K=1,I
      sum1(J) = sum1(J)+v2(K,J)   
    END DO
    sum1(J) = sum1(J)/FLOAT(I)
  END DO
  DO L=1,N
    sum2 = sum2 + sum1(L)*factor
    sum3 = sum3 + v2(I,L)*factor
  END DO
  sum2 = sum2/FLOAT(N)
  sum3 = sum3/FLOAT(N)
  auxsum3(I) = sum3
  WRITE(1,*)time(I),sum2,sum3
  DO L=1,N
    sum1(L) = 0.0d0
  END DO
  sum2    = 0.0d0
  sum3    = 0.0d0
END DO
CLOSE(1,STATUS='KEEP')

OPEN(3,FILE='tempxt.gpt',STATUS='UNKNOWN')
WRITE(3,*)'#'
WRITE(3,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(3,*)'#LUIZ A. RIBEIRO JR'
WRITE(3,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(3,*)'#'
WRITE(3,*)
WRITE(3,*)'#' 
WRITE(3,*)'#THERMALIZATION CURVES FOR SSH.F PROGRAM'
WRITE(3,*)'#'
WRITE(3,*)
WRITE(3,*)'set term postscript color enhanced font "Times,20"'
WRITE(3,*)"set out 'thermalization.eps'" 
WRITE(3,*)'set encoding iso_8859_1'
WRITE(3,*)'set ylabel "T[K]" font "Times-Bold,25"'
WRITE(3,*)'set xtics 100'
WRITE(3,*)'set xlabel "Time[fs]" font "Times-Bold,25"'
WRITE(3,*)'set style line 1  lt -1 lc -1 lw 2'
WRITE(3,*)'set style line 2  lt -1 lc  3 lw 2'
WRITE(3,*)'set yrange[-2:',INT(MAXVAL(auxsum3))+5,']'
WRITE(3,*)"plot 'tempxt.dat' us 1:2 w l ls 1 t '','tempxt.dat' us 1:3 w l ls 2 t ''"
WRITE(3,*)
WRITE(3,*)'set out'
CLOSE(3,STATUS='KEEP')

CALL SYSTEM('gnuplot tempxt.gpt')
CALL SYSTEM('ps2pdf thermalization.eps')

RETURN
END
!***********************************************************************************************************
!***********************************************************************************************************

!***********************************************************************************************************
!***********************************************************************************************************
SUBROUTINE occupation(steps,N,dt,iwritemode)
IMPLICIT NONE
INTEGER      :: I,J,K,L
INTEGER      :: steps,N,iwritemode
REAL(KIND=8) :: occU(steps,N),occD(steps,N),occUD(steps,N)
REAL(KIND=8) :: dt,time(steps)

OPEN(19,STATUS='OLD')
DO I=1,steps
  READ(19,*)
  READ(19,*)
  READ(19,*) (occU(I,J),J=1,N)
END DO
CLOSE(19,STATUS='KEEP')

OPEN(20,STATUS='OLD')
DO I=1,steps
  READ(20,*)
  READ(20,*)
  READ(20,*) (occD(I,J),J=1,N)
END DO
CLOSE(20,STATUS='KEEP')

OPEN(21,STATUS='OLD')
DO I=1,steps
  READ(21,*)
  READ(21,*)
  READ(21,*) (occUD(I,J),J=1,N)
END DO
CLOSE(21,STATUS='KEEP')

DO I=1,steps
  time(I) = FLOAT(I-1)*4.0d0*dt*iwritemode
END DO

OPEN(1,FILE='occU.dat',STATUS='UNKNOWN')
DO J=1,N
  WRITE(1,100)'#',J
  DO I=1,steps
     WRITE(1,*)time(I),occU(I,J)
  END DO
  WRITE(1,*)
  WRITE(1,*)
END DO
CLOSE(1,STATUS='KEEP')

OPEN(2,FILE='occD.dat',STATUS='UNKNOWN')
DO J=1,N
  WRITE(2,100)'#',J
  DO I=1,steps
     WRITE(2,*)time(I),occD(I,J)
  END DO
  WRITE(2,*)
  WRITE(2,*)
END DO
CLOSE(2,STATUS='KEEP')

OPEN(3,FILE='occUD.dat',STATUS='UNKNOWN')
DO J=1,N
  WRITE(3,100)'#',J
  DO I=1,steps
     WRITE(3,*)time(I),occUD(I,J)
  END DO
  WRITE(3,*)
  WRITE(3,*)
END DO
CLOSE(3,STATUS='KEEP')

OPEN(7,FILE='occ3D.dat',STATUS='UNKNOWN')
DO J=1,N
   IF(J .EQ. 99 .OR. J .EQ. 100 .OR. J .EQ. 101 .OR. J .EQ. 102)THEN
      DO I=1,steps
         WRITE(7,*)time(I),FLOAT(J),occUD(I,J)
      END DO
      WRITE(7,*)
   END IF
END DO
CLOSE(7,STATUS='KEEP')


!
!FOR SEE THE DESIRED LEVEL YOU MUST DO INDEX - 1 IN THE GNUPLOT FILE  
!

OPEN(4,FILE='occU.gpt',STATUS='UNKNOWN')
WRITE(4,*)'#'
WRITE(4,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(4,*)'#LUIZ A. RIBEIRO JR'
WRITE(4,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(4,*)'#'
WRITE(4,*)
WRITE(4,*)'#' 
WRITE(4,*)'#OCCUPATION NUMBER TIME EVOLUTION FOR UP ELECTRONS'
WRITE(4,*)'#'
WRITE(4,*)
WRITE(4,*)'set term postscript color enhanced font "Times,20"'
WRITE(4,*)"set out 'occU.eps'" 
WRITE(4,*)'set encoding iso_8859_1'
WRITE(4,*)'set ylabel "occupation number" font "Times-Bold,25"'
WRITE(4,*)'set xtics 100'
WRITE(4,*)'set xlabel "Time[fs]" font "Times-Bold,25"'
WRITE(4,*)'set style line 1  lt -1 lc -1 lw 2'
WRITE(4,*)'set style line 2  lt -1 lc  3 lw 2'
WRITE(4,*)'set yrange[-0.2:1.2]'
WRITE(4,*)"plot 'occU.dat' index 101 us 1:2 w l ls 1 t '' "
WRITE(4,*)
WRITE(4,*)'set out'
CLOSE(4,STATUS='KEEP')

OPEN(5,FILE='occD.gpt',STATUS='UNKNOWN')
WRITE(5,*)'#'
WRITE(5,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(5,*)'#LUIZ A. RIBEIRO JR'
WRITE(5,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(5,*)'#'
WRITE(5,*)
WRITE(5,*)'#' 
WRITE(5,*)'#OCCUPATION NUMBER TIME EVOLUTION FOR DOWN ELECTRONS'
WRITE(5,*)'#'
WRITE(5,*)
WRITE(5,*)'set term postscript color enhanced font "Times,20"'
WRITE(5,*)"set out 'occD.eps'" 
WRITE(5,*)'set encoding iso_8859_1'
WRITE(5,*)'set ylabel "occupation number" font "Times-Bold,25"'
WRITE(5,*)'set xtics 100'
WRITE(5,*)'set xlabel "Time[fs]" font "Times-Bold,25"'
WRITE(5,*)'set style line 1  lt -1 lc -1 lw 2'
WRITE(5,*)'set style line 2  lt -1 lc  3 lw 2'
WRITE(5,*)'set yrange[-0.2:1.2]'
WRITE(5,*)"plot 'occD.dat' index 100 us 1:2 w l ls 1 t '' "
WRITE(5,*)
WRITE(5,*)'set out'
CLOSE(5,STATUS='KEEP')

OPEN(6,FILE='occUD.gpt',STATUS='UNKNOWN')
WRITE(6,*)'#'
WRITE(6,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(6,*)'#LUIZ A. RIBEIRO JR'
WRITE(6,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(6,*)'#'
WRITE(6,*)
WRITE(6,*)'#' 
WRITE(6,*)'#OCCUPATION NUMBER TIME EVOLUTION FOR UP + DOWN ELECTRONS'
WRITE(6,*)'#'
WRITE(6,*)
WRITE(6,*)'set term postscript color enhanced font "Times,20"'
WRITE(6,*)"set out 'occUD.eps'" 
WRITE(6,*)'set encoding iso_8859_1'
WRITE(6,*)'set ylabel "occupation number" font "Times-Bold,25"'
WRITE(6,*)'set xtics 100'
WRITE(6,*)'set xlabel "Time[fs]" font "Times-Bold,25"'
WRITE(6,*)'set style line 1  lt -1 lc -1 lw 2'
WRITE(6,*)'set style line 2  lt -1 lc  3 lw 2'
WRITE(6,*)'set style line 3  lt -1 lc  2 lw 2'
WRITE(6,*)'set style line 4  lt -1 lc  5 lw 2'
WRITE(6,*)'set yrange[-0.2:2.2]'
WRITE(6,*)"plot 'occUD.dat' index 98 us 1:2 w l ls 1 t '',\"
WRITE(6,*)"'occUD.dat' index 99 us 1:2 w l ls 2 t '',\"
WRITE(6,*)"'occUD.dat' index 100 us 1:2 w l ls 3 t '',\"
WRITE(6,*)"'occUD.dat' index 101 us 1:2 w l ls 4 t ''"
WRITE(6,*)
WRITE(6,*)'set out'
CLOSE(6,STATUS='KEEP')

OPEN(8,FILE='occ3D.gpt',STATUS='UNKNOWN')
WRITE(8,*)'#'
WRITE(8,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(8,*)'#LUIZ A. RIBEIRO JR'
WRITE(8,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(8,*)'#'
WRITE(8,*)
WRITE(8,*)'#' 
WRITE(8,*)'#3D OCCUPATION NUMBER FOR SSH.F PROGRAM'
WRITE(8,*)'#'
WRITE(8,*)
WRITE(8,*)'set term postscript color enhanced font "Times,20"'
WRITE(8,*)"set out 'occ3D.eps'" 
WRITE(8,*)'set encoding iso_8859_1'
!WRITE(8,*)'set ylabel "Energy Level" offset 1,1,1 font "Times-Bold,25"'
!WRITE(8,*)'set xlabel "Time[fs]" offset 4,1,1 font "Times-Bold,25"'
WRITE(8,*)'set pm3d'
WRITE(8,*)'set samples 100'
WRITE(8,*)'set isosamples 100'
WRITE(8,*)'set xtics 100'
WRITE(8,*)'set ytics 1'
WRITE(8,*)'set ztics 1'
WRITE(8,*)'set cbtics 0.5'
WRITE(8,*)'set cbtics scale 0'
!WRITE(8,*)'set zlabel "ON" offset -2,2,8 font "Times-Bold,25"'
WRITE(8,*)'set xyplane 0.0'
WRITE(8,*)'set view 57,168'

!HOT SCALE PALETTE
!WRITE(8,*)"set palette model RGB defined (0'white',0.04'yellow',0.06'red',0.12'black')"

!WRITE(8,*)'set xrange[',INT(MINVAL(time)),':',INT(MAXVAL(time))+20,']'
WRITE(8,*)"set palette model RGB defined (0'blue',.35'cyan',.5'green',.65'yellow',1'red')"
!WRITE(8,*)'set yrange[0.8:',N+10,']'
!WRITE(8,*)'set yrange[98:101]'
WRITE(8,*)"splot 'occ3D.dat' w l lt -1 t '' "
WRITE(8,*)
WRITE(8,*)'set out'
CLOSE(8,STATUS='KEEP')


CALL SYSTEM('gnuplot occU.gpt')
CALL SYSTEM('ps2pdf occU.eps')

CALL SYSTEM('gnuplot occD.gpt')
CALL SYSTEM('ps2pdf occD.eps')

CALL SYSTEM('gnuplot occUD.gpt')
CALL SYSTEM('ps2pdf occUD.eps')

CALL SYSTEM('gnuplot occ3D.gpt')
CALL SYSTEM('ps2pdf occ3D.eps')


100 FORMAT(T4,A,I5)

RETURN
END
!***********************************************************************************************************
!***********************************************************************************************************

!***********************************************************************************************************
!***********************************************************************************************************
SUBROUTINE bondlength(steps,N,dt,iwritemode)
IMPLICIT NONE
INTEGER      :: I,J
INTEGER      :: steps,N,iwritemode
REAL(KIND=8) :: y(N),auxy(steps,N),time(steps)
REAL(KIND=8) :: dt

DO I=1,steps
  time(I) = FLOAT(I-1)*4.0d0*dt*iwritemode
END DO

DO I=1,steps
   DO J=1,N
      auxy(I,J)=0.0d0
   END DO
END DO

REWIND(10)
OPEN(10,STATUS='OLD')
DO J=1,steps
   READ(10,*)
   READ(10,*)
   READ(10,*)(y(I),I=1,N)
   auxy(J,1) = -(y(N)-2.0d0*y(1)+y(2))/4.0d0
   auxy(J,N) = (-1.0d0)**N*(y(N-1)-2.0d0*y(N)+y(1))/4.0d0
   DO I=2,N-1
      auxy(J,I) = (-1.0d0)**I*(y(I-1)-2.0d0*y(I)+y(I+1))/4.0d0
   END DO
END DO
CLOSE(10,STATUS='KEEP')

OPEN(UNIT=1,FILE='bl2D.dat',STATUS='UNKNOWN')
DO I=1,steps
   WRITE(1,110)'#',I
   DO J=1,N
      WRITE(1,100)J,auxy(I,J)
   END DO
   WRITE(1,*)
   WRITE(1,*)
END DO
CLOSE(1,STATUS='KEEP')

OPEN(2,FILE='bl3D.dat',STATUS='UNKNOWN')
DO J=1,N
   DO I=1,steps
      WRITE(2,*)time(I),FLOAT(J),auxy(I,J)
   END DO
   WRITE(2,*)
END DO
CLOSE(2,STATUS='KEEP')

OPEN(3,FILE='bl2D.gpt',STATUS='UNKNOWN')
WRITE(3,*)'#'
WRITE(3,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(3,*)'#LUIZ A. RIBEIRO JR'
WRITE(3,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(3,*)'#'
WRITE(3,*)
WRITE(3,*)'#' 
WRITE(3,*)'#2D BOND-LENTH CURVES FOR SSH.F PROGRAM'
WRITE(3,*)'#'
WRITE(3,*)
WRITE(3,*)'set term postscript color enhanced font "Times,20"'
WRITE(3,*)"set out 'bl2D.eps'" 
WRITE(3,*)'set encoding iso_8859_1'
WRITE(3,*)'set ylabel "~y{.5-}_n[\305]" font "Times-Bold,25"'
WRITE(3,*)'set xlabel "Site" font "Times-Bold,25"'
WRITE(3,*)'set xrange [0:',N+2,']'
WRITE(3,*)'set style line 1  lt -1 lc -1 lw 2'
WRITE(3,*)"plot 'bl2D.dat' index 1 us 1:2 w l ls 1 t ''"
WRITE(3,*)
WRITE(3,*)'set out'
CLOSE(3,STATUS='KEEP')

OPEN(4,FILE='bl3D.gpt',STATUS='UNKNOWN')
WRITE(4,*)'#'
WRITE(4,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(4,*)'#LUIZ A. RIBEIRO JR'
WRITE(4,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(4,*)'#'
WRITE(4,*)
WRITE(4,*)'#' 
WRITE(4,*)'#3D BOND-LENTH CURVES FOR SSH.F PROGRAM'
WRITE(4,*)'#'
WRITE(4,*)
WRITE(4,*)'set term postscript color enhanced font "Times,20"'
WRITE(4,*)"set out 'bl3D.eps'" 
WRITE(4,*)'set encoding iso_8859_1'
WRITE(4,*)'set ylabel "Site" font "Times-Bold,25"'
WRITE(4,*)'set xlabel "Time[fs]" font "Times-Bold,25"'
!WRITE(4,*)'set dgrid3d 300,300 gauss 2'
WRITE(4,*)'set pm3d'
WRITE(4,*)'set xtics 100'
WRITE(4,*)'set ytics 50'
WRITE(4,*)'set view map'
WRITE(4,*)'set cbtics 0.02'
WRITE(4,*)'set cbtics scale 0'
WRITE(4,*)'unset surface'
!WRITE(4,*)'set label "~y{.5-}_n[\305]" font "Times-Bold,25" at ',INT(MAXVAL(time))+10 ,',',N+20
WRITE(4,*)'set label "~y{.5-}_n[\305]" font "Times-Bold,25" at 410,',N+20

!GRAY SCALE PALETTE
!WRITE(4,*)'set palette defined ( 0 0 0 0, 1 1 1 1 )'

!HOT SCALE PALETTE
WRITE(4,*)"set palette model RGB defined (0'white',0.04'yellow',0.06'red',0.12'black')"

!WRITE(4,*)'set xrange[',INT(MINVAL(time)),':',INT(MAXVAL(time)),']'
WRITE(4,*)'set xrange[0:400]'
WRITE(4,*)'set yrange[1:',N,']'
WRITE(4,*)"splot 'bl3D.dat' w l pal t '' "
WRITE(4,*)
WRITE(4,*)'set out'
CLOSE(4,STATUS='KEEP')

OPEN(5,FILE='blsurface1.gpt',STATUS='UNKNOWN')
WRITE(5,*)'#'
WRITE(5,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(5,*)'#LUIZ A. RIBEIRO JR'
WRITE(5,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(5,*)'#'
WRITE(5,*)
WRITE(5,*)'#' 
WRITE(5,*)'#3D BOND-LENTH CURVES FOR SSH.F PROGRAM'
WRITE(5,*)'#'
WRITE(5,*)
WRITE(5,*)'set term postscript color enhanced font "Times,20"'
WRITE(5,*)"set out 'blsurface1.eps'" 
WRITE(5,*)'set encoding iso_8859_1'
WRITE(5,*)'set ylabel "Site" offset 1,1,1 font "Times-Bold,25"'
WRITE(5,*)'set xlabel "Time[fs]" offset 4,1,1 font "Times-Bold,25"'
WRITE(5,*)'set pm3d'
WRITE(5,*)'set samples 100'
WRITE(5,*)'set isosamples 100'
WRITE(5,*)'set xtics 100'
WRITE(5,*)'set ytics 50'
WRITE(5,*)'set ztics 0.04'
WRITE(5,*)'set cbtics 0.02'
WRITE(5,*)'set cbtics scale 0'
WRITE(5,*)'set zlabel "~y{.5-}_n[\305]" offset -2,2,8 font "Times-Bold,25"'
WRITE(5,*)'set xyplane 0'
WRITE(5,*)'set view 25,312'

!HOT SCALE PALETTE
WRITE(5,*)"set palette model RGB defined (0'white',0.04'yellow',0.06'red',0.12'black')"

WRITE(5,*)'set xrange[',INT(MINVAL(time)),':',INT(MAXVAL(time))+20,']'
WRITE(5,*)'set yrange[0.8:',N+10,']'
WRITE(5,*)"splot 'bl3D.dat' w l lt -1 t '' "
WRITE(5,*)
WRITE(5,*)'set out'
CLOSE(5,STATUS='KEEP')

OPEN(6,FILE='blsurface2.gpt',STATUS='UNKNOWN')
WRITE(6,*)'#'
WRITE(6,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(6,*)'#LUIZ A. RIBEIRO JR'
WRITE(6,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(6,*)'#'
WRITE(6,*)
WRITE(6,*)'#' 
WRITE(6,*)'#3D BOND-LENTH CURVES FOR SSH.F PROGRAM'
WRITE(6,*)'#'
WRITE(6,*)
WRITE(6,*)'set term postscript color enhanced font "Times,20"'
WRITE(6,*)"set out 'blsurface2.eps'" 
WRITE(6,*)'set encoding iso_8859_1'
WRITE(6,*)'set ylabel "Site" offset 1,1,1 font "Times-Bold,25"'
WRITE(6,*)'set xlabel "Time[fs]" offset 4,1,1 font "Times-Bold,25"'
WRITE(6,*)'set pm3d'
WRITE(6,*)'set samples 100'
WRITE(6,*)'set isosamples 100'
WRITE(6,*)'set xtics 100'
WRITE(6,*)'set ytics 50'
WRITE(6,*)'set ztics 0.04'
WRITE(6,*)'set cbtics 0.02'
WRITE(6,*)'set cbtics scale 0'
WRITE(6,*)'set zlabel "~y{.5-}_n[\305]" offset -2,2,8 font "Times-Bold,25"'
WRITE(6,*)'set xyplane 0'
WRITE(6,*)'set view 25,312'
WRITE(6,*)'set hidden3d trianglepattern 4'

!HOT SCALE PALETTE
WRITE(6,*)"set palette model RGB defined (0'white',0.04'yellow',0.06'red',0.12'black')"

WRITE(6,*)'set xrange[',INT(MINVAL(time)),':',INT(MAXVAL(time))+20,']'
WRITE(6,*)'set yrange[0.8:',N+10,']'
WRITE(6,*)"splot 'bl3D.dat' w l lt -1 t '' "
WRITE(6,*)
WRITE(6,*)'set out'
CLOSE(6,STATUS='KEEP')

OPEN(7,FILE='blsurface3.gpt',STATUS='UNKNOWN')
WRITE(7,*)'#'
WRITE(7,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(7,*)'#LUIZ A. RIBEIRO JR'
WRITE(7,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(7,*)'#'
WRITE(7,*)
WRITE(7,*)'#' 
WRITE(7,*)'#3D BOND-LENTH CURVES FOR SSH.F PROGRAM'
WRITE(7,*)'#'
WRITE(7,*)
WRITE(7,*)'set term postscript color enhanced font "Times,20"'
WRITE(7,*)"set out 'blsurface3.eps'" 
WRITE(7,*)'set encoding iso_8859_1'
WRITE(7,*)'set ylabel "Site" offset 1,1,1 font "Times-Bold,25"'
WRITE(7,*)'set xlabel "Time[fs]" offset 4,1,1 font "Times-Bold,25"'
WRITE(7,*)'set pm3d at b'
WRITE(7,*)'set samples 100'
WRITE(7,*)'set isosamples 100'
WRITE(7,*)'set xtics 100'
WRITE(7,*)'set ytics 50'
WRITE(7,*)'set ztics 0.04'
WRITE(7,*)'set cbtics 0.02'
WRITE(7,*)'set cbtics scale 0'
WRITE(7,*)'set zlabel "~y{.5-}_n[\305]" offset -3,2,8 font "Times-Bold,25"'
WRITE(7,*)'set xyplane 0.2'
WRITE(7,*)'set view 25,312'
WRITE(7,*)'set tics out'

!HOT SCALE PALETTE
WRITE(7,*)"set palette model RGB defined (0'white',0.04'yellow',0.06'red',0.12'black')"

WRITE(7,*)'set xrange[',INT(MINVAL(time)),':',INT(MAXVAL(time)),']'
WRITE(7,*)'set yrange[0.8:',N,']'
WRITE(7,*)"splot 'bl3D.dat' w l lt -1 t '' "
WRITE(7,*)
WRITE(7,*)'set out'
CLOSE(7,STATUS='KEEP')

CALL SYSTEM('gnuplot bl2D.gpt')
CALL SYSTEM('ps2pdf bl2D.eps')

CALL SYSTEM('gnuplot bl3D.gpt')
CALL SYSTEM('ps2pdf bl3D.eps')

!CALL SYSTEM('gnuplot blsurface1.gpt')
!CALL SYSTEM('ps2pdf blsurface1.eps')

!CALL SYSTEM('gnuplot blsurface2.gpt')
!CALL SYSTEM('ps2pdf blsurface2.eps')

!CALL SYSTEM('gnuplot blsurface3.gpt')
!CALL SYSTEM('ps2pdf blsurface3.eps')

100 FORMAT(I4,2X,20(e17.10,1x))
110 FORMAT(T4,A,I5)

RETURN
END
!***********************************************************************************************************
!***********************************************************************************************************

!***********************************************************************************************************
!***********************************************************************************************************
SUBROUTINE chargedensity(steps,N,dt,iwritemode)
IMPLICIT NONE
INTEGER      :: I,J
INTEGER      :: steps,N,iwritemode
REAL(KIND=8) :: cd(N),auxcd(steps,N),time(steps)
REAL(KIND=8) :: dt

DO I=1,steps
  time(I) = FLOAT(I-1)*4.0d0*dt*iwritemode
END DO

DO I=1,steps
   DO J=1,N
      auxcd(I,J)=0.0d0
   END DO
END DO

OPEN(14,STATUS='OLD')
DO J=1,steps
   READ(14,*)
   READ(14,*)
   READ(14,*)(cd(I),I=1,N)
   auxcd(J,1) = (cd(N)+2.0d0*cd(1)+cd(2))/4.0d0
   auxcd(J,N) = (cd(N-1)+2.0d0*cd(N)+cd(1))/4.0d0
   DO I=2,N-1
      auxcd(J,I)=(cd(I-1)+2.0d0*cd(I)+cd(I+1))/4.0d0
   END DO
END DO
CLOSE(14,STATUS='KEEP')

OPEN(UNIT=1,FILE='cd2D.dat',STATUS='UNKNOWN')
DO I=1,steps
   WRITE(1,110)'#',I
   DO J=1,N
      WRITE(1,100)J,1.0d0-auxcd(I,J)
   END DO
   WRITE(1,*)
   WRITE(1,*)
END DO
CLOSE(1,STATUS='KEEP')

OPEN(2,FILE='cd3D.dat',STATUS='UNKNOWN')
DO J=1,N
   DO I=1,steps
      WRITE(2,*)time(I),FLOAT(J),1.0d0-auxcd(I,J)
   END DO
   WRITE(2,*)
END DO
CLOSE(2,STATUS='KEEP')

OPEN(3,FILE='cd2D.gpt',STATUS='UNKNOWN')
WRITE(3,*)'#'
WRITE(3,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(3,*)'#LUIZ A. RIBEIRO JR'
WRITE(3,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(3,*)'#'
WRITE(3,*)
WRITE(3,*)'#' 
WRITE(3,*)'#2D CHARGE DENSITY CURVES FOR SSH.F PROGRAM'
WRITE(3,*)'#'
WRITE(3,*)
WRITE(3,*)'set term postscript color enhanced font "Times,20"'
WRITE(3,*)"set out 'cd2D.eps'" 
WRITE(3,*)'set encoding iso_8859_1'
WRITE(3,*)'set ylabel "{/Symbol ~r{.5-}}_n[e]" font "Times-Bold,25"'
WRITE(3,*)'set xlabel "Site" font "Times-Bold,25"'
WRITE(3,*)'set xrange [0:',N+2,']'
WRITE(3,*)'set style line 1  lt -1 lc -1 lw 2'
WRITE(3,*)"plot 'cd2D.dat' index 1 us 1:2 w l ls 1 t ''"
WRITE(3,*)
WRITE(3,*)'set out'
CLOSE(3,STATUS='KEEP')

OPEN(4,FILE='cd3D.gpt',STATUS='UNKNOWN')
WRITE(4,*)'#'
WRITE(4,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(4,*)'#LUIZ A. RIBEIRO JR'
WRITE(4,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(4,*)'#'
WRITE(4,*)
WRITE(4,*)'#' 
WRITE(4,*)'#3D CHARGE DENSITY CURVES FOR SSH.F PROGRAM'
WRITE(4,*)'#'
WRITE(4,*)
WRITE(4,*)'set term postscript color enhanced font "Times,20"'
WRITE(4,*)"set out 'cd3D.eps'" 
WRITE(4,*)'set encoding iso_8859_1'
WRITE(4,*)'set ylabel "Site" font "Times-Bold,25"'
WRITE(4,*)'set xlabel "Time[fs]" font "Times-Bold,25"'
!WRITE(4,*)'set dgrid3d 300,300 gauss 2'
WRITE(4,*)'set pm3d'
WRITE(4,*)'set view map'
WRITE(4,*)'unset surface'
WRITE(4,*)'set xtics 100'
WRITE(4,*)'set ytics 50'
WRITE(4,*)'set view map'
WRITE(4,*)'set cbtics 0.02'
WRITE(4,*)'set cbtics scale 0'
WRITE(4,120)'set label "{/Symbol ~r{.5-}}_n[e]" font "Times-Bold,25" at ',INT(MAXVAL(time))+10 ,',',N+20
WRITE(4,*)"set palette model RGB defined (0'blue',.35'cyan',.5'green',.65'yellow',1'red')"
WRITE(4,*)'set xrange[',INT(MINVAL(time)),':',INT(MAXVAL(time)),']'
WRITE(4,*)'set yrange[1:',N,']'
WRITE(4,*)"splot 'cd3D.dat' w l pal t '' "
WRITE(4,*)
WRITE(4,*)'set out'
CLOSE(4,STATUS='KEEP')

CALL SYSTEM('gnuplot cd2D.gpt')
CALL SYSTEM('ps2pdf cd2D.eps')

CALL SYSTEM('gnuplot cd3D.gpt')
CALL SYSTEM('ps2pdf cd3D.eps')

100 FORMAT(I4,2X,20(e17.10,1x))
110 FORMAT(T4,A,I5)
120 FORMAT(A,I4,A,I4)

RETURN
END
!***********************************************************************************************************
!***********************************************************************************************************

!***********************************************************************************************************
!***********************************************************************************************************
SUBROUTINE spindensity(steps,N,dt,iwritemode)
IMPLICIT NONE
INTEGER      :: I,J
INTEGER      :: steps,N,iwritemode
REAL(KIND=8) :: sd(N),auxsd(steps,N),time(steps)
REAL(KIND=8) :: dt

DO I=1,steps
  time(I) = FLOAT(I-1)*4.0d0*dt*iwritemode
END DO

DO I=1,steps
   DO J=1,N
      auxsd(I,J)=0.0d0
   END DO
END DO

OPEN(15,STATUS='OLD')
DO J=1,steps
   READ(15,*)
   READ(15,*)
   READ(15,*)(sd(I),I=1,N)
   auxsd(J,1) = (sd(N)+2.0d0*sd(1)+sd(2))/4.0d0
   auxsd(J,N) = (sd(N-1)+2.0d0*sd(N)+sd(1))/4.0d0
   DO I=2,N-1
      auxsd(J,I)=(sd(I-1)+2.0d0*sd(I)+sd(I+1))/4.0d0
   END DO
END DO
CLOSE(15,STATUS='KEEP')

OPEN(UNIT=1,FILE='sd2D.dat',STATUS='UNKNOWN')
DO I=1,steps
   WRITE(1,110)'#',I
   DO J=1,N
      WRITE(1,100)J,auxsd(I,J)
   END DO
   WRITE(1,*)
   WRITE(1,*)
END DO
CLOSE(1,STATUS='KEEP')

OPEN(2,FILE='sd3D.dat',STATUS='UNKNOWN')
DO J=1,N
   DO I=1,steps
      WRITE(2,*)time(I),FLOAT(J),auxsd(I,J)
   END DO
   WRITE(2,*)
END DO
CLOSE(2,STATUS='KEEP')

OPEN(3,FILE='sd2D.gpt',STATUS='UNKNOWN')
WRITE(3,*)'#'
WRITE(3,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(3,*)'#LUIZ A. RIBEIRO JR'
WRITE(3,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(3,*)'#'
WRITE(3,*)
WRITE(3,*)'#' 
WRITE(3,*)'#2D SPIN DENSITY CURVES FOR SSH.F PROGRAM'
WRITE(3,*)'#'
WRITE(3,*)
WRITE(3,*)'set term postscript color enhanced font "Times,20"'
WRITE(3,*)"set out 'sd2D.eps'" 
WRITE(3,*)'set encoding iso_8859_1'
WRITE(3,*)'set ylabel "~s{.5-}_n[e]" font "Times-Bold,25"'
WRITE(3,*)'set xlabel "Site" font "Times-Bold,25"'
WRITE(3,*)'set xrange [0:',N+2,']'
WRITE(3,*)'set style line 1  lt -1 lc -1 lw 2'
WRITE(3,*)"plot 'sd2D.dat' index 1 us 1:2 w l ls 1 t ''"
WRITE(3,*)
WRITE(3,*)'set out'
CLOSE(3,STATUS='KEEP')

OPEN(4,FILE='sd3D.gpt',STATUS='UNKNOWN')
WRITE(4,*)'#'
WRITE(4,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(4,*)'#LUIZ A. RIBEIRO JR'
WRITE(4,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(4,*)'#'
WRITE(4,*)
WRITE(4,*)'#' 
WRITE(4,*)'#3D SPIN DENSITY CURVES FOR SSH.F PROGRAM'
WRITE(4,*)'#'
WRITE(4,*)
WRITE(4,*)'set term postscript color enhanced font "Times,20"'
WRITE(4,*)"set out 'sd3D.eps'" 
WRITE(4,*)'set encoding iso_8859_1'
WRITE(4,*)'set ylabel "Site" font "Times-Bold,25"'
WRITE(4,*)'set xlabel "Time[fs]" font "Times-Bold,25"'
WRITE(4,*)'set dgrid3d 300,300 gauss 2'
WRITE(4,*)'set pm3d'
WRITE(4,*)'set view map'
WRITE(4,*)'unset surface'
WRITE(4,*)'set xtics 100'
WRITE(4,*)'set ytics 50'
WRITE(4,*)'set view map'
WRITE(4,*)'set cbtics 0.02'
WRITE(4,*)'set cbtics scale 0'
WRITE(4,120)'set label "~s{.5-}_n[e]" font "Times-Bold,25" at ',INT(MAXVAL(time))+10 ,',',N+20

!GRAY SCALE PALETTE POSITIVE
!WRITE(4,*)'set palette defined ( 0 0 0 0, 1 1 1 1 )'

!GRAY SCALE PALETTE NEGATIVE
WRITE(4,*)'set palette gray negative'
!GRAY SCALE PALETTE OCEAN
WRITE(4,*)'set palette rgb 23,28,3'

WRITE(4,*)'set xrange[',INT(MINVAL(time)),':',INT(MAXVAL(time)),']'
WRITE(4,*)'set yrange[1:',N,']'
WRITE(4,*)"splot 'sd3D.dat' w l pal t '' "
WRITE(4,*)
WRITE(4,*)'set out'
CLOSE(4,STATUS='KEEP')

CALL SYSTEM('gnuplot sd2D.gpt')
CALL SYSTEM('ps2pdf sd2D.eps')

!CALL SYSTEM('gnuplot sd3D.gpt')
!CALL SYSTEM('ps2pdf sd3D.eps')

100 FORMAT(I4,2X,20(e17.10,1x))
110 FORMAT(T4,A,I5)
120 FORMAT(A,I4,A,I4)

RETURN
END
!***********************************************************************************************************
!***********************************************************************************************************

!***********************************************************************************************************
!***********************************************************************************************************
SUBROUTINE electricfield(IFIELDOFF,IFIELDON,dt)
IMPLICIT NONE
INTEGER      :: I,IFIELDOFF,IFIELDON,auxIFIELDON
REAL(KIND=8) :: efield(IFIELDOFF-IFIELDON),dt,time(IFIELDOFF-IFIELDON)

IF(IFIELDON .EQ. 0)THEN
  auxIFIELDON = 1
ELSE
  auxIFIELDON = IFIELDON
END IF 

DO I=auxIFIELDON,IFIELDOFF
  IF(IFIELDON .EQ. 0)THEN
     time(I) = FLOAT(I-1)*4.0d0*dt
  ELSE
     time(I) = FLOAT(I)*4.0d0*dt
  END IF   
END DO

DO I=auxIFIELDON,IFIELDOFF
   READ(18,*) efield(I)
END DO

OPEN(1,FILE='efield.dat',STATUS='UNKNOWN')
WRITE(1,*) 0.0d0,0.0d0
DO I=auxIFIELDON,IFIELDOFF
   IF(MINVAL(efield) .LT. 0.0d0)THEN
      WRITE(1,*) time(I),-1.0d0*(efield(I)*1.3d0)/0.01d0
   ELSE
      WRITE(1,*) time(I),(efield(I)*1.3d0)/0.01d0
   END IF
END DO
CLOSE(1,STATUS='KEEP')

OPEN(2,FILE='efield.gpt',STATUS='UNKNOWN')
WRITE(2,*)'#'
WRITE(2,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(2,*)'#LUIZ A. RIBEIRO JR'
WRITE(2,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(2,*)'#'
WRITE(2,*)
WRITE(2,*)'#' 
WRITE(2,*)'#ELECTRIC FIELD BEHAVIOUR FOR SSH.F PROGRAM'
WRITE(2,*)'#'
WRITE(2,*)
WRITE(2,*)'set term postscript color enhanced font "Times,20"'
WRITE(2,*)"set out 'efield.eps'" 
WRITE(2,*)'set encoding iso_8859_1'
WRITE(2,*)'set ylabel "E[mV/\305]" font "Times-Bold,25"'
WRITE(2,*)'set xlabel "Time[fs]" font "Times-Bold,25"'
WRITE(2,*)'set style line 1  lt -1 lc -1 lw 2'
WRITE(2,*)'set xtics 100'
WRITE(2,*)"plot 'efield.dat' us 1:2 w l ls 1 t ''"
WRITE(2,*)
WRITE(2,*)'set out'
CLOSE(2,STATUS='KEEP')

CALL SYSTEM('gnuplot efield.gpt')
CALL SYSTEM('ps2pdf efield.eps')

RETURN
END
!***********************************************************************************************************
!***********************************************************************************************************

!***********************************************************************************************************
!***********************************************************************************************************
SUBROUTINE energylevels(steps,N,dt,iwritemode)
IMPLICIT NONE
INTEGER      :: I,J,N,iwritemode,steps
REAL(KIND=8) :: time(steps),el(steps,N),dt

DO I=1,steps
  time(I) = FLOAT(I-1)*4.0d0*dt*iwritemode
END DO

OPEN(12,STATUS='OLD')
DO I=1,steps
   READ(12,*)
   READ(12,*)
   READ(12,*)(el(I,J),J=1,N)
END DO
CLOSE(12,STATUS='KEEP')

OPEN(1,FILE='elevelsup.dat',STATUS='UNKNOWN')
DO J=N/2-20,N/2+21
   WRITE(1,110)'#',J
   DO I=1,steps
      WRITE(1,*)time(I),2.5d0*el(I,J)
   END DO
   WRITE(1,*)
   WRITE(1,*)
END DO   
CLOSE(1,STATUS='KEEP')

OPEN(2,FILE='elevelsup.gpt',STATUS='UNKNOWN')
WRITE(2,*)'#'
WRITE(2,*)'#INSITTUTE OF PHYSICS - UNIVERSITY OF BRASILIA'
WRITE(2,*)'#LUIZ A. RIBEIRO JR'
WRITE(2,*)'#BRASILIA, OCTOBER 12 2012'
WRITE(2,*)'#'
WRITE(2,*)
WRITE(2,*)'#' 
WRITE(2,*)'#ENERGY LEVELS UP BEHAVIOUR FOR SSH.F PROGRAM'
WRITE(2,*)'#'
WRITE(2,*)
WRITE(2,*)'set term postscript color enhanced font "Times,20"'
WRITE(2,*)"set out 'elevelsup.eps'" 
WRITE(2,*)'set encoding iso_8859_1'
WRITE(2,*)'set ylabel "Energy[eV]" font "Times-Bold,25"'
WRITE(2,*)'set xlabel "Time[fs]" font "Times-Bold,25"'
WRITE(2,*)'set style line 1  lt -1 lw 2'
WRITE(2,*)'set style line 2  lt  3 lw 2'
WRITE(2,*)'set xtics 100'

!SHOW LINES INSIDE THE GAP WITH DIFFERENT FORMAT (PAY ATTENTION IN THE SYMMETRY)
!THE BELLOW CONFIGURATION IS FOR TWO DEFECTS INSIDE THE GAP WITH N=200
WRITE(2,*)"plot 'elevelsup.dat' index 1:18  us 1:2 w l ls 1 t '',\"
WRITE(2,*)"     'elevelsup.dat' index 23:41 us 1:2 w l ls 1 t '',\"
WRITE(2,*)"     'elevelsup.dat' index 19:22 us 1:2 w l ls 2 t ''"

WRITE(2,*)
WRITE(2,*)'set out'
CLOSE(2,STATUS='KEEP')

CALL SYSTEM('gnuplot elevelsup.gpt')
CALL SYSTEM('ps2pdf elevelsup.eps')

RETURN

110 FORMAT(T4,A,I5)

END
!***********************************************************************************************************
!***********************************************************************************************************

!***********************************************************************************************************
!***********************************************************************************************************
SUBROUTINE velocities(steps,N,dt,iwritemode)
IMPLICIT NONE
PARAMETER    PI=4.0d0*DATAN(1.0d0)
INTEGER      :: I,J,N,steps,iwritemode
REAL(KIND=8) :: A,B,dt,THETA
REAL(KIND=8) :: C(N),CS(N),position(steps),velocity(steps)

REWIND(14)
DO J=1,steps
   READ(14,*)
   READ(14,*)
   READ(14,*)(C(I),I=1,N)
   
   CS(1)=1.0d0-(C(N)+2.0d0*C(1)+C(2))/4.0d0
   
   DO I=2,N-1
      CS(I)=1.0d0-(C(I-1)+2.0d0*C(I)+C(I+1))/4.0d0
   END DO
   
   CS(N)=1.0d0-(C(N-1)+2.0d0*C(N)+C(1))/4.0d0
   A=0.0d0
   B=0.0d0
   
   DO I=1,N 
      A = A + DCOS((2.0d0*I+1.0d0)*PI/N)*CS(I)
      B = B + DSIN((2.0d0*I+1.0d0)*PI/N)*CS(I)
   END DO
   
   THETA=DATAN(B/A)
   IF(B.GT.0.0D0.AND.A.LT.0.0D0) THETA=THETA+PI
   IF(B.LT.0.0D0.AND.A.LT.0.0D0) THETA=THETA+PI
   IF(B.LT.0.0D0.AND.A.GT.0.0D0) THETA=THETA+PI*2.0d0
   position(J)=N*0.5d0*THETA/PI
END DO

OPEN(1,FILE='position.dat',STATUS='UNKNOWN')
DO I=1,steps
   WRITE(1,*) FLOAT(I-1)*dt*4.0d0*iwritemode,position(I)*1.22d0
END DO
CLOSE(1,STATUS='KEEP')

DO I=1,steps-1
   velocity(I)=((position(I+1)-position(I))*1.22d0)/((FLOAT(I+1)-FLOAT(I))*dt*4.0d0*iwritemode)
END DO
   
OPEN(2,FILE='velocity.dat',STATUS='UNKNOWN')
DO I=1,steps-1
   WRITE(2,*) FLOAT(I)*dt*4.0d0*iwritemode,velocity(I)
END DO
CLOSE(2,STATUS='KEEP')

RETURN   
END
!***********************************************************************************************************
!***********************************************************************************************************
