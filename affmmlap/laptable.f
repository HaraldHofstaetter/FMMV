	SUBROUTINE VWTS(X,W,N)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION W(N),X(N)
c
c
c this routine returns a set of gaussian nodes and weights for
c integrating the functions j0(r*x)*exp(z*x)dx over the range x=0 to
c x=infinity.
c they work only for values of z within the range [1,4] and r within
c the range [0,4*sqrt(2)].
c
c
c input arguments:
c	n - number of weights and nodes in the quadrature.  this must
c	  be an integer in the range [2,39].
c
c output arguments:
c	w - weights
c	x - nodes
c
c
c approximate accuracies of the quadrature:
c
c	n	maximum error
c
c        2     0.15318e+00
c        3     0.76505e-01
c        4     0.32149e-01
c        5     0.15630e-01
c        6     0.75110e-02
c        7     0.35030e-02
c        8     0.16243e-02
c        9     0.72230e-03
c       10     0.33074e-03
c       11     0.15035e-03
c       12     0.70952e-04
c       13     0.31751e-04
c       14     0.14589e-04
c       15     0.64300e-05
c       16     0.29477e-05
c       17     0.13222e-05
c       18     0.61488e-06
c       19     0.27435e-06
c       20     0.12534e-06
c       21     0.55324e-07
c       22     0.25257e-07
c       23     0.11293e-07
c       24     0.52063e-08
c       25     0.23256e-08
c       26     0.10580e-08
c       27     0.46835e-09
c       28     0.21286e-09
c       29     0.95164e-10
c       30     0.43599e-10
c       31     0.19516e-10
c       32     0.88491e-11
c       33     0.39313e-11
c       34     0.17821e-11
c       35     0.79603e-12
c       36     0.36460e-12
c       37     0.16331e-12
c       38     0.73497e-13
c       39     0.31530e-13
c
c
	IF(N.LT.2 .OR. N.GT.39) STOP 'n out of bounds in vwts()'

	IF(N.EQ.  9) THEN
	  X(1)=0.99273996739714473469540223504736787E-01
	  X(2)=0.47725674637049431137114652301534079E+00
	  X(3)=0.10553366138218296388373573790886439E+01
	  X(4)=0.17675934335400844688024335482623428E+01
	  X(5)=0.25734262935147067530294862081063911E+01
	  X(6)=0.34482433920158257478760788217186928E+01
	  X(7)=0.43768098355472631055818055756390095E+01
	  X(8)=0.53489575720546005399569367000367492E+01
	  X(9)=0.63576578531337464283978988532908261E+01
	  W(1)=0.24776441819008371281185532097879332E+00
	  W(2)=0.49188566500464336872511239562300034E+00
	  W(3)=0.65378749137677805158830324216978624E+00
	  W(4)=0.76433038408784093054038066838984378E+00
	  W(5)=0.84376180565628111640563702167128213E+00
	  W(6)=0.90445883985098263213586733400006779E+00
	  W(7)=0.95378613136833456653818075210438110E+00
	  W(8)=0.99670261613218547047665651916759089E+00
	  W(9)=0.10429422730252668749528766056755558E+01
	ELSEIF(N.EQ. 18) THEN
	  X(1)=0.52788527661177607475107009804560221E-01
	  X(2)=0.26949859838931256028615734976483509E+00
	  X(3)=0.63220353174689392083962502510985360E+00
	  X(4)=0.11130756427760852833586113774799742E+01
	  X(5)=0.16893949614021379623807206371566281E+01
	  X(6)=0.23437620046953044905535534780938178E+01
	  X(7)=0.30626998290780611533534738555317745E+01
	  X(8)=0.38356294126529686394633245072327554E+01
	  X(9)=0.46542473432156272750148673367220908E+01
	  X(10)=0.55120938659358147404532246582675725E+01
	  X(11)=0.64042126837727888499784967279992998E+01
	  X(12)=0.73268800190617540124549122992902994E+01
	  X(13)=0.82774009925823861522076185792684555E+01
	  X(14)=0.92539718060248947750778825138695538E+01
	  X(15)=0.10255602723746401139237605093512684E+02
	  X(16)=0.11282088297877740146191172243561596E+02
	  X(17)=0.12334067909676926788620221486780792E+02
	  X(18)=0.13414920240172401477707353478763252E+02
	  W(1)=0.13438265914335215112096477696468355E+00
	  W(2)=0.29457752727395436487256574764614925E+00
	  W(3)=0.42607819361148618897416895379137713E+00
	  W(4)=0.53189220776549905878027857397682965E+00
	  W(5)=0.61787306245538586857435348065337166E+00
	  W(6)=0.68863156078905074508611505734734237E+00
	  W(7)=0.74749099381426187260757387775811367E+00
	  W(8)=0.79699192718599998208617307682288811E+00
	  W(9)=0.83917454386997591964103548889397644E+00
	  W(10)=0.87570092283745315508980411323136650E+00
	  W(11)=0.90792943590067498593754180546966381E+00
	  W(12)=0.93698393742461816291466902839601971E+00
	  W(13)=0.96382546688788062194674921556725167E+00
	  W(14)=0.98932985769673820186653756536543369E+00
	  W(15)=0.10143828459791703888726033255807124E+01
	  W(16)=0.10400365437416452252250564924906939E+01
	  W(17)=0.10681548926956736522697610780596733E+01
	  W(18)=0.11090758097553685690428437737864442E+01
        ELSE
          PRINT *, 'wrong accuracy requirement, stop'
          STOP
	ENDIF
c
        RETURN
	END
c
	SUBROUTINE NUMTHETAHALF(NUMTETS,NLAMS)
	INTEGER *4 NUMTETS(NLAMS),NLAMS
c
c this routine returns the number of fourier modes needed in the
c phi integral for each of the discrete lambda values given
c by norman's quadratures. (see weights.f).
c
c input arguments:
c	nlams - number of nodes in the lambda quadrature.  this must
c	  be an integer in the range [2,39].
c
c output arguments:
c	numtets(i) - number of fourier modes needed for phi
c                    integral with lambda_i.
c
c approximate accuracies of the quadrature:
c
c	nlams	maximum error
c
c        2     0.15318e+00
c        3     0.76505e-01
c        4     0.32149e-01
c        5     0.15630e-01
c        6     0.75110e-02
c        7     0.35030e-02
c        8     0.16243e-02
c        9     0.72230e-03
c       10     0.33074e-03
c       11     0.15035e-03
c       12     0.70952e-04
c       13     0.31751e-04
c       14     0.14589e-04
c       15     0.64300e-05
c       16     0.29477e-05
c       17     0.13222e-05
c       18     0.61488e-06
c       19     0.27435e-06
c       20     0.12534e-06
c       21     0.55324e-07
c       22     0.25257e-07
c       23     0.11293e-07
c       24     0.52063e-08
c       25     0.23256e-08
c       26     0.10580e-08
c       27     0.46835e-09
c       28     0.21286e-09
c       29     0.95164e-10
c       30     0.43599e-10
c       31     0.19516e-10
c       32     0.88491e-11
c       33     0.39313e-11
c       34     0.17821e-11
c       35     0.79603e-12
c       36     0.36460e-12
c       37     0.16331e-12
c       38     0.73497e-13
c       39     0.31530e-13
c
      IF(NLAMS.EQ. 9) THEN
         NUMTETS( 1) = 2
         NUMTETS( 2) = 4
         NUMTETS( 3) = 4
         NUMTETS( 4) = 6
         NUMTETS( 5) = 6
         NUMTETS( 6) = 4
         NUMTETS( 7) = 6
         NUMTETS( 8) = 4
         NUMTETS( 9) = 2
ccc         numtets( 1) = 2
ccc         numtets( 2) = 3
ccc         numtets( 3) = 4
ccc         numtets( 4) = 5
ccc         numtets( 5) = 5
ccc         numtets( 6) = 4
ccc         numtets( 7) = 5
ccc         numtets( 8) = 4
ccc         numtets( 9) = 1
      ELSEIF (NLAMS.EQ.18) THEN
         NUMTETS( 1) = 4
         NUMTETS( 2) = 6
         NUMTETS( 3) = 6
         NUMTETS( 4) = 8
         NUMTETS( 5) = 8
         NUMTETS( 6) = 8
         NUMTETS( 7) =10
         NUMTETS( 8) =10
         NUMTETS( 9) =10
         NUMTETS(10) =10
         NUMTETS(11) =12
         NUMTETS(12) =12
         NUMTETS(13) =12
         NUMTETS(14) =12
         NUMTETS(15) =12
         NUMTETS(16) =12
         NUMTETS(17) = 8
         NUMTETS(18) = 2
      ELSE
        PRINT *, 'wrong accuracy requirement, stop'
        STOP
      ENDIF
c
      RETURN
      END
c
	SUBROUTINE NUMTHETAFOUR(NUMTETS,NLAMS)
	INTEGER *4 NUMTETS(NLAMS),NLAMS
c
c this routine returns the number of fourier modes needed in the
c phi integral for each of the discrete lambda values given
c by norman's quadratures. (see weights.f).
c
c input arguments:
c	nlams - number of nodes in the lambda quadrature.  this must
c	  be an integer in the range [2,39].
c
c output arguments:
c	numtets(i) - number of fourier modes needed for phi
c                    integral with lambda_i.
c
c approximate accuracies of the quadrature:
c
c	nlams	maximum error
c
c        2     0.15318e+00
c        3     0.76505e-01
c        4     0.32149e-01
c        5     0.15630e-01
c        6     0.75110e-02
c        7     0.35030e-02
c        8     0.16243e-02
c        9     0.72230e-03
c       10     0.33074e-03
c       11     0.15035e-03
c       12     0.70952e-04
c       13     0.31751e-04
c       14     0.14589e-04
c       15     0.64300e-05
c       16     0.29477e-05
c       17     0.13222e-05
c       18     0.61488e-06
c       19     0.27435e-06
c       20     0.12534e-06
c       21     0.55324e-07
c       22     0.25257e-07
c       23     0.11293e-07
c       24     0.52063e-08
c       25     0.23256e-08
c       26     0.10580e-08
c       27     0.46835e-09
c       28     0.21286e-09
c       29     0.95164e-10
c       30     0.43599e-10
c       31     0.19516e-10
c       32     0.88491e-11
c       33     0.39313e-11
c       34     0.17821e-11
c       35     0.79603e-12
c       36     0.36460e-12
c       37     0.16331e-12
c       38     0.73497e-13
c       39     0.31530e-13
c
      IF(NLAMS.EQ. 9) THEN
         NUMTETS( 1) =  4
         NUMTETS( 2) =  8
         NUMTETS( 3) = 12
         NUMTETS( 4) = 16
         NUMTETS( 5) = 20
         NUMTETS( 6) = 20
         NUMTETS( 7) = 24
         NUMTETS( 8) =  8
         NUMTETS( 9) =  2
      ELSEIF(NLAMS.EQ.18) THEN
         NUMTETS( 1) =  6
         NUMTETS( 2) =  8
         NUMTETS( 3) = 12
         NUMTETS( 4) = 16
         NUMTETS( 5) = 20
         NUMTETS( 6) = 26
         NUMTETS( 7) = 30
         NUMTETS( 8) = 34
         NUMTETS( 9) = 38
         NUMTETS(10) = 44
         NUMTETS(11) = 48
         NUMTETS(12) = 52
         NUMTETS(13) = 56
         NUMTETS(14) = 60
         NUMTETS(15) = 60
         NUMTETS(16) = 52
         NUMTETS(17) =  4
         NUMTETS(18) =  2
      ELSE
        PRINT *, 'wrong accuracy requirement, stop'
        STOP
      ENDIF
c
      RETURN
      END

