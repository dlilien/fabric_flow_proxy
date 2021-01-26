      FUNCTION VSIA2000M01MY(Model, nodenumber, dumy) RESULT(yvel)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: dumy,yvel
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac
            REAL(kind=dp) :: n, H, acc

            H = 2000._dp
            n = 4._dp
            acc = 0.1_dp

            y=Model % Nodes % y (nodenumber)

            
            frac = (h - y) / h
            yvel = -acc * (1.0_dp - frac * ((n + 2.0_dp) / (n + 1.0_dp) - 1.0_dp / (1.0_dp + n) * frac ** (n + 1)))

            RETURN 
        END


      FUNCTION USIA2000M125MY(Model, nodenumber, dumy) RESULT(xvel)
          ! Profile of horizontal velocity assuming no sliding
          ! Taken from Charlie's 1983 modeling paper
          ! Basically we assume n, shape function, surf vel, and thickness
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: xvel
            REAL(kind=dp) :: us, surf, bed
            REAL(kind=dp) :: dumy
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac
            REAL(kind=dp) :: n

            n = 4.0_dp

            us = 12.0_dp
            surf = 2000.0_dp
            bed = 0.0_dp
 
            y=Model % Nodes % y (nodenumber)

            frac = (surf - y) / (surf - bed)
            xvel = us * (1.0_dp - frac**(n+1))
            RETURN 
        END


      FUNCTION POSSPIN1e0(Model, nodenumber, dumy) RESULT(spin)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: spin
            REAL(kind=dp) :: us, surf, bed
            REAL(kind=dp) :: dumy
            INTEGER :: nodenumber
            REAL(kind=dp) :: x, frac
            REAL(kind=dp) :: n, maxdudy, width, center

            center = 0.
            width = 5000.
            maxdudy = 1.0
            x=Model % Nodes % x (nodenumber)
            spin = maxdudy * exp(-0.5*((x - center)/width)**2.0)
            ! Need to zero it outside
            IF (x.GE.(2.0*width + center)) THEN
                spin = 0.0_dp
            END IF
            IF (x.LE.(-2.0*width + center)) THEN
                spin = 0.0_dp
            END IF

            RETURN 
        END


      FUNCTION POSSPIN1e1(Model, nodenumber, dumy) RESULT(spin)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: spin
            REAL(kind=dp) :: us, surf, bed
            REAL(kind=dp) :: dumy
            INTEGER :: nodenumber
            REAL(kind=dp) :: x, frac
            REAL(kind=dp) :: n, maxdudy, width, center

            center = 0.
            width = 5000.
            maxdudy = 1.0e-1
            x=Model % Nodes % x (nodenumber)
            spin = maxdudy * exp(-0.5*((x - center)/width)**2.0)
            ! Need to zero it outside
            IF (x.GE.(2.0*width + center)) THEN
                spin = 0.0_dp
            END IF
            IF (x.LE.(-2.0*width + center)) THEN
                spin = 0.0_dp
            END IF

            RETURN 
        END


      FUNCTION POSSPIN1e2(Model, nodenumber, dumy) RESULT(spin)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: spin
            REAL(kind=dp) :: us, surf, bed
            REAL(kind=dp) :: dumy
            INTEGER :: nodenumber
            REAL(kind=dp) :: x, frac
            REAL(kind=dp) :: n, maxdudy, width, center

            center = 0.
            width = 5000.
            maxdudy = 1.0e-2
            x=Model % Nodes % x (nodenumber)
            spin = maxdudy * exp(-0.5*((x - center)/width)**2.0)
            ! Need to zero it outside
            IF (x.GE.(2.0*width + center)) THEN
                spin = 0.0_dp
            END IF
            IF (x.LE.(-2.0*width + center)) THEN
                spin = 0.0_dp
            END IF

            RETURN 
        END


      FUNCTION POSSPIN1e3(Model, nodenumber, dumy) RESULT(spin)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: spin
            REAL(kind=dp) :: us, surf, bed
            REAL(kind=dp) :: dumy
            INTEGER :: nodenumber
            REAL(kind=dp) :: x, frac
            REAL(kind=dp) :: n, maxdudy, width, center

            center = 0.
            width = 5000.
            maxdudy = 1e-3
            x=Model % Nodes % x (nodenumber)
            spin = maxdudy * exp(-0.5*((x - center)/width)**2.0)
            ! Need to zero it outside
            IF (x.GE.(2.0*width + center)) THEN
                spin = 0.0_dp
            END IF
            IF (x.LE.(-2.0*width + center)) THEN
                spin = 0.0_dp
            END IF

            RETURN 
        END
        

      FUNCTION OOP13(Model, nodenumber, inputs) RESULT(theoop)
            ! This computes
            ! $\dot{varepsilon}^{(xy)}=\frac{1}{2}\left(\frac{\partial u_s}{\partial y}\frac{u(z)}{u_s} - \frac{u(z)}{r}\right)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: theoop, dudy, vel, surfvel, r
            REAL(kind=dp), dimension(4) :: inputs
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac

            dudy = inputs(1)
            vel = inputs(2)
            surfvel = inputs(3)
            r = inputs(4)

            theoop = (1./2.) * (dudy * vel / surfvel - vel / r)
            RETURN 
        END


        FUNCTION ftw(model, nodenumber, dumy) RESULT(theftw)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: theftw
            REAL(kind=dp):: dumy
            INTEGER :: nodenumber
            REAL(kind=dp) :: x

            x=Model % Nodes % x (nodenumber)
            theftw=1.0_dp
            IF (x.LT.0) THEN
                theftw = 1.0_dp - x / 10000.
            END IF
            RETURN 
        END


        FUNCTION ftwlc(model, nodenumber, dumy) RESULT(theftw)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: theftw
            REAL(kind=dp):: dumy
            INTEGER :: nodenumber
            REAL(kind=dp) :: x

            x=Model % Nodes % x (nodenumber)
            theftw=1.0_dp
            IF (x.LT.0) THEN
                theftw = 1.0_dp - x / 33333.
            END IF
            RETURN 
        END


        FUNCTION hiacc(model, nodenumber, dumy) RESULT(theftw)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: theftw
            REAL(kind=dp):: dumy
            INTEGER :: nodenumber

            theftw=1.0_dp
            RETURN 
        END


        FUNCTION ftw_uc(model, nodenumber, dumy) RESULT(theftw)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: theftw
            REAL(kind=dp):: dumy
            INTEGER :: nodenumber
            REAL(kind=dp) :: x

            x=Model % Nodes % x (nodenumber)
            theftw = 6.0_dp - x / 20000.
            RETURN 
        END


        FUNCTION nuc(model, nodenumber, dumy) RESULT(theftw)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: theftw
            REAL(kind=dp):: dumy
            INTEGER :: nodenumber
            REAL(kind=dp) :: x

            x=Model % Nodes % x (nodenumber)
            theftw = 3.0_dp - x / 50000.
            RETURN 
        END


        FUNCTION nuclc(model, nodenumber, dumy) RESULT(theftw)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: theftw
            REAL(kind=dp):: dumy
            INTEGER :: nodenumber
            REAL(kind=dp) :: x

            x=Model % Nodes % x (nodenumber)
            theftw = 2.0_dp - x / 66666.
            RETURN 
        END


      FUNCTION VSSA2000M01MY(Model, nodenumber, inputs) RESULT(yvel)
            ! Just use the fact that strain is constant
            ! so vertical vel is linear in depth
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: yvel, smb, surf, bed
            REAL(kind=dp) :: inputs
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac

            smb = 0.1_dp
            surf = 2000.0_dp
            bed = 0.0_dp

            y=Model % Nodes % y (nodenumber)
            frac = (y - bed) / (surf - bed)

            yvel = -frac * smb
            RETURN 
        END


      FUNCTION VSSA2000M03MY(Model, nodenumber, inputs) RESULT(yvel)
            ! Just use the fact that strain is constant
            ! so vertical vel is linear in depth
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: yvel, smb, surf, bed
            REAL(kind=dp) :: inputs
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac

            smb = 0.3_dp
            surf = 2000.0_dp
            bed = 0.0_dp

            y=Model % Nodes % y (nodenumber)
            frac = (y - bed) / (surf - bed)

            yvel = -frac * smb
            RETURN 
        END


      FUNCTION VSSA2000M02MY(Model, nodenumber, inputs) RESULT(yvel)
            ! Just use the fact that strain is constant
            ! so vertical vel is linear in depth
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: yvel, smb, surf, bed
            REAL(kind=dp) :: inputs
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac

            smb = 0.2_dp
            surf = 2000.0_dp
            bed = 0.0_dp

            y=Model % Nodes % y (nodenumber)
            frac = (y - bed) / (surf - bed)

            yvel = -frac * smb
            RETURN 
        END


      FUNCTION VSSA2000M017MY(Model, nodenumber, inputs) RESULT(yvel)
            ! Just use the fact that strain is constant
            ! so vertical vel is linear in depth
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: yvel, smb, surf, bed
            REAL(kind=dp) :: inputs
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac

            smb = 0.1714_dp
            surf = 2000.0_dp
            bed = 0.0_dp

            y=Model % Nodes % y (nodenumber)
            frac = (y - bed) / (surf - bed)

            yvel = -frac * smb
            RETURN 
        END


      FUNCTION VSSA2000M017MY01bm(Model, nodenumber, inputs) RESULT(yvel)
            ! Just use the fact that strain is constant
            ! so vertical vel is linear in depth
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: yvel, smb, surf, bed, bm
            REAL(kind=dp) :: inputs
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac

            bm = 0.01_dp
            smb = 0.1714_dp
            surf = 2000.0_dp
            bed = 0.0_dp

            y=Model % Nodes % y (nodenumber)
            frac = (y - bed) / (surf - bed)

            yvel = -bm-frac * (smb - bm)
            RETURN 
        END


      FUNCTION VSSA2000M017MY001bm(Model, nodenumber, inputs) RESULT(yvel)
            ! Just use the fact that strain is constant
            ! so vertical vel is linear in depth
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: yvel, smb, surf, bed, bm
            REAL(kind=dp) :: inputs
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac

            bm = 0.001_dp
            smb = 0.1714_dp
            surf = 2000.0_dp
            bed = 0.0_dp

            y=Model % Nodes % y (nodenumber)
            frac = (y - bed) / (surf - bed)

            yvel = -bm-frac * (smb - bm)
            RETURN 
        END


      FUNCTION VSSA2000M017MY02bm(Model, nodenumber, inputs) RESULT(yvel)
            ! Just use the fact that strain is constant
            ! so vertical vel is linear in depth
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: yvel, smb, surf, bed, bm
            REAL(kind=dp) :: inputs
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac

            bm = 0.02_dp
            smb = 0.1714_dp
            surf = 2000.0_dp
            bed = 0.0_dp

            y=Model % Nodes % y (nodenumber)
            frac = (y - bed) / (surf - bed)

            yvel = -bm-frac * (smb - bm)
            RETURN 
        END



      FUNCTION VSSA2000M06MY(Model, nodenumber, inputs) RESULT(yvel)
            ! Just use the fact that strain is constant
            ! so vertical vel is linear in depth
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: yvel, smb, surf, bed
            REAL(kind=dp) :: inputs
            INTEGER :: nodenumber
            REAL(kind=dp) :: y, frac

            smb = 0.6_dp
            surf = 2000.0_dp
            bed = 0.0_dp

            y=Model % Nodes % y (nodenumber)
            frac = (y - bed) / (surf - bed)

            yvel = -frac * smb
            RETURN 
        END
