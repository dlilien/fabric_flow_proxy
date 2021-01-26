        FUNCTION TempIni( Model, nodenumber, dumy) RESULT(Temp)
            USE types
            USE read_routines
            USE DefUtils

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: dumy,Temp
            INTEGER :: nodenumber
            REAL,ALLOCATABLE :: dem(:,:),xx(:),yy(:)
            REAL(kind=dp) :: x,y
            INTEGER :: nx,ny
            INTEGER :: i,j
            CHARACTER(len=16) :: glacier
            CHARACTER(len=MAX_NAME_LEN) :: filin

            LOGICAL :: Firsttime=.true.

            SAVE dem,xx,yy,nx,ny
            SAVE Firsttime

            WRITE(Message,'(A)') 'Called'
               CALL INFO('TempIni',Message,Level=11)
            filin = 'temp.xyz'

            IF (Firsttime) THEN
                    Firsttime=.False.
                    call get_twod_grid(filin, xx, yy, dem)
                    WRITE(Message,'(A)') 'Called get_twod_grid'
                        CALL INFO('TempIni',Message,Level=9)
                    nx = SIZE(xx)
                    ny = SIZE(yy)
                    WRITE(Message,'(A)') 'First time complete'
                        CALL INFO('TempIni',Message,Level=3)
            END IF

            ! position current point
            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            Temp=LID(dem,xx,yy,nx,ny,x,y,263.15)
            RETURN 
        END
