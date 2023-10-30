!PG RaFSIP PARAMETERS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!This MODULE holds the subroutines which are used to initialize all  +
!built random forest regressors.                                     +
!This MODULE CONTAINS the following routines:                        +
!  *forestbrhm                                                       +
!  *forestbr                                                         +
!  *forestall                                                        +
!  *forestbrds                                                       +
!  *forestbrwarm                                                     +
!Each subroutine opens, reads and stores the parameters of all 4     +
!random forest regressors. The initial .txt files are first          +
!converted into binary files so that the processing is faster.       +
!                                                                    +
!This module also includes the three subroutines that make all the   +
!random forest predictions needed in the microphysics routine.       +
!These are the following:                                            +
!  *runforest                                                        +
!  *runforestriv                                                     +
!  *runforestmulti                                                   +     
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 MODULE module_random_forests

      use micro_mg_utils, only: r8
      use spmd_utils,     only: masterproc
      use phys_control,   only: use_simple_phys
      use cam_abortutils, only: endrun

      IMPLICIT NONE
      
      PUBLIC :: sec_ice_readnl

      PUBLIC  :: forestbrhm,forestbr,forestall,forestbrds,forestbrwarm,runforest,runforestriv,runforestmulti

      !!MDIM DEFINES THE NUMBER OF FEATURES/INPUTS TO THE RaFSIP PARAMETERIZATION
      INTEGER, PARAMETER, PUBLIC :: MDIM5=5
      INTEGER, PARAMETER, PUBLIC :: MDIM6=6
      INTEGER, PARAMETER, PUBLIC :: JBT=10  !!The number of trees in each random forest regressor

      !!The maximum number of nodes across trees
      INTEGER, PARAMETER, PUBLIC :: MAX_NODES1=7705  !forestBRHM
      INTEGER, PARAMETER, PUBLIC :: MAX_NODES2=8219  !forestBR
      INTEGER, PARAMETER, PUBLIC :: MAX_NODES3=7833  !forestALL
      INTEGER, PARAMETER, PUBLIC :: MAX_NODES4=7093  !forestBRDS
      INTEGER, PARAMETER, PUBLIC :: MAX_NODES5=8593  !forestBRwarm

      !!Thresh = threshold value at each internal node
      !!Outi = prediction for a given node
      REAL(r8), DIMENSION(JBT,MAX_NODES1), PUBLIC    :: THRESH1,OUT11,OUT12,OUT13
      REAL(r8), DIMENSION(JBT,MAX_NODES2), PUBLIC    :: THRESH2,OUT21
      REAL(r8), DIMENSION(JBT,MAX_NODES3), PUBLIC    :: THRESH3,OUT31,OUT32,OUT33,OUT34,OUT35
      REAL(r8), DIMENSION(JBT,MAX_NODES4), PUBLIC    :: THRESH4,OUT41,OUT42,OUT43
      REAL(r8), DIMENSION(JBT,MAX_NODES5), PUBLIC    :: THRESH5,OUT51

      !!Splitfeat = feature used for splitting the node
      !!Leftchild = left child of node 
      !!Rightchild = right child of node
      INTEGER, DIMENSION(JBT,MAX_NODES1), PUBLIC :: SPLITFEAT1,LEFTCHILD1,RIGHTCHILD1
      INTEGER, DIMENSION(JBT,MAX_NODES2), PUBLIC :: SPLITFEAT2,LEFTCHILD2,RIGHTCHILD2
      INTEGER, DIMENSION(JBT,MAX_NODES3), PUBLIC :: SPLITFEAT3,LEFTCHILD3,RIGHTCHILD3
      INTEGER, DIMENSION(JBT,MAX_NODES4), PUBLIC :: SPLITFEAT4,LEFTCHILD4,RIGHTCHILD4
      INTEGER, DIMENSION(JBT,MAX_NODES5), PUBLIC :: SPLITFEAT5,LEFTCHILD5,RIGHTCHILD5

      !!The exact number of nodes across in consecutive trees of the forest
      INTEGER, DIMENSION(JBT) :: NRNODES1,NRNODES2,NRNODES3,NRNODES4,NRNODES5

      LOGICAL, PUBLIC :: FIRST_RAFSIP = .TRUE.
      
      character(len=256), public :: forestfileALL,forestfileBRDS
      character(len=256), public :: forestfileBRHM,forestfileBR
      character(len=256), public :: forestfileBRwarm
 
  CONTAINS


!---------------------------------------------------------------------------------------------------------------


 subroutine sec_ice_readnl(nlfile)
       ! Read files needed for random forest tables of seconary ice formation
    
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input


   ! Local variables
   integer :: unitn, ierr, i
   character(len=2) :: suffix
   character(len=1), pointer   :: ctype(:)
   character(len=*), parameter :: subname = 'sec_ice_readnl'


   namelist /sec_ice_nl/ forestfileALL,    &
                             forestfileBRDS,   &
                             forestfileBRHM,   &
                             forestfileBR,     &
                             forestfileBRwarm

   if (use_simple_phys) return


   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'sec_ice_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, sec_ice_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD

   call mpibcast (forestfileALL,   len(forestfileALL),    mpichar, 0, mpicom)
   call mpibcast (forestfileBRDS,  len(forestfileBRDS),   mpichar, 0, mpicom)
   call mpibcast (forestfileBRHM,  len(forestfileBRHM),   mpichar, 0, mpicom)
   call mpibcast (forestfileBR,    len(forestfileBR),     mpichar, 0, mpicom)
   call mpibcast (forestfileBRwarm,len(forestfileBRwarm), mpichar, 0, mpicom)

#endif
end subroutine sec_ice_readnl

      SUBROUTINE forestbrhm(jbt,max_nodes1,leftchild1,rightchild1,splitfeat1, &
                               nrnodes1,thresh1,out11,out12,out13)
      use units,           only: getunit, freeunit
      IMPLICIT NONE

      INTEGER,intent(in) :: jbt, max_nodes1

      REAL (r8),DIMENSION(jbt,max_nodes1),intent(inout) :: thresh1,out11,out12,out13
      INTEGER,DIMENSION(jbt,max_nodes1),intent(inout) :: splitfeat1,leftchild1,rightchild1
      INTEGER,DIMENSION(jbt),intent(inout) :: nrnodes1

      INTEGER :: jb,n

   integer :: unitn, ierr, i
!      unitn = 137
      open( 137, file=trim(forestfileBRHM), form='formatted',status='old' )
         !Open the ASCII file
!         OPEN(unit=137,file="forestBRHM.txt",status="old",action="read")
!         OPEN(unit=137,file=forestfileBRHM,status="old",action="read")
         DO jb=1,jbt
            read (137,*) nrnodes1(jb)
            read (137,*) (leftchild1(jb,n),rightchild1(jb,n),out11(jb,n),out12(jb,n),out13(jb,n), &
                    & thresh1(jb,n),splitfeat1(jb,n), n=1,nrnodes1(jb))
         ENDDO
         CLOSE(137)
!      call freeunit(unitn)

      END subroutine forestbrhm
!---------------------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------------------
      SUBROUTINE forestbr(jbt,max_nodes2,leftchild2,rightchild2,splitfeat2, &
                               nrnodes2,thresh2,out21)

      IMPLICIT NONE

      INTEGER,intent(in) :: jbt, max_nodes2

      REAL(r8),DIMENSION(jbt,max_nodes2),intent(inout) :: thresh2,out21
      INTEGER,DIMENSION(jbt,max_nodes2),intent(inout) :: splitfeat2,leftchild2,rightchild2
      INTEGER,DIMENSION(jbt),intent(inout) :: nrnodes2

      INTEGER :: jb,n

!         OPEN(unit=138,file="forestBR.txt",status="old",action="read")
         OPEN(unit=138,file=forestfileBR,status="old",action="read")
         DO jb=1,jbt
            read (138,*) nrnodes2(jb)
            read (138,*) (leftchild2(jb,n),rightchild2(jb,n),out21(jb,n), &
                    & thresh2(jb,n),splitfeat2(jb,n), n=1,nrnodes2(jb))
         ENDDO
         CLOSE(138)

      END subroutine forestbr
!---------------------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------------------
      SUBROUTINE forestall(jbt,max_nodes3,leftchild3,rightchild3,splitfeat3, &
                               nrnodes3,thresh3,out31,out32,out33,out34,out35)

      IMPLICIT NONE

      INTEGER,intent(in) :: jbt, max_nodes3

      REAL(r8),DIMENSION(jbt,max_nodes3),intent(inout) :: thresh3,out31,out32,out33,out34,out35
      INTEGER,DIMENSION(jbt,max_nodes3),intent(inout) :: splitfeat3,leftchild3,rightchild3
      INTEGER,DIMENSION(jbt),intent(inout) :: nrnodes3

      INTEGER :: jb,n

!         OPEN(unit=139,file="forestALL.txt",status="old",action="read")
         OPEN(unit=139,file=forestfileALL,status="old",action="read")
         DO jb=1,jbt
            read (139,*) nrnodes3(jb)
            read (139,*) (leftchild3(jb,n),rightchild3(jb,n),out31(jb,n),out32(jb,n),out33(jb,n), &
                    & out34(jb,n),out35(jb,n),thresh3(jb,n),splitfeat3(jb,n), n=1,nrnodes3(jb))
         ENDDO
         CLOSE(139)

      END subroutine forestall
!---------------------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------------------
      SUBROUTINE forestbrds(jbt,max_nodes4,leftchild4,rightchild4,splitfeat4, &
                               nrnodes4,thresh4,out41,out42,out43)

      IMPLICIT NONE

      INTEGER,intent(in) :: jbt, max_nodes4

      REAL(r8),DIMENSION(jbt,max_nodes4),intent(inout) :: thresh4,out41,out42,out43
      INTEGER,DIMENSION(jbt,max_nodes4),intent(inout) :: splitfeat4,leftchild4,rightchild4
      INTEGER,DIMENSION(jbt),intent(inout) :: nrnodes4

      INTEGER :: jb,n

!         OPEN(unit=140,file="forestBRDS.txt",status="old",action="read")
         OPEN(unit=140,file=forestfileBRDS,status="old",action="read")
         DO jb=1,jbt
            read (140,*) nrnodes4(jb)
            read (140,*) (leftchild4(jb,n),rightchild4(jb,n),out41(jb,n),out42(jb,n),out43(jb,n), &
                    & thresh4(jb,n),splitfeat4(jb,n), n=1,nrnodes4(jb))
         ENDDO
         CLOSE(140)

      END subroutine forestbrds
!---------------------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------------------
      SUBROUTINE forestbrwarm(jbt,max_nodes5,leftchild5,rightchild5,splitfeat5, &
                               nrnodes5,thresh5,out51)

      IMPLICIT NONE

      INTEGER,intent(in) :: jbt, max_nodes5

      REAL(r8),DIMENSION(jbt,max_nodes5),intent(inout) :: thresh5,out51
      INTEGER,DIMENSION(jbt,max_nodes5),intent(inout) :: splitfeat5,leftchild5,rightchild5
      INTEGER,DIMENSION(jbt),intent(inout) :: nrnodes5

      INTEGER :: jb,n

!         OPEN(unit=141,file="forestBRwarm.txt",status="old",action="read")
         OPEN(unit=141,file=forestfileBRwarm,status="old",action="read")
         DO jb=1,jbt
            read (141,*) nrnodes5(jb)
            read (141,*) (leftchild5(jb,n),rightchild5(jb,n),out51(jb,n), &
                    & thresh5(jb,n),splitfeat5(jb,n), n=1,nrnodes5(jb))
         ENDDO
         CLOSE(141)

      END subroutine forestbrwarm
!---------------------------------------------------------------------------------------------------------------



!======================================================================+
!   THREE SUBROUTINES CALLED BY THE RaFSIP PARAMETERIZATION            !
!======================================================================+

     !This subroutine is called only when the requirements for the
     !activation of the forestBR model are met (i.e., -25<T<-8 in
     !the absence of rainwater). In this case only the BR process
     !can contribute to the ice production, and hence the RF gives
     !only one prediction: log(IEFBR)
     SUBROUTINE runforest(mdim,max_nodes,jbt,features,ypred1,leftchild,rightchild, &
        & splitfeat,thresh,out1)

      IMPLICIT NONE
      integer,intent(in) :: jbt,mdim,max_nodes
      integer,dimension(jbt,max_nodes),intent(in) :: splitfeat,leftchild,rightchild
      real(r8),dimension(jbt,max_nodes),intent(in)    :: out1,thresh
      real(r8),dimension(mdim),intent(in) :: features
      real(r8), intent(out)   :: ypred1
      integer :: jb,inode,next_node

!      PRINT*, "runforest_feature",features

      ! Initialize variables
      ypred1 = 0._r8

      ! START DOWN FOREST TO CALCULATE THE PREDICTED VALUES
      ! loop over trees in forest
      DO jb=1,jbt
         
         ! set current node to root node         
         inode = 1

         ! loop as long as we reach a leaf node
         do while (leftchild(jb,inode) .ne. rightchild(jb,inode))
             if (features(splitfeat(jb,inode)).le.thresh(jb,inode)) then
                next_node = leftchild(jb,inode)
             else
                next_node = rightchild(jb,inode)
             endif

             inode = next_node

          enddo  !do while

          YPRED1 = YPRED1 + out1(jb,inode)

      ENDDO  !tree loop


      YPRED1 = YPRED1/jbt  !YPRED1=log10(IEFBR)


      End subroutine runforest


!-------------------------------------------------------------------------------------+
     !This subroutine is called when the requirements for either the forestBRDS
     !or the forestBRHM are met (i.e., -25<T<-8 in the presence of raindrops or 
     !-8<T<-3 without raindrops). In this case the RF gives in total 3
     !predictions: log(IEFBR), log(IEFDS) or log(IEFHM), log(QIRSIP) or log(QICSIP) 
      SUBROUTINE runforestriv(mdim,max_nodes,jbt,features,ypred1,ypred2,ypred3, &
            & leftchild,rightchild,splitfeat,thresh,out1,out2,out3)

      IMPLICIT NONE
      integer,intent(in) :: jbt,mdim,max_nodes
      integer,dimension(jbt,max_nodes),intent(in) :: splitfeat,leftchild,rightchild
      real(r8),dimension(jbt,max_nodes),intent(in)    :: out1,out2,out3,thresh
      real(r8),dimension(mdim),intent(in) :: features
      real(r8), intent(out)   :: ypred1,ypred2,ypred3
      integer :: jb,inode,next_node

!      PRINT*, "runforestriv_feature",features

      ! Initialize variables
      ypred1 = 0._r8
      ypred2 = 0._r8
      ypred3 = 0._r8

      ! START DOWN FOREST TO CALCULATE THE PREDICTED VALUES
      ! loop over trees in forest
      DO jb=1,jbt

          ! set current node to root node
          inode = 1

          ! loop as long as we reach a leaf node
          do while (leftchild(jb,inode) .ne. rightchild(jb,inode))
             if (features(splitfeat(jb,inode)).le.thresh(jb,inode)) then
                next_node = leftchild(jb,inode)
             else 
                next_node = rightchild(jb,inode)
             endif

             inode = next_node
       
          enddo  !do while

          YPRED1 = YPRED1 + out1(jb,inode)
          YPRED2 = YPRED2 + out2(jb,inode)  
          YPRED3 = YPRED3 + out3(jb,inode)

      ENDDO  !tree loop


      YPRED1 = YPRED1/jbt  !YPRED1=log10(IEFBR)
      YPRED2 = YPRED2/jbt  !YPRED2=log10(IEFDS) or log10(IEFHM)
      YPRED3 = YPRED3/jbt  !YPRED3=log10(QIRSIP) or log10(QICSIP)


      End subroutine runforestriv


!-------------------------------------------------------------------------------------+
     !This subroutine is called when the requirements for the forestALL are met 
     !(i.e., -8<T<-3 in the presence of raindrops). In this case the RF gives 5
     !predictions: log(IEFBR), log(IEFHM), log(IEFDS), log(QICSIP), log(QIRSIP)
      SUBROUTINE runforestmulti(mdim,max_nodes,jbt,features,ypred1,ypred2,ypred3,ypred4,ypred5, &
        & leftchild,rightchild,splitfeat,thresh,out1,out2,out3,out4,out5)

      IMPLICIT NONE
      integer,intent(in) :: jbt,mdim,max_nodes
      integer,dimension(jbt,max_nodes),intent(in) :: splitfeat,leftchild,rightchild
      real(r8),dimension(jbt,max_nodes),intent(in)    :: out1,out2,out3,out4,out5,thresh
      real(r8),dimension(mdim),intent(in) :: features
      real(r8), intent(out)   :: ypred1,ypred2,ypred3,ypred4,ypred5
      integer :: jb,inode,next_node

!      PRINT*, "runforest_multi",features

      ! Initialize variables
      ypred1 = 0._r8
      ypred2 = 0._r8
      ypred3 = 0._r8
      ypred4 = 0._r8
      ypred5 = 0._r8

      ! START DOWN FOREST TO CALCULATE THE PREDICTED VALUES
      ! loop over trees in forest
      DO jb=1,jbt

!         PRINT*, "forestALL_nrnodes", nrnodes3(jb)

          ! set current node to root node
          inode = 1

          ! loop as long as we reach a leaf node
          do while (leftchild(jb,inode) .ne. rightchild(jb,inode))
            if (features(splitfeat(jb,inode)).le.thresh(jb,inode)) then
                next_node = leftchild(jb,inode)
             else
                next_node = rightchild(jb,inode)
             endif

             inode = next_node

          enddo  !do while

          YPRED1 = YPRED1 + out1(jb,inode)
          YPRED2 = YPRED2 + out2(jb,inode)
          YPRED3 = YPRED3 + out3(jb,inode)
          YPRED4 = YPRED4 + out4(jb,inode)
          YPRED5 = YPRED5 + out5(jb,inode)

      ENDDO  !tree loop


      YPRED1 = YPRED1/jbt  !YPRED1=log10(IEFBR)
      YPRED2 = YPRED2/jbt  !YPRED2=log10(IEFHM)
      YPRED3 = YPRED3/jbt  !YPRED3=log10(IEFDS)
      YPRED4 = YPRED4/jbt  !YPRED4=log10(QICSIP)
      YPRED5 = YPRED5/jbt  !YPRED5=log10(QIRSIP)


      End subroutine runforestmulti

!+---+-----------------------------------------------------------------+



 END module module_random_forests
