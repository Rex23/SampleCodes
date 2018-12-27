MODULE Global
IMPLICIT NONE

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Sub_StrainsXZ1_R,Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,Sub_StrainsYZ2_R
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: sigma_11_Eff,sigma_22_Eff,sigma_12_Eff,rho_11_Eff,rho_22_Eff,rho_12_Eff
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ElementType

DOUBLE PRECISION :: MacroStrainsXX,MacroStrainsYY,MacroStrainsZZ,MacroStrainsXY1,MacroStrainsXY2,MacroStrainsXZ1,MacroStrainsXZ2,&
MacroStrainsYZ1,MacroStrainsYZ2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: Nodes
INTEGER, ALLOCATABLE, DIMENSION(:) :: NodeIndex
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)  :: CurrentNodes 
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: Connectivities
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: Boundaries
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: BoundaryIndex
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: CC
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) :: Sub_CC
INTEGER, ALLOCATABLE, DIMENSION (:) :: PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: Vol
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) ::  DNDX,DNDY,DNDZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: KComponent
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:)  :: KIC
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:)   :: KII
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: K1GRR
INTEGER, ALLOCATABLE, DIMENSION (:) :: K1GRCDD
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: K1GRCII
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: F1Globall
INTEGER :: NumberofDE
INTEGER :: Highest_GRCD
INTEGER :: TestIndex1,TestIndex2,BeginofModeIndex,EndofModeIndex,Anisotropic
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: RotationMatrix
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: RotationMatrix2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: E_Boundaries
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: E_BoundaryIndex
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: E_Connectivities
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: E_sub_sigmasigma_in
DOUBLE PRECISION :: MacroEx_Knot,MacroEy_Knot,MacroEz_Knot
INTEGER :: E_Highest_GRCD
INTEGER  :: i,j
INTEGER :: E_ZeroTimeStep
DOUBLE PRECISION :: E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ2_Eff
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ElementType2
INTEGER, ALLOCATABLE, DIMENSION(:) :: ElementTypeIndex
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: E_rotated_sigma_temp,E_sigma_temp,E_sigma_ori
INTEGER :: i1,i2,i3,i4
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM main
USE Global
IMPLICIT NONE
INTEGER :: k,ii,ii1,ii2
INTEGER :: M_Highest_GRCD,ME_Highest_GRCD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL ReadControlFlags ()

ALLOCATE(Nodes(NumberofNodes,4))
ALLOCATE(NodeIndex(NumberofNodes))
CALL ReadNodes (NodeIndex,Nodes)

ALLOCATE(Connectivities(NumberofElements,17))
CALL ReadConnectivities (Connectivities)

ALLOCATE(Boundaries(NumberofBoundaryNodes,TotalTimeSteps+2))
ALLOCATE(BoundaryIndex(NumberofBoundaryNodes,2))
CALL ReadBoundaryConditions (BoundaryIndex,Boundaries)

ALLOCATE(CC(NumberofMaterialModes,6,6))
ALLOCATE(Sub_CC(NumberofGaussPoints,NumberofElements,6,6))
CALL ReadMechanicalProperties (CC)

ALLOCATE(PositiveXNodes(NumberofPositiveXNodes))
ALLOCATE(PositiveYNodes(NumberofPositiveYNodes))
ALLOCATE(PositiveZNodes(NumberofPositiveZNodes))
ALLOCATE(NegativeXNodes(NumberofNegativeXNodes))
ALLOCATE(NegativeYNodes(NumberofNegativeYNodes))
ALLOCATE(NegativeZNodes(NumberofNegativeZNodes))
CALL SixFaces(Nodes,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes)

ALLOCATE(Vol(NumberofGaussPoints,NumberofElements))
ALLOCATE(KComponent(NumberofElements,24,24))
ALLOCATE(DNDX(NumberofGaussPoints,NumberofElements,8+3))
ALLOCATE(DNDY(NumberofGaussPoints,NumberofElements,8+3))
ALLOCATE(DNDZ(NumberofGaussPoints,NumberofElements,8+3))
ALLOCATE(RotationMatrix(NumberofMaterialModes,3,3))
ALLOCATE(RotationMatrix2(NumberofMaterialModes,3,3))
ALLOCATE(KIC(NumberofElements,9,24))
ALLOCATE(KII(NumberofElements,9,9))
CALL Components(Nodes,Connectivities,CC,Sub_CC,Vol,KComponent,DNDX,DNDY,DNDZ,TestIndex1,TestIndex2,BeginofModeIndex,EndofModeIndex,&
Anisotropic,RotationMatrix,RotationMatrix2,KIC,KII)
NumberofDE=3*(NumberofPositiveXNodes+NumberofPositiveYNodes+NumberofPositiveZNodes)
ALLOCATE(K1GRR(3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec))
ALLOCATE(K1GRCDD(3*NumberofNodes+NumberofDE))
ALLOCATE(K1GRCII(3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec))
ALLOCATE(F1Globall(3*NumberofNodes+NumberofDE))
CALL KPreprations(Connectivities,KComponent,NumberofDE,K1GRR,K1GRCDD,K1GRCII,F1Globall)

ALLOCATE (Sub_StrainsXX_R(NumberofElements,8))
ALLOCATE (Sub_StrainsYY_R(NumberofElements,8))
ALLOCATE (Sub_StrainsZZ_R(NumberofElements,8))
ALLOCATE (Sub_StrainsXY1_R(NumberofElements,8))
ALLOCATE (Sub_StrainsXY2_R(NumberofElements,8))
ALLOCATE (Sub_StrainsXZ1_R(NumberofElements,8))
ALLOCATE (Sub_StrainsXZ2_R(NumberofElements,8))
ALLOCATE (Sub_StrainsYZ1_R(NumberofElements,8))
ALLOCATE (Sub_StrainsYZ2_R(NumberofElements,8))
ALLOCATE (CurrentNodes(NumberofNodes,4))
!CALL ReadStrains (Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,Sub_StrainsXY2,&
!Sub_StrainsXZ1,Sub_StrainsXZ2,Sub_StrainsYZ1,Sub_StrainsYZ2)
ALLOCATE (ElementType(NumberofElements,9))
ALLOCATE (ElementType2(NumberofRVETypes,NumberofElements*NumberofGaussPoints,2))
ALLOCATE (ElementTypeIndex(NumberofRVETypes))
CALL ReadElementType(ElementType,ElementType2,ElementTypeIndex)

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
CALL E_ReadControlFlags()
ALLOCATE(E_Boundaries(E_NumberofBoundaryNodes,TotalTimeSteps+1))
ALLOCATE(E_BoundaryIndex(E_NumberofBoundaryNodes,1))
CALL E_ReadBoundaryConditions(E_BoundaryIndex,E_Boundaries)
ALLOCATE(E_sub_sigmasigma_in(TotalTimeSteps+1,E_NumberofGaussPoints,NumberofElements,3,3))
ALLOCATE(E_rotated_sigma_temp(3,3))
ALLOCATE(E_sigma_temp(3,3))
ALLOCATE(E_sigma_ori(3,3))
!E_sub_sigmasigma=0.0D0
!DO j=1,E_NumberofGaussPoints
!DO i=1,NumberofElements
!E_sub_sigmasigma(j,i,1,1)=1.0E5
!E_sub_sigmasigma(j,i,2,2)=1.0E5
!E_sub_sigmasigma(j,i,3,3)=1.0E5
!ENDDO
!ENDDO
IF (TestType==2 .OR. TestType==3) THEN
ALLOCATE(E_Connectivities(NumberofElements,17))
CALL E_ReadElectroPropConnec (E_sub_sigmasigma_in,E_Connectivities)
!ELSEIF (TestType==3) THEN
!ALLOCATE(E_Connectivities(NumberofElements,17))
!E_Connectivities=Connectivities
ELSE
ENDIF
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!

DO k=1,TotalTimeSteps

CALL Macro_Mechanical(k,Nodes,NodeIndex,Connectivities,DNDX,DNDY,DNDZ,Sub_CC,CC,Vol,Boundaries,BoundaryIndex,NumberofDE,CurrentNodes,&
KIC,KII,K1GRR,K1GRCDD,K1GRCII,F1Globall,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,&
NegativeYNodes,NegativeZNodes,TestIndex1,TestIndex2,BeginofModeIndex,EndofModeIndex,Anisotropic,RotationMatrix,Sub_StrainsXX_R,&
Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXZ1_R,Sub_StrainsYZ1_R,Sub_StrainsXY2_R,Sub_StrainsXZ2_R,Sub_StrainsYZ2_R,Highest_GRCD)

IF (TestType==3 .OR.TestType==4) THEN
!E_sub_sigmasigma_in=0.0D0
DO ii=1,NumberofRVETypes
IF (ii==3) THEN
DO ii1=1,ElementTypeIndex(ii)
IF (k==1) THEN
!E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=2.4875621890547300000D-9
!E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=2.4875621890547300000D-9
!E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=2.4875621890547300000D-9
!E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=2.4875621890547300000D-9
!E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=2.4875621890547300000D-9
!E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=2.4875621890547300000D-9

E_sigma_ori(1,1)=E_sub_sigmasigma_in(1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)
E_sigma_ori(2,2)=E_sub_sigmasigma_in(1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)
E_sigma_ori(3,3)=E_sub_sigmasigma_in(1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)
E_sigma_ori(1,2)=E_sub_sigmasigma_in(1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,3)
E_sigma_ori(1,3)=E_sub_sigmasigma_in(1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,3)
E_sigma_ori(2,3)=E_sub_sigmasigma_in(1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,2)
E_sigma_ori(2,1)=E_sigma_ori(1,2)
E_sigma_ori(3,1)=E_sigma_ori(1,3)
E_sigma_ori(3,2)=E_sigma_ori(2,3)

E_sigma_temp=E_sigma_ori

E_rotated_sigma_temp=0.0D0
DO i1=1,3
DO i2=1,3
DO i3=1,3
DO i4=1,3
E_rotated_sigma_temp(i1,i2)=E_rotated_sigma_temp(i1,i2)+RotationMatrix2(Connectivities(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)+9)-BeginofModeIndex+1,i1,i3)*&
RotationMatrix2(Connectivities(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)+9)-BeginofModeIndex+1,i2,i4)*E_sigma_temp(i3,i4)
ENDDO
ENDDO
ENDDO
ENDDO
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=E_rotated_sigma_temp(1,1)
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,2)=E_rotated_sigma_temp(1,2)
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,3)=E_rotated_sigma_temp(1,3)
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=E_rotated_sigma_temp(2,2)
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,3)=E_rotated_sigma_temp(2,3)
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=E_rotated_sigma_temp(3,3)
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,1)=E_rotated_sigma_temp(2,1)
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,1)=E_rotated_sigma_temp(3,1)
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,2)=E_rotated_sigma_temp(3,2)

!0.0286199D0
!IF ( (1.0D-5+0.0286199D0*Sub_StrainsXX_R(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2))) > 0.0D0 ) THEN
!E_sigma_temp(1,1)=E_sigma_ori(1,1)*1.0D-5/(1.0D-5+0.0286199D0*Sub_StrainsXX_R(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2))) !*1.0D-5/(1.0D-5+8.60615D-5*Sub_StrainsXX_R(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)))
!E_sigma_temp(1,1)=E_sigma_ori(1,1)*1.0D-5/(1.0D-5+8.60615D-5*Sub_StrainsXX_R(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)))
E_sigma_temp(1,1)=E_sigma_ori(1,1)
!ELSE
!E_sigma_temp(1,1)=E_sigma_ori(1,1)*1.0D-5/1.0D-6
!ENDIF
E_sigma_temp(2,2)=E_sigma_ori(2,2)
E_sigma_temp(3,3)=E_sigma_ori(3,3)
E_sigma_temp(1,2)=E_sigma_ori(1,2)
E_sigma_temp(1,3)=E_sigma_ori(1,3)
E_sigma_temp(2,3)=E_sigma_ori(2,3)
E_sigma_temp(2,1)=E_sigma_temp(1,2)
E_sigma_temp(3,1)=E_sigma_temp(1,3)
E_sigma_temp(3,2)=E_sigma_temp(2,3)


E_rotated_sigma_temp=0.0D0
DO i1=1,3
DO i2=1,3
DO i3=1,3
DO i4=1,3
E_rotated_sigma_temp(i1,i2)=E_rotated_sigma_temp(i1,i2)+RotationMatrix2(Connectivities(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)+9)-BeginofModeIndex+1,i1,i3)*&
RotationMatrix2(Connectivities(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)+9)-BeginofModeIndex+1,i2,i4)*E_sigma_temp(i3,i4)
ENDDO
ENDDO
ENDDO
ENDDO
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=E_rotated_sigma_temp(1,1)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,2)=E_rotated_sigma_temp(1,2)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,3)=E_rotated_sigma_temp(1,3)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=E_rotated_sigma_temp(2,2)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,3)=E_rotated_sigma_temp(2,3)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=E_rotated_sigma_temp(3,3)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,1)=E_rotated_sigma_temp(2,1)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,1)=E_rotated_sigma_temp(3,1)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,2)=E_rotated_sigma_temp(3,2)

ELSE
!E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=2.4875621890547300000D-9
!E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=2.4875621890547300000D-9
!E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=2.4875621890547300000D-9
!0.0286199D0
!IF ( (1.0D-5+0.0286199D0*Sub_StrainsXX_R(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2))) > 0.0D0 ) THEN
!E_sigma_temp(1,1)=E_sigma_ori(1,1)*1.0D-5/(1.0D-5+0.0286199D0*Sub_StrainsXX_R(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)))  !*1.0D-5/(1.0D-5+8.60615D-5*Sub_StrainsXX_R(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)))
!E_sigma_temp(1,1)=E_sigma_ori(1,1)*1.0D-5/(1.0D-5+8.60615D-5*Sub_StrainsXX_R(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)))
E_sigma_temp(1,1)=E_sigma_ori(1,1)
!ELSE
!E_sigma_temp(1,1)=E_sigma_ori(1,1)*1.0D-5/1.0D-6
!ENDIF
E_sigma_temp(2,2)=E_sigma_ori(2,2)
E_sigma_temp(3,3)=E_sigma_ori(3,3)
E_sigma_temp(1,2)=E_sigma_ori(1,2)
E_sigma_temp(1,3)=E_sigma_ori(1,3)
E_sigma_temp(2,3)=E_sigma_ori(2,3)
E_sigma_temp(2,1)=E_sigma_temp(1,2)
E_sigma_temp(3,1)=E_sigma_temp(1,3)
E_sigma_temp(3,2)=E_sigma_temp(2,3)

E_rotated_sigma_temp=0.0D0
DO i1=1,3
DO i2=1,3
DO i3=1,3
DO i4=1,3
E_rotated_sigma_temp(i1,i2)=E_rotated_sigma_temp(i1,i2)+RotationMatrix2(Connectivities(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)+9)-BeginofModeIndex+1,i1,i3)*&
RotationMatrix2(Connectivities(ElementType2(ii,ii1,1),ElementType2(ii,ii1,2)+9)-BeginofModeIndex+1,i2,i4)*E_sigma_temp(i3,i4)
ENDDO
ENDDO
ENDDO
ENDDO
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=E_rotated_sigma_temp(1,1)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,2)=E_rotated_sigma_temp(1,2)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,3)=E_rotated_sigma_temp(1,3)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=E_rotated_sigma_temp(2,2)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,3)=E_rotated_sigma_temp(2,3)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=E_rotated_sigma_temp(3,3)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,1)=E_rotated_sigma_temp(2,1)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,1)=E_rotated_sigma_temp(3,1)
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,2)=E_rotated_sigma_temp(3,2)
ENDIF

ENDDO
ELSEIF (ii==1) THEN
CALL Macro_Each_RVE_Type(k,ii,ElementType,ElementType2,ElementTypeIndex,Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R,&
Sub_StrainsXZ1_R,Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,Sub_StrainsYZ2_R,M_Highest_GRCD,ME_Highest_GRCD,E_sub_sigmasigma_in,RotationMatrix2,BeginofModeIndex,&
EndofModeIndex,Connectivities)
ELSEIF (ii==2) THEN
DO ii1=1,ElementTypeIndex(ii)
IF (k==1) THEN
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=1.0D-15 !2.48756218905473D-13
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=1.0D-15 !2.48756218905473D-13
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=1.0D-15 !2.48756218905473D-13
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=1.0D-15 !2.48756218905473D-13
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=1.0D-15 !2.48756218905473D-13
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=1.0D-15 !2.48756218905473D-13
ELSE
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=1.0D-15 !2.48756218905473D-13
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=1.0D-15 !2.48756218905473D-13
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=1.0D-15 !2.48756218905473D-13
ENDIF
ENDDO
ELSEIF(ii==4) THEN
DO ii1=1,ElementTypeIndex(ii)
IF (k==1) THEN
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=1.0D-15 !1.0D-15
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=1.0D-15 !1.0D-15
E_sub_sigmasigma_in(k,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=1.0D-15 !1.0D-15
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=1.0D-15 !1.0D-15
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=1.0D-15 !1.0D-15
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=1.0D-15 !1.0D-15
ELSE
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),1,1)=1.0D-15 !1.0D-15
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),2,2)=1.0D-15 !1.0D-15
E_sub_sigmasigma_in(k+1,ElementType2(ii,ii1,2),ElementType2(ii,ii1,1),3,3)=1.0D-15 !1.0D-15
ENDIF
ENDDO
ENDIF
ENDDO !For ii
ENDIF !For IF (TestType==3 .OR.TestType==4) THEN

IF (TestType==2 .OR. TestType==3) THEN

IF (k==1) THEN
E_ZeroTimeStep=1

DO ii=1,6
IF (ii==1) THEN ! For sigma_xx and sigma_zz (the rule of mixtures)
MacroEx_Knot=1.0D0
MacroEy_Knot=0.0D0
MacroEz_Knot=0.0D0
ELSEIF (ii==3) THEN ! For sigma_yy
MacroEx_Knot=0.0D0
MacroEy_Knot=1.0D0
MacroEz_Knot=0.0D0
ELSEIF (ii==4) THEN ! For sigma_xy
MacroEx_Knot=1.0D0
MacroEy_Knot=1.0D0
MacroEz_Knot=0.0D0
ELSEIF (ii==5) THEN ! For sigma_xz
MacroEx_Knot=1.0D0
MacroEy_Knot=0.0D0
MacroEz_Knot=1.0D0
ELSEIF (ii==6) THEN ! For sigma_yz
MacroEx_Knot=0.0D0
MacroEy_Knot=1.0D0
MacroEz_Knot=1.0D0
ELSEIF (ii==2) THEN ! For sigma_zz2
MacroEx_Knot=0.0D0
MacroEy_Knot=0.0D0
MacroEz_Knot=1.0D0
ENDIF

CALL Macro_Electrostatic(E_ZeroTimeStep,k,ii,MacroEx_Knot,MacroEy_Knot,MacroEz_Knot,E_sub_sigmasigma_in,E_Boundaries,E_BoundaryIndex,Nodes,NodeIndex,CurrentNodes,&
Connectivities,E_Connectivities,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes,E_Highest_GRCD,&
E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ2_Eff)

ENDDO !For ii

ENDIF

E_ZeroTimeStep=0

DO ii=1,6
IF (ii==1) THEN ! For sigma_xx and sigma_zz (the rule of mixtures)
MacroEx_Knot=1.0D0
MacroEy_Knot=0.0D0
MacroEz_Knot=0.0D0
ELSEIF (ii==3) THEN ! For sigma_yy
MacroEx_Knot=0.0D0
MacroEy_Knot=1.0D0
MacroEz_Knot=0.0D0
ELSEIF (ii==4) THEN ! For sigma_xy
MacroEx_Knot=1.0D0
MacroEy_Knot=1.0D0
MacroEz_Knot=0.0D0
ELSEIF (ii==5) THEN ! For sigma_xz
MacroEx_Knot=1.0D0
MacroEy_Knot=0.0D0
MacroEz_Knot=1.0D0
ELSEIF (ii==6) THEN ! For sigma_yz
MacroEx_Knot=0.0D0
MacroEy_Knot=1.0D0
MacroEz_Knot=1.0D0
ELSEIF (ii==2) THEN ! For sigma_zz2
MacroEx_Knot=0.0D0
MacroEy_Knot=0.0D0
MacroEz_Knot=1.0D0
ENDIF

CALL Macro_Electrostatic(E_ZeroTimeStep,k,ii,MacroEx_Knot,MacroEy_Knot,MacroEz_Knot,E_sub_sigmasigma_in,E_Boundaries,E_BoundaryIndex,Nodes,NodeIndex,CurrentNodes,&
Connectivities,E_Connectivities,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes,E_Highest_GRCD,&
E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ2_Eff)

ENDDO !For ii


ENDIF !IF (TestType==2 .OR. TestType==3) THEN

IF (TestType==1) THEN
WRITE(*,*) "Highest Column Dimension  ",Highest_GRCD
ELSEIF (TestType==2) THEN
WRITE(*,*) "Highest Column Dimension  ",Highest_GRCD,E_Highest_GRCD
ELSEIF (TestType==3) THEN
WRITE(*,*) "Highest Column Dimension  ",Highest_GRCD,E_Highest_GRCD,M_Highest_GRCD,ME_Highest_GRCD
ELSEIF (TestType==4) THEN
WRITE(*,*) "Highest Column Dimension  ",Highest_GRCD,M_Highest_GRCD,ME_Highest_GRCD
ENDIF

ENDDO !For k

DO ii=TestIndex1,TestIndex2
CALL TecplotGenerator (ii)
ENDDO

IF (TestType==2 .OR. TestType==3) THEN
CALL E_TecplotGenerator (1)
CALL E_TecplotGenerator (2)
CALL E_TecplotGenerator (3)
CALL E_TecplotGenerator (4)
CALL E_TecplotGenerator (5)
ENDIF


ENDPROGRAM main

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EffectiveProperties(k,m,MacroStrainsXX,MacroStrainsYY,MacroStrainsZZ,MacroStrainsXY1,MacroStrainsXZ1,MacroStrainsYZ1,MacroStrainsXY2,&
MacroStrainsXZ2,MacroStrainsYZ2,Vol,Vol_Ele,Nodes,Connectivities,CC,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,&
NegativeYNodes,NegativeZNodes,StrainsXX,StrainsYY,StrainsZZ,StrainsXY,StrainsXZ,StrainsYZ,StressesXX,StressesYY,StressesZZ,&
StressesXY,StressesXZ,StressesYZ,Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,Sub_StrainsXZ1,Sub_StrainsYZ1,&
Sub_StrainsXY2,Sub_StrainsXZ2,Sub_StrainsYZ2,Sub_StressesXX,Sub_StressesYY,Sub_StressesZZ,Sub_StressesXY,Sub_StressesXZ,Sub_StressesYZ,&
Kappa12,Mu12,Mu23,E33,Nu32,E33_Real,C11,C12,C33)
IMPLICIT NONE
INTEGER :: i,j,k1,OUT1

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(IN) :: k,m
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveXNodes) :: PositiveXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveYNodes) :: PositiveYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveZNodes) :: PositiveZNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeXNodes) :: NegativeXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeYNodes) :: NegativeYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeZNodes) :: NegativeZNodes

DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsXX
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsYY
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsXY
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsXZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsYZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesXX
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesYY
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesXY
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesXZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesYZ

DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements,8) :: Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,Sub_StrainsXZ1,Sub_StrainsYZ1,&
Sub_StrainsXY2,Sub_StrainsXZ2,Sub_StrainsYZ2,Sub_StressesXX,Sub_StressesYY,Sub_StressesZZ,Sub_StressesXY,Sub_StressesXZ,Sub_StressesYZ

DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofMaterialModes,6,6) :: CC
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements) :: Vol_Ele
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofGaussPoints,NumberofElements) :: Vol
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: Nodes

DOUBLE PRECISION :: Lx,Ly,Lz
DOUBLE PRECISION, INTENT(INOUT) :: Kappa12,Mu12,Mu23,E33,Nu32,E33_Real,C11,C12,C33

DOUBLE PRECISION, DIMENSION(NumberofElements,NumberofGaussPoints) :: Kappa12_Input,Nu32_Input,E33_Input

DOUBLE PRECISION :: Vol_Solid

DOUBLE PRECISION :: Vol_Ave_SigmaEpsilon

DOUBLE PRECISION, INTENT(IN) :: MacroStrainsXX,MacroStrainsYY,MacroStrainsZZ,MacroStrainsXY1,MacroStrainsXZ1,MacroStrainsYZ1,MacroStrainsXY2,&
MacroStrainsXZ2,MacroStrainsYZ2

OUT1=20

!Kappa12=0.0D0
!Mu12=0.0D0
!Mu23=0.0D0
E33=0.0D0
!Nu32=0.0D0
Vol_Solid=0.0D0
Vol_Ave_SigmaEpsilon=0.0D0

Lx=Nodes(PositiveXNodes(1),2)-Nodes(NegativeXNodes(1),2)
Ly=Nodes(PositiveYNodes(1),3)-Nodes(NegativeYNodes(1),3)
Lz=Nodes(PositiveZNodes(1),4)-Nodes(NegativeZNodes(1),4)

DO k1=1,NumberofElements
DO i=1,NumberofGaussPoints
Kappa12_Input(k1,i)=(CC(Connectivities(k1,9+i),1,1)+CC(Connectivities(k1,9+i),1,2))/2.0D0
Nu32_Input(k1,i)=CC(Connectivities(k1,9+i),1,3)/2.0D0/Kappa12_Input(k1,i)
E33_Input(k1,i)=CC(Connectivities(k1,9+i),3,3)-4.0D0*Nu32_Input(k1,i)*Nu32_Input(k1,i)*Kappa12_Input(k1,i)
ENDDO
ENDDO

DO k1=1,NumberofElements
DO i=1,NumberofGaussPoints
E33=E33+E33_Input(k1,i)*Vol(i,k1)/(Lx*Ly*Lz)
ENDDO
ENDDO

DO k1=1,NumberofElements
DO i=1,NumberofGaussPoints

Vol_Ave_SigmaEpsilon=Vol_Ave_SigmaEpsilon+0.5D0*(Sub_StrainsXX(k1,i)*Sub_StressesXX(k1,i)+Sub_StrainsYY(k1,i)*Sub_StressesYY(k1,i)+&
Sub_StrainsZZ(k1,i)*Sub_StressesZZ(k1,i)+2.0D0*(Sub_StrainsXY1(k1,i)+Sub_StrainsXY2(k1,i))*Sub_StressesXY(k1,i)+&
2.0D0*(Sub_StrainsXZ1(k1,i)+Sub_StrainsXZ2(k1,i))*Sub_StressesXZ(k1,i)+2.0D0*(Sub_StrainsYZ1(k1,i)+Sub_StrainsYZ2(k1,i))*Sub_StressesYZ(k1,i) )*Vol(i,k1)

Vol_Solid=Vol_Solid+Vol(i,k1)

ENDDO
ENDDO

Vol_Ave_SigmaEpsilon=Vol_Ave_SigmaEpsilon/(Lx*Ly*Lz)

OPEN(OUT1,FILE="3D_Output_Micro_F/3D_Output/EffectiveProperties.txt",STATUS="UNKNOWN")

IF (m==1) THEN

Mu12=Vol_Ave_SigmaEpsilon/((MacroStrainsXY1+MacroStrainsXY2)*(MacroStrainsXY1+MacroStrainsXY2))/2.0D0

WRITE(OUT1,24) "Mu12",Mu12,"E33_VA",E33

24 FORMAT(A6,2X,E25.15,2X,A9,2X,E25.15,/)

ELSEIF (m==2) THEN

Mu23=Vol_Ave_SigmaEpsilon/((MacroStrainsYZ1+MacroStrainsYZ2)*(MacroStrainsYZ1+MacroStrainsYZ2))/2.0D0

WRITE(OUT1,25) "Mu12",Mu12,"Mu23",Mu23,"E33_VA",E33

25 FORMAT(A6,2X,E25.15,2X,A6,2X,E25.15,A9,2X,E25.15,/)

ELSEIF (m==3) THEN

C11=Vol_Ave_SigmaEpsilon/0.5D0/(MacroStrainsXX*MacroStrainsXX)

Kappa12=Vol_Ave_SigmaEpsilon/0.5D0/(MacroStrainsXX*MacroStrainsXX)-Mu12

WRITE(OUT1,26) "Mu12",Mu12,"Mu23",Mu23,"Kappa12",Kappa12,"E33_VA",E33

26 FORMAT(A6,2X,E25.15,2X,A6,2X,E25.15,2X,A9,2X,E25.15,2X,A9,2X,E25.15,/)

ELSEIF (m==4) THEN

C33=Vol_Ave_SigmaEpsilon/0.5D0/(MacroStrainsZZ*MacroStrainsZZ)

!Nu32=SQRT( ( Vol_Ave_SigmaEpsilon/0.5D0/(MacroStrainsZZ*MacroStrainsZZ)-E33 )/4.0D0/Kappa12 )

ELSEIF (m==5) THEN

C12=(Vol_Ave_SigmaEpsilon/0.5D0-C11*MacroStrainsXX*MacroStrainsXX-C33*MacroStrainsZZ*MacroStrainsZZ)/2.0D0/&
(MacroStrainsXX*MacroStrainsZZ)

Nu32=C12/2.0D0/Kappa12

E33_Real=C33-4.0D0*Nu32*Nu32*Kappa12

WRITE(OUT1,27) "Mu12",Mu12,"Mu23",Mu23,"Kappa12",Kappa12,"Nu32",Nu32,"E33",E33_Real,"E33_VA",E33

27 FORMAT(A6,2X,E25.15,2X,A6,2X,E25.15,2X,A9,2X,E25.15,2X,A6,2X,E25.15,2X,A6,2X,E25.15,A9,2X,E25.15/)

ENDIF

CLOSE(OUT1)


ENDSUBROUTINE EffectiveProperties

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INCLUDE 'mkl_pardiso.f90'
SUBROUTINE HarwellBoeing_Symmetric_Reduced(k,Dimen1,DUMMYY,DUMMYY2,A1,JA1,JA2,IA,F1Global,x)
USE mkl_pardiso
IMPLICIT NONE
INTEGER, INTENT(IN) :: Dimen1,DUMMYY,DUMMYY2
INTEGER, INTENT(IN) :: k
DOUBLE PRECISION, INTENT(IN), DIMENSION (Dimen1) :: F1Global

INTEGER :: i,j

INTEGER, INTENT(IN), DIMENSION(DUMMYY2) ::IA
INTEGER, INTENT(IN), DIMENSION(DUMMYY) :: JA1
INTEGER, INTENT(IN), DIMENSION(DUMMYY) :: JA2
DOUBLE PRECISION, INTENT(IN), DIMENSION(DUMMYY) :: A1

!INTEGER, PARAMETER :: dp = KIND(1.0D0)
!.. Internal solver memory pointer 
TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
!.. All other variables
INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
INTEGER error1
INTEGER, ALLOCATABLE :: iparm( : )
DOUBLE PRECISION, ALLOCATABLE :: b( : )
DOUBLE PRECISION,  INTENT(OUT), Dimension(Dimen1) :: x
INTEGER :: idum(1)
DOUBLE PRECISION :: ddum(1)
!.. Fill all arrays containing matrix data.
nrhs = 1 
maxfct = 1 
mnum = 1
ALLOCATE( b ( Dimen1 ) )
!..
!.. Set up PARDISO control parameter
!..
ALLOCATE( iparm ( 64 ) )

do i = 1, 64
   iparm(i) = 0
end do 

iparm(1) = 1 ! no solver default
iparm(2) = 2 ! fill-in reordering from METIS
iparm(3) = 1
iparm(4) = 0 ! no iterative-direct algorithm
iparm(5) = 0 ! no user fill-in reducing permutation
iparm(6) = 0 ! =0 solution on the first n compoments of x
iparm(8) = 20 ! numbers of iterative refinement steps
iparm(10) = 16 ! perturbe the pivot elements with 1E-13
iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
iparm(14) = 0 ! Output: number of perturbed pivots
iparm(18) = -1 ! Output: number of nonzeros in the factor LU
iparm(19) = -1 ! Output: Mflops for LU factorization
iparm(20) = 0 ! Output: Numbers of CG Iterations
iparm(52) = 1

error  = 1 ! initialize error flag
msglvl = 1 ! print statistical information
mtype  = -2 ! symmetric, indefinite

!.. Initiliaze the internal solver memory pointer. This is only
! necessary for the FIRST call of the PARDISO solver.

ALLOCATE ( pt ( 64 ) )
do i = 1, 64
   pt( i )%DUMMY =  0 
end do

!.. Reordering and Symbolic Factorization, This step also allocates
! all memory that is necessary for the factorization

phase = 11 ! only reordering and symbolic factorization

CALL pardiso (pt, maxfct, mnum, mtype, phase, DIMEN1, A1, IA, JA1, &
idum, nrhs, iparm, msglvl, ddum, ddum, error)

IF (error /= 0) THEN
   WRITE(*,*) '1 The following ERROR was detected: ', error
!GOTO 1000
ENDIF

!.. Factorization.
phase = 22 ! only factorization
CALL pardiso (pt, maxfct, mnum, mtype, phase, DIMEN1, A1, IA, JA1, &
idum, nrhs, iparm, msglvl, ddum, ddum, error)
!WRITE(*,*) 'Factorization completed ... '

IF (error /= 0) THEN
   WRITE(*,*) '2 The following ERROR was detected: ', error
!GOTO 1000
ENDIF

!.. Back substitution and iterative refinement
!iparm(8) = 2 ! max numbers of iterative refinement steps
phase = 33 ! only factorization
do i = 1, Dimen1
b(i) = F1Global(i)
end do

CALL pardiso (pt, maxfct, mnum, mtype, phase, DIMEN1, A1, IA, JA1, &
idum, nrhs, iparm, msglvl, b, x, error)
IF (error /= 0) THEN
   WRITE(*,*) '3 The following ERROR was detected: ', error
!   GOTO 1000
ENDIF

!WRITE(*,*) Dimen1,DUMMYY,DUMMYY2

WRITE(*,*) 'The solution of the system is '
DO i = 1, Dimen1
   WRITE(*,*) ' x(',i,') = ', x(i)
END DO
      
!1000 CONTINUE
!.. Termination and release of memory
phase = -1 ! release internal memory
CALL pardiso (pt, maxfct, mnum, mtype, phase, DIMEN1, ddum, idum, idum, &
idum, nrhs, iparm, msglvl, ddum, ddum, error1)

IF ( ALLOCATED( iparm ) )   DEALLOCATE( iparm )

!IF (error1 /= 0) THEN
!   WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
!   STOP 1
!ENDIF
!
!IF ( error /= 0 ) STOP 1
!STOP 0
!!END
!     
END SUBROUTINE HarwellBoeing_Symmetric_Reduced


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_Matrix_Reduced(k,E_NumberofDE,E_One_Step_DEs,E_BoundaryIndex,E_Boundaries,E_K1GRR,E_K1GRCDD,E_K1GRCII,E_F1Globall,Connectivities,&
E_Dis1,E_Highest_GRCD)
IMPLICIT NONE
INTEGER :: i,j,m,n,i1,j1,k1,k2,k3
INTEGER, INTENT(IN) :: k

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

INTEGER, INTENT(IN) :: E_NumberofDE
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec) :: E_K1GRR
INTEGER, INTENT(IN), DIMENSION (NumberofNodes+E_NumberofDE) :: E_K1GRCDD
INTEGER, INTENT(IN), DIMENSION (NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec) :: E_K1GRCII
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes+E_NumberofDE) :: E_F1Globall
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: E_K1GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: E_K1GRCD
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: E_K1GRCI
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: E_F1Global

INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: E_K2GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: E_K2GRRD

DOUBLE PRECISION, INTENT(IN), DIMENSION(E_NumberofDE,6) :: E_One_Step_DEs
DOUBLE PRECISION, INTENT(IN), DIMENSION (E_NumberofBoundaryNodes,TotalTimeSteps+1) :: E_Boundaries
INTEGER, INTENT(IN), DIMENSION(E_NumberofBoundaryNodes,1) :: E_BoundaryIndex
DOUBLE PRECISION, DIMENSION (NumberofNodes+E_NumberofDE) :: E_Dis
DOUBLE PRECISION, INTENT(OUT),DIMENSION (NumberofNodes+E_NumberofDE) :: E_Dis1
INTEGER :: E_Dimen1
DOUBLE PRECISION, DIMENSION(NumberofNodes) :: E_x
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
DOUBLE PRECISION :: ElapsedTime,Temp
INTEGER :: E_TEMP_GRCD
INTEGER, INTENT(OUT) :: E_Highest_GRCD
INTEGER :: E_Temp_INT
DOUBLE PRECISION :: E_Temp
INTEGER :: E_DUMMY,E_DUMMY2
INTEGER, ALLOCATABLE, DIMENSION(:) :: E_IA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: E_IA
INTEGER, ALLOCATABLE, DIMENSION(:) :: E_JA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: E_JA2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E_A1

ALLOCATE(E_K1GR(NumberofNodes,E_MaxDimenofInverseConnec))
ALLOCATE(E_K1GRCD(NumberofNodes))
ALLOCATE(E_K1GRCI(NumberofNodes,E_MaxDimenofInverseConnec))
ALLOCATE(E_F1Global(NumberofNodes))

!E_K1GR=E_K1GRR
!E_K1GRCD=E_K1GRCDD
!E_K1GRCI=E_K1GRCII
!E_F1Global=E_F1Globall

DO j=1,E_MaxDimenofInverseConnec
DO i=1,NumberofNodes
E_K1GR(i,j)=E_K1GRR(i,j)
E_K1GRCI(i,j)=E_K1GRCII(i,j)
ENDDO
ENDDO

DO j=1,NumberofNodes
E_K1GRCD(j)=E_K1GRCDD(j)
E_F1Global(j)=E_F1Globall(j)
ENDDO

ALLOCATE(E_K2GR(NumberofNodes,E_MaxDimenofInverseConnec,2))
ALLOCATE(E_K2GRRD(NumberofNodes))
E_K2GR=0
E_K2GRRD=0

DO i=1,NumberofNodes
DO j=1,E_K1GRCD(i)
E_K2GRRD(E_K1GRCI(i,j))=E_K2GRRD(E_K1GRCI(i,j))+1
E_K2GR(E_K1GRCI(i,j),E_K2GRRD(E_K1GRCI(i,j)),1)=i
E_K2GR(E_K1GRCI(i,j),E_K2GRRD(E_K1GRCI(i,j)),2)=j
ENDDO
ENDDO

E_Dis=0.0D0

DO i=1,E_NumberofBoundaryNodes
E_Dis(1*(E_BoundaryIndex(i,1)-1)+1)=E_Boundaries(i,k+1)
ENDDO


DO i=1,E_NumberofBoundaryNodes

DO m=1,E_K2GRRD(1*(E_BoundaryIndex(i,1)-1)+1)
E_F1Global( E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,1) )=E_F1Global( E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,1) )&
-E_K1GR( E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,1), E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,2) )*E_Dis(1*(E_BoundaryIndex(i,1)-1)+1)
ENDDO

E_F1Global(1*(E_BoundaryIndex(i,1)-1)+1)=E_Dis(1*(E_BoundaryIndex(i,1)-1)+1)

DO m=1,E_K2GRRD(1*(E_BoundaryIndex(i,1)-1)+1)
E_K1GR( E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,1), E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,2) )=0.0D0
ENDDO

DO m=1,E_K1GRCD(1*(E_BoundaryIndex(i,1)-1)+1)
IF ( 1*(E_BoundaryIndex(i,1)-1)+1 /= E_K1GRCI(1*(E_BoundaryIndex(i,1)-1)+1,m) ) THEN
E_K1GR( 1*(E_BoundaryIndex(i,1)-1)+1, m )=0.0D0
ELSE
E_K1GR( 1*(E_BoundaryIndex(i,1)-1)+1, m )=1.0D0
ENDIF
ENDDO

ENDDO

DEALLOCATE(E_K2GR)
DEALLOCATE(E_K2GRRD)

E_TEMP_GRCD=E_K1GRCD(1)
E_Highest_GRCD=E_TEMP_GRCD
DO i=2,1*NumberofNodes
IF ( E_K1GRCD(i)>E_TEMP_GRCD ) THEN
E_TEMP_GRCD=E_K1GRCD(i)
E_Highest_GRCD=E_TEMP_GRCD
ENDIF
ENDDO

E_DUMMY=0

DO i=1,1*NumberofNodes
DO j=1,E_K1GRCD(i)
IF (E_K1GR(i,j) /= 0.0D0) THEN
E_DUMMY=E_DUMMY+1
ELSEIF (E_K1GR(i,j) == 0.0D0) THEN
GOTO 3
ENDIF
3 ENDDO
ENDDO

ALLOCATE(E_IA1(E_DUMMY))
ALLOCATE(E_JA1(E_DUMMY))
ALLOCATE(E_JA2(E_DUMMY))
ALLOCATE(E_A1(E_DUMMY))

E_DUMMY=0

DO i=1,1*NumberofNodes
DO j=1,E_K1GRCD(i)
IF (E_K1GR(i,j) /= 0.0D0) THEN
E_DUMMY=E_DUMMY+1
E_A1(E_DUMMY)=E_K1GR(i,j)
E_JA1(E_DUMMY)=E_K1GRCI(i,j)
E_JA2(E_DUMMY)=i
ELSEIF (E_K1GR(i,j) == 0.0D0) THEN
GOTO 4
ENDIF
4 ENDDO
ENDDO

!OPEN(33,FILE="3D_Output_Micro_F/3D_Output_Microscale_1/Check2D_KGR_ME.txt",STATUS="UNKNOWN")
!DO i=1,NumberofNodes+E_NumberofDE
!WRITE(33,34) (E_K1GR(i,k1),k1=1,E_MaxDimenofInverseConnec)
!ENDDO
!34 FORMAT(<40>E15.8,/)
!WRITE(33,*) "E_Highest_GRCD  ",E_Highest_GRCD
!CLOSE(33)

DEALLOCATE(E_K1GR)
DEALLOCATE(E_K1GRCD)
DEALLOCATE(E_K1GRCI)

E_IA1(1)=1

E_DUMMY2=1

DO i=2,E_DUMMY

IF (E_JA2(i)>E_JA2(i-1) .AND. i>1) THEN

E_DUMMY2=E_DUMMY2+1

E_IA1(E_DUMMY2)=i

ENDIF

ENDDO

E_DUMMY2=E_DUMMY2+1

E_IA1(E_DUMMY2)=E_DUMMY+1

ALLOCATE(E_IA(E_DUMMY2))

DO i=1,E_DUMMY2
E_IA(i)=E_IA1(i)
ENDDO

DEALLOCATE(E_IA1)

!OPEN(27,FILE="3D_Output_Micro_F/3D_Output_Microscale_1/A1JA1JA2_ME.txt",STATUS="UNKNOWN")
!DO i=1,E_DUMMY
!WRITE(27,29) E_A1(i),E_JA1(i),E_JA2(i)
!ENDDO
!CLOSE(27)
!
!29 FORMAT(E15.8,2X,I,2X,I,/)
!
!OPEN(39,FILE="3D_Output_Micro_F/3D_Output_Microscale_1/IA_ME.txt",STATUS="UNKNOWN")
!DO i=1,E_DUMMY2
!WRITE(39,33) E_IA(i)
!ENDDO
!CLOSE(39)
!
!33 FORMAT(I,/)

E_Dimen1=1*NumberofNodes

CALL HarwellBoeing_Reduced(k,E_Dimen1,E_DUMMY,E_DUMMY2,E_A1,E_JA1,E_JA2,E_IA,E_F1Global,E_x)

DEALLOCATE(E_F1Global)
DEALLOCATE(E_IA)
DEALLOCATE(E_JA1)
DEALLOCATE(E_JA2)
DEALLOCATE(E_A1)

DO j=1,NumberofNodes
E_Dis1(j)=E_x(1*(j-1)+1)
!WRITE(*,*) E_Dis1(j)
ENDDO

RETURN
END SUBROUTINE E_Matrix_Reduced

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Matrix_Reduced(k,NumberofDE,One_Step_DEs,BoundaryIndex,Boundaries,K1GRR,K1GRCDD,K1GRCII,F1Globall,Connectivities,&
Dis1,Highest_GRCD)
IMPLICIT NONE
INTEGER :: i,j,m,n,i1,j1,k1,k2,k3
INTEGER, INTENT(IN) :: k

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(IN) :: NumberofDE
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec) :: K1GRR
INTEGER, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE) :: K1GRCDD
INTEGER, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec) :: K1GRCII
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE) :: F1Globall
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: K1GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: K1GRCD
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: K1GRCI
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: F1Global

INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: K2GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: K2GRRD

DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofDE,8) :: One_Step_DEs
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofBoundaryNodes,TotalTimeSteps+2) :: Boundaries
INTEGER, INTENT(IN), DIMENSION(NumberofBoundaryNodes,2) :: BoundaryIndex
DOUBLE PRECISION, DIMENSION (3*NumberofNodes+NumberofDE) :: Dis
DOUBLE PRECISION, INTENT(OUT),DIMENSION (3*NumberofNodes+NumberofDE) :: Dis1
INTEGER :: Dimen1
DOUBLE PRECISION, DIMENSION(3*NumberofNodes) :: x
INTEGER :: OUT1=40,OUT2=41
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
!DOUBLE PRECISION :: ElapsedTime,Temp
INTEGER :: TEMP_GRCD
INTEGER, INTENT(OUT) :: Highest_GRCD
INTEGER :: Temp_INT
DOUBLE PRECISION :: Temp
INTEGER :: DUMMY,DUMMY2
INTEGER, ALLOCATABLE, DIMENSION(:) :: IA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: IA
INTEGER, ALLOCATABLE, DIMENSION(:) :: JA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: JA2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: A1

ALLOCATE(K1GR(3*NumberofNodes,MaxDimenofInverseConnec))
ALLOCATE(K1GRCD(3*NumberofNodes))
ALLOCATE(K1GRCI(3*NumberofNodes,MaxDimenofInverseConnec))
ALLOCATE(F1Global(3*NumberofNodes))

!K1GR=K1GRR(1:3*NumberofNodes,MaxDimenofInverseConnec)
!K1GRCD=K1GRCDD(1:3*NumberofNodes)
!K1GRCI=K1GRCII(1:3*NumberofNodes,MaxDimenofInverseConnec)
!F1Global=F1Globall(1:3*NumberofNodes)

DO j=1,MaxDimenofInverseConnec
DO i=1,3*NumberofNodes
K1GR(i,j)=K1GRR(i,j)
K1GRCI(i,j)=K1GRCII(i,j)
ENDDO
ENDDO

DO i=1,3*NumberofNodes
K1GRCD(i)=K1GRCDD(i)
F1Global(i)=F1Globall(i)
ENDDO

ALLOCATE(K2GR(3*NumberofNodes,MaxDimenofInverseConnec,2))
ALLOCATE(K2GRRD(3*NumberofNodes))
K2GR=0
K2GRRD=0

DO i=1,3*NumberofNodes
DO j=1,K1GRCD(i)
K2GRRD(K1GRCI(i,j))=K2GRRD(K1GRCI(i,j))+1
K2GR(K1GRCI(i,j),K2GRRD(K1GRCI(i,j)),1)=i
K2GR(K1GRCI(i,j),K2GRRD(K1GRCI(i,j)),2)=j
ENDDO
ENDDO

Dis=0.0D0

DO i=1,NumberofBoundaryNodes

IF (BoundaryIndex(i,2)==1) THEN

Dis(3*(BoundaryIndex(i,1)-1)+1)=Boundaries(i,k+2)

ELSEIF (BoundaryIndex(i,2)==2) THEN

Dis(3*(BoundaryIndex(i,1)-1)+2)=Boundaries(i,k+2)

ELSEIF (BoundaryIndex(i,2)==3) THEN

Dis(3*(BoundaryIndex(i,1)-1)+3)=Boundaries(i,k+2)

ENDIF

ENDDO


DO i=1,NumberofBoundaryNodes

IF (BoundaryIndex(i,2)==1) THEN

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+1)
F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+1,m,1) )=F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+1,m,1) )&
-K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+1,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+1,m,2) )*Dis(3*(BoundaryIndex(i,1)-1)+1)

ENDDO

F1Global(3*(BoundaryIndex(i,1)-1)+1)=Dis(3*(BoundaryIndex(i,1)-1)+1)

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+1)
K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+1,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+1,m,2) )=0.0D0
ENDDO

DO m=1,K1GRCD(3*(BoundaryIndex(i,1)-1)+1)
IF ( 3*(BoundaryIndex(i,1)-1)+1 /= K1GRCI(3*(BoundaryIndex(i,1)-1)+1,m) ) THEN
K1GR( 3*(BoundaryIndex(i,1)-1)+1, m )=0.0D0
ELSE
K1GR( 3*(BoundaryIndex(i,1)-1)+1, m )=1.0D0
ENDIF
ENDDO


ELSEIF (BoundaryIndex(i,2)==2) THEN

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+2)
F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+2,m,1) )=F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+2,m,1) )&
-K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+2,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+2,m,2) )*Dis(3*(BoundaryIndex(i,1)-1)+2)

ENDDO

F1Global(3*(BoundaryIndex(i,1)-1)+2)=Dis(3*(BoundaryIndex(i,1)-1)+2)

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+2)
K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+2,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+2,m,2) )=0.0D0
ENDDO

DO m=1,K1GRCD(3*(BoundaryIndex(i,1)-1)+2)
IF ( 3*(BoundaryIndex(i,1)-1)+2 /= K1GRCI(3*(BoundaryIndex(i,1)-1)+2,m) ) THEN
K1GR( 3*(BoundaryIndex(i,1)-1)+2, m )=0.0D0
ELSE
K1GR( 3*(BoundaryIndex(i,1)-1)+2, m )=1.0D0
ENDIF
ENDDO

ELSEIF (BoundaryIndex(i,2)==3) THEN

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+3)
F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+3,m,1) )=F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+3,m,1) )&
-K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+3,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+3,m,2) )*Dis(3*(BoundaryIndex(i,1)-1)+3)

ENDDO

F1Global(3*(BoundaryIndex(i,1)-1)+3)=Dis(3*(BoundaryIndex(i,1)-1)+3)

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+3)
K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+3,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+3,m,2) )=0.0D0
ENDDO

DO m=1,K1GRCD(3*(BoundaryIndex(i,1)-1)+3)
IF ( 3*(BoundaryIndex(i,1)-1)+3 /= K1GRCI(3*(BoundaryIndex(i,1)-1)+3,m) ) THEN
K1GR( 3*(BoundaryIndex(i,1)-1)+3, m )=0.0D0
ELSE
K1GR( 3*(BoundaryIndex(i,1)-1)+3, m )=1.0D0
ENDIF
ENDDO

ENDIF

ENDDO

DEALLOCATE(K2GR)
DEALLOCATE(K2GRRD)

TEMP_GRCD=K1GRCD(1)
Highest_GRCD=TEMP_GRCD
DO i=2,3*NumberofNodes
IF ( K1GRCD(i)>TEMP_GRCD ) THEN
TEMP_GRCD=K1GRCD(i)
Highest_GRCD=TEMP_GRCD
ENDIF
ENDDO

!OPEN(33,FILE="3D_Output_Micro_F/Check3.txt",STATUS="UNKNOWN")
!DO i=1,1000  !,MaxDimenofInverseConnec
!WRITE(33,34) (K1GR(i,k1),k1=1,80)
!ENDDO
!34 FORMAT(<80>E15.8,/)
!CLOSE(33)

!CALL CPU_TIME (ElapsedTime)
!WRITE(*,*) "Highest_GRCD  ",Highest_GRCD

DUMMY=0

DO i=1,3*NumberofNodes
DO j=1,K1GRCD(i)
IF (K1GR(i,j) /= 0.0D0) THEN
DUMMY=DUMMY+1
ELSEIF (K1GR(i,j) == 0.0D0) THEN
GOTO 3
ENDIF
3 ENDDO
ENDDO

ALLOCATE(IA1(DUMMY))
ALLOCATE(JA1(DUMMY))
ALLOCATE(JA2(DUMMY))
ALLOCATE(A1(DUMMY))

DUMMY=0

DO i=1,3*NumberofNodes
DO j=1,K1GRCD(i)
IF (K1GR(i,j) /= 0.0D0) THEN
DUMMY=DUMMY+1
A1(DUMMY)=K1GR(i,j)
JA1(DUMMY)=K1GRCI(i,j)
JA2(DUMMY)=i
ELSEIF (K1GR(i,j) == 0.0D0) THEN
GOTO 4
ENDIF
4 ENDDO
ENDDO

DEALLOCATE(K1GR)
DEALLOCATE(K1GRCD)
DEALLOCATE(K1GRCI)

!CALL CPU_TIME (ElapsedTime)
!WRITE(20,*) "DUMMY  ",DUMMY,"  Elapsed Time=",ElapsedTime/60.0D0

IA1(1)=1

DUMMY2=1

DO i=2,DUMMY

IF (JA2(i)>JA2(i-1) .AND. i>1) THEN

DUMMY2=DUMMY2+1

IA1(DUMMY2)=i

ENDIF

ENDDO

DUMMY2=DUMMY2+1

IA1(DUMMY2)=DUMMY+1

ALLOCATE(IA(DUMMY2))

DO i=1,DUMMY2
IA(i)=IA1(i)
ENDDO

DEALLOCATE(IA1)

Dimen1=3*NumberofNodes

CALL HarwellBoeing_Reduced(k,Dimen1,DUMMY,DUMMY2,A1,JA1,JA2,IA,F1Global,x)

!OPEN(27,FILE="3D_Output_Micro_F/A1JA1JA2.txt",STATUS="UNKNOWN")
!DO i=1,6300
!WRITE(27,29) A1(i),JA1(i),JA2(i)
!ENDDO
!CLOSE(27)
!
!29 FORMAT(E15.8,2X,I,2X,I,/)
!
!OPEN(39,FILE="3D_Output_Micro_F/IA.txt",STATUS="UNKNOWN")
!DO i=1,6300
!WRITE(39,33) IA(i)
!ENDDO
!CLOSE(39)
!
!33 FORMAT(I,/)

DEALLOCATE(F1Global)
DEALLOCATE(IA)
DEALLOCATE(JA1)
DEALLOCATE(JA2)
DEALLOCATE(A1)

DO j=1,NumberofNodes
Dis1(j)=x(3*(j-1)+1)
Dis1(NumberofNodes+j)=x(3*(j-1)+2)
Dis1(2*NumberofNodes+j)=x(3*(j-1)+3)
ENDDO

!DO j=1,3*NumberofNodes
!WRITE(*,*) Dis1(j)
!ENDDO

RETURN
ENDSUBROUTINE Matrix_Reduced

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_RotatedConductivities(ii1,ii2,Connectivities,BeginofModeIndex,EndofModeIndex,RotationMatrix2,ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,&
ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff,ME_Rotated_Sigma_Eff)
IMPLICIT NONE
INTEGER :: i1,i2,i3,i4,i5,i6,i7,i8

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

INTEGER, INTENT(IN) :: ii1,ii2
INTEGER, INTENT(IN) :: BeginofModeIndex,EndofModeIndex
INTEGER, INTENT(IN), DIMENSION(NumberofElements,17) :: Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofMaterialModes,3,3) :: RotationMatrix2
DOUBLE PRECISION, INTENT(IN) :: ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff
DOUBLE PRECISION, DIMENSION(3,3) :: ME_Rotated_Sigma_Eff_Temp
DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,3) :: ME_Rotated_Sigma_Eff

ME_Rotated_Sigma_Eff_Temp(1,1)=ME_Sigma_ZZ_Eff
ME_Rotated_Sigma_Eff_Temp(2,2)=ME_Sigma_XX_Eff
ME_Rotated_Sigma_Eff_Temp(3,3)=ME_Sigma_YY_Eff
ME_Rotated_Sigma_Eff_Temp(1,2)=ME_Sigma_XZ_Eff
ME_Rotated_Sigma_Eff_Temp(1,3)=ME_Sigma_YZ_Eff
ME_Rotated_Sigma_Eff_Temp(2,3)=ME_Sigma_XY_Eff
ME_Rotated_Sigma_Eff_Temp(2,1)=ME_Rotated_Sigma_Eff_Temp(1,2)
ME_Rotated_Sigma_Eff_Temp(3,1)=ME_Rotated_Sigma_Eff_Temp(1,3)
ME_Rotated_Sigma_Eff_Temp(3,2)=ME_Rotated_Sigma_Eff_Temp(2,3)

ME_Rotated_Sigma_Eff=0.0D0

IF (Connectivities(ii1,ii2+9)>=BeginofModeIndex .AND. Connectivities(ii1,ii2+9)<=EndofModeIndex) THEN
DO i1=1,3
DO i2=1,3
DO i3=1,3
DO i4=1,3
ME_Rotated_Sigma_Eff(i1,i2)=ME_Rotated_Sigma_Eff(i1,i2)+RotationMatrix2(Connectivities(ii1,ii2+9)-BeginofModeIndex+1,i1,i3)*&
RotationMatrix2(Connectivities(ii1,ii2+9)-BeginofModeIndex+1,i2,i4)*ME_Rotated_Sigma_Eff_Temp(i3,i4)
ENDDO
ENDDO
ENDDO
ENDDO
ENDIF

END SUBROUTINE ME_RotatedConductivities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_EffectiveProperties(ZeroTimeStep,k,ii,ii3,ii1,ii2,iiii,Ex_Knot,Ey_Knot,Ez_Knot,M_Nodes,ME_Vol,ME_Sub_EXX,ME_Sub_EYY,ME_Sub_EZZ,&
ME_Sub_JXX,ME_Sub_JYY,ME_Sub_JZZ,M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes,&
ME_sub_sigmasigma_updated,ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff,Effec_Electro_Prop_Microscale_E,&
StrainsZZ)
IMPLICIT NONE

INTEGER :: i,j,i1,i2,i3,i4,k1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

INTEGER, INTENT(IN) :: ZeroTimeStep,k,ii,ii3,ii1,ii2,iiii
INTEGER :: ME_OUT7
DOUBLE PRECISION, INTENT(IN) ::  Ex_Knot,Ey_Knot,Ez_Knot
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements) :: ME_Vol
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_EXX
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_EYY
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_EZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_JXX
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_JYY
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_JZZ

INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveXNodes) :: M_PositiveXNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveYNodes) :: M_PositiveYNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveZNodes) :: M_PositiveZNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeXNodes) :: M_NegativeXNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeYNodes) :: M_NegativeYNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeZNodes) :: M_NegativeZNodes

DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofNodes,3) :: M_Nodes
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_sigmasigma_updated
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_sigmasigma_updated2
DOUBLE PRECISION :: M_Length_X,M_Length_Y,M_Length_Z,ME_TotalEnergy
DOUBLE PRECISION, INTENT(INOUT) :: ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff
CHARACTER*90, INTENT(IN) :: Effec_Electro_Prop_Microscale_E
DOUBLE PRECISION :: StrainsZZ

ME_OUT7=51+100*(ii-1)

ME_sub_sigmasigma_updated2=ME_sub_sigmasigma_updated

!DO i=1,ME_NumberofGaussPoints
!DO k1=1,M_NumberofElements
!DO i1=1,3
!IF ( ME_sub_sigmasigma_updated2(i,k1,i1,i1) == 1.0D-8) THEN
!ME_sub_sigmasigma_updated2(i,k1,i1,i1)=1.0D5
!ENDIF
!ENDDO
!ENDDO
!ENDDO

M_Length_X=M_Nodes(M_PositiveXNodes(1),2)-M_Nodes(M_NegativeXNodes(1),2)
M_Length_Y=M_Nodes(M_PositiveYNodes(1),3)-M_Nodes(M_NegativeYNodes(1),3)
M_Length_Z=M_Nodes(M_PositiveZNodes(1),4)-M_Nodes(M_NegativeZNodes(1),4)
ME_TotalEnergy=0.0D0
DO i=1,ME_NumberofGaussPoints
DO j=1,M_NumberofElements
ME_TotalEnergy=ME_TotalEnergy+0.5D0*( ME_Vol(i,j)*ME_Sub_EXX(j,i)*ME_Sub_JXX(j,i)+ME_Vol(i,j)*ME_Sub_EYY(j,i)*ME_Sub_JYY(j,i)+&
ME_Vol(i,j)*ME_Sub_EZZ(j,i)*ME_Sub_JZZ(j,i) )
ENDDO
ENDDO

IF (iiii==1) THEN ! For sigma_xx and sigma_zz (the rule of mixtures)
ME_Sigma_XX_Eff=ME_TotalEnergy/M_Length_X/M_Length_Y/M_Length_Z/0.5D0/Ex_Knot/Ex_Knot
!ME_Sigma_XX_Eff=1.0D-1
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!IF (ME_Sigma_XX_Eff <= 1.0D-1) THEN
!ME_Sigma_XX_Eff=1.0D-6
!ENDIF
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ME_Sigma_ZZ_Eff=0.0D0
DO i=1,ME_NumberofGaussPoints
DO j=1,M_NumberofElements
ME_Sigma_ZZ_Eff=ME_Sigma_ZZ_Eff+ME_sub_sigmasigma_updated2(i,j,3,3)*ME_Vol(i,j)/M_Length_X/M_Length_Y/M_Length_Z
ENDDO
ENDDO
ELSEIF (iiii==2) THEN ! For sigma_yy
ME_Sigma_YY_Eff=ME_TotalEnergy/M_Length_X/M_Length_Y/M_Length_Z/0.5D0/Ey_Knot/Ey_Knot
!ME_Sigma_YY_Eff=1.0D-1
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!IF (ME_Sigma_YY_Eff <= 1.0D-1) THEN
!ME_Sigma_YY_Eff=1.0D-6
!ENDIF
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ELSEIF (iiii==3) THEN ! For sigma_xy
ME_Sigma_XY_Eff=(ME_TotalEnergy/M_Length_X/M_Length_Y/M_Length_Z-0.5D0*Ex_Knot*Ex_Knot*ME_Sigma_XX_Eff-0.5D0*Ey_Knot*Ey_Knot*&
ME_Sigma_YY_Eff)/Ex_Knot/Ey_Knot
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ME_Sigma_XY_Eff=0.0D0
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ELSEIF (iiii==4) THEN ! For sigma_xz
ME_Sigma_XZ_Eff=(ME_TotalEnergy/M_Length_X/M_Length_Y/M_Length_Z-0.5D0*Ex_Knot*Ex_Knot*ME_Sigma_XX_Eff-0.5D0*Ez_Knot*Ez_Knot*&
ME_Sigma_ZZ_Eff)/Ex_Knot/Ez_Knot
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ME_Sigma_XZ_Eff=0.0D0
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ELSEIF (iiii==5) THEN ! For sigma_yz
ME_Sigma_YZ_Eff=(ME_TotalEnergy/M_Length_X/M_Length_Y/M_Length_Z-0.5D0*Ey_Knot*Ey_Knot*ME_Sigma_YY_Eff-0.5D0*Ez_Knot*Ez_Knot*&
ME_Sigma_ZZ_Eff)/Ey_Knot/Ez_Knot
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ME_Sigma_YZ_Eff=0.0D0
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!IF(ZeroTimeStep==1 .AND. k==1) THEN
!ME_Sigma_ZZ_Eff=1.0/1.0D7  !1.0D0/2.0D14 !1.0D0/(3.95148D14)
!ELSE
!ME_Sigma_ZZ_Eff=1.0D0/2.0D14 !1.0D0/(-2.90687D18*StrainsZZ**3+2.64456D17*StrainsZZ**2-1.07798D16*StrainsZZ+3.95148D14)  !Temporily changed
!ENDIF
ENDIF

IF (ii1==ElementShown .AND. ii2==SubElementShown .AND. ZeroTimeStep==1 .AND. k==1 .AND. iiii==5) THEN
OPEN(ME_OUT7,FILE=Effec_Electro_Prop_Microscale_E,STATUS="UNKNOWN")
WRITE(ME_OUT7,20) "Element_Index","Sub_Element_Index","ME_Sigma_XX_Eff","ME_Sigma_YY_Eff","ME_Sigma_ZZ_Eff","ME_Sigma_XY_Eff","ME_Sigma_XZ_Eff","ME_Sigma_YZ_Eff"
WRITE(ME_OUT7,10) ii1,ii2,ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff
CLOSE(ME_OUT7)
ELSEIF (ii1==ElementShown .AND. ii2==SubElementShown .AND. ZeroTimeStep==0 .AND. k==1 .AND. iiii==5) THEN
OPEN(ME_OUT7,FILE=Effec_Electro_Prop_Microscale_E,STATUS="UNKNOWN",POSITION="APPEND")
WRITE(ME_OUT7,10) ii1,ii2,ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff
CLOSE(ME_OUT7)
ELSEIF (ii1==ElementShown .AND. ii2==SubElementShown .AND. k/=1 .AND. iiii==5) THEN
OPEN(ME_OUT7,FILE=Effec_Electro_Prop_Microscale_E,STATUS="UNKNOWN",POSITION="APPEND")
WRITE(ME_OUT7,10) ii1,ii2,ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff
CLOSE(ME_OUT7)
ENDIF

20 FORMAT(A25,2X,A25,2X,A25,2X,A25,2X,A25,2X,A25,2X,A25,2X,A25)
10 FORMAT(I25,2X,I25,2X,E25.15,2X,E25.15,2X,E25.15,2X,E25.15,2X,E25.15,2X,E25.15)

END SUBROUTINE ME_EffectiveProperties

SUBROUTINE E_EffectiveProperties(E_ZeroTimeStep,k,ii,MacroEx_Knot,MacroEy_Knot,MacroEz_Knot,Nodes,E_Vol,E_Sub_EXX,E_Sub_EYY,E_Sub_EZZ,E_Sub_JXX,E_Sub_JYY,&
E_Sub_JZZ,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes,E_sub_sigmasigma,&
E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ2_Eff)
IMPLICIT NONE

INTEGER :: i,j,i1,i2,i3,i4
INTEGER :: E_OUT3=18 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

INTEGER, INTENT(IN) :: E_ZeroTimeStep,k,ii
DOUBLE PRECISION, INTENT(IN) ::  MacroEx_Knot,MacroEy_Knot,MacroEz_Knot
DOUBLE PRECISION, INTENT(IN), DIMENSION(E_NumberofGaussPoints,NumberofElements) :: E_Vol
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,E_NumberofGaussPoints) :: E_Sub_EXX
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,E_NumberofGaussPoints) :: E_Sub_EYY
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,E_NumberofGaussPoints) :: E_Sub_EZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,E_NumberofGaussPoints) :: E_Sub_JXX
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,E_NumberofGaussPoints) :: E_Sub_JYY
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,E_NumberofGaussPoints) :: E_Sub_JZZ

INTEGER, INTENT(IN), DIMENSION(NumberofPositiveXNodes) :: PositiveXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveYNodes) :: PositiveYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveZNodes) :: PositiveZNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeXNodes) :: NegativeXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeYNodes) :: NegativeYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeZNodes) :: NegativeZNodes

DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofNodes,3) :: Nodes
DOUBLE PRECISION, INTENT(IN), DIMENSION(E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_sigmasigma
DOUBLE PRECISION :: Length_X,Length_Y,Length_Z,E_TotalEnergy
DOUBLE PRECISION, INTENT(INOUT) :: E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ2_Eff

Length_X=Nodes(PositiveXNodes(1),2)-Nodes(NegativeXNodes(1),2)
Length_Y=Nodes(PositiveYNodes(1),3)-Nodes(NegativeYNodes(1),3)
Length_Z=Nodes(PositiveZNodes(1),4)-Nodes(NegativeZNodes(1),4)
E_TotalEnergy=0.0D0
DO i=1,E_NumberofGaussPoints
DO j=1,NumberofElements
E_TotalEnergy=E_TotalEnergy+0.5D0*( E_Vol(i,j)*E_Sub_EXX(j,i)*E_Sub_JXX(j,i)+E_Vol(i,j)*E_Sub_EYY(j,i)*E_Sub_JYY(j,i)+&
E_Vol(i,j)*E_Sub_EZZ(j,i)*E_Sub_JZZ(j,i) )
ENDDO
ENDDO

IF (ii==1) THEN ! For sigma_xx and sigma_zz (the rule of mixtures)
E_Sigma_XX_Eff=E_TotalEnergy/Length_X/Length_Y/Length_Z/0.5D0/MacroEx_Knot/MacroEx_Knot
E_Sigma_ZZ_Eff=0.0D0
DO i=1,E_NumberofGaussPoints
DO j=1,NumberofElements
E_Sigma_ZZ_Eff=E_Sigma_ZZ_Eff+E_sub_sigmasigma(i,j,3,3)*E_Vol(i,j)/Length_X/Length_Y/Length_Z
ENDDO
ENDDO
ELSEIF (ii==3) THEN ! For sigma_yy
E_Sigma_YY_Eff=E_TotalEnergy/Length_X/Length_Y/Length_Z/0.5D0/MacroEy_Knot/MacroEy_Knot
ELSEIF (ii==4) THEN ! For sigma_xy
E_Sigma_XY_Eff=(E_TotalEnergy/Length_X/Length_Y/Length_Z-0.5D0*MacroEx_Knot*MacroEx_Knot*E_Sigma_XX_Eff-0.5D0*MacroEy_Knot*MacroEy_Knot*&
E_Sigma_YY_Eff)/MacroEx_Knot/MacroEy_Knot
ELSEIF (ii==5) THEN ! For sigma_xz
E_Sigma_XZ_Eff=(E_TotalEnergy/Length_X/Length_Y/Length_Z-0.5D0*MacroEx_Knot*MacroEx_Knot*E_Sigma_XX_Eff-0.5D0*MacroEz_Knot*MacroEz_Knot*&
E_Sigma_ZZ2_Eff)/MacroEx_Knot/MacroEz_Knot
ELSEIF (ii==6) THEN ! For sigma_yz
E_Sigma_YZ_Eff=(E_TotalEnergy/Length_X/Length_Y/Length_Z-0.5D0*MacroEy_Knot*MacroEy_Knot*E_Sigma_YY_Eff-0.5D0*MacroEz_Knot*MacroEz_Knot*&
E_Sigma_ZZ2_Eff)/MacroEy_Knot/MacroEz_Knot
ELSEIF (ii==2) THEN ! For sigma_zz
E_Sigma_ZZ2_Eff=E_TotalEnergy/Length_X/Length_Y/Length_Z/0.5D0/MacroEz_Knot/MacroEz_Knot
ENDIF

IF (E_ZeroTimeStep==1 .AND. k==1 .AND. ii==6) THEN
OPEN(E_OUT3,FILE="3D_Output_Micro_F/3D_Output/Effec_Electro_Prop_E.txt",STATUS="UNKNOWN")
WRITE(E_OUT3,20) "Time_Step","E_Sigma_XX_Eff","E_Sigma_YY_Eff","E_Sigma_ZZ2_Eff","E_Sigma_XY_Eff","E_Sigma_XZ_Eff","E_Sigma_YZ_Eff","E_Sigma_ZZ_Eff"
WRITE(E_OUT3,10) k-1,E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ2_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ_Eff
CLOSE(E_OUT3)
ELSEIF (E_ZeroTimeStep==0 .AND. k==1 .AND. ii==6) THEN
OPEN(E_OUT3,FILE="3D_Output_Micro_F/3D_Output/Effec_Electro_Prop_E.txt",STATUS="UNKNOWN",POSITION="APPEND")
WRITE(E_OUT3,10) k,E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ2_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ_Eff
CLOSE(E_OUT3)
ELSEIF (k/=1 .AND. ii==6) THEN
OPEN(E_OUT3,FILE="3D_Output_Micro_F/3D_Output/Effec_Electro_Prop_E.txt",STATUS="UNKNOWN",POSITION="APPEND")
WRITE(E_OUT3,10) k,E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ2_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ_Eff
CLOSE(E_OUT3)
ENDIF

20 FORMAT(A25,2X,A25,2X,A25,2X,A25,2X,A25,2X,A25,2X,A25,2X,A25)
10 FORMAT(I,2X,E25.15,2X,E25.15,2X,E25.15,2X,E25.15,2X,E25.15,2X,E25.15,2X,E25.15)

END SUBROUTINE E_EffectiveProperties

SUBROUTINE E_Anisotropic(Nodes,E_Connectivities,E_BeginofModeIndex,E_EndofModeIndex,E_XYofAxis,E_MaximumRadius,E_sub_sigmasigma,E_NN)
IMPLICIT NONE
INTEGER :: i,j,i1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: Nodes
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: E_Connectivities
INTEGER, INTENT(IN) :: E_BeginofModeIndex,E_EndofModeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION(11,2) :: E_XYofAxis
DOUBLE PRECISION, INTENT(IN) :: E_MaximumRadius
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_sigmasigma
DOUBLE PRECISION, DIMENSION(E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_Ori_sigmasigma

DOUBLE PRECISION, DIMENSION(E_EndofModeIndex-E_BeginofModeIndex+1,3) :: E_Centroid

DOUBLE PRECISION, INTENT(IN), DIMENSION(8,E_NumberofGaussPoints) :: E_NN

DOUBLE PRECISION, DIMENSION(E_EndofModeIndex-E_BeginofModeIndex+1) :: E_COSTHETA
DOUBLE PRECISION, DIMENSION(E_EndofModeIndex-E_BeginofModeIndex+1) :: E_SINTHETA

DOUBLE PRECISION, DIMENSION(7) :: E_r

DO i=1,NumberofElements

DO j=1,E_NumberofGaussPoints

IF ( E_Connectivities(i,j+9) >= E_BeginofModeIndex .AND. E_Connectivities(i,j+9) <= E_EndofModeIndex) THEN

E_Centroid(E_Connectivities(i,j+9)-E_BeginofModeIndex+1,1)=Nodes(E_Connectivities(i,2),2)*E_NN(1,j)+Nodes(E_Connectivities(i,3),2)*E_NN(2,j)+&
Nodes(E_Connectivities(i,4),2)*E_NN(3,j)+Nodes(E_Connectivities(i,5),2)*E_NN(4,j)+Nodes(E_Connectivities(i,6),2)*E_NN(5,j)+&
Nodes(E_Connectivities(i,7),2)*E_NN(6,j)+Nodes(E_Connectivities(i,8),2)*E_NN(7,j)+Nodes(E_Connectivities(i,9),2)*E_NN(8,j)

E_Centroid(E_Connectivities(i,j+9)-E_BeginofModeIndex+1,2)=Nodes(E_Connectivities(i,2),3)*E_NN(1,j)+Nodes(E_Connectivities(i,3),3)*E_NN(2,j)+&
Nodes(E_Connectivities(i,4),3)*E_NN(3,j)+Nodes(E_Connectivities(i,5),3)*E_NN(4,j)+Nodes(E_Connectivities(i,6),3)*E_NN(5,j)+&
Nodes(E_Connectivities(i,7),3)*E_NN(6,j)+Nodes(E_Connectivities(i,8),3)*E_NN(7,j)+Nodes(E_Connectivities(i,9),3)*E_NN(8,j)

E_Centroid(E_Connectivities(i,j+9)-E_BeginofModeIndex+1,3)=Nodes(E_Connectivities(i,2),4)*E_NN(1,j)+Nodes(E_Connectivities(i,3),4)*E_NN(2,j)+&
Nodes(E_Connectivities(i,4),4)*E_NN(3,j)+Nodes(E_Connectivities(i,5),4)*E_NN(4,j)+Nodes(E_Connectivities(i,6),4)*E_NN(5,j)+&
Nodes(E_Connectivities(i,7),4)*E_NN(6,j)+Nodes(E_Connectivities(i,8),4)*E_NN(7,j)+Nodes(E_Connectivities(i,9),4)*E_NN(8,j)

ENDIF

ENDDO

ENDDO

DO i=1,E_EndofModeIndex-E_BeginofModeIndex+1

E_r=0.0D0

!DO i1=1,E_NumberofGaussPoints

DO j=1,7

E_r(j)=SQRT((E_Centroid(i,1)-E_XYofAxis(j,1))**2+(E_Centroid(i,2)-E_XYofAxis(j,2))**2)

IF (E_r(j) < E_MaximumRadius) THEN

E_COSTHETA(i)= (E_Centroid(i,1)-E_XYofAxis(j,1)) /SQRT( (E_Centroid(i,1)-E_XYofAxis(j,1))**2+ (E_Centroid(i,2)-E_XYofAxis(j,2))**2)
E_SINTHETA(i)= (E_Centroid(i,2)-E_XYofAxis(j,2)) /SQRT( (E_Centroid(i,1)-E_XYofAxis(j,1))**2+ (E_Centroid(i,2)-E_XYofAxis(j,2))**2)

GOTO 1

ELSE

ENDIF

1 ENDDO

!ENDDO

ENDDO


E_sub_Ori_sigmasigma=E_sub_sigmasigma


DO i1=1,E_NumberofGaussPoints

DO i=1,NumberofElements

IF ( E_Connectivities(i,i1+9) >= E_BeginofModeIndex .AND. E_Connectivities(i,i1+9) <= E_EndofModeIndex) THEN

E_sub_sigmasigma(i1,i,1,1)=E_sub_Ori_sigmasigma(i1,i,1,1)*E_COSTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)**2-2.0D0*E_sub_Ori_sigmasigma(i1,i,1,2)*&
E_COSTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)*E_SINTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)+&
E_sub_Ori_sigmasigma(i1,i,2,2)*E_SINTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)**2

E_sub_sigmasigma(i1,i,1,2)=E_sub_Ori_sigmasigma(i1,i,1,2)*E_COSTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)**2+E_sub_Ori_sigmasigma(i1,i,1,1)*&
E_COSTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)*E_SINTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)-&
E_sub_Ori_sigmasigma(i1,i,2,2)*E_COSTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)*&
E_SINTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)-E_sub_Ori_sigmasigma(i1,i,1,2)*E_SINTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)**2

E_sub_sigmasigma(i1,i,1,3)=E_sub_Ori_sigmasigma(i1,i,1,3)*E_COSTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)-&
E_sub_Ori_sigmasigma(i1,i,2,3)*E_SINTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)

E_sub_sigmasigma(i1,i,2,1)=E_sub_sigmasigma(i1,i,1,2)

E_sub_sigmasigma(i1,i,2,2)=E_sub_Ori_sigmasigma(i1,i,2,2)*E_COSTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)**2+2.0D0*E_sub_Ori_sigmasigma(i1,i,1,2)*&
E_COSTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)*E_SINTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)+&
E_sub_Ori_sigmasigma(i1,i,1,1)*E_SINTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)**2

E_sub_sigmasigma(i1,i,2,3)=E_sub_Ori_sigmasigma(i1,i,2,3)*E_COSTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)+&
E_sub_Ori_sigmasigma(i1,i,1,3)*E_SINTHETA(E_Connectivities(i,i1+9)-E_BeginofModeIndex+1)

E_sub_sigmasigma(i1,i,3,1)=E_sub_sigmasigma(i1,i,1,3)

E_sub_sigmasigma(i1,i,3,2)=E_sub_sigmasigma(i1,i,2,3)

E_sub_sigmasigma(i1,i,3,3)=E_sub_Ori_sigmasigma(i1,i,3,3)

ENDIF

ENDDO

ENDDO

ENDSUBROUTINE E_Anisotropic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_ReadElectroPropConnec (E_sub_sigmasigma_in,E_Connectivities)
IMPLICIT NONE
INTEGER :: E_IN3=18
INTEGER :: E_IN4=20
INTEGER :: i,j,i1
INTEGER :: E_ModeIndex

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

INTEGER, INTENT(OUT), DIMENSION(NumberofElements,17) :: E_Connectivities 
DOUBLE PRECISION, DIMENSION(E_NumberofMaterialModes) :: E_sigma11,E_sigma12,E_sigma13,E_sigma22,E_sigma23,E_sigma33
DOUBLE PRECISION, INTENT(OUT), DIMENSION(TotalTimeSteps+1,E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_sigmasigma_in

DOUBLE PRECISION, DIMENSION(E_NumberofMaterialModes) :: E_rho11,E_rho12,E_rho13,E_rho22,E_rho23,E_rho33
DOUBLE PRECISION, DIMENSION(E_NumberofMaterialModes,3,3) :: E_rho,E_sigma

E_rho=0.0D0
E_sigma=0.0D0

E_sub_sigmasigma_in=0.0D0

OPEN(E_IN3,file = "3D_Output_Micro_F/3D_Input/Electrical_material_properties_E.txt",status ='unknown')
DO i=1,E_NumberofMaterialModes
READ(E_IN3,*) E_ModeIndex,E_rho11(i),E_rho12(i),E_rho13(i),E_rho22(i),E_rho23(i),E_rho33(i)
E_rho(i,1,1)=E_rho11(i)
E_rho(i,1,2)=E_rho12(i)
E_rho(i,1,3)=E_rho13(i)
E_rho(i,2,2)=E_rho22(i)
E_rho(i,2,3)=E_rho23(i)
E_rho(i,3,3)=E_rho33(i)
E_rho(i,2,1)=E_rho(i,1,2)
E_rho(i,3,1)=E_rho(i,1,3)
E_rho(i,3,2)=E_rho(i,2,3)
CALL INV(i,3,E_rho(i,:,:),E_sigma(i,:,:))
E_sigma11(i)=E_sigma(i,1,1)
E_sigma12(i)=E_sigma(i,1,2)
E_sigma13(i)=E_sigma(i,1,3)
E_sigma22(i)=E_sigma(i,2,2)
E_sigma23(i)=E_sigma(i,2,3)
E_sigma33(i)=E_sigma(i,3,3)
ENDDO
CLOSE(E_IN3)

OPEN(E_IN4,file = "3D_Output_Micro_F/3D_Input/Connectivities_E_New.txt",status ='unknown')
READ (E_IN4,*) ((E_Connectivities(i,j),j=1,17),i=1,NumberofElements)
CLOSE(E_IN4)

DO j=1,E_NumberofGaussPoints
DO i=1,NumberofElements
E_sub_sigmasigma_in(:,j,i,1,1)=E_sigma11(E_Connectivities(i,j+9))
E_sub_sigmasigma_in(:,j,i,1,2)=E_sigma12(E_Connectivities(i,j+9))
E_sub_sigmasigma_in(:,j,i,1,3)=E_sigma13(E_Connectivities(i,j+9))
E_sub_sigmasigma_in(:,j,i,2,2)=E_sigma22(E_Connectivities(i,j+9))
E_sub_sigmasigma_in(:,j,i,2,3)=E_sigma23(E_Connectivities(i,j+9))
E_sub_sigmasigma_in(:,j,i,3,3)=E_sigma33(E_Connectivities(i,j+9))
E_sub_sigmasigma_in(:,j,i,2,1)=E_sub_sigmasigma_in(:,j,i,1,2)
E_sub_sigmasigma_in(:,j,i,3,1)=E_sub_sigmasigma_in(:,j,i,1,3)
E_sub_sigmasigma_in(:,j,i,3,2)=E_sub_sigmasigma_in(:,j,i,2,3)
ENDDO
ENDDO

RETURN

END SUBROUTINE E_ReadElectroPropConnec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_Tunneling_Effect(ZeroTimeStep,k,Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R,Sub_StrainsXZ1_R,&
Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,Sub_StrainsYZ2_R,M_NumberofDE,M_Nodes,ME_Connectivities,ME_sub_rhorho,M_NumberofElementsPerCNT,ME_ElementsPerCNT,&
ME_TunnelingIndex,ME_TunnelingIndexDimension,ME_Width1,ME_Width2,ME_CNTIndex,ME_Critical_Distance,ME_Coordinates2,M_Dis1,M_CurrentNodes,M_NN,M_Vol,&
ME_sub_Tunneling_rhorho)
IMPLICIT NONE
INTEGER :: i,j,k1,i1,i2,i3,i4

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec
!The electrical tunneling control flags in the microscale
INTEGER :: ME_Tunneling
DOUBLE PRECISION :: ME_Vf
DOUBLE PRECISION :: ME_RofCNT
DOUBLE PRECISION :: ME_Lambda
DOUBLE PRECISION :: ME_rho
INTEGER :: ME_BeginningPolymer
INTEGER :: ME_EndingPolymer
INTEGER :: ME_NumberofCNTs
INTEGER :: ME_NumberofCNTs2
INTEGER :: ME_MaxNumberofElementsPerCNT
COMMON/Micro_Electrical_More_Inputs/ME_Tunneling,ME_BeginningPolymer,ME_EndingPolymer,ME_NumberofCNTs,ME_NumberofCNTs2,&
ME_MaxNumberofElementsPerCNT,ME_Vf,ME_RofCNT,ME_Lambda,ME_rho

INTEGER, INTENT(IN) :: ZeroTimeStep,k,M_NumberofDE
DOUBLE PRECISION, INTENT(IN) :: Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R,Sub_StrainsXZ1_R,&
Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,Sub_StrainsYZ2_R    !Changed
DOUBLE PRECISION :: Sub_StrainsXX_R2,Sub_StrainsYY_R2,Sub_StrainsZZ_R2,Sub_StrainsXY1_R2,Sub_StrainsXY2_R2,Sub_StrainsXZ1_R2,&
Sub_StrainsXZ2_R2,Sub_StrainsYZ1_R2,Sub_StrainsYZ2_R2 !Changed
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofNodes,4) :: M_Nodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofElements,10) :: ME_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_rhorho !Changed
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,4,3) :: ME_TunnelingIndex !Changed
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: ME_TunnelingIndex3  !Changed
INTEGER, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements) :: ME_TunnelingIndexDimension
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs2) :: ME_CNTIndex
DOUBLE PRECISION, INTENT(IN) :: ME_Critical_Distance
DOUBLE PRECISION, INTENT(IN), DIMENSION(9,ME_NumberofCNTs2,2) :: ME_Coordinates2
DOUBLE PRECISION, INTENT(IN) :: ME_Width1,ME_Width2
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs) :: M_NumberofElementsPerCNT
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs,ME_MaxNumberofElementsPerCNT) :: ME_ElementsPerCNT
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofGaussPoints,M_NumberofElements) :: M_Vol !Changed
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofGaussPoints,8) :: M_NN !Changed

DOUBLE PRECISION, INTENT(IN), DIMENSION(3*M_NumberofNodes+M_NumberofDE) :: M_Dis1    !Changed
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofNodes,4) :: M_CurrentNodes         !Changed
DOUBLE PRECISION, DIMENSION(ME_NumberofCNTs,3) :: ME_DisPerCNT                        !Changed
DOUBLE PRECISION, DIMENSION(9,ME_NumberofCNTs2,2) :: ME_Current_LocationPerCNT
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,2) :: ME_Original_Location !Changed
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,2) :: ME_Current_Location  !Changed
DOUBLE PRECISION, INTENT(OUT), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_Tunneling_rhorho !Changed
DOUBLE PRECISION :: Vol_Temp,Temp1,Temp2,Temp3,Temp4,Temp5,Temp6,Resistivity_Temp,Temp_X,Temp_Y,Temp_Distance
DOUBLE PRECISION :: ME_Width11,ME_Width12,ME_Width21,ME_Width22
DOUBLE PRECISION :: Average_Distance

IF (k==1 .AND. ZeroTimeStep==1) THEN
Sub_StrainsXX_R2=0.0D0
Sub_StrainsYY_R2=0.0D0
Sub_StrainsZZ_R2=0.0D0
Sub_StrainsXY1_R2=0.0D0
Sub_StrainsXY2_R2=0.0D0
Sub_StrainsXZ1_R2=0.0D0
Sub_StrainsXZ2_R2=0.0D0
Sub_StrainsYZ1_R2=0.0D0
Sub_StrainsYZ2_R2=0.0D0
ELSE

!StrainsZZ=Sub_StrainsXX_R(ii1,ii2)
!StrainsXX=Sub_StrainsYY_R(ii1,ii2)
!StrainsYY=Sub_StrainsZZ_R(ii1,ii2)
!StrainsXZ2=Sub_StrainsXY1_R(ii1,ii2)
!StrainsXZ1=Sub_StrainsXY2_R(ii1,ii2)
!StrainsYZ2=Sub_StrainsXZ1_R(ii1,ii2)
!StrainsYZ1=Sub_StrainsXZ2_R(ii1,ii2)
!StrainsXY1=Sub_StrainsYZ1_R(ii1,ii2)
!StrainsXY2=Sub_StrainsYZ2_R(ii1,ii2)

Sub_StrainsZZ_R2=Sub_StrainsXX_R
Sub_StrainsXX_R2=Sub_StrainsYY_R
Sub_StrainsYY_R2=Sub_StrainsZZ_R
Sub_StrainsXZ2_R2=Sub_StrainsXY1_R
Sub_StrainsXZ1_R2=Sub_StrainsXY2_R
Sub_StrainsYZ2_R2=Sub_StrainsXZ1_R
Sub_StrainsYZ1_R2=Sub_StrainsXZ2_R
Sub_StrainsXY1_R2=Sub_StrainsYZ1_R
Sub_StrainsXY2_R2=Sub_StrainsYZ2_R
ENDIF

ME_Width11=Sub_StrainsXX_R2*ME_Width1
ME_Width12=2.0D0*Sub_StrainsXY2_R2*ME_Width1
ME_Width21=2.0D0*Sub_StrainsXY1_R2*ME_Width2
ME_Width22=Sub_StrainsYY_R2*ME_Width2

!ME_Width11=0.0D0 !Sub_StrainsYY_R2*ME_Width1
!ME_Width12=2.0D0*Sub_StrainsYZ2_R2*ME_Width1
!ME_Width21=2.0D0*Sub_StrainsYZ1_R2*ME_Width2
!ME_Width22=Sub_StrainsZZ_R2*ME_Width2

!WRITE(*,*) "11,12,21,22",Sub_StrainsYY_R2,Sub_StrainsYZ2_R2,Sub_StrainsYZ1_R2,Sub_StrainsZZ_R2

ALLOCATE(ME_TunnelingIndex3(ME_NumberofGaussPoints,M_NumberofElements,4,3))
ME_TunnelingIndex3=ME_TunnelingIndex

ME_Current_LocationPerCNT=0.0D0
ME_DisPerCNT=0.0D0
ME_sub_Tunneling_rhorho=0.0D0



DO i=1,ME_NumberofCNTs

Vol_Temp=0.0D0

DO i1=1,ME_NumberofGaussPoints

DO j=1,M_NumberofElementsPerCNT(i)

DO k1=1,ME_NumberofCNTs2

Temp_X=M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),2)*M_NN(i1,1)+M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),2)*M_NN(i1,2)+&
M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),2)*M_NN(i1,3)+M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),2)*M_NN(i1,4)+&
M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),2)*M_NN(i1,5)+M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),2)*M_NN(i1,6)+&
M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),2)*M_NN(i1,7)+M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),2)*M_NN(i1,8)

Temp_Y=M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),3)*M_NN(i1,1)+M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),3)*M_NN(i1,2)+&
M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),3)*M_NN(i1,3)+M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),3)*M_NN(i1,4)+&
M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),3)*M_NN(i1,5)+M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),3)*M_NN(i1,6)+&
M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),3)*M_NN(i1,7)+M_Nodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),3)*M_NN(i1,8)

Temp_Distance= SQRT((Temp_X-ME_Coordinates2(1,k1,1))*(Temp_X-ME_Coordinates2(1,k1,1))+(Temp_Y-ME_Coordinates2(1,k1,2))*(Temp_Y-ME_Coordinates2(1,k1,2)))

IF (  Temp_Distance < ME_RofCNT ) THEN

!WRITE(*,*) Temp_Distance

IF ( ABS(ME_Coordinates2(1,k1,1)-(-ME_Width1/2.0D0))<ME_RofCNT .AND. ABS(ME_Coordinates2(1,k1,2)-(-ME_Width2/2.0D0))>ME_RofCNT .AND.&
ABS(ME_Coordinates2(1,k1,2)-(ME_Width2/2.0D0))>ME_RofCNT ) THEN ! Left

ME_DisPerCNT(i,1)=ME_DisPerCNT(i,1)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),2)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),2)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),2)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),2)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),2)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),2)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),2)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),2)*M_NN(i1,8)+ME_Width1*M_Mag+ME_Width1*Sub_StrainsXX_R2*M_Mag )*&
M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag

ME_DisPerCNT(i,2)=ME_DisPerCNT(i,2)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),3)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),3)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),3)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),3)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),3)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),3)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),3)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),3)*M_NN(i1,8)+ME_Width1*Sub_StrainsXY2_R2*2.0D0*M_Mag )* &
M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag

Vol_Temp=Vol_Temp+M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0

GOTO 11

ELSEIF ( ABS(ME_Coordinates2(1,k1,1)-(-ME_Width1/2.0D0))<ME_RofCNT .AND. ABS(ME_Coordinates2(1,k1,2)-(-ME_Width2/2.0D0))<ME_RofCNT ) THEN !Bottom Left

ME_DisPerCNT(i,1)=ME_DisPerCNT(i,1)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),2)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),2)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),2)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),2)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),2)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),2)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),2)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),2)*M_NN(i1,8)+ME_Width1*M_Mag+&
ME_Width1*Sub_StrainsXX_R2*M_Mag )*M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag

ME_DisPerCNT(i,2)=ME_DisPerCNT(i,2)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),3)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),3)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),3)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),3)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),3)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),3)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),3)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),3)*M_NN(i1,8)+ME_Width1*Sub_StrainsXY2_R2*2.0D0*M_Mag )* &
M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag

Vol_Temp=Vol_Temp+M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0

GOTO 11

ELSEIF ( ABS(ME_Coordinates2(1,k1,1)-(-ME_Width1/2.0D0))<ME_RofCNT .AND. ABS(ME_Coordinates2(1,k1,2)-(ME_Width2/2.0D0))<ME_RofCNT ) THEN  !Top Left

ME_DisPerCNT(i,1)=ME_DisPerCNT(i,1)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),2)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),2)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),2)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),2)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),2)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),2)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),2)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),2)*M_NN(i1,8)+ME_Width1*M_Mag+ME_Width1*Sub_StrainsXX_R2*M_Mag&
-2.0D0*ME_Width2*Sub_StrainsXY1_R2*M_Mag )*M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag

ME_DisPerCNT(i,2)=ME_DisPerCNT(i,2)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),3)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),3)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),3)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),3)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),3)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),3)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),3)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),3)*M_NN(i1,8)-ME_Width2*M_Mag&
+2.0D0*ME_Width1*Sub_StrainsXY2_R2*M_Mag-ME_Width2*Sub_StrainsYY_R2*M_Mag )* &
M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag

Vol_Temp=Vol_Temp+M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0

GOTO 11

ELSEIF ( ABS(ME_Coordinates2(1,k1,2)-(ME_Width2/2.0D0))<ME_RofCNT .AND. ABS(ME_Coordinates2(1,k1,1)-(-ME_Width1/2.0D0))>ME_RofCNT .AND.&
ABS(ME_Coordinates2(1,k1,1)-(ME_Width1/2.0D0))>ME_RofCNT ) THEN ! Top

ME_DisPerCNT(i,1)=ME_DisPerCNT(i,1)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),2)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),2)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),2)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),2)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),2)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),2)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),2)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),2)*M_NN(i1,8)-2.0D0*ME_Width2*Sub_StrainsXY1_R2*M_Mag )*&
M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag

ME_DisPerCNT(i,2)=ME_DisPerCNT(i,2)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),3)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),3)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),3)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),3)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),3)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),3)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),3)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),3)*M_NN(i1,8)-ME_Width2*M_Mag-ME_Width2*Sub_StrainsYY_R2*M_Mag )* &
M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag

Vol_Temp=Vol_Temp+M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0

GOTO 11

ELSEIF ( ABS(ME_Coordinates2(1,k1,2)-(ME_Width2/2.0D0))<ME_RofCNT .AND. ABS(ME_Coordinates2(1,k1,1)-(ME_Width1/2.0D0))<ME_RofCNT ) THEN !Top Right
ME_DisPerCNT(i,1)=ME_DisPerCNT(i,1)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),2)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),2)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),2)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),2)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),2)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),2)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),2)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),2)*M_NN(i1,8)-2.0D0*ME_Width2*Sub_StrainsXY1_R2*M_Mag )*&
M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag
ME_DisPerCNT(i,2)=ME_DisPerCNT(i,2)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),3)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),3)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),3)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),3)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),3)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),3)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),3)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),3)*M_NN(i1,8)-ME_Width2*M_Mag-ME_Width2*Sub_StrainsYY_R2*M_Mag )* &
M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag
Vol_Temp=Vol_Temp+M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0
GOTO 11

ELSE
ME_DisPerCNT(i,1)=ME_DisPerCNT(i,1)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),2)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),2)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),2)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),2)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),2)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),2)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),2)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),2)*M_NN(i1,8) )*&
M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag

ME_DisPerCNT(i,2)=ME_DisPerCNT(i,2)+( M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),2),3)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),3),3)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),4),3)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),5),3)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),6),3)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),7),3)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),8),3)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(ME_ElementsPerCNT(i,j),9),3)*M_NN(i1,8) )* &
M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0/M_Mag
Vol_Temp=Vol_Temp+M_Vol( i1,ME_Connectivities(ME_ElementsPerCNT(i,j),1) ) /2.0D0
GOTO 11
ENDIF

ENDIF
11 ENDDO

ENDDO
ENDDO !For i1

ME_DisPerCNT(i,1)=ME_DisPerCNT(i,1)/Vol_Temp
ME_DisPerCNT(i,2)=ME_DisPerCNT(i,2)/Vol_Temp
ENDDO



DO j=1,9
DO i=1,ME_NumberofCNTs2
IF (k==1 .AND. ZeroTimeStep==1) THEN
ME_Current_LocationPerCNT(j,i,1)=ME_Coordinates2(j,i,1)
ME_Current_LocationPerCNT(j,i,2)=ME_Coordinates2(j,i,2)
ELSE

IF (j==1) THEN

IF ( ABS(ME_Coordinates2(1,i,1)-(-ME_Width1/2.0D0))<ME_RofCNT .AND. ABS(ME_Coordinates2(1,i,2)-(-ME_Width2/2.0D0))>ME_RofCNT .AND.&
ABS(ME_Coordinates2(1,i,2)-(ME_Width2/2.0D0))>ME_RofCNT ) THEN ! Left
ME_Current_LocationPerCNT(j,i,1)=ME_DisPerCNT(ME_CNTIndex(i),1)-ME_Width1-ME_Width1*Sub_StrainsXX_R2*M_Mag
ME_Current_LocationPerCNT(j,i,2)=ME_DisPerCNT(ME_CNTIndex(i),2)-2.0D0*ME_Width1*Sub_StrainsXY2_R2*M_Mag
ELSEIF ( ABS(ME_Coordinates2(1,i,1)-(-ME_Width1/2.0D0))<ME_RofCNT .AND. ABS(ME_Coordinates2(1,i,2)-(-ME_Width2/2.0D0))<ME_RofCNT ) THEN !Bottom Left
ME_Current_LocationPerCNT(j,i,1)=ME_DisPerCNT(ME_CNTIndex(i),1)-ME_Width1-ME_Width1*Sub_StrainsXX_R2*M_Mag
ME_Current_LocationPerCNT(j,i,2)=ME_DisPerCNT(ME_CNTIndex(i),2)-2.0D0*ME_Width1*Sub_StrainsXY2_R2*M_Mag
ELSEIF ( ABS(ME_Coordinates2(1,i,1)-(-ME_Width1/2.0D0))<ME_RofCNT .AND. ABS(ME_Coordinates2(1,i,2)-(ME_Width2/2.0D0))<ME_RofCNT ) THEN  !Top Left
ME_Current_LocationPerCNT(j,i,1)=ME_DisPerCNT(ME_CNTIndex(i),1)-ME_Width1-ME_Width1*Sub_StrainsXX_R2*M_Mag+2.0D0*ME_Width2*Sub_StrainsXY1_R2*M_Mag
ME_Current_LocationPerCNT(j,i,2)=ME_DisPerCNT(ME_CNTIndex(i),2)+ME_Width2-2.0D0*ME_Width1*Sub_StrainsXY2_R2*M_Mag+ME_Width2*Sub_StrainsYY_R2*M_Mag
ELSEIF ( ABS(ME_Coordinates2(1,i,2)-(ME_Width2/2.0D0))<ME_RofCNT .AND. ABS(ME_Coordinates2(1,i,1)-(-ME_Width1/2.0D0))>ME_RofCNT .AND.&
ABS(ME_Coordinates2(1,i,1)-(ME_Width1/2.0D0))>ME_RofCNT ) THEN ! Top
ME_Current_LocationPerCNT(j,i,1)=ME_DisPerCNT(ME_CNTIndex(i),1)+2.0D0*ME_Width2*Sub_StrainsXY1_R2*M_Mag
ME_Current_LocationPerCNT(j,i,2)=ME_DisPerCNT(ME_CNTIndex(i),2)+ME_Width2+ME_Width2*Sub_StrainsYY_R2*M_Mag
ELSEIF ( ABS(ME_Coordinates2(1,i,2)-(ME_Width2/2.0D0))<ME_RofCNT .AND. ABS(ME_Coordinates2(1,i,1)-(ME_Width1/2.0D0))<ME_RofCNT ) THEN !Top Right
ME_Current_LocationPerCNT(j,i,1)=ME_DisPerCNT(ME_CNTIndex(i),1)+2.0D0*ME_Width2*Sub_StrainsXY1_R2*M_Mag
ME_Current_LocationPerCNT(j,i,2)=ME_DisPerCNT(ME_CNTIndex(i),2)+ME_Width2+ME_Width2*Sub_StrainsYY_R2*M_Mag
ELSE
ME_Current_LocationPerCNT(j,i,1)=ME_DisPerCNT(ME_CNTIndex(i),1)
ME_Current_LocationPerCNT(j,i,2)=ME_DisPerCNT(ME_CNTIndex(i),2)
ENDIF

ELSEIF (j==2) THEN
ME_Current_LocationPerCNT(j,i,1)=-ME_Width1-ME_Width11+ME_Current_LocationPerCNT(1,i,1)
ME_Current_LocationPerCNT(j,i,2)=-ME_Width12+ME_Current_LocationPerCNT(1,i,2)
ELSEIF (j==3) THEN
ME_Current_LocationPerCNT(j,i,1)=-ME_Width1-ME_Width11+ME_Width21+ME_Current_LocationPerCNT(1,i,1)
ME_Current_LocationPerCNT(j,i,2)=-ME_Width12+ME_Width2+ME_Width22+ME_Current_LocationPerCNT(1,i,2)
ELSEIF (j==4) THEN
ME_Current_LocationPerCNT(j,i,1)=ME_Width21+ME_Current_LocationPerCNT(1,i,1)
ME_Current_LocationPerCNT(j,i,2)=ME_Width22+ME_Width2+ME_Current_LocationPerCNT(1,i,2)
ELSEIF (j==5) THEN
ME_Current_LocationPerCNT(j,i,1)=ME_Width21+ME_Width1+ME_Width11+ME_Current_LocationPerCNT(1,i,1)
ME_Current_LocationPerCNT(j,i,2)=ME_Width22+ME_Width2+ME_Width12+ME_Current_LocationPerCNT(1,i,2)
ELSEIF (j==6) THEN
ME_Current_LocationPerCNT(j,i,1)=ME_Width1+ME_Width11+ME_Current_LocationPerCNT(1,i,1)
ME_Current_LocationPerCNT(j,i,2)=ME_Width12+ME_Current_LocationPerCNT(1,i,2)
ELSEIF (j==7) THEN
ME_Current_LocationPerCNT(j,i,1)=ME_Width1+ME_Width11-ME_Width21+ME_Current_LocationPerCNT(1,i,1)
ME_Current_LocationPerCNT(j,i,2)=ME_Width12-ME_Width2-ME_Width22+ME_Current_LocationPerCNT(1,i,2)
ELSEIF (j==8) THEN
ME_Current_LocationPerCNT(j,i,1)=-ME_Width21+ME_Current_LocationPerCNT(1,i,1)
ME_Current_LocationPerCNT(j,i,2)=-ME_Width22-ME_Width2+ME_Current_LocationPerCNT(1,i,2)
ELSEIF (j==9) THEN
ME_Current_LocationPerCNT(j,i,1)=-ME_Width21-ME_Width1-ME_Width11+ME_Current_LocationPerCNT(1,i,1)
ME_Current_LocationPerCNT(j,i,2)=-ME_Width22-ME_Width2-ME_Width12+ME_Current_LocationPerCNT(1,i,2)
ENDIF

ENDIF
ENDDO
ENDDO

WRITE(*,*) "The CNT Location information"
DO i=1,ME_NumberofCNTs2
WRITE(*,*) "CNT",ME_CNTIndex(i),"CNT2",i,ME_Current_LocationPerCNT(1,i,1),ME_Current_LocationPerCNT(1,i,2),"Next"
ENDDO

!DOUBLE PRECISION, DIMENSION(9,ME_NumberofCNTs2,2) :: 
!WRITE(*,*) ME_Current_LocationPerCNT(1,1,1),ME_Current_LocationPerCNT(1,1,2)
!WRITE(*,*) ME_Current_LocationPerCNT(1,2,1),ME_Current_LocationPerCNT(1,2,2)
!WRITE(*,*) ME_Current_LocationPerCNT(1,3,1),ME_Current_LocationPerCNT(1,3,2)
!WRITE(*,*) ME_Current_LocationPerCNT(1,4,1),ME_Current_LocationPerCNT(1,4,2)

DO i1=1,ME_NumberofGaussPoints
DO i=1,M_NumberofElements
ME_Original_Location(i1,i,1)=(M_Nodes(ME_Connectivities(i,2),2)*M_NN(i1,1)+M_Nodes(ME_Connectivities(i,3),2)*M_NN(i1,2)+&
M_Nodes(ME_Connectivities(i,4),2)*M_NN(i1,3)+M_Nodes(ME_Connectivities(i,5),2)*M_NN(i1,4)+M_Nodes(ME_Connectivities(i,6),2)*M_NN(i1,5)+&
M_Nodes(ME_Connectivities(i,7),2)*M_NN(i1,6)+M_Nodes(ME_Connectivities(i,8),2)*M_NN(i1,7)+M_Nodes(ME_Connectivities(i,9),2)*M_NN(i1,8))
ME_Original_Location(i1,i,2)=(M_Nodes(ME_Connectivities(i,2),3)*M_NN(i1,1)+M_Nodes(ME_Connectivities(i,3),3)*M_NN(i1,2)+&
M_Nodes(ME_Connectivities(i,4),3)*M_NN(i1,3)+M_Nodes(ME_Connectivities(i,5),3)*M_NN(i1,4)+M_Nodes(ME_Connectivities(i,6),3)*M_NN(i1,5)+&
M_Nodes(ME_Connectivities(i,7),3)*M_NN(i1,6)+M_Nodes(ME_Connectivities(i,8),3)*M_NN(i1,7)+M_Nodes(ME_Connectivities(i,9),3)*M_NN(i1,8))

IF (k==1 .AND. ZeroTimeStep==1) THEN
ME_Current_Location(i1,i,1)=ME_Original_Location(i1,i,1)
ME_Current_Location(i1,i,2)=ME_Original_Location(i1,i,2)
ELSE
ME_Current_Location(i1,i,1)=( M_CurrentNodes(ME_Connectivities(i,2),2)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(i,3),2)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(i,4),2)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(i,5),2)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(i,6),2)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(i,7),2)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(i,8),2)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(i,9),2)*M_NN(i1,8))/M_Mag
ME_Current_Location(i1,i,2)=( M_CurrentNodes(ME_Connectivities(i,2),3)*M_NN(i1,1)+&
M_CurrentNodes(ME_Connectivities(i,3),3)*M_NN(i1,2)+M_CurrentNodes(ME_Connectivities(i,4),3)*M_NN(i1,3)+&
M_CurrentNodes(ME_Connectivities(i,5),3)*M_NN(i1,4)+M_CurrentNodes(ME_Connectivities(i,6),3)*M_NN(i1,5)+&
M_CurrentNodes(ME_Connectivities(i,7),3)*M_NN(i1,6)+M_CurrentNodes(ME_Connectivities(i,8),3)*M_NN(i1,7)+&
M_CurrentNodes(ME_Connectivities(i,9),3)*M_NN(i1,8))/M_Mag
ENDIF
ENDDO
ENDDO


DO i3=1,ME_NumberofGaussPoints
DO i=1,M_NumberofElements
DO j=1,4
IF ( i >= ME_BeginningPolymer .AND. i<=ME_EndingPolymer) THEN
i1=ME_TunnelingIndex3(i3,i,j,1)
i2=ME_TunnelingIndex3(i3,i,j,2)
ME_TunnelingIndex3(i3,i,j,3)=SQRT( (ME_Current_LocationPerCNT(i1,i2,1)-ME_Current_Location(i3,i,1))*(ME_Current_LocationPerCNT(i1,i2,1)-ME_Current_Location(i3,i,1))+&
(ME_Current_LocationPerCNT(i1,i2,2)-ME_Current_Location(i3,i,2))*(ME_Current_LocationPerCNT(i1,i2,2)-ME_Current_Location(i3,i,2)) )-ME_RofCNT
ENDIF
ENDDO
ENDDO
ENDDO

DO i3=1,ME_NumberofGaussPoints
DO i=1,M_NumberofElements

IF ( i >= ME_BeginningPolymer .AND. i<=ME_EndingPolymer) THEN

DO i1=1,3
DO i2=i1+1,4

IF ( ME_TunnelingIndex3(i3,i,i1,3) <= ME_TunnelingIndex3(i3,i,i2,3) ) THEN
GOTO 7
ELSE
Temp1=ME_TunnelingIndex3(i3,i,i1,1)
Temp2=ME_TunnelingIndex3(i3,i,i1,2)
Temp3=ME_TunnelingIndex3(i3,i,i1,3)
ME_TunnelingIndex3(i3,i,i1,1)=ME_TunnelingIndex3(i3,i,i2,1)
ME_TunnelingIndex3(i3,i,i1,2)=ME_TunnelingIndex3(i3,i,i2,2)
ME_TunnelingIndex3(i3,i,i1,3)=ME_TunnelingIndex3(i3,i,i2,3)
ME_TunnelingIndex3(i3,i,i2,1)=Temp1
ME_TunnelingIndex3(i3,i,i2,2)=Temp2
ME_TunnelingIndex3(i3,i,i2,3)=Temp3
7 ENDIF

ENDDO
ENDDO
ENDIF

ENDDO
ENDDO !For i3

DO i3=1,ME_NumberofGaussPoints
DO i=1,M_NumberofElements

IF ( i >= ME_BeginningPolymer .AND. i<=ME_EndingPolymer) THEN

Temp1=ME_TunnelingIndex3(i3,i,1,3)+ME_TunnelingIndex3(i3,i,2,3)

IF ( Temp1< ME_Critical_Distance ) THEN
Temp4=0.0000316651D0
Temp5=10.24610D0
Temp6=2.71828182845904523536028747135D0
Resistivity_Temp=Temp4*(Temp6**(Temp5*Temp1*SQRT(ME_Lambda)))/SQRT(ME_Lambda)
ME_sub_Tunneling_rhorho(i3,i,1,1)=Resistivity_Temp
ME_sub_Tunneling_rhorho(i3,i,2,2)=Resistivity_Temp
ME_sub_Tunneling_rhorho(i3,i,3,3)=Resistivity_Temp

ELSE
ME_sub_Tunneling_rhorho(i3,i,1,1)=ME_sub_rhorho(i3,i,1,1)
ME_sub_Tunneling_rhorho(i3,i,2,2)=ME_sub_rhorho(i3,i,2,2)
ME_sub_Tunneling_rhorho(i3,i,3,3)=ME_sub_rhorho(i3,i,3,3)
ENDIF

!WRITE(*,*) ME_Tunneling_rhorho(i,1,1),ME_Tunneling_rhorho(i,2,2),ME_Tunneling_rhorho(i,3,3)

!IF (i==3128) THEN
!WRITE(*,*) ME_sub_Tunneling_rhorho(i3,i,1,1)
!ENDIF

ENDIF

ENDDO
ENDDO !For i3

!IF (k==1 .AND. ZeroTimeStep==1) THEN
!OPEN(30,file = '3D_Output_Micro_F/3D_Output_Microscale_1/Tunneling_Index_3_ME.txt',status ='unknown')
!DO i=1,7270 !1,M_NumberofElements !13081,13081 !11298,11298
!WRITE(30,10) i,INT(ME_TunnelingIndex3(1,i,1,1)),INT(ME_TunnelingIndex3(1,i,2,1)),INT(ME_TunnelingIndex3(1,i,3,1)),INT(ME_TunnelingIndex3(1,i,4,1))
!WRITE(30,10) i,INT(ME_TunnelingIndex3(1,i,1,2)),INT(ME_TunnelingIndex3(1,i,2,2)),INT(ME_TunnelingIndex3(1,i,3,2)),INT(ME_TunnelingIndex3(1,i,4,2))
!WRITE(30,20) i,ME_TunnelingIndex3(1,i,1,3),ME_TunnelingIndex3(1,i,2,3),ME_TunnelingIndex3(1,i,3,3),ME_TunnelingIndex3(1,i,4,3)
!ENDDO
!ELSE
!OPEN(30,file = '3D_Output_Micro_F/3D_Output_Microscale_1/Tunneling_Index_3_ME.txt',status ='unknown',position='append')
!DO i=1,7270   !1,M_NumberofElements  !13081,13081 !11298,11298
!WRITE(30,10) i,INT(ME_TunnelingIndex3(1,i,1,1)),INT(ME_TunnelingIndex3(1,i,2,1)),INT(ME_TunnelingIndex3(1,i,3,1)),INT(ME_TunnelingIndex3(1,i,4,1))
!WRITE(30,10) i,INT(ME_TunnelingIndex3(1,i,1,2)),INT(ME_TunnelingIndex3(1,i,2,2)),INT(ME_TunnelingIndex3(1,i,3,2)),INT(ME_TunnelingIndex3(1,i,4,2))
!WRITE(30,20) i,ME_TunnelingIndex3(1,i,1,3),ME_TunnelingIndex3(1,i,2,3),ME_TunnelingIndex3(1,i,3,3),ME_TunnelingIndex3(1,i,4,3)
!ENDDO
!ENDIF
!10 FORMAT(I,2X,I,2X,I,2X,I,2X,I)
!20 FORMAT(I,2X,E30.15,2X,E30.15,2X,E30.15,2X,E30.15)
!CLOSE(30)

DEALLOCATE(ME_TunnelingIndex3)

ENDSUBROUTINE ME_Tunneling_Effect

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_Prep_Tunneling(M_Nodes,ME_Connectivities,M_NN,ME_Coordinates,M_NumberofElementsPerCNT,ME_ElementsPerCNT,ME_TunnelingIndex,&
ME_TunnelingIndexDimension,ME_Width1,ME_Width2,ME_CNTIndex,ME_Critical_Distance,ME_Coordinates2)
IMPLICIT NONE
INTEGER :: i,j,k,i1,i2,i3,i4

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec
!The electrical tunneling control flags in the microscale
INTEGER :: ME_Tunneling
DOUBLE PRECISION :: ME_Vf
DOUBLE PRECISION :: ME_RofCNT
DOUBLE PRECISION :: ME_Lambda
DOUBLE PRECISION :: ME_rho
INTEGER :: ME_BeginningPolymer
INTEGER :: ME_EndingPolymer
INTEGER :: ME_NumberofCNTs
INTEGER :: ME_NumberofCNTs2
INTEGER :: ME_MaxNumberofElementsPerCNT
COMMON/Micro_Electrical_More_Inputs/ME_Tunneling,ME_BeginningPolymer,ME_EndingPolymer,ME_NumberofCNTs,ME_NumberofCNTs2,&
ME_MaxNumberofElementsPerCNT,ME_Vf,ME_RofCNT,ME_Lambda,ME_rho

DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes,4) :: M_Nodes
INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,10) :: ME_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofGaussPoints,8) :: M_NN      !Changed
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs) :: M_NumberofElementsPerCNT
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofCNTs2,2) :: ME_Coordinates
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs,ME_MaxNumberofElementsPerCNT) :: ME_ElementsPerCNT
DOUBLE PRECISION, INTENT(OUT), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,4,3) :: ME_TunnelingIndex    !Changed
INTEGER, INTENT(OUT), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements) :: ME_TunnelingIndexDimension  !Changed
INTEGER, INTENT(OUT), DIMENSION(ME_NumberofCNTs2) :: ME_CNTIndex
DOUBLE PRECISION, INTENT(OUT) :: ME_Critical_Distance
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: ME_TunnelingIndex2  !Changed
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ME_TunnelingIndexDimension2       !Changed
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3) :: ME_IntegrationPoints  !Changed
DOUBLE PRECISION, INTENT(OUT), DIMENSION(9,ME_NumberofCNTs2,2) :: ME_Coordinates2
DOUBLE PRECISION, INTENT(OUT) :: ME_Width1,ME_Width2
DOUBLE PRECISION :: ME_Temp_Distance
DOUBLE PRECISION :: Temp1,Temp2,Temp4,Temp5,Temp6

ALLOCATE(ME_TunnelingIndex2(ME_NumberofGaussPoints,M_NumberofElements,9*ME_NumberofCNTs2,3))
ALLOCATE(ME_TunnelingIndexDimension2(ME_NumberofGaussPoints,M_NumberofElements))
ME_TunnelingIndex2=0.0D0
ME_TunnelingIndexDimension2=0
ME_TunnelingIndex=0.0D0
ME_TunnelingIndexDimension=0

IF (ME_Tunneling==2) THEN
ME_Width1=SQRT(ME_NumberofCNTs*ME_RofCNT*ME_RofCNT*3.141592653589793/ME_Vf)
ME_Width2=SQRT(ME_NumberofCNTs*ME_RofCNT*ME_RofCNT*3.141592653589793/ME_Vf)
ELSEIF (ME_Tunneling==3) THEN
ME_Width1=SQRT( ME_NumberofCNTs*ME_RofCNT*ME_RofCNT*3.141592653589793/ME_Vf/2.0D0/SQRT(3.0D0) )*SQRT(3.0D0)
ME_Width2=SQRT( ME_NumberofCNTs*ME_RofCNT*ME_RofCNT*3.141592653589793/ME_Vf/2.0D0/SQRT(3.0D0) )*2.0D0
ENDIF

Temp1=0.0000316651D0
Temp2=10.24610D0
ME_Critical_Distance= LOG(ME_rho*SQRT(ME_Lambda)/Temp1)/SQRT(ME_Lambda)/Temp2

DO j=1,ME_NumberofGaussPoints
DO i=1,M_NumberofElements
ME_IntegrationPoints(j,i,1)=M_Nodes(ME_Connectivities(i,2),2)*M_NN(j,1)+M_Nodes(ME_Connectivities(i,3),2)*M_NN(j,2)+&
M_Nodes(ME_Connectivities(i,4),2)*M_NN(j,3)+M_Nodes(ME_Connectivities(i,5),2)*M_NN(j,4)+M_Nodes(ME_Connectivities(i,6),2)*M_NN(j,5)+&
M_Nodes(ME_Connectivities(i,7),2)*M_NN(j,6)+M_Nodes(ME_Connectivities(i,8),2)*M_NN(j,7)+M_Nodes(ME_Connectivities(i,9),2)*M_NN(j,8)

ME_IntegrationPoints(j,i,2)=M_Nodes(ME_Connectivities(i,2),3)*M_NN(j,1)+M_Nodes(ME_Connectivities(i,3),3)*M_NN(j,2)+&
M_Nodes(ME_Connectivities(i,4),3)*M_NN(j,3)+M_Nodes(ME_Connectivities(i,5),3)*M_NN(j,4)+M_Nodes(ME_Connectivities(i,6),3)*M_NN(j,5)+&
M_Nodes(ME_Connectivities(i,7),3)*M_NN(j,6)+M_Nodes(ME_Connectivities(i,8),3)*M_NN(j,7)+M_Nodes(ME_Connectivities(i,9),3)*M_NN(j,8)

ME_IntegrationPoints(j,i,3)=M_Nodes(ME_Connectivities(i,2),4)*M_NN(j,1)+M_Nodes(ME_Connectivities(i,3),4)*M_NN(j,2)+&
M_Nodes(ME_Connectivities(i,4),4)*M_NN(j,3)+M_Nodes(ME_Connectivities(i,5),4)*M_NN(j,4)+M_Nodes(ME_Connectivities(i,6),4)*M_NN(j,5)+&
M_Nodes(ME_Connectivities(i,7),4)*M_NN(j,6)+M_Nodes(ME_Connectivities(i,8),4)*M_NN(j,7)+M_Nodes(ME_Connectivities(i,9),4)*M_NN(j,8)
ENDDO
ENDDO


DO i=1,ME_NumberofCNTs2

ME_Coordinates2(1,i,1)=ME_Coordinates(i,1)
ME_Coordinates2(1,i,2)=ME_Coordinates(i,2)

ME_Coordinates2(2,i,1)=ME_Coordinates(i,1)-ME_Width1
ME_Coordinates2(2,i,2)=ME_Coordinates(i,2)

ME_Coordinates2(3,i,1)=ME_Coordinates(i,1)-ME_Width1
ME_Coordinates2(3,i,2)=ME_Coordinates(i,2)+ME_Width2

ME_Coordinates2(4,i,1)=ME_Coordinates(i,1)
ME_Coordinates2(4,i,2)=ME_Coordinates(i,2)+ME_Width2

ME_Coordinates2(5,i,1)=ME_Coordinates(i,1)+ME_Width1
ME_Coordinates2(5,i,2)=ME_Coordinates(i,2)+ME_Width2

ME_Coordinates2(6,i,1)=ME_Coordinates(i,1)+ME_Width1
ME_Coordinates2(6,i,2)=ME_Coordinates(i,2)

ME_Coordinates2(7,i,1)=ME_Coordinates(i,1)+ME_Width1
ME_Coordinates2(7,i,2)=ME_Coordinates(i,2)-ME_Width2

ME_Coordinates2(8,i,1)=ME_Coordinates(i,1)
ME_Coordinates2(8,i,2)=ME_Coordinates(i,2)-ME_Width2

ME_Coordinates2(9,i,1)=ME_Coordinates(i,1)-ME_Width1
ME_Coordinates2(9,i,2)=ME_Coordinates(i,2)-ME_Width2
ENDDO

DO i1=1,ME_NumberofGaussPoints
DO i=1,M_NumberofElements
DO j=1,9
DO k=1,ME_NumberofCNTs2

IF ( i<ME_BeginningPolymer .OR. i>ME_EndingPolymer) THEN
GOTO 5
ENDIF

ME_Temp_Distance=SQRT( ( ME_IntegrationPoints(i1,i,1)-ME_Coordinates2(j,k,1) )*( ME_IntegrationPoints(i1,i,1)-ME_Coordinates2(j,k,1) )+&
( ME_IntegrationPoints(i1,i,2)-ME_Coordinates2(j,k,2) )*( ME_IntegrationPoints(i1,i,2)-ME_Coordinates2(j,k,2) ) )-ME_RofCNT

IF (ME_Temp_Distance > ME_Critical_Distance) THEN
GOTO 1
ELSE

IF (j==2 .AND. ABS(ME_Coordinates2(j,k,1)-(-ME_Width1/2.0D0)) < ME_RofCNT ) THEN
GOTO 1
ELSEIF (j==3) THEN
IF (ABS(ME_Coordinates2(j,k,1)-(-ME_Width1/2.0D0)) < ME_RofCNT .OR. ABS(ME_Coordinates2(j,k,2)-(ME_Width2/2.0D0)) < ME_RofCNT) THEN
GOTO 1
ENDIF
ELSEIF (j==4 .AND. ABS(ME_Coordinates2(j,k,2)-(ME_Width2/2.0D0)) < ME_RofCNT ) THEN
GOTO 1
ELSEIF (j==5) THEN
IF (ABS(ME_Coordinates2(j,k,1)-(ME_Width1/2.0D0)) < ME_RofCNT .OR. ABS(ME_Coordinates2(j,k,2)-(ME_Width2/2.0D0)) < ME_RofCNT) THEN
GOTO 1
ENDIF
ELSEIF (j==6 .AND. ABS(ME_Coordinates2(j,k,1)-(ME_Width1/2.0D0)) < ME_RofCNT ) THEN
GOTO 1
ELSEIF (j==7) THEN
IF (ABS(ME_Coordinates2(j,k,1)-(ME_Width1/2.0D0)) < ME_RofCNT .OR. ABS(ME_Coordinates2(j,k,2)-(-ME_Width2/2.0D0)) < ME_RofCNT) THEN
GOTO 1
ENDIF
ELSEIF (j==8 .AND. ABS(ME_Coordinates2(j,k,2)-(-ME_Width2/2.0D0)) < ME_RofCNT ) THEN
GOTO 1
ELSEIF (j==9) THEN
IF (ABS(ME_Coordinates2(j,k,1)-(-ME_Width1/2.0D0)) < ME_RofCNT .OR. ABS(ME_Coordinates2(j,k,2)-(-ME_Width2/2.0D0)) < ME_RofCNT) THEN
GOTO 1
ENDIF
ELSE
ME_TunnelingIndexDimension2(i1,i)=ME_TunnelingIndexDimension2(i1,i)+1
ME_TunnelingIndex2(i1,i,ME_TunnelingIndexDimension2(i1,i),1)=j
ME_TunnelingIndex2(i1,i,ME_TunnelingIndexDimension2(i1,i),2)=k
ME_TunnelingIndex2(i1,i,ME_TunnelingIndexDimension2(i1,i),3)=ME_Temp_Distance

ENDIF

1 ENDIF

ENDDO
5 ENDDO
ENDDO
ENDDO

DO i1=1,ME_NumberofGaussPoints
DO i=1,M_NumberofElements

IF ( ME_TunnelingIndexDimension2(i1,i) >1 ) THEN

DO j=1,ME_TunnelingIndexDimension2(i1,i)-1
IF (j > 4) THEN
GOTO 2
ELSE
ME_TunnelingIndex(i1,i,j,1)=ME_TunnelingIndex2(i1,i,j,1)
ME_TunnelingIndex(i1,i,j,2)=ME_TunnelingIndex2(i1,i,j,2)
ME_TunnelingIndex(i1,i,j,3)=ME_TunnelingIndex2(i1,i,j,3)
DO k=j+1,ME_TunnelingIndexDimension2(i1,i)
IF ( ME_TunnelingIndex2(i1,i,k,3) <= ME_TunnelingIndex2(i1,i,j,3) ) THEN
ME_TunnelingIndex(i1,i,j,1)=ME_TunnelingIndex2(i1,i,k,1)
ME_TunnelingIndex(i1,i,j,2)=ME_TunnelingIndex2(i1,i,k,2)
ME_TunnelingIndex(i1,i,j,3)=ME_TunnelingIndex2(i1,i,k,3)
Temp4=ME_TunnelingIndex2(i1,i,k,1)
Temp5=ME_TunnelingIndex2(i1,i,k,2)
Temp6=ME_TunnelingIndex2(i1,i,k,3)
ME_TunnelingIndex2(i1,i,k,1)=ME_TunnelingIndex2(i1,i,j,1)
ME_TunnelingIndex2(i1,i,k,2)=ME_TunnelingIndex2(i1,i,j,2)
ME_TunnelingIndex2(i1,i,k,3)=ME_TunnelingIndex2(i1,i,j,3)
ME_TunnelingIndex2(i1,i,j,1)=Temp4
ME_TunnelingIndex2(i1,i,j,2)=Temp5
ME_TunnelingIndex2(i1,i,j,3)=Temp6

ENDIF
ENDDO

IF ( ME_TunnelingIndexDimension2(i1,i) <= 4) THEN
ME_TunnelingIndex(i1,i,ME_TunnelingIndexDimension2(i1,i),1)=ME_TunnelingIndex2(i1,i,ME_TunnelingIndexDimension2(i1,i),1)
ME_TunnelingIndex(i1,i,ME_TunnelingIndexDimension2(i1,i),2)=ME_TunnelingIndex2(i1,i,ME_TunnelingIndexDimension2(i1,i),2)
ME_TunnelingIndex(i1,i,ME_TunnelingIndexDimension2(i1,i),3)=ME_TunnelingIndex2(i1,i,ME_TunnelingIndexDimension2(i1,i),3)
ENDIF

ENDIF !For IF (j > 4) THEN 
ENDDO

2 ENDIF !IF ( ME_TunnelingIndexDimension(i1,i) > 1 ) THEN

ENDDO

ENDDO




DO i1=1,ME_NumberofGaussPoints
DO i=1,M_NumberofElements
DO j=1,4

IF (ME_TunnelingIndex(i1,i,j,3) /=0.0D0 ) THEN
ME_TunnelingIndexDimension(i1,i)=ME_TunnelingIndexDimension(i1,i)+1
ELSE
GOTO 3
ENDIF
3 ENDDO
ENDDO
ENDDO

DEALLOCATE(ME_TunnelingIndex2)
DEALLOCATE(ME_TunnelingIndexDimension2)


DO i1=1,ME_NumberofGaussPoints
DO i=1,ME_NumberofCNTs2

DO j=1,ME_NumberofCNTs
DO k=1,M_NumberofElementsPerCNT(j)
ME_Temp_Distance= SQRT( ( ME_IntegrationPoints(i1,ME_ElementsPerCNT(j,k),1)-ME_Coordinates(i,1) )*( ME_IntegrationPoints(i1,ME_ElementsPerCNT(j,k),1)-ME_Coordinates(i,1) )+&
( ME_IntegrationPoints(i1,ME_ElementsPerCNT(j,k),2)-ME_Coordinates(i,2) )*( ME_IntegrationPoints(i1,ME_ElementsPerCNT(j,k),2)-ME_Coordinates(i,2) ) )
IF ( ME_Temp_Distance<ME_RofCNT) THEN
ME_CNTIndex(i)=j
GOTO 6
ENDIF
ENDDO
6 ENDDO

ENDDO
ENDDO

!WRITE(*,*) "ME_Critical_Distance",ME_Critical_Distance,"ME_Width1",ME_Width1,"ME_Width2",ME_Width2

!OPEN(20,file = '3D_Output_Micro_F/3D_Output_Microscale_1/Tunneling_Index_ME.txt',status ='unknown')
!DO i=1,M_NumberofElements
!WRITE(20,10) i,INT(ME_TunnelingIndex(1,i,1,1)),INT(ME_TunnelingIndex(1,i,2,1)),INT(ME_TunnelingIndex(1,i,3,1)),INT(ME_TunnelingIndex(1,i,4,1))
!WRITE(20,10) i,INT(ME_TunnelingIndex(1,i,1,2)),INT(ME_TunnelingIndex(1,i,2,2)),INT(ME_TunnelingIndex(1,i,3,2)),INT(ME_TunnelingIndex(1,i,4,2))
!WRITE(20,20) i,ME_TunnelingIndex(1,i,1,3),ME_TunnelingIndex(1,i,2,3),ME_TunnelingIndex(1,i,3,3),ME_TunnelingIndex(1,i,4,3)
!ENDDO
!10 FORMAT(I,2X,I,2X,I,2X,I,2X,I)
!20 FORMAT(I,2X,E30.15,2X,E30.15,2X,E30.15,2X,E30.15)
!CLOSE(20)

RETURN
END SUBROUTINE ME_Prep_Tunneling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_MoreInputs2 (ii,More_Inputs_Microscale_2_E,M_NumberofElementsPerCNT,ME_Coordinates,ME_ElementsPerCNT)
IMPLICIT NONE
INTEGER, INTENT(IN) :: ii
CHARACTER*90, INTENT(IN) :: More_Inputs_Microscale_2_E
INTEGER :: ME_IN6
INTEGER :: i,j

!The electrical tunneling control flags in the microscale
INTEGER :: ME_Tunneling
DOUBLE PRECISION :: ME_Vf
DOUBLE PRECISION :: ME_RofCNT
DOUBLE PRECISION :: ME_Lambda
DOUBLE PRECISION :: ME_rho
INTEGER :: ME_BeginningPolymer
INTEGER :: ME_EndingPolymer
INTEGER :: ME_NumberofCNTs
INTEGER :: ME_NumberofCNTs2
INTEGER :: ME_MaxNumberofElementsPerCNT
COMMON/Micro_Electrical_More_Inputs/ME_Tunneling,ME_BeginningPolymer,ME_EndingPolymer,ME_NumberofCNTs,ME_NumberofCNTs2,&
ME_MaxNumberofElementsPerCNT,ME_Vf,ME_RofCNT,ME_Lambda,ME_rho

INTEGER, INTENT(OUT), DIMENSION(ME_NumberofCNTs) :: M_NumberofElementsPerCNT
DOUBLE PRECISION, INTENT(OUT), DIMENSION(ME_NumberofCNTs2,2) :: ME_Coordinates
INTEGER, INTENT(OUT), DIMENSION(ME_NumberofCNTs,ME_MaxNumberofElementsPerCNT) :: ME_ElementsPerCNT

ME_IN6=16+100*(ii-1)

OPEN(ME_IN6,file = More_Inputs_Microscale_2_E,status ='unknown')

DO i=1,ME_NumberofCNTs
READ(ME_IN6,*) M_NumberofElementsPerCNT(i)
ENDDO

DO i=1,ME_NumberofCNTs2
READ(ME_IN6,*) ME_Coordinates(i,1),ME_Coordinates(i,2)
ENDDO

ME_ElementsPerCNT=0

DO i=1,ME_NumberofCNTs
READ(ME_IN6,*) (ME_ElementsPerCNT(i,j),j=1,M_NumberofElementsPerCNT(i))
ENDDO

ClOSE(ME_IN6)

RETURN

END SUBROUTINE ME_MoreInputs2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_MoreInputs (ii,More_Inputs_Microscale_E) 
IMPLICIT NONE
INTEGER, INTENT(IN) :: ii
INTEGER :: ME_IN5
INTEGER :: i
CHARACTER*90, INTENT(IN) :: More_Inputs_Microscale_E

!The electrical tunneling control flags in the microscale
INTEGER :: ME_Tunneling
DOUBLE PRECISION :: ME_Vf
DOUBLE PRECISION :: ME_RofCNT
DOUBLE PRECISION :: ME_Lambda
DOUBLE PRECISION :: ME_rho
INTEGER :: ME_BeginningPolymer
INTEGER :: ME_EndingPolymer
INTEGER :: ME_NumberofCNTs
INTEGER :: ME_NumberofCNTs2
INTEGER :: ME_MaxNumberofElementsPerCNT
COMMON/Micro_Electrical_More_Inputs/ME_Tunneling,ME_BeginningPolymer,ME_EndingPolymer,ME_NumberofCNTs,ME_NumberofCNTs2,&
ME_MaxNumberofElementsPerCNT,ME_Vf,ME_RofCNT,ME_Lambda,ME_rho

DOUBLE PRECISION, DIMENSION (10) :: ME_More_Inputs

ME_IN5=15+100*(ii-1)

OPEN(ME_IN5,file =More_Inputs_Microscale_E,status ='unknown')

READ (ME_IN5,*) ME_More_Inputs(1)
READ (ME_IN5,*) ME_More_Inputs(2)
READ (ME_IN5,*) ME_More_Inputs(3)
READ (ME_IN5,*) ME_More_Inputs(4)
READ (ME_IN5,*) ME_More_Inputs(5)
READ (ME_IN5,*) ME_More_Inputs(6)
READ (ME_IN5,*) ME_More_Inputs(7)
READ (ME_IN5,*) ME_More_Inputs(8)
READ (ME_IN5,*) ME_More_Inputs(9)
READ (ME_IN5,*) ME_More_Inputs(10)

ME_Tunneling=ME_More_Inputs(1)
ME_Vf=ME_More_Inputs(2)
ME_RofCNT=ME_More_Inputs(3)
ME_Lambda=ME_More_Inputs(4)
ME_rho=ME_More_Inputs(5)
ME_BeginningPolymer=ME_More_Inputs(6)
ME_EndingPolymer=ME_More_Inputs(7)
ME_NumberofCNTs=ME_More_Inputs(8)
ME_NumberofCNTs2=ME_More_Inputs(9)
ME_MaxNumberofElementsPerCNT=ME_More_Inputs(10)

ClOSE(ME_IN5)

RETURN

END SUBROUTINE ME_MoreInputs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_TecplotGenerator (ii,iiii)
IMPLICIT NONE

INTEGER, INTENT(IN) :: ii,iiii
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Table1
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Table2
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: Table3
INTEGER :: IN1, IN2, IN3, IN4, OUT
INTEGER :: NumberofNodes, NumberofElements, TotalTimeSteps
INTEGER :: m,I,J
      
ALLOCATE(Table1(1000000,4))
ALLOCATE(Table2(1000000,18))
ALLOCATE(Table3(1000000,8))

IN1=51
IN2=61
IN3=71
IN4=81
OUT=91

IF (ii==1) THEN

IF (iiii==1) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_11.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_11.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale_E_11.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale_E_11.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_1/XiangRentecplot_Microscale_E_11.dat',status = 'unknown')
ELSEIF (iiii==2) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_22.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_22.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale_E_22.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale_E_22.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_1/XiangRentecplot_Microscale_E_22.dat',status = 'unknown')
ELSEIF (iiii==3) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_12.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_12.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale_E_12.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale_E_12.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_1/XiangRentecplot_Microscale_E_12.dat',status = 'unknown')
ELSEIF (iiii==4) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_13.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_13.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale_E_13.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale_E_13.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_1/XiangRentecplot_Microscale_E_13.dat',status = 'unknown')
ELSEIF (iiii==5) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_23.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_23.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale_E_23.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale_E_23.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_1/XiangRentecplot_Microscale_E_23.dat',status = 'unknown')
ENDIF

ELSEIF (ii==2) THEN

IF (iiii==1) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_11.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_11.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale_E_11.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale_E_11.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_2/XiangRentecplot_Microscale_E_11.dat',status = 'unknown')
ELSEIF (iiii==2) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_22.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_22.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale_E_22.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale_E_22.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_2/XiangRentecplot_Microscale_E_22.dat',status = 'unknown')
ELSEIF (iiii==3) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_12.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_12.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale_E_12.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale_E_12.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_2/XiangRentecplot_Microscale_E_12.dat',status = 'unknown')
ELSEIF (iiii==4) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_13.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_13.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale_E_13.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale_E_13.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_2/XiangRentecplot_Microscale_E_13.dat',status = 'unknown')
ELSEIF (iiii==5) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_23.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_23.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale_E_23.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale_E_23.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_2/XiangRentecplot_Microscale_E_23.dat',status = 'unknown')
ENDIF

ELSEIF (ii==3) THEN

IF (iiii==1) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_11.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_11.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale_E_11.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale_E_11.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_3/XiangRentecplot_Microscale_E_11.dat',status = 'unknown')
ELSEIF (iiii==2) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_22.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_22.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale_E_22.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale_E_22.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_3/XiangRentecplot_Microscale_E_22.dat',status = 'unknown')
ELSEIF (iiii==3) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_12.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_12.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale_E_12.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale_E_12.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_3/XiangRentecplot_Microscale_E_12.dat',status = 'unknown')
ELSEIF (iiii==4) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_13.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_13.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale_E_13.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale_E_13.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_3/XiangRentecplot_Microscale_E_13.dat',status = 'unknown')
ELSEIF (iiii==5) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_23.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_23.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale_E_23.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale_E_23.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_3/XiangRentecplot_Microscale_E_23.dat',status = 'unknown')
ENDIF

ELSEIF (ii==4) THEN

IF (iiii==1) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_11.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_11.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale_E_11.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale_E_11.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_4/XiangRentecplot_Microscale_E_11.dat',status = 'unknown')
ELSEIF (iiii==2) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_22.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_22.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale_E_22.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale_E_22.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_4/XiangRentecplot_Microscale_E_22.dat',status = 'unknown')
ELSEIF (iiii==3) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_12.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_12.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale_E_12.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale_E_12.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_4/XiangRentecplot_Microscale_E_12.dat',status = 'unknown')
ELSEIF (iiii==4) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_13.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_13.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale_E_13.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale_E_13.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_4/XiangRentecplot_Microscale_E_13.dat',status = 'unknown')
ELSEIF (iiii==5) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_23.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_23.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale_E_23.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale_E_23.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_4/XiangRentecplot_Microscale_E_23.dat',status = 'unknown')
ENDIF

ELSEIF (ii==5) THEN

IF (iiii==1) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_11.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_11.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale_E_11.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale_E_11.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_5/XiangRentecplot_Microscale_E_11.dat',status = 'unknown')
ELSEIF (iiii==2) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_22.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_22.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale_E_22.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale_E_22.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_5/XiangRentecplot_Microscale_E_22.dat',status = 'unknown')
ELSEIF (iiii==3) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_12.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_12.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale_E_12.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale_E_12.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_5/XiangRentecplot_Microscale_E_12.dat',status = 'unknown')
ELSEIF (iiii==4) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_13.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_13.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale_E_13.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale_E_13.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_5/XiangRentecplot_Microscale_E_13.dat',status = 'unknown')
ELSEIF (iiii==5) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_23.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_23.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale_E_23.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale_E_23.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_5/XiangRentecplot_Microscale_E_23.dat',status = 'unknown')
ENDIF

ENDIF
      
READ(IN3,*) NumberofNodes,NumberofElements,TotalTimeSteps

WRITE(OUT,101) 'TITLE = "XR Plot"','VARIABLES = x1, x2, x3, "Phi", "E_11","E_22","E_33","J_11", "J_22","J_33","kappa_11",&
"kappa_12","kappa_13","kappa_22","kappa_23","kappa_33","rho_11","rho_12","rho_13","rho_22","rho_23","rho_33"'

DO m=1,TotalTimeSteps+1
  
READ(IN1,*) ((Table1(I,J),J=1,4),I=1,NumberofNodes)
READ(IN2,*) ((Table2(I,J),J=1,18),I=1,NumberofElements)
READ(IN4,*) ((Table3(I,J),J=1,8),I=1,NumberofElements)
REWIND(IN4)

IF (m==1) THEN

WRITE(OUT,109) NumberofNodes,NumberofElements,1.0D0*(m-1)

WRITE(OUT,*) ((Table1(I,J),I=1,NumberofNodes),J=1,4)
WRITE(OUT,*) ((Table2(I,J),I=1,NumberofElements),J=1,18)
WRITE(OUT,130) ((Table3(I,J),J=1,8),I=1,NumberofElements)

ELSE

WRITE(OUT,110) NumberofNodes,NumberofElements,1.0D0*(m-1)

WRITE(OUT,*) ((Table1(I,J),I=1,NumberofNodes),J=1,4)
WRITE(OUT,*) ((Table2(I,J),I=1,NumberofElements),J=1,18)
WRITE(OUT,130) ((Table3(I,J),J=1,8),I=1,NumberofElements)

ENDIF

ENDDO

101  FORMAT (2A/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

109   Format ('ZONE T= "Undeformed Mesh XRFEM", datapacking=block,&
N =',I,', E =',I,', zonetype= FEBRICK, varlocation=&
([1,2,3,4]=nodal,[5-22]=cellcentered), SOLUTIONTIME=',E/)

110   Format ('ZONE T= "Deformed Mesh XRFEM", datapacking=block,&
N =',I,', E =',I,', zonetype= FEBRICK, varlocation=&
([1,2,3,4]=nodal,[5-22]=cellcentered), SOLUTIONTIME=',E/)

!120   Format (<NumberofNodes>E)
!121   Format (<NumberofElements>E)
130   Format (8I)       

DEALLOCATE(Table1)
DEALLOCATE(Table2)
DEALLOCATE(Table3)

CLOSE (IN1)
CLOSE (IN2)
CLOSE (IN3)
CLOSE (IN4)
CLOSE (OUT)


END SUBROUTINE ME_TecplotGenerator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_FeedForTecplot (ZeroTimeStep,k,ii,iiii,ME_Connectivities,M_Nodes,M_NodeIndex,ME_Dis1,ME_EXX,ME_EYY,ME_EZZ,ME_JXX,ME_JYY,ME_JZZ,&
M_CurrentNodes,ME_sigmasigma,ME_sub_sigmasigma_updated)
IMPLICIT NONE

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

INTEGER, INTENT(IN) :: ZeroTimeStep,k,ii,iiii
INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,10) :: ME_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes,4) :: M_Nodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNodes) :: M_NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes) :: ME_Dis1
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: ME_EXX
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: ME_EYY
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: ME_EZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: ME_JXX
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: ME_JYY
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: ME_JZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes,4) :: M_CurrentNodes
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofMaterialModes,3,3) :: ME_sigmasigma
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_sigmasigma_updated
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_rhorho_updated
DOUBLE PRECISION, DIMENSION(M_NumberofElements,3,3) :: ME_sigmasigma2
DOUBLE PRECISION, DIMENSION(M_NumberofElements,3,3) :: ME_rhorho2
INTEGER :: i,k1,i1,j1,OUT1,OUT2,OUT3,OUT4

ME_sigmasigma2=0.0D0
ME_rhorho2=0.0D0

DO k1=1,M_NumberofElements
DO i=1,8
CALL INV(k1,3,ME_sub_sigmasigma_updated(i,k1,:,:),ME_sub_rhorho_updated(i,k1,:,:))
ENDDO
ENDDO

DO k1=1,M_NumberofElements
DO i=1,8
DO i1=1,3
DO j1=1,3

ME_sigmasigma2(k1,i1,j1)=ME_sigmasigma2(k1,i1,j1)+ME_sub_sigmasigma_updated(i,k1,i1,j1)
ME_rhorho2(k1,i1,j1)=ME_rhorho2(k1,i1,j1)+ME_sub_rhorho_updated(i,k1,i1,j1)

ENDDO
ENDDO
ENDDO
ENDDO

DO k1=1,M_NumberofElements
DO i1=1,3
DO j1=1,3

ME_sigmasigma2(k1,i1,j1)=ME_sigmasigma2(k1,i1,j1)/8.0D0
ME_rhorho2(k1,i1,j1)=ME_rhorho2(k1,i1,j1)/8.0D0

ENDDO
ENDDO
ENDDO

OUT1=40
OUT2=41
OUT3=42
OUT4=43

IF (ii==1) THEN

IF (iiii==1) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_11.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_11.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_11.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_11.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==2) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_22.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_22.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_22.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_22.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==3) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_12.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_12.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_12.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_12.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==4) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_13.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_13.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_13.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_13.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==5) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_23.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_23.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale_E_23.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale_E_23.txt',status ='unknown',position='append')
ENDIF
ENDIF

ELSEIF (ii==2) THEN

IF (iiii==1) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_11.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_11.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_11.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_11.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==2) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_22.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_22.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_22.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_22.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==3) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_12.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_12.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_12.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_12.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==4) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_13.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_13.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_13.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_13.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==5) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_23.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_23.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale_E_23.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale_E_23.txt',status ='unknown',position='append')
ENDIF
ENDIF

ELSEIF (ii==3) THEN

IF (iiii==1) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_11.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_11.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_11.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_11.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==2) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_22.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_22.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_22.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_22.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==3) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_12.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_12.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_12.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_12.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==4) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_13.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_13.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_13.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_13.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==5) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_23.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_23.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale_E_23.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale_E_23.txt',status ='unknown',position='append')
ENDIF
ENDIF

ELSEIF (ii==4) THEN

IF (iiii==1) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_11.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_11.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_11.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_11.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==2) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_22.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_22.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_22.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_22.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==3) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_12.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_12.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_12.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_12.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==4) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_13.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_13.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_13.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_13.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==5) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_23.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_23.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale_E_23.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale_E_23.txt',status ='unknown',position='append')
ENDIF
ENDIF

ELSEIF (ii==5) THEN

IF (iiii==1) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_11.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_11.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_11.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_11.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==2) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_22.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_22.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_22.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_22.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==3) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_12.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_12.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_12.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_12.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==4) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_13.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_13.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_13.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_13.txt',status ='unknown',position='append')
ENDIF
ELSEIF (iiii==5) THEN
IF (k==1 .AND. ZeroTimeStep==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_23.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_23.txt',status ='unknown')
ELSE
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale_E_23.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale_E_23.txt',status ='unknown',position='append')
ENDIF
ENDIF

ENDIF

!For OUT1
IF (k==1 .AND. ZeroTimeStep==1) THEN

DO i=1,M_NumberofNodes
WRITE(OUT1,10) M_Nodes(i,2), M_Nodes(i,3), M_Nodes(i,4),ME_Dis1(i)
ENDDO

!DO i=1,M_NumberofNodes
!
!WRITE(OUT1,10) M_CurrentNodes(i,2), M_CurrentNodes(i,3), M_CurrentNodes(i,4), ME_Dis1(i)
!
!ENDDO

ELSE

DO i=1,M_NumberofNodes

WRITE(OUT1,10) M_CurrentNodes(i,2), M_CurrentNodes(i,3), M_CurrentNodes(i,4), ME_Dis1(i)


ENDDO

ENDIF

!For OUT2
IF (k==1 .AND. ZeroTimeStep==1) THEN

!DO i=1,M_NumberofElements
!
!WRITE(OUT2,20) 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,ME_sigmasigma2(i,1,1),&
!ME_sigmasigma2(i,1,2),ME_sigmasigma2(i,1,3),ME_sigmasigma2(i,2,2),&
!ME_sigmasigma2(i,2,3),ME_sigmasigma2(i,3,3)
!
!ENDDO

DO i=1,M_NumberofElements

WRITE(OUT2,20) ME_EXX(i),ME_EYY(i),ME_EZZ(i),ME_JXX(i),ME_JYY(i),ME_JZZ(i),ME_sigmasigma2(i,1,1),&
ME_sigmasigma2(i,1,2),ME_sigmasigma2(i,1,3),ME_sigmasigma2(i,2,2),ME_sigmasigma2(i,2,3),&
ME_sigmasigma2(i,3,3),ME_rhorho2(i,1,1),ME_rhorho2(i,1,2),ME_rhorho2(i,1,3),ME_rhorho2(i,2,2),&
ME_rhorho2(i,2,3),ME_rhorho2(i,3,3)

ENDDO

ELSE

DO i=1,M_NumberofElements

WRITE(OUT2,20) ME_EXX(i),ME_EYY(i),ME_EZZ(i),ME_JXX(i),ME_JYY(i),ME_JZZ(i),ME_sigmasigma2(i,1,1),&
ME_sigmasigma2(i,1,2),ME_sigmasigma2(i,1,3),ME_sigmasigma2(i,2,2),ME_sigmasigma2(i,2,3),&
ME_sigmasigma2(i,3,3),ME_rhorho2(i,1,1),ME_rhorho2(i,1,2),ME_rhorho2(i,1,3),ME_rhorho2(i,2,2),&
ME_rhorho2(i,2,3),ME_rhorho2(i,3,3)

ENDDO

ENDIF

IF (k==1 .AND. ZeroTimeStep==1) THEN


IF (ii==1) THEN

IF (iiii==1) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale_E_11.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale_E_11.txt',status ='unknown')
ELSEIF (iiii==2) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale_E_22.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale_E_22.txt',status ='unknown')
ELSEIF (iiii==3) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale_E_12.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale_E_12.txt',status ='unknown')
ELSEIF (iiii==4) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale_E_13.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale_E_13.txt',status ='unknown')
ELSEIF (iiii==5) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale_E_23.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale_E_23.txt',status ='unknown')
ENDIF

ELSEIF (ii==2) THEN

IF (iiii==1) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale_E_11.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale_E_11.txt',status ='unknown')
ELSEIF (iiii==2) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale_E_22.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale_E_22.txt',status ='unknown')
ELSEIF (iiii==3) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale_E_12.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale_E_12.txt',status ='unknown')
ELSEIF (iiii==4) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale_E_13.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale_E_13.txt',status ='unknown')
ELSEIF (iiii==5) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale_E_23.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale_E_23.txt',status ='unknown')
ENDIF

ELSEIF (ii==3) THEN

IF (iiii==1) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale_E_11.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale_E_11.txt',status ='unknown')
ELSEIF (iiii==2) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale_E_22.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale_E_22.txt',status ='unknown')
ELSEIF (iiii==3) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale_E_12.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale_E_12.txt',status ='unknown')
ELSEIF (iiii==4) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale_E_13.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale_E_13.txt',status ='unknown')
ELSEIF (iiii==5) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale_E_23.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale_E_23.txt',status ='unknown')
ENDIF

ELSEIF (ii==4) THEN

IF (iiii==1) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale_E_11.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale_E_11.txt',status ='unknown')
ELSEIF (iiii==2) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale_E_22.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale_E_22.txt',status ='unknown')
ELSEIF (iiii==3) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale_E_12.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale_E_12.txt',status ='unknown')
ELSEIF (iiii==4) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale_E_13.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale_E_13.txt',status ='unknown')
ELSEIF (iiii==5) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale_E_23.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale_E_23.txt',status ='unknown')
ENDIF

ELSEIF (ii==5) THEN

IF (iiii==1) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale_E_11.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale_E_11.txt',status ='unknown')
ELSEIF (iiii==2) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale_E_22.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale_E_22.txt',status ='unknown')
ELSEIF (iiii==3) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale_E_12.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale_E_12.txt',status ='unknown')
ELSEIF (iiii==4) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale_E_13.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale_E_13.txt',status ='unknown')
ELSEIF (iiii==5) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale_E_23.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale_E_23.txt',status ='unknown')
ENDIF

ENDIF

WRITE(OUT3,30) M_NumberofNodes, M_NumberofElements, TotalTimeSteps


DO i=1,M_NumberofElements

WRITE(OUT4,40) ME_Connectivities(i,2), ME_Connectivities(i,3), ME_Connectivities(i,4), ME_Connectivities(i,5),&
ME_Connectivities(i,6), ME_Connectivities(i,7), ME_Connectivities(i,8), ME_Connectivities(i,9)

ENDDO

ENDIF


10 FORMAT(4E/)
20 FORMAT(18E/)
30 FORMAT(I,I,I,E,I,I)
40 FORMAT(I,I,I,I,I,I,I,I)

CLOSE(OUT1)
CLOSE(OUT2)
CLOSE(OUT3)
CLOSE(OUT4)


END SUBROUTINE ME_FeedForTecplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_EandJ (ME_Vol,ME_Vol_Ele,ME_sigmasigma,ME_sub_sigmasigma_updated,ME_Connectivities,M_Nodes,M_NodeIndex,ME_Dis1,ME_EXX,ME_EYY,ME_EZZ,ME_JXX,&
ME_JYY,ME_JZZ,ME_DNDX,ME_DNDY,ME_DNDZ,ME_Sub_EXX,ME_Sub_EYY,ME_Sub_EZZ,ME_Sub_JXX,ME_Sub_JYY,ME_Sub_JZZ)
IMPLICIT NONE
INTEGER :: i,j,k

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,10) :: ME_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes,4) :: M_Nodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNodes) :: M_NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes) :: ME_Dis1
DOUBLE PRECISION, INTENT(IN), DIMENSION (ME_NumberofGaussPoints,M_NumberofElements) :: ME_Vol
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: ME_Vol_Ele
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofMaterialModes,3,3) :: ME_sigmasigma
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_sigmasigma_updated
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: ME_EXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: ME_EYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: ME_EZZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: ME_JXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: ME_JYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: ME_JZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,8) :: ME_DNDX,ME_DNDY,ME_DNDZ

DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_EXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_EYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_EZZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_JXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_JYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_JZZ


ME_EXX=0.0D0
ME_EYY=0.0D0
ME_EZZ=0.0D0
ME_JXX=0.0D0
ME_JYY=0.0D0
ME_JZZ=0.0D0

ME_Sub_EXX=0.0D0
ME_Sub_EYY=0.0D0
ME_Sub_EZZ=0.0D0
ME_Sub_JXX=0.0D0
ME_Sub_JYY=0.0D0
ME_Sub_JZZ=0.0D0
ME_Vol_Ele=0.0D0

DO k=1,M_NumberofElements
DO j=1,ME_NumberofGaussPoints

DO i=1,8
ME_Sub_EXX(k,j)=ME_Sub_EXX(k,j)+ME_Dis1(ME_Connectivities(k,i+1))*(-ME_DNDX(j,k,i))
ME_Sub_EYY(k,j)=ME_Sub_EYY(k,j)+ME_Dis1(ME_Connectivities(k,i+1))*(-ME_DNDY(j,k,i))
ME_Sub_EZZ(k,j)=ME_Sub_EZZ(k,j)+ME_Dis1(ME_Connectivities(k,i+1))*(-ME_DNDZ(j,k,i))
ENDDO

ME_Sub_JXX(k,j)=ME_sub_sigmasigma_updated(j,k,1,1)*ME_Sub_EXX(k,j)+ME_sub_sigmasigma_updated(j,k,1,2)*ME_Sub_EYY(k,j)+&
ME_sub_sigmasigma_updated(j,k,1,3)*ME_Sub_EZZ(k,j)
ME_Sub_JYY(k,j)=ME_sub_sigmasigma_updated(j,k,1,2)*ME_Sub_EXX(k,j)+ME_sub_sigmasigma_updated(j,k,2,2)*ME_Sub_EYY(k,j)+&
ME_sub_sigmasigma_updated(j,k,2,3)*ME_Sub_EZZ(k,j)
ME_Sub_JZZ(k,j)=ME_sub_sigmasigma_updated(j,k,1,3)*ME_Sub_EXX(k,j)+ME_sub_sigmasigma_updated(j,k,2,3)*ME_Sub_EYY(k,j)+&
ME_sub_sigmasigma_updated(j,k,3,3)*ME_Sub_EZZ(k,j)

ENDDO
ENDDO



DO k=1,M_NumberofElements
DO j=1,ME_NumberofGaussPoints

ME_EXX(k)=ME_EXX(k)+ME_Sub_EXX(k,j)*ME_Vol(j,k)
ME_EYY(k)=ME_EYY(k)+ME_Sub_EYY(k,j)*ME_Vol(j,k)
ME_EZZ(k)=ME_EZZ(k)+ME_Sub_EZZ(k,j)*ME_Vol(j,k)
ME_JXX(k)=ME_JXX(k)+ME_Sub_JXX(k,j)*ME_Vol(j,k)
ME_JYY(k)=ME_JYY(k)+ME_Sub_JYY(k,j)*ME_Vol(j,k)
ME_JZZ(k)=ME_JZZ(k)+ME_Sub_JZZ(k,j)*ME_Vol(j,k)
ME_Vol_Ele(k)=ME_Vol_Ele(k)+ME_Vol(j,k)

ENDDO

ME_EXX(k)=ME_EXX(k)/ME_Vol_Ele(k)
ME_EYY(k)=ME_EYY(k)/ME_Vol_Ele(k)
ME_EZZ(k)=ME_EZZ(k)/ME_Vol_Ele(k)
ME_JXX(k)=ME_JXX(k)/ME_Vol_Ele(k)
ME_JYY(k)=ME_JYY(k)/ME_Vol_Ele(k)
ME_JZZ(k)=ME_JZZ(k)/ME_Vol_Ele(k)

ENDDO

END SUBROUTINE ME_EandJ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_TecplotGenerator (ii)
IMPLICIT NONE

INTEGER, INTENT(IN) :: ii
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: M_Table1
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: M_Table2
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: M_Table3
INTEGER :: M_IN1, M_IN2, M_IN3, M_IN4, M_OUT
INTEGER :: M_NumberofNodes, M_NumberofElements, TotalTimeSteps
INTEGER :: m,I,J
      
ALLOCATE(M_Table1(2000000,6))
ALLOCATE(M_Table2(2000000,14))
ALLOCATE(M_Table3(2000000,8))

M_IN1=51
M_IN2=61
M_IN3=71
M_IN4=81
M_OUT=91

IF (ii==1) THEN      
OPEN (M_IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale.txt',status = 'old')
OPEN (M_IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale.txt',status = 'old')
OPEN (M_IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale.txt',status = 'old')
OPEN (M_IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale.txt',status = 'old')
OPEN (M_OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_1/XiangRentecplot_Microscale.dat',status = 'unknown')
ELSEIF (ii==2) THEN
OPEN (M_IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale.txt',status = 'old')
OPEN (M_IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale.txt',status = 'old')
OPEN (M_IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale.txt',status = 'old')
OPEN (M_IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale.txt',status = 'old')
OPEN (M_OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_2/XiangRentecplot_Microscale.dat',status = 'unknown')
ELSEIF (ii==3) THEN
OPEN (M_IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale.txt',status = 'old')
OPEN (M_IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale.txt',status = 'old')
OPEN (M_IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale.txt',status = 'old')
OPEN (M_IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale.txt',status = 'old')
OPEN (M_OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_3/XiangRentecplot_Microscale.dat',status = 'unknown')
ELSEIF (ii==4) THEN
OPEN (M_IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale.txt',status = 'old')
OPEN (M_IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale.txt',status = 'old')
OPEN (M_IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale.txt',status = 'old')
OPEN (M_IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale.txt',status = 'old')
OPEN (M_OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_4/XiangRentecplot_Microscale.dat',status = 'unknown')
ELSEIF (ii==5) THEN
OPEN (M_IN1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale.txt',status = 'old')
OPEN (M_IN2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale.txt',status = 'old')
OPEN (M_IN3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale.txt',status = 'old')
OPEN (M_IN4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale.txt',status = 'old')
OPEN (M_OUT,file = '3D_Output_Micro_F/3D_Output_Microscale_5/XiangRentecplot_Microscale.dat',status = 'unknown')
ENDIF

READ(M_IN3,*) M_NumberofNodes,M_NumberofElements,TotalTimeSteps

WRITE(M_OUT,101) 'TITLE = "XR Plot"','VARIABLES = x1, x2, x3, "u1", "u2", "u3", "epsilon_11","epsilon_22","epsilon_33","epsilon_12",&
"epsilon_13","epsilon_23","epsilon_v","sigma_11", "sigma_22","sigma_33","sigma_12","sigma_13","sigma_23","sigma_v"'

DO m=1,TotalTimeSteps+1
  
READ(M_IN1,*) ((M_Table1(I,J),J=1,6),I=1,M_NumberofNodes)
READ(M_IN2,*) ((M_Table2(I,J),J=1,14),I=1,M_NumberofElements)
READ(M_IN4,*) ((M_Table3(I,J),J=1,8),I=1,M_NumberofElements)
REWIND(M_IN4)

IF (m==1) THEN

WRITE(M_OUT,109) M_NumberofNodes,M_NumberofElements,1.0D0*(m-1)

WRITE(M_OUT,*) ((M_Table1(I,J),I=1,M_NumberofNodes),J=1,6)
WRITE(M_OUT,*) ((M_Table2(I,J),I=1,M_NumberofElements),J=1,14)
WRITE(M_OUT,130) ((M_Table3(I,J),J=1,8),I=1,M_NumberofElements)

ELSE

WRITE(M_OUT,110) M_NumberofNodes,M_NumberofElements,1.0D0*(m-1)

WRITE(M_OUT,*) ((M_Table1(I,J),I=1,M_NumberofNodes),J=1,6)
WRITE(M_OUT,*) ((M_Table2(I,J),I=1,M_NumberofElements),J=1,14)
WRITE(M_OUT,130) ((M_Table3(I,J),J=1,8),I=1,M_NumberofElements)

ENDIF

ENDDO

101  FORMAT (2A/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

109   Format ('ZONE T= "Undeformed Mesh XRFEM", datapacking=block,&
N =',I,', E =',I,', zonetype= FEBRICK, varlocation=&
([1,2,3,4,5,6]=nodal,[7-20]=cellcentered), SOLUTIONTIME=',E/)

110   Format ('ZONE T= "Deformed Mesh XRFEM", datapacking=block,&
N =',I,', E =',I,', zonetype= FEBRICK, varlocation=&
([1,2,3,4,5,6]=nodal,[7-20]=cellcentered), SOLUTIONTIME=',E/)

!120   Format (<NumberofNodes>E)
!121   Format (<NumberofElements>E)
130   Format (8I)       

DEALLOCATE(M_Table1)
DEALLOCATE(M_Table2)
DEALLOCATE(M_Table3)

CLOSE (M_IN1)
CLOSE (M_IN2)
CLOSE (M_IN3)
CLOSE (M_IN4)
CLOSE (M_OUT)

END SUBROUTINE M_TecplotGenerator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_FeedForTecplot (k,ii,M_Connectivities,M_Nodes,M_NodeIndex,M_Dis1,M_StrainsXX,M_StrainsYY,M_StrainsZZ,M_StrainsXY,M_StrainsXZ,M_StrainsYZ,&
M_StrainsVonMises,M_StressesXX,M_StressesYY,M_StressesZZ,M_StressesXY,M_StressesXZ,M_StressesYZ,M_StressesVonMises,M_CurrentNodes)
IMPLICIT NONE

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

INTEGER, INTENT(IN) :: k,ii
INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,6) :: M_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes,4) :: M_Nodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNodes) :: M_NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*M_NumberofNodes) :: M_Dis1
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StrainsXX
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StrainsYY
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StrainsZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StrainsXY
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StrainsXZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StrainsYZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StrainsVonMises
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StressesXX
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StressesYY
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StressesZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StressesXY
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StressesXZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StressesYZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements) :: M_StressesVonMises
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes,4) :: M_CurrentNodes
INTEGER :: i,M_OUT1,M_OUT2,M_OUT3,M_OUT4

M_OUT1=40
M_OUT2=41
M_OUT3=42
M_OUT4=43

IF (ii==1) THEN
IF (k==1) THEN
OPEN(M_OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale.txt',status ='unknown')
OPEN(M_OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(M_OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_1_Microscale.txt',status ='unknown',position='append')
OPEN(M_OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_2_Microscale.txt',status ='unknown',position='append')
ENDIF
ELSEIF (ii==2) THEN
IF (k==1) THEN
OPEN(M_OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale.txt',status ='unknown')
OPEN(M_OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(M_OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_1_Microscale.txt',status ='unknown',position='append')
OPEN(M_OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_2_Microscale.txt',status ='unknown',position='append')
ENDIF
ELSEIF (ii==3) THEN
IF (k==1) THEN
OPEN(M_OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale.txt',status ='unknown')
OPEN(M_OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(M_OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_1_Microscale.txt',status ='unknown',position='append')
OPEN(M_OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_2_Microscale.txt',status ='unknown',position='append')
ENDIF
ELSEIF (ii==4) THEN
IF (k==1) THEN
OPEN(M_OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale.txt',status ='unknown')
OPEN(M_OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(M_OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_1_Microscale.txt',status ='unknown',position='append')
OPEN(M_OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_2_Microscale.txt',status ='unknown',position='append')
ENDIF
ELSEIF (ii==5) THEN
IF (k==1) THEN
OPEN(M_OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale.txt',status ='unknown')
OPEN(M_OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(M_OUT1,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_1_Microscale.txt',status ='unknown',position='append')
OPEN(M_OUT2,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_2_Microscale.txt',status ='unknown',position='append')
ENDIF
ENDIF



!For M_OUT1
IF(k==1) THEN

DO i=1,M_NumberofNodes

WRITE(M_OUT1,10) M_Nodes(i,2), M_Nodes(i,3), M_Nodes(i,4),0.0D0, 0.0D0, 0.0D0


ENDDO

DO i=1,M_NumberofNodes

WRITE(M_OUT1,10) M_CurrentNodes(i,2), M_CurrentNodes(i,3), M_CurrentNodes(i,4), M_Dis1(i), M_Dis1(M_NumberofNodes+i), M_Dis1(2*M_NumberofNodes+i)

ENDDO

ELSEIF (k /= 1) THEN

DO i=1,M_NumberofNodes

WRITE(M_OUT1,10) M_CurrentNodes(i,2), M_CurrentNodes(i,3), M_CurrentNodes(i,4), M_Dis1(i), M_Dis1(M_NumberofNodes+i), M_Dis1(2*M_NumberofNodes+i)

ENDDO

ENDIF

!For M_OUT2
IF(k==1) THEN

DO i=1,M_NumberofElements

WRITE(M_OUT2,20) 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0


ENDDO

DO i=1,M_NumberofElements

WRITE(M_OUT2,20) M_StrainsXX(i),M_StrainsYY(i),M_StrainsZZ(i),M_StrainsXY(i),M_StrainsXZ(i),M_StrainsYZ(i),M_StrainsVonMises(i),&
M_StressesXX(i),M_StressesYY(i),M_StressesZZ(i),M_StressesXY(i),M_StressesXZ(i),M_StressesYZ(i),M_StressesVonMises(i)

ENDDO

ELSEIF (k /= 1) THEN

DO i=1,M_NumberofElements

WRITE(M_OUT2,20) M_StrainsXX(i),M_StrainsYY(i),M_StrainsZZ(i),M_StrainsXY(i),M_StrainsXZ(i),M_StrainsYZ(i),M_StrainsVonMises(i),&
M_StressesXX(i),M_StressesYY(i),M_StressesZZ(i),M_StressesXY(i),M_StressesXZ(i),M_StressesYZ(i),M_StressesVonMises(i)

ENDDO

ENDIF

IF (k==1) THEN

IF (ii==1) THEN
OPEN(M_OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_3_Microscale.txt',status ='unknown')
OPEN(M_OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_1/The_data_for_tecplot_4_Microscale.txt',status ='unknown')
ELSEIF (ii==2) THEN
OPEN(M_OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_3_Microscale.txt',status ='unknown')
OPEN(M_OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_2/The_data_for_tecplot_4_Microscale.txt',status ='unknown')
ELSEIF (ii==3) THEN
OPEN(M_OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_3_Microscale.txt',status ='unknown')
OPEN(M_OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_3/The_data_for_tecplot_4_Microscale.txt',status ='unknown')
ELSEIF (ii==4) THEN
OPEN(M_OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_3_Microscale.txt',status ='unknown')
OPEN(M_OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_4/The_data_for_tecplot_4_Microscale.txt',status ='unknown')
ELSEIF (ii==5) THEN
OPEN(M_OUT3,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_3_Microscale.txt',status ='unknown')
OPEN(M_OUT4,file = '3D_Output_Micro_F/3D_Output_Microscale_5/The_data_for_tecplot_4_Microscale.txt',status ='unknown')
ENDIF

WRITE(M_OUT3,30) M_NumberofNodes, M_NumberofElements, TotalTimeSteps


DO i=1,M_NumberofElements

WRITE(M_OUT4,40) M_Connectivities(i,2), M_Connectivities(i,3), M_Connectivities(i,4), M_Connectivities(i,5),&
M_Connectivities(i,6), M_Connectivities(i,7), M_Connectivities(i,8), M_Connectivities(i,9)    !,M_Connectivities(i,2)

ENDDO

ENDIF


10 FORMAT(6E/)
20 FORMAT(14E/)
30 FORMAT(I,I,I)
40 FORMAT(I,I,I,I,I,I,I,I)

CLOSE(M_OUT1)
CLOSE(M_OUT2)
CLOSE(M_OUT3)
CLOSE(M_OUT4)


END SUBROUTINE M_FeedForTecplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TecplotGenerator (ii)
IMPLICIT NONE

INTEGER, INTENT(IN) :: ii
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Table1
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Table2
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: Table3
INTEGER :: IN1, IN2, IN3, IN4, OUT
INTEGER :: NumberofNodes, NumberofElements, TotalTimeSteps, NumberofVariablesatNodes, TotalNumberofVariables
DOUBLE PRECISION :: DeltaT
INTEGER :: m,I,J
      
ALLOCATE(Table1(2000000,6))
ALLOCATE(Table2(2000000,14))
ALLOCATE(Table3(2000000,8))

IN1=51
IN2=61
IN3=71
IN4=81
OUT=91

IF (ii==1) THEN      
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_Mu12.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_Mu12.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_Mu12.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_Mu12.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot_Mu12.dat',status = 'unknown')
ELSEIF (ii==2) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_Mu23.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_Mu23.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_Mu23.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_Mu23.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot_Mu23.dat',status = 'unknown')
ELSEIF (ii==3) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_C11K12.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_C11K12.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_C11K12.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_C11K12.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot_C11K12.dat',status = 'unknown')
ELSEIF (ii==4) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_C33.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_C33.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_C33.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_C33.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot_C33Nu32.dat',status = 'unknown')
ELSEIF (ii==5) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_C12Nu32E33.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_C12Nu32E33.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_C12Nu32E33.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_C12Nu32E33.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot_C12Nu32E33.dat',status = 'unknown')
ELSEIF (ii==6) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot.dat',status = 'unknown')
ENDIF

READ(IN3,*) NumberofNodes,NumberofElements,TotalTimeSteps

WRITE(OUT,101) 'TITLE = "XR Plot"','VARIABLES = x1, x2, x3, "u1", "u2", "u3", "epsilon_11","epsilon_22","epsilon_33","epsilon_12",&
"epsilon_13","epsilon_23","epsilon_v","sigma_11", "sigma_22","sigma_33","sigma_12","sigma_13","sigma_23","sigma_v"'

DO m=1,TotalTimeSteps+1
  
READ(IN1,*) ((Table1(I,J),J=1,6),I=1,NumberofNodes)
READ(IN2,*) ((Table2(I,J),J=1,14),I=1,NumberofElements)
READ(IN4,*) ((Table3(I,J),J=1,8),I=1,NumberofElements)
REWIND(IN4)

IF (m==1) THEN

WRITE(OUT,109) NumberofNodes,NumberofElements,1.0D0*(m-1)

WRITE(OUT,*) ((Table1(I,J),I=1,NumberofNodes),J=1,6)
WRITE(OUT,*) ((Table2(I,J),I=1,NumberofElements),J=1,14)
WRITE(OUT,130) ((Table3(I,J),J=1,8),I=1,NumberofElements)

ELSE

WRITE(OUT,110) NumberofNodes,NumberofElements,1.0D0*(m-1)

WRITE(OUT,*) ((Table1(I,J),I=1,NumberofNodes),J=1,6)
WRITE(OUT,*) ((Table2(I,J),I=1,NumberofElements),J=1,14)
WRITE(OUT,130) ((Table3(I,J),J=1,8),I=1,NumberofElements)

ENDIF

ENDDO

101  FORMAT (2A/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

109   Format ('ZONE T= "Undeformed Mesh XRFEM", datapacking=block,&
N =',I,', E =',I,', zonetype= FEBRICK, varlocation=&
([1,2,3,4,5,6]=nodal,[7-20]=cellcentered), SOLUTIONTIME=',E/)

110   Format ('ZONE T= "Deformed Mesh XRFEM", datapacking=block,&
N =',I,', E =',I,', zonetype= FEBRICK, varlocation=&
([1,2,3,4,5,6]=nodal,[7-20]=cellcentered), SOLUTIONTIME=',E/)

!120   Format (<NumberofNodes>E)
!121   Format (<NumberofElements>E)
130   Format (8I)       

DEALLOCATE(Table1)
DEALLOCATE(Table2)
DEALLOCATE(Table3)

CLOSE (IN1)
CLOSE (IN2)
CLOSE (IN3)
CLOSE (IN4)
CLOSE (OUT)

END SUBROUTINE TecplotGenerator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FeedForTecplot (k,m,Connectivities,Nodes,NodeIndex,Dis1,StrainsXX,StrainsYY,StrainsZZ,StrainsXY,StrainsXZ,StrainsYZ,StrainsVonMises,StressesXX,&
StressesYY,StressesZZ,StressesXY,StressesXZ,StressesYZ,StressesVonMises,CurrentNodes)
IMPLICIT NONE

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(IN) :: k,m
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: Nodes
INTEGER, INTENT(IN), DIMENSION(NumberofNodes) :: NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*NumberofNodes) :: Dis1
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsXX
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsYY
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsXY
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsXZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsYZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StrainsVonMises
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesXX
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesYY
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesXY
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesXZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesYZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: StressesVonMises
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: CurrentNodes
INTEGER :: i,OUT1,OUT2,OUT3,OUT4

OUT1=40
OUT2=41
OUT3=42
OUT4=43

IF (m==1) THEN
IF (k==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_Mu12.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_Mu12.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_Mu12.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_Mu12.txt',status ='unknown',position='append')
ENDIF
ELSEIF (m==2) THEN
IF (k==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_Mu23.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_Mu23.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_Mu23.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_Mu23.txt',status ='unknown',position='append')
ENDIF
ELSEIF (m==3) THEN
IF (k==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_C11K12.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_C11K12.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_C11K12.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_C11K12.txt',status ='unknown',position='append')
ENDIF
ELSEIF (m==4) THEN
IF (k==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_C33.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_C33.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_C33.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_C33.txt',status ='unknown',position='append')
ENDIF
ELSEIF (m==5) THEN
IF (k==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_C12Nu32E33.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_C12Nu32E33.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_C12Nu32E33.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_C12Nu32E33.txt',status ='unknown',position='append')
ENDIF
ELSEIF (m==6) THEN
IF (k==1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1.txt',status ='unknown')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2.txt',status ='unknown')
ELSEIF (k /= 1) THEN
OPEN(OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1.txt',status ='unknown',position='append')
OPEN(OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2.txt',status ='unknown',position='append')
ENDIF
ENDIF

!For OUT1
IF(k==1) THEN

DO i=1,NumberofNodes

WRITE(OUT1,10) Nodes(i,2), Nodes(i,3), Nodes(i,4),0.0D0, 0.0D0, 0.0D0


ENDDO

DO i=1,NumberofNodes

WRITE(OUT1,10) CurrentNodes(i,2), CurrentNodes(i,3), CurrentNodes(i,4), Dis1(i), Dis1(NumberofNodes+i), Dis1(2*NumberofNodes+i)

ENDDO

ELSEIF (k /= 1) THEN

DO i=1,NumberofNodes

WRITE(OUT1,10) CurrentNodes(i,2), CurrentNodes(i,3), CurrentNodes(i,4), Dis1(i), Dis1(NumberofNodes+i), Dis1(2*NumberofNodes+i)

ENDDO

ENDIF

!For OUT2
IF(k==1) THEN

DO i=1,NumberofElements

WRITE(OUT2,20) 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0


ENDDO

DO i=1,NumberofElements

WRITE(OUT2,20) StrainsXX(i),StrainsYY(i),StrainsZZ(i),StrainsXY(i),StrainsXZ(i),StrainsYZ(i),StrainsVonMises(i),StressesXX(i),StressesYY(i),StressesZZ(i),&
StressesXY(i),StressesXZ(i),StressesYZ(i),StressesVonMises(i)

ENDDO

ELSEIF (k /= 1) THEN

DO i=1,NumberofElements

WRITE(OUT2,20) StrainsXX(i),StrainsYY(i),StrainsZZ(i),StrainsXY(i),StrainsXZ(i),StrainsYZ(i),StrainsVonMises(i),StressesXX(i),StressesYY(i),StressesZZ(i),&
StressesXY(i),StressesXZ(i),StressesYZ(i),StressesVonMises(i)

ENDDO

ENDIF

IF (k==1) THEN
IF (m==1) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_Mu12.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_Mu12.txt',status ='unknown')
ELSEIF (m==2) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_Mu23.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_Mu23.txt',status ='unknown')
ELSEIF (m==3) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_C11K12.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_C11K12.txt',status ='unknown')
ELSEIF (m==4) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_C33.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_C33.txt',status ='unknown')
ELSEIF (m==5) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_C12Nu32E33.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_C12Nu32E33.txt',status ='unknown')
ELSEIF (m==6) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3.txt',status ='unknown')
OPEN(OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4.txt',status ='unknown')
ENDIF

WRITE(OUT3,30) NumberofNodes, NumberofElements, TotalTimeSteps


DO i=1,NumberofElements

WRITE(OUT4,40) Connectivities(i,2), Connectivities(i,3), Connectivities(i,4), Connectivities(i,5),&
Connectivities(i,6), Connectivities(i,7), Connectivities(i,8), Connectivities(i,9)    !,Connectivities(i,2)

ENDDO

ENDIF


10 FORMAT(6E/)
20 FORMAT(14E/)
30 FORMAT(I,I,I)
40 FORMAT(I,I,I,I,I,I,I,I)

CLOSE(OUT1)
CLOSE(OUT2)
CLOSE(OUT3)
CLOSE(OUT4)


END SUBROUTINE FeedForTecplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_TecplotGenerator (ii)
IMPLICIT NONE

INTEGER, INTENT(IN) :: ii
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Table1
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: Table2
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: Table3
INTEGER :: IN1, IN2, IN3, IN4, OUT
INTEGER :: NumberofNodes, NumberofElements, TotalTimeSteps
INTEGER :: m,I,J
      
ALLOCATE(Table1(1000000,4))
ALLOCATE(Table2(1000000,18))
ALLOCATE(Table3(1000000,8))

IN1=51
IN2=61
IN3=71
IN4=81
OUT=91

IF (ii==1) THEN      
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_11.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_11.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_E_11.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_E_11.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot_E_11.dat',status = 'unknown')
ELSEIF (ii==2) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_22.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_22.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_E_22.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_E_22.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot_E_22.dat',status = 'unknown')
ELSEIF (ii==3) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_12.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_12.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_E_12.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_E_12.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot_E_12.dat',status = 'unknown')
ELSEIF (ii==4) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_13.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_13.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_E_13.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_E_13.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot_E_13.dat',status = 'unknown')
ELSEIF (ii==5) THEN
OPEN (IN1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_23.txt',status = 'old')
OPEN (IN2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_23.txt',status = 'old')
OPEN (IN3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_E_23.txt',status = 'old')
OPEN (IN4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_E_23.txt',status = 'old')
OPEN (OUT,file = '3D_Output_Micro_F/3D_Output/XiangRentecplot_E_23.dat',status = 'unknown')
ENDIF

READ(IN3,*) NumberofNodes,NumberofElements,TotalTimeSteps

WRITE(OUT,101) 'TITLE = "XR Plot"','VARIABLES = x1, x2, x3, "Phi", "E_11","E_22","E_33","J_11", "J_22","J_33","kappa_11",&
"kappa_12","kappa_13","kappa_22","kappa_23","kappa_33","rho_11","rho_12","rho_13","rho_22","rho_23","rho_33"'

DO m=1,TotalTimeSteps+1
  
READ(IN1,*) ((Table1(I,J),J=1,4),I=1,NumberofNodes)
READ(IN2,*) ((Table2(I,J),J=1,18),I=1,NumberofElements)
READ(IN4,*) ((Table3(I,J),J=1,8),I=1,NumberofElements)
REWIND(IN4)

IF (m==1) THEN

WRITE(OUT,109) NumberofNodes,NumberofElements,1.0D0*(m-1)

WRITE(OUT,*) ((Table1(I,J),I=1,NumberofNodes),J=1,4)
WRITE(OUT,*) ((Table2(I,J),I=1,NumberofElements),J=1,18)
WRITE(OUT,130) ((Table3(I,J),J=1,8),I=1,NumberofElements)

ELSE

WRITE(OUT,110) NumberofNodes,NumberofElements,1.0D0*(m-1)

WRITE(OUT,*) ((Table1(I,J),I=1,NumberofNodes),J=1,4)
WRITE(OUT,*) ((Table2(I,J),I=1,NumberofElements),J=1,18)
WRITE(OUT,130) ((Table3(I,J),J=1,8),I=1,NumberofElements)

ENDIF

ENDDO

101  FORMAT (2A/)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

109   Format ('ZONE T= "Undeformed Mesh XRFEM", datapacking=block,&
N =',I,', E =',I,', zonetype= FEBRICK, varlocation=&
([1,2,3,4]=nodal,[5-22]=cellcentered), SOLUTIONTIME=',E/)

110   Format ('ZONE T= "Deformed Mesh XRFEM", datapacking=block,&
N =',I,', E =',I,', zonetype= FEBRICK, varlocation=&
([1,2,3,4]=nodal,[5-22]=cellcentered), SOLUTIONTIME=',E/)

!120   Format (<NumberofNodes>E)
!121   Format (<NumberofElements>E)
130   Format (8I)       

DEALLOCATE(Table1)
DEALLOCATE(Table2)
DEALLOCATE(Table3)

CLOSE (IN1)
CLOSE (IN2)
CLOSE (IN3)
CLOSE (IN4)
CLOSE (OUT)


END SUBROUTINE E_TecplotGenerator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_FeedForTecplot (E_ZeroTimeStep,k,ii,Connectivities,Nodes,NodeIndex,CurrentNodes,E_Dis1,E_EXX,E_EYY,E_EZZ,E_JXX,E_JYY,E_JZZ,E_sub_sigmasigma,&
E_sub_rhorho)
IMPLICIT NONE

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

INTEGER, INTENT(IN) :: E_ZeroTimeStep,k,ii
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: Nodes
INTEGER, INTENT(IN), DIMENSION(NumberofNodes) :: NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes) :: E_Dis1
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: E_EXX
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: E_EYY
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: E_EZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: E_JXX
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: E_JYY
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements) :: E_JZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: CurrentNodes
DOUBLE PRECISION, INTENT(IN), DIMENSION(E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_sigmasigma
DOUBLE PRECISION, INTENT(IN), DIMENSION(E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_rhorho
DOUBLE PRECISION, DIMENSION(NumberofElements,3,3) :: E_sigmasigma2
DOUBLE PRECISION, DIMENSION(NumberofElements,3,3) :: E_rhorho2
INTEGER :: i,k1,i1,j1,E_OUT1,E_OUT2,E_OUT3,E_OUT4

E_sigmasigma2=0.0D0
E_rhorho2=0.0D0

DO k1=1,NumberofElements
DO i=1,8
DO i1=1,3
DO j1=1,3

E_sigmasigma2(k1,i1,j1)=E_sigmasigma2(k1,i1,j1)+E_sub_sigmasigma(i,k1,i1,j1)
E_rhorho2(k1,i1,j1)=E_rhorho2(k1,i1,j1)+E_sub_rhorho(i,k1,i1,j1)

ENDDO
ENDDO
ENDDO
ENDDO

DO k1=1,NumberofElements
DO i1=1,3
DO j1=1,3

E_sigmasigma2(k1,i1,j1)=E_sigmasigma2(k1,i1,j1)/8.0D0
E_rhorho2(k1,i1,j1)=E_rhorho2(k1,i1,j1)/8.0D0

ENDDO
ENDDO
ENDDO

E_OUT1=16
E_OUT2=17
E_OUT3=18
E_OUT4=19

IF (k==1 .AND. E_ZeroTimeStep==1) THEN
IF (ii==1) THEN
OPEN(E_OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_11.txt',status ='unknown')
OPEN(E_OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_11.txt',status ='unknown')
ELSEIF (ii==2) THEN
OPEN(E_OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_22.txt',status ='unknown')
OPEN(E_OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_22.txt',status ='unknown')
ELSEIF (ii==3) THEN
OPEN(E_OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_12.txt',status ='unknown')
OPEN(E_OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_12.txt',status ='unknown')
ELSEIF (ii==4) THEN
OPEN(E_OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_13.txt',status ='unknown')
OPEN(E_OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_13.txt',status ='unknown')
ELSE
OPEN(E_OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_23.txt',status ='unknown')
OPEN(E_OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_23.txt',status ='unknown')
ENDIF

ELSE
IF (ii==1) THEN
OPEN(E_OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_11.txt',status ='unknown',position='append')
OPEN(E_OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_11.txt',status ='unknown',position='append')
ELSEIF (ii==2) THEN
OPEN(E_OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_22.txt',status ='unknown',position='append')
OPEN(E_OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_22.txt',status ='unknown',position='append')
ELSEIF (ii==3) THEN
OPEN(E_OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_12.txt',status ='unknown',position='append')
OPEN(E_OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_12.txt',status ='unknown',position='append')
ELSEIF (ii==4) THEN
OPEN(E_OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_13.txt',status ='unknown',position='append')
OPEN(E_OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_13.txt',status ='unknown',position='append')
ELSE
OPEN(E_OUT1,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_1_E_23.txt',status ='unknown',position='append')
OPEN(E_OUT2,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_2_E_23.txt',status ='unknown',position='append')
ENDIF


ENDIF


!For E_OUT1
IF(k==1 .AND. E_ZeroTimeStep==1) THEN

DO i=1,NumberofNodes
WRITE(E_OUT1,10) Nodes(i,2), Nodes(i,3), Nodes(i,4),E_Dis1(i)
ENDDO

!DO i=1,NumberofNodes
!
!WRITE(E_OUT1,10) CurrentNodes(i,2), CurrentNodes(i,3), CurrentNodes(i,4), E_Dis1(i)
!
!ENDDO

ELSE

DO i=1,NumberofNodes

WRITE(E_OUT1,10) CurrentNodes(i,2), CurrentNodes(i,3), CurrentNodes(i,4), E_Dis1(i)


ENDDO

ENDIF

!For E_OUT2
IF (k==1 .AND. E_ZeroTimeStep==1) THEN

!DO i=1,NumberofElements

!WRITE(E_OUT2,20) 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0,E_sigmasigma2(i,1,1),&
!E_sigmasigma2(i,1,2),E_sigmasigma2(i,1,3),E_sigmasigma2(i,2,2),&
!E_sigmasigma2(i,2,3),E_sigmasigma2(i,3,3)
!
!ENDDO

DO i=1,NumberofElements

WRITE(E_OUT2,20) E_EXX(i),E_EYY(i),E_EZZ(i),E_JXX(i),E_JYY(i),E_JZZ(i),E_sigmasigma2(i,1,1),&
E_sigmasigma2(i,1,2),E_sigmasigma2(i,1,3),E_sigmasigma2(i,2,2),E_sigmasigma2(i,2,3),E_sigmasigma2(i,3,3),&
E_rhorho2(i,1,1),E_rhorho2(i,1,2),E_rhorho2(i,1,3),E_rhorho2(i,2,2),E_rhorho2(i,2,3),E_rhorho2(i,3,3)

ENDDO

ELSE

DO i=1,NumberofElements

WRITE(E_OUT2,20)  E_EXX(i),E_EYY(i),E_EZZ(i),E_JXX(i),E_JYY(i),E_JZZ(i),E_sigmasigma2(i,1,1),&
E_sigmasigma2(i,1,2),E_sigmasigma2(i,1,3),E_sigmasigma2(i,2,2),E_sigmasigma2(i,2,3),E_sigmasigma2(i,3,3),&
E_rhorho2(i,1,1),E_rhorho2(i,1,2),E_rhorho2(i,1,3),E_rhorho2(i,2,2),E_rhorho2(i,2,3),E_rhorho2(i,3,3)


ENDDO

ENDIF

IF (k==1 .AND. E_ZeroTimeStep==1) THEN

IF (ii==1) THEN
OPEN(E_OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_E_11.txt',status ='unknown')
OPEN(E_OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_E_11.txt',status ='unknown')
ELSEIF (ii==2) THEN
OPEN(E_OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_E_22.txt',status ='unknown')
OPEN(E_OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_E_22.txt',status ='unknown')
ELSEIF (ii==3) THEN
OPEN(E_OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_E_12.txt',status ='unknown')
OPEN(E_OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_E_12.txt',status ='unknown')
ELSEIF (ii==4) THEN
OPEN(E_OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_E_13.txt',status ='unknown')
OPEN(E_OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_E_13.txt',status ='unknown')
ELSE
OPEN(E_OUT3,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_3_E_23.txt',status ='unknown')
OPEN(E_OUT4,file = '3D_Output_Micro_F/3D_Output/The_data_for_tecplot_4_E_23.txt',status ='unknown')
ENDIF

WRITE(E_OUT3,30) NumberofNodes, NumberofElements, TotalTimeSteps


DO i=1,NumberofElements

WRITE(E_OUT4,40) Connectivities(i,2), Connectivities(i,3), Connectivities(i,4), Connectivities(i,5),&
Connectivities(i,6), Connectivities(i,7), Connectivities(i,8), Connectivities(i,9)

ENDDO

ENDIF


10 FORMAT(4E/)
20 FORMAT(18E/)
30 FORMAT(I,I,I)
40 FORMAT(I,I,I,I,I,I,I,I)

CLOSE(E_OUT1)
CLOSE(E_OUT2)
CLOSE(E_OUT3)
CLOSE(E_OUT4)


END SUBROUTINE E_FeedForTecplot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_EandJ (E_Vol,E_Vol_Ele,E_sub_sigmasigma,Connectivities,Nodes,NodeIndex,E_Dis1,E_EXX,E_EYY,E_EZZ,E_JXX,E_JYY,E_JZZ,E_DNDX,E_DNDY,&
E_DNDZ,E_Sub_EXX,E_Sub_EYY,E_Sub_EZZ,E_Sub_JXX,E_Sub_JYY,E_Sub_JZZ)
IMPLICIT NONE
INTEGER :: i,j,k

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: Nodes
INTEGER, INTENT(IN), DIMENSION(NumberofNodes) :: NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes) :: E_Dis1
DOUBLE PRECISION, INTENT(IN), DIMENSION (E_NumberofGaussPoints,NumberofElements) :: E_Vol
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: E_Vol_Ele
DOUBLE PRECISION, INTENT(IN), DIMENSION(E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_sigmasigma
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: E_EXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: E_EYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: E_EZZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: E_JXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: E_JYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: E_JZZ
DOUBLE PRECISION, INTENT(IN), DIMENSION(E_NumberofGaussPoints,NumberofElements,8) :: E_DNDX,E_DNDY,E_DNDZ

DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,E_NumberofGaussPoints) :: E_Sub_EXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,E_NumberofGaussPoints) :: E_Sub_EYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,E_NumberofGaussPoints) :: E_Sub_EZZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,E_NumberofGaussPoints) :: E_Sub_JXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,E_NumberofGaussPoints) :: E_Sub_JYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,E_NumberofGaussPoints) :: E_Sub_JZZ

E_EXX=0.0D0
E_EYY=0.0D0
E_EZZ=0.0D0
E_JXX=0.0D0
E_JYY=0.0D0
E_JZZ=0.0D0

E_Sub_EXX=0.0D0
E_Sub_EYY=0.0D0
E_Sub_EZZ=0.0D0
E_Sub_JXX=0.0D0
E_Sub_JYY=0.0D0
E_Sub_JZZ=0.0D0
E_Vol_Ele=0.0D0

DO k=1,NumberofElements
DO j=1,E_NumberofGaussPoints

DO i=1,8
E_Sub_EXX(k,j)=E_Sub_EXX(k,j)+E_Dis1(Connectivities(k,i+1))*(-E_DNDX(j,k,i))
E_Sub_EYY(k,j)=E_Sub_EYY(k,j)+E_Dis1(Connectivities(k,i+1))*(-E_DNDY(j,k,i))
E_Sub_EZZ(k,j)=E_Sub_EZZ(k,j)+E_Dis1(Connectivities(k,i+1))*(-E_DNDZ(j,k,i))
ENDDO

E_Sub_JXX(k,j)=E_sub_sigmasigma(j,k,1,1)*E_Sub_EXX(k,j)+E_sub_sigmasigma(j,k,1,2)*E_Sub_EYY(k,j)+&
E_sub_sigmasigma(j,k,1,3)*E_Sub_EZZ(k,j)
E_Sub_JYY(k,j)=E_sub_sigmasigma(j,k,2,1)*E_Sub_EXX(k,j)+E_sub_sigmasigma(j,k,2,2)*E_Sub_EYY(k,j)+&
E_sub_sigmasigma(j,k,2,3)*E_Sub_EZZ(k,j)
E_Sub_JZZ(k,j)=E_sub_sigmasigma(j,k,3,1)*E_Sub_EXX(k,j)+E_sub_sigmasigma(j,k,3,2)*E_Sub_EYY(k,j)+&
E_sub_sigmasigma(j,k,3,3)*E_Sub_EZZ(k,j)

ENDDO
ENDDO


DO k=1,NumberofElements
DO j=1,E_NumberofGaussPoints

E_EXX(k)=E_EXX(k)+E_Sub_EXX(k,j)*E_Vol(j,k)
E_EYY(k)=E_EYY(k)+E_Sub_EYY(k,j)*E_Vol(j,k)
E_EZZ(k)=E_EZZ(k)+E_Sub_EZZ(k,j)*E_Vol(j,k)
E_JXX(k)=E_JXX(k)+E_Sub_JXX(k,j)*E_Vol(j,k)
E_JYY(k)=E_JYY(k)+E_Sub_JYY(k,j)*E_Vol(j,k)
E_JZZ(k)=E_JZZ(k)+E_Sub_JZZ(k,j)*E_Vol(j,k)
E_Vol_Ele(k)=E_Vol_Ele(k)+E_Vol(j,k)

ENDDO

E_EXX(k)=E_EXX(k)/E_Vol_Ele(k)
E_EYY(k)=E_EYY(k)/E_Vol_Ele(k)
E_EZZ(k)=E_EZZ(k)/E_Vol_Ele(k)
E_JXX(k)=E_JXX(k)/E_Vol_Ele(k)
E_JYY(k)=E_JYY(k)/E_Vol_Ele(k)
E_JZZ(k)=E_JZZ(k)/E_Vol_Ele(k)

ENDDO


END SUBROUTINE E_EandJ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_Matrix_DE_Reduced(k,E_NumberofDE,E_One_Step_DEs,E_BoundaryIndex,E_Boundaries,E_K1GRR,E_K1GRCDD,E_K1GRCII,E_F1Globall,Connectivities,&
E_Dis1,E_Highest_GRCD)
IMPLICIT NONE
INTEGER :: i,j,m,n,i1,j1,k1,k2,k3
INTEGER, INTENT(IN) :: k

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

INTEGER, INTENT(IN) :: E_NumberofDE
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec) :: E_K1GRR
INTEGER, INTENT(IN), DIMENSION (NumberofNodes+E_NumberofDE) :: E_K1GRCDD
INTEGER, INTENT(IN), DIMENSION (NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec) :: E_K1GRCII
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes+E_NumberofDE) :: E_F1Globall
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: E_K1GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: E_K1GRCD
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: E_K1GRCI
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: E_F1Global

INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: E_K2GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: E_K2GRRD

DOUBLE PRECISION, INTENT(IN), DIMENSION(E_NumberofDE,6) :: E_One_Step_DEs
DOUBLE PRECISION, INTENT(IN), DIMENSION (E_NumberofBoundaryNodes,TotalTimeSteps+1) :: E_Boundaries
INTEGER, INTENT(IN), DIMENSION(E_NumberofBoundaryNodes,1) :: E_BoundaryIndex
DOUBLE PRECISION, DIMENSION (NumberofNodes+E_NumberofDE) :: E_Dis
DOUBLE PRECISION, INTENT(OUT),DIMENSION (NumberofNodes+E_NumberofDE) :: E_Dis1
INTEGER :: E_Dimen1
DOUBLE PRECISION, DIMENSION(NumberofNodes+E_NumberofDE) :: E_x
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
DOUBLE PRECISION :: ElapsedTime,Temp
INTEGER :: E_TEMP_GRCD
INTEGER, INTENT(OUT) :: E_Highest_GRCD
INTEGER :: E_Temp_INT
DOUBLE PRECISION :: E_Temp
INTEGER :: E_DUMMY,E_DUMMY2
INTEGER, ALLOCATABLE, DIMENSION(:) :: E_IA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: E_IA
INTEGER, ALLOCATABLE, DIMENSION(:) :: E_JA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: E_JA2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E_A1

ALLOCATE(E_K1GR(NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec))
ALLOCATE(E_K1GRCD(NumberofNodes+E_NumberofDE))
ALLOCATE(E_K1GRCI(NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec))
ALLOCATE(E_F1Global(NumberofNodes+E_NumberofDE))

E_K1GR=E_K1GRR
E_K1GRCD=E_K1GRCDD
E_K1GRCI=E_K1GRCII
E_F1Global=E_F1Globall

!!!Apply P.B.C.s (Top Right and Bottom Right of F)

DO i=1,E_NumberofDE
DO j=1,int(E_One_Step_DEs(i,1))

E_K1GRCD( 1*(int(E_One_Step_DEs(i,2*j))-1)+1 ) = E_K1GRCD( 1*(int(E_One_Step_DEs(i,2*j))-1)+1 )+1

E_K1GR( 1*(int(E_One_Step_DEs(i,2*j))-1)+1,E_K1GRCD( 1*(int(E_One_Step_DEs(i,2*j))-1)+1 ) ) = E_One_Step_DEs(i,2*j+1) 

E_K1GRCI( 1*(int(E_One_Step_DEs(i,2*j))-1)+1, E_K1GRCD( 1*(int(E_One_Step_DEs(i,2*j))-1)+1 ) ) = NumberofNodes+i

E_F1Global(1*NumberofNodes+i)=E_One_Step_DEs(i,2*int(E_One_Step_DEs(i,1))+1+1)

ENDDO
ENDDO

!!!Apply P.B.C.s (Bottom left)
DO i=1,E_NumberofDE
DO j=1,int(E_One_Step_DEs(i,1))
E_K1GRCD( NumberofNodes+i ) = E_K1GRCD( NumberofNodes+i )+1
E_K1GR( NumberofNodes+i,E_K1GRCD( NumberofNodes+i ) ) = E_One_Step_DEs(i,2*j+1) 
E_K1GRCI( NumberofNodes+i, E_K1GRCD( NumberofNodes+i ) ) = 1*(int(E_One_Step_DEs(i,2*j))-1)+1
ENDDO
ENDDO

ALLOCATE(E_K2GR(NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec,2))
ALLOCATE(E_K2GRRD(NumberofNodes+E_NumberofDE))
E_K2GR=0
E_K2GRRD=0

DO i=1,NumberofNodes+E_NumberofDE
DO j=1,E_K1GRCD(i)
E_K2GRRD(E_K1GRCI(i,j))=E_K2GRRD(E_K1GRCI(i,j))+1
E_K2GR(E_K1GRCI(i,j),E_K2GRRD(E_K1GRCI(i,j)),1)=i
E_K2GR(E_K1GRCI(i,j),E_K2GRRD(E_K1GRCI(i,j)),2)=j
ENDDO
ENDDO

E_Dis=0.0D0

DO i=1,E_NumberofBoundaryNodes
E_Dis(1*(E_BoundaryIndex(i,1)-1)+1)=E_Boundaries(i,k+1)
ENDDO


DO i=1,E_NumberofBoundaryNodes

DO m=1,E_K2GRRD(1*(E_BoundaryIndex(i,1)-1)+1)
E_F1Global( E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,1) )=E_F1Global( E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,1) )&
-E_K1GR( E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,1), E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,2) )*E_Dis(1*(E_BoundaryIndex(i,1)-1)+1)
ENDDO

E_F1Global(1*(E_BoundaryIndex(i,1)-1)+1)=E_Dis(1*(E_BoundaryIndex(i,1)-1)+1)

DO m=1,E_K2GRRD(1*(E_BoundaryIndex(i,1)-1)+1)
E_K1GR( E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,1), E_K2GR(1*(E_BoundaryIndex(i,1)-1)+1,m,2) )=0.0D0
ENDDO

DO m=1,E_K1GRCD(1*(E_BoundaryIndex(i,1)-1)+1)
IF ( 1*(E_BoundaryIndex(i,1)-1)+1 /= E_K1GRCI(1*(E_BoundaryIndex(i,1)-1)+1,m) ) THEN
E_K1GR( 1*(E_BoundaryIndex(i,1)-1)+1, m )=0.0D0
ELSE
E_K1GR( 1*(E_BoundaryIndex(i,1)-1)+1, m )=1.0D0
ENDIF
ENDDO

ENDDO

DEALLOCATE(E_K2GR)
DEALLOCATE(E_K2GRRD)

E_TEMP_GRCD=E_K1GRCD(1)
E_Highest_GRCD=E_TEMP_GRCD
DO i=2,1*NumberofNodes+E_NumberofDE
IF ( E_K1GRCD(i)>E_TEMP_GRCD ) THEN
E_TEMP_GRCD=E_K1GRCD(i)
E_Highest_GRCD=E_TEMP_GRCD
ENDIF
ENDDO

E_DUMMY=0

DO i=1,1*NumberofNodes+E_NumberofDE
DO j=1,E_K1GRCD(i)
IF (E_K1GR(i,j) /= 0.0D0) THEN
E_DUMMY=E_DUMMY+1
ELSEIF (E_K1GR(i,j) == 0.0D0) THEN
GOTO 3
ENDIF
3 ENDDO
ENDDO

ALLOCATE(E_IA1(E_DUMMY))
ALLOCATE(E_JA1(E_DUMMY))
ALLOCATE(E_JA2(E_DUMMY))
ALLOCATE(E_A1(E_DUMMY))

E_DUMMY=0

DO i=1,1*NumberofNodes+E_NumberofDE

!DO m=1,E_K1GRCD(i)                 !
!IF ( E_K1GRCI(i,m) >=i ) THEN     !
!E_Temp_INT=m                      !
!GOTO 21                          !
!ENDIF                            !
!ENDDO                            !
!
!IF ( E_K1GRCI(i,E_K1GRCD(i)) < i ) THEN  !
!E_DUMMY=E_DUMMY+1                 !
!E_A1(E_DUMMY)=0.0D0               !
!E_JA1(E_DUMMY)=i                  !
!E_JA2(E_DUMMY)=i                  !
!ENDIF                            !
!
!21 DO j=E_Temp_INT,E_K1GRCD(i) ! 1 is changed to E_Temp_INT
!E_DUMMY=E_DUMMY+1              ! 
!E_A1(E_DUMMY)=E_K1GR(i,j)      !
!E_JA1(E_DUMMY)=E_K1GRCI(i,j)   !
!E_JA2(E_DUMMY)=i               !
!ENDDO                         !

!!!!!!!!!!!!!!Original!!!!!!!!!!!!!!!!!!!!!
DO j=1,E_K1GRCD(i)
IF (E_K1GR(i,j) /= 0.0D0) THEN
E_DUMMY=E_DUMMY+1
E_A1(E_DUMMY)=E_K1GR(i,j)
E_JA1(E_DUMMY)=E_K1GRCI(i,j)
E_JA2(E_DUMMY)=i
ELSEIF (E_K1GR(i,j) == 0.0D0) THEN
GOTO 4
ENDIF
4 ENDDO

ENDDO

!OPEN(33,FILE="3D_Output_Micro_F/3D_Output/Check2D_KGR_E.txt",STATUS="UNKNOWN")
!DO i=1,NumberofNodes+E_NumberofDE
!WRITE(33,34) (E_K1GR(i,k1),k1=1,E_MaxDimenofInverseConnec)
!ENDDO
!34 FORMAT(<40>E15.8,/)
!WRITE(33,*) "E_Highest_GRCD  ",E_Highest_GRCD
!CLOSE(33)

DEALLOCATE(E_K1GR)
DEALLOCATE(E_K1GRCD)
DEALLOCATE(E_K1GRCI)

E_IA1(1)=1

E_DUMMY2=1

DO i=2,E_DUMMY

IF (E_JA2(i)>E_JA2(i-1) .AND. i>1) THEN

E_DUMMY2=E_DUMMY2+1

E_IA1(E_DUMMY2)=i

ENDIF

ENDDO

E_DUMMY2=E_DUMMY2+1

E_IA1(E_DUMMY2)=E_DUMMY+1

ALLOCATE(E_IA(E_DUMMY2))

DO i=1,E_DUMMY2
E_IA(i)=E_IA1(i)
ENDDO

DEALLOCATE(E_IA1)

!OPEN(27,FILE="3D_Output_Micro_F/3D_Output/A1JA1JA2_E.txt",STATUS="UNKNOWN")
!DO i=1,E_DUMMY
!WRITE(27,29) E_A1(i),E_JA1(i),E_JA2(i)
!ENDDO
!CLOSE(27)
!
!29 FORMAT(E15.8,2X,I,2X,I,/)
!
!OPEN(39,FILE="3D_Output_Micro_F/3D_Output/IA_E.txt",STATUS="UNKNOWN")
!DO i=1,E_DUMMY2
!WRITE(39,33) E_IA(i)
!ENDDO
!CLOSE(39)
!
!33 FORMAT(I,/)

E_Dimen1=1*NumberofNodes+E_NumberofDE

!CALL HarwellBoeing_Symmetric_Reduced(k,E_Dimen1,E_DUMMY,E_DUMMY2,E_A1,E_JA1,E_JA2,E_IA,E_F1Global,E_x)

CALL HarwellBoeing_Reduced(k,E_Dimen1,E_DUMMY,E_DUMMY2,E_A1,E_JA1,E_JA2,E_IA,E_F1Global,E_x)

DEALLOCATE(E_F1Global)
DEALLOCATE(E_IA)
DEALLOCATE(E_JA1)
DEALLOCATE(E_JA2)
DEALLOCATE(E_A1)

DO j=1,NumberofNodes
E_Dis1(j)=E_x(1*(j-1)+1)
!WRITE(*,*) E_Dis1(j)
ENDDO

RETURN
END SUBROUTINE E_Matrix_DE_Reduced

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_OneStepDEs(Nodes,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes,&
MacroEx_Knot,MacroEy_Knot,MacroEz_Knot,E_NumberofDE,E_One_Step_DEs)
IMPLICIT NONE
INTEGER :: i,j
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofNodes,4) :: Nodes
INTEGER, INTENT(IN) :: E_NumberofDE

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

INTEGER, INTENT(IN), DIMENSION(NumberofPositiveXNodes) :: PositiveXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveYNodes) :: PositiveYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveZNodes) :: PositiveZNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeXNodes) :: NegativeXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeYNodes) :: NegativeYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeZNodes) :: NegativeZNodes

DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NumberofPositiveXNodes+NumberofPositiveYNodes+NumberofPositiveZNodes,6) :: E_One_Step_DEs

DOUBLE PRECISION, INTENT(IN) :: MacroEx_Knot,MacroEy_Knot,MacroEz_Knot

INTEGER :: PosiTemp,NegaTemp

DO j=1,NumberofPositiveXNodes

PosiTemp=PositiveXNodes(j)
NegaTemp=NegativeXNodes(j)

IF ( 1*(PosiTemp-1)+1 < 1*(NegaTemp-1)+1 ) THEN
E_One_Step_DEs(j,1)=2
E_One_Step_DEs(j,2)=PositiveXNodes(j)
E_One_Step_DEs(j,3)=1.0D0
E_One_Step_DEs(j,4)=NegativeXNodes(j)
E_One_Step_DEs(j,5)=-1.0D0
E_One_Step_DEs(j,6)=MacroEx_Knot*(Nodes(PositiveXNodes(j),2)-Nodes(NegativeXNodes(j),2))
ELSE
E_One_Step_DEs(j,1)=2
E_One_Step_DEs(j,2)=NegativeXNodes(j)
E_One_Step_DEs(j,3)=-1.0D0
E_One_Step_DEs(j,4)=PositiveXNodes(j)
E_One_Step_DEs(j,5)=1.0D0
E_One_Step_DEs(j,6)=MacroEx_Knot*(Nodes(PositiveXNodes(j),2)-Nodes(NegativeXNodes(j),2))
ENDIF

ENDDO

DO j=NumberofPositiveXNodes+1,NumberofPositiveXNodes+NumberofPositiveYNodes

PosiTemp=PositiveYNodes(j-NumberofPositiveXNodes)
NegaTemp=NegativeYNodes(j-NumberofPositiveXNodes)

IF ( 1*(PosiTemp-1)+1 < 1*(NegaTemp-1)+1 ) THEN
E_One_Step_DEs(j,1)=2
E_One_Step_DEs(j,2)=PositiveYNodes(j-NumberofPositiveXNodes)
E_One_Step_DEs(j,3)=1.0D0
E_One_Step_DEs(j,4)=NegativeYNodes(j-NumberofPositiveXNodes)
E_One_Step_DEs(j,5)=-1.0D0
E_One_Step_DEs(j,6)=MacroEy_Knot*(Nodes(PositiveYNodes(j-NumberofPositiveXNodes),3)-&
Nodes(NegativeYNodes(j-NumberofPositiveXNodes),3))
ELSE
E_One_Step_DEs(j,1)=2
E_One_Step_DEs(j,2)=NegativeYNodes(j-NumberofPositiveXNodes)
E_One_Step_DEs(j,3)=-1.0D0
E_One_Step_DEs(j,4)=PositiveYNodes(j-NumberofPositiveXNodes)
E_One_Step_DEs(j,5)=1.0D0
E_One_Step_DEs(j,6)=MacroEy_Knot*(Nodes(PositiveYNodes(j-NumberofPositiveXNodes),3)-&
Nodes(NegativeYNodes(j-NumberofPositiveXNodes),3))
ENDIF

ENDDO

DO j=NumberofPositiveXNodes+NumberofPositiveYNodes+1,NumberofPositiveXNodes+NumberofPositiveYNodes+NumberofPositiveZNodes
PosiTemp=PositiveZNodes(j-NumberofPositiveXNodes-NumberofPositiveYNodes)
NegaTemp=NegativeZNodes(j-NumberofPositiveXNodes-NumberofPositiveYNodes)

IF ( 1*(PosiTemp-1)+1 < 1*(NegaTemp-1)+1 ) THEN
E_One_Step_DEs(j,1)=2
E_One_Step_DEs(j,2)=PositiveZNodes(j-NumberofPositiveXNodes-NumberofPositiveYNodes)
E_One_Step_DEs(j,3)=1.0D0
E_One_Step_DEs(j,4)=NegativeZNodes(j-NumberofPositiveXNodes-NumberofPositiveYNodes)
E_One_Step_DEs(j,5)=-1.0D0
E_One_Step_DEs(j,6)=MacroEz_Knot*(Nodes(PositiveZNodes(j-NumberofPositiveXNodes-NumberofPositiveYNodes),4)-&
Nodes(NegativeZNodes(j-NumberofPositiveXNodes-NumberofPositiveYNodes),4))
ELSE
E_One_Step_DEs(j,1)=2
E_One_Step_DEs(j,2)=NegativeZNodes(j-NumberofPositiveXNodes-NumberofPositiveYNodes)
E_One_Step_DEs(j,3)=-1.0D0
E_One_Step_DEs(j,4)=PositiveZNodes(j-NumberofPositiveXNodes-NumberofPositiveYNodes)
E_One_Step_DEs(j,5)=1.0D0
E_One_Step_DEs(j,6)=MacroEz_Knot*(Nodes(PositiveZNodes(j-NumberofPositiveXNodes-NumberofPositiveYNodes),4)-&
Nodes(NegativeZNodes(j-NumberofPositiveXNodes-NumberofPositiveYNodes),4))
ENDIF

ENDDO

RETURN
ENDSUBROUTINE E_OneStepDEs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_KPreprations(Connectivities,E_KComponent,E_NumberofDE,E_K1GRR,E_K1GRCDD,E_K1GRCII,E_F1Globall)
IMPLICIT NONE
INTEGER :: i,j,m,n,i1,j1,k1,k2,k3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements,8,8) :: E_KComponent
INTEGER, INTENT(IN) :: E_NumberofDE
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec) :: E_K1GRR
INTEGER, INTENT(OUT), DIMENSION (NumberofNodes+E_NumberofDE) :: E_K1GRCDD
INTEGER, INTENT(OUT), DIMENSION (NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec) :: E_K1GRCII
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofNodes+E_NumberofDE) :: E_F1Globall

E_K1GRR=0.0D0
E_K1GRCDD=0
E_K1GRCII=0

E_F1Globall=0.0D0

DO j1=1,1
DO i1=1,1
DO j=1,8
DO i=1,8
DO k1=1,NumberofElements

IF ( E_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j) /= 0.0D0 ) THEN

IF ( E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)==0 ) THEN

E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,1)=E_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)

E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)=E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)+1

E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,1)=1*(Connectivities(k1,j+1)-1)+j1

GOTO 1

ELSEIF ( E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1) == 1 ) THEN

IF ( 1*(Connectivities(k1,j+1)-1)+j1 == E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,1) ) THEN
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,1)=E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,1)+E_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ELSEIF ( 1*(Connectivities(k1,j+1)-1)+j1 > E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,1) ) THEN
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,2)=E_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)=E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)+1
E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,2)=1*(Connectivities(k1,j+1)-1)+j1
GOTO 1
ELSEIF ( 1*(Connectivities(k1,j+1)-1)+j1 < E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,1) ) THEN
E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)=E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)+1
DO k3=2,E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k3)=E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k3-1)
E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k3)=E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,1)=E_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,1)=1*(Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF

ELSEIF ( E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1) > 1 ) THEN

IF ( 1*(Connectivities(k1,j+1)-1)+j1 < E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,1) ) THEN
E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)=E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)+1
DO k3=E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1),2,-1
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k3)=E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k3-1)
E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k3)=E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,1)=E_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,1)=1*(Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF

DO k2=1,E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)-1

IF ( 1*(Connectivities(k1,j+1)-1)+j1 == E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k2) ) THEN
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k2)=E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k2)+E_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ENDIF

IF (1*(Connectivities(k1,j+1)-1)+j1 > E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k2) .AND. &
1*(Connectivities(k1,j+1)-1)+j1 < E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k2+1) ) THEN
E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)=E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)+1
DO k3=E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1),k2+2,-1
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k3)=E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k3-1)
E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k3)=E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k2+1)=E_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k2+1)=1*(Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF  

ENDDO  !k2

DO k2=E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1),E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)

IF ( 1*(Connectivities(k1,j+1)-1)+j1 == E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k2) ) THEN
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k2)=E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k2)+E_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ENDIF

IF ( 1*(Connectivities(k1,j+1)-1)+j1 > E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k2)  ) THEN
E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)=E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)+1
DO k3=k2+1,E_K1GRCDD(1*(Connectivities(k1,i+1)-1)+i1)
E_K1GRR(1*(Connectivities(k1,i+1)-1)+i1,k3)=E_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
E_K1GRCII(1*(Connectivities(k1,i+1)-1)+i1,k3)=1*(Connectivities(k1,j+1)-1)+j1
ENDDO
GOTO 1
ENDIF  

ENDDO !For k2

1 ENDIF  ! For ( K1GRCD(3*(Connectivities(k1,i+1)-1)+i1)==0 )
 
ENDIF  ! For ( KComponent(k1,8*(i1-1)+i,8*(j1-1)+j) /= 0.0D0 )

ENDDO
ENDDO
ENDDO
ENDDO
ENDDO

!DO i=1,20
!DO j=1,20
!IF (E_K1GRR(i,j) /=0.0D0) THEN
!WRITE(*,*) i,j,E_K1GRR(i,j)
!ENDIF
!ENDDO
!ENDDO

END SUBROUTINE E_KPreprations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_Components_Hexahedron(Nodes,Connectivities,E_Connectivities,E_sub_sigmasigma,E_Vol,E_KComponent,E_DNDX,E_DNDY,E_DNDZ)
IMPLICIT NONE
INTEGER :: i,j,k,i1,i2,i3,i4
INTEGER :: E_IN5=21

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: Nodes
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: E_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION(E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_sigmasigma
DOUBLE PRECISION, INTENT(OUT), DIMENSION (E_NumberofGaussPoints,NumberofElements) :: E_Vol
DOUBLE PRECISION, DIMENSION (E_NumberofGaussPoints,NumberofElements,3,8) :: E_BB
DOUBLE PRECISION, DIMENSION (E_NumberofGaussPoints,NumberofElements,8,3) :: E_BBT
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,8,8) :: E_KComponent

DOUBLE PRECISION, DIMENSION(NumberofElements,8) :: E_x,E_y,E_z

DOUBLE PRECISION, DIMENSION(8) :: E_Xi,E_Eta,E_Mu

DOUBLE PRECISION, DIMENSION(E_NumberofGaussPoints) :: E_Xi1,E_Eta1,E_Mu1,E_W

DOUBLE PRECISION, DIMENSION(E_NumberofGaussPoints,NumberofElements,3,3) :: E_JJ,E_InverseJJ

DOUBLE PRECISION, DIMENSION(E_NumberofGaussPoints,8) :: E_DNDXi,E_DNDEta,E_DNDMu
DOUBLE PRECISION, DIMENSION(8,E_NumberofGaussPoints) :: E_NN
DOUBLE PRECISION, INTENT(OUT), DIMENSION(E_NumberofGaussPoints,NumberofElements,8) :: E_DNDX,E_DNDY,E_DNDZ
DOUBLE PRECISION, DIMENSION(E_NumberofGaussPoints,NumberofElements,8) :: E_QX,E_QY,E_QZ

DOUBLE PRECISION :: E_Delta_A

DOUBLE PRECISION :: E_a11,E_a12,E_a13,E_a21,E_a22,E_a23,E_a31,E_a32,E_a33

DOUBLE PRECISION, DIMENSION(11,2) :: E_XYofAxis

DOUBLE PRECISION :: E_MaximumRadius

INTEGER :: E_BeginofModeIndex,E_EndofModeIndex

DOUBLE PRECISION :: E_Vf,E_R

INTEGER :: E_AnisotropicorElectric

E_JJ=0.0D0
E_InverseJJ=0.0D0

E_BB=0.0D0
E_BBT=0.0D0
E_KComponent=0.0D0

E_Xi=(/1.0D0,+1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,-1.0D0,-1.0D0/)
E_Eta=(/-1.0D0,1.0D0,1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,-1.0D0/)
E_Mu=(/-1.0D0,-1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,1.0D0,1.0D0/)

DO i=1,E_NumberofGaussPoints
E_Xi1(i)=E_Xi(i)/SQRT(3.0D0)
E_Eta1(i)=E_Eta(i)/SQRT(3.0D0)
E_Mu1(i)=E_Mu(i)/SQRT(3.0D0)
E_W(i)=1.0D0
ENDDO

DO k=1,NumberofElements
DO j=1,8
E_x(k,j)=Nodes(Connectivities(k,j+1),2)
E_y(k,j)=Nodes(Connectivities(k,j+1),3)
E_z(k,j)=Nodes(Connectivities(k,j+1),4)
ENDDO
ENDDO

DO i=1,E_NumberofGaussPoints

E_DNDXi(i,4)=-1.0D0/8.0D0*(1.0D0-E_Eta1(i))*(1.0D0-E_Mu1(i))
E_DNDEta(i,4)=-1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0-E_Mu1(i))
E_DNDMu(i,4)=-1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0-E_Eta1(i))

E_DNDXi(i,1)=1.0D0/8.0D0*(1.0D0-E_Eta1(i))*(1.0D0-E_Mu1(i))
E_DNDEta(i,1)=-1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0-E_Mu1(i))
E_DNDMu(i,1)=-1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0-E_Eta1(i))

E_DNDXi(i,2)=1.0D0/8.0D0*(1.0D0+E_Eta1(i))*(1.0D0-E_Mu1(i))
E_DNDEta(i,2)=1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0-E_Mu1(i))
E_DNDMu(i,2)=-1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0+E_Eta1(i))

E_DNDXi(i,3)=-1.0D0/8.0D0*(1.0D0+E_Eta1(i))*(1.0D0-E_Mu1(i))
E_DNDEta(i,3)=1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0-E_Mu1(i))
E_DNDMu(i,3)=-1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0+E_Eta1(i))

E_DNDXi(i,8)=-1.0D0/8.0D0*(1.0D0-E_Eta1(i))*(1.0D0+E_Mu1(i))
E_DNDEta(i,8)=-1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0+E_Mu1(i))
E_DNDMu(i,8)=1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0-E_Eta1(i))

E_DNDXi(i,5)=1.0D0/8.0D0*(1.0D0-E_Eta1(i))*(1.0D0+E_Mu1(i))
E_DNDEta(i,5)=-1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0+E_Mu1(i))
E_DNDMu(i,5)=1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0-E_Eta1(i))

E_DNDXi(i,6)=1.0D0/8.0D0*(1.0D0+E_Eta1(i))*(1.0D0+E_Mu1(i))
E_DNDEta(i,6)=1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0+E_Mu1(i))
E_DNDMu(i,6)=1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0+E_Eta1(i))

E_DNDXi(i,7)=-1.0D0/8.0D0*(1.0D0+E_Eta1(i))*(1.0D0+E_Mu1(i))
E_DNDEta(i,7)=1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0+E_Mu1(i))
E_DNDMu(i,7)=1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0+E_Eta1(i))

E_W(i)=1.0D0

ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (TestType==2) THEN

OPEN(E_IN5,file = '3D_Output_Micro_F/3D_Input/More_Inputs_E.txt',status ='unknown')
READ (E_IN5,*) E_Vf
READ (E_IN5,*) E_R
READ (E_IN5,*) E_BeginofModeIndex
READ (E_IN5,*) E_EndofModeIndex
READ (E_IN5,*) E_MaximumRadius !For the subroutine anisotropic
READ (E_IN5,*) E_AnisotropicorElectric
ClOSE(E_IN5)

E_XYofAxis(1,1)=0.0D0
E_XYofAxis(1,2)=-SQRT(4.0D0*3.141592653589793*E_R*E_R/2.0D0/SQRT(3.0D0)/E_Vf)
E_XYofAxis(2,1)=SQRT(4.0D0*3.141592653589793*E_R*E_R/2.0D0/SQRT(3.0D0)/E_Vf)*SQRT(3.0D0)/2.0D0
E_XYofAxis(2,2)=-0.5D0*SQRT(4.0D0*3.141592653589793*E_R*E_R/2.0D0/SQRT(3.0D0)/E_Vf)
E_XYofAxis(3,1)=E_XYofAxis(2,1)
E_XYofAxis(3,2)=-E_XYofAxis(2,2)
E_XYofAxis(4,1)=E_XYofAxis(1,1)
E_XYofAxis(4,2)=-E_XYofAxis(1,2)
E_XYofAxis(5,1)=-E_XYofAxis(3,1)
E_XYofAxis(5,2)=E_XYofAxis(3,2)
E_XYofAxis(6,1)=-E_XYofAxis(2,1)
E_XYofAxis(6,2)=E_XYofAxis(2,2)
E_XYofAxis(7,1)=0.0D0
E_XYofAxis(7,2)=0.0D0
E_XYofAxis(8,1)=E_XYofAxis(2,1)
E_XYofAxis(8,2)=-1.5D0*SQRT(4.0D0*3.141592653589793*E_R*E_R/2.0D0/SQRT(3.0D0)/E_Vf)
E_XYofAxis(9,1)=E_XYofAxis(8,1)
E_XYofAxis(9,2)=-E_XYofAxis(8,2)
E_XYofAxis(10,1)=-E_XYofAxis(9,1)
E_XYofAxis(10,2)=E_XYofAxis(9,2)
E_XYofAxis(11,1)=-E_XYofAxis(8,1)
E_XYofAxis(11,2)=E_XYofAxis(8,2)

DO i=1,E_NumberofGaussPoints
E_NN(1,i)=1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0-E_Eta1(i))*(1.0D0-E_Mu1(i))
E_NN(2,i)=1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0+E_Eta1(i))*(1.0D0-E_Mu1(i))
E_NN(3,i)=1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0+E_Eta1(i))*(1.0D0-E_Mu1(i))
E_NN(4,i)=1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0-E_Eta1(i))*(1.0D0-E_Mu1(i))
E_NN(5,i)=1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0-E_Eta1(i))*(1.0D0+E_Mu1(i))
E_NN(6,i)=1.0D0/8.0D0*(1.0D0+E_Xi1(i))*(1.0D0+E_Eta1(i))*(1.0D0+E_Mu1(i))
E_NN(7,i)=1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0+E_Eta1(i))*(1.0D0+E_Mu1(i))
E_NN(8,i)=1.0D0/8.0D0*(1.0D0-E_Xi1(i))*(1.0D0-E_Eta1(i))*(1.0D0+E_Mu1(i))
ENDDO

IF (E_AnisotropicorElectric == 1) THEN
CALL E_Anisotropic(Nodes,E_Connectivities,E_BeginofModeIndex,E_EndofModeIndex,E_XYofAxis,E_MaximumRadius,E_sub_sigmasigma,E_NN)
ELSE
ENDIF

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


DO i=1,E_NumberofGaussPoints
DO k=1,NumberofElements

DO j=1,8
E_JJ(i,k,1,1)=E_JJ(i,k,1,1)+E_x(k,j)*E_DNDXi(i,j)
E_JJ(i,k,1,2)=E_JJ(i,k,1,2)+E_y(k,j)*E_DNDXi(i,j)
E_JJ(i,k,1,3)=E_JJ(i,k,1,3)+E_z(k,j)*E_DNDXi(i,j)

E_JJ(i,k,2,1)=E_JJ(i,k,2,1)+E_x(k,j)*E_DNDEta(i,j)
E_JJ(i,k,2,2)=E_JJ(i,k,2,2)+E_y(k,j)*E_DNDEta(i,j)
E_JJ(i,k,2,3)=E_JJ(i,k,2,3)+E_z(k,j)*E_DNDEta(i,j)

E_JJ(i,k,3,1)=E_JJ(i,k,3,1)+E_x(k,j)*E_DNDMu(i,j)
E_JJ(i,k,3,2)=E_JJ(i,k,3,2)+E_y(k,j)*E_DNDMu(i,j)
E_JJ(i,k,3,3)=E_JJ(i,k,3,3)+E_z(k,j)*E_DNDMu(i,j)

ENDDO

ENDDO
ENDDO

DO i=1,E_NumberofGaussPoints
DO k=1,NumberofElements

E_a11=E_JJ(i,k,1,1)
E_a12=E_JJ(i,k,1,2)
E_a13=E_JJ(i,k,1,3)
E_a21=E_JJ(i,k,2,1)
E_a22=E_JJ(i,k,2,2)
E_a23=E_JJ(i,k,2,3)
E_a31=E_JJ(i,k,3,1)
E_a32=E_JJ(i,k,3,2)
E_a33=E_JJ(i,k,3,3)

E_Vol(i,k)=-E_a13*E_a22*E_a31 + E_a12*E_a23*E_a31 + E_a13*E_a21*E_a32&
- E_a11*E_a23*E_a32 - E_a12*E_a21*E_a33 +E_a11*E_a22*E_a33

E_InverseJJ(i,k,1,1)=(-E_a23*E_a32 + E_a22*E_a33)/E_Vol(i,k)
E_InverseJJ(i,k,1,2)=(E_a13*E_a32 - E_a12*E_a33)/E_Vol(i,k)
E_InverseJJ(i,k,1,3)=(-E_a13*E_a22 + E_a12*E_a23)/E_Vol(i,k)
E_InverseJJ(i,k,2,1)=(E_a23*E_a31 - E_a21*E_a33)/E_Vol(i,k)
E_InverseJJ(i,k,2,2)=(-E_a13*E_a31 + E_a11*E_a33)/E_Vol(i,k)
E_InverseJJ(i,k,2,3)=(E_a13*E_a21 - E_a11*E_a23)/E_Vol(i,k)
E_InverseJJ(i,k,3,1)=(-E_a22*E_a31 + E_a21*E_a32)/E_Vol(i,k)
E_InverseJJ(i,k,3,2)=(E_a12*E_a31 - E_a11*E_a32)/E_Vol(i,k)
E_InverseJJ(i,k,3,3)=(-E_a12*E_a21 + E_a11*E_a22)/E_Vol(i,k)

ENDDO
ENDDO

DO i=1,E_NumberofGaussPoints
DO k=1,NumberofElements
DO j=1,8

E_DNDX(i,k,j)=E_InverseJJ(i,k,1,1)*E_DNDXi(i,j)+E_InverseJJ(i,k,1,2)*E_DNDEta(i,j)+E_InverseJJ(i,k,1,3)*E_DNDMu(i,j)
E_DNDY(i,k,j)=E_InverseJJ(i,k,2,1)*E_DNDXi(i,j)+E_InverseJJ(i,k,2,2)*E_DNDEta(i,j)+E_InverseJJ(i,k,2,3)*E_DNDMu(i,j)
E_DNDZ(i,k,j)=E_InverseJJ(i,k,3,1)*E_DNDXi(i,j)+E_InverseJJ(i,k,3,2)*E_DNDEta(i,j)+E_InverseJJ(i,k,3,3)*E_DNDMu(i,j)

ENDDO
ENDDO
ENDDO

DO i=1,E_NumberofGaussPoints
DO k=1,NumberofElements
DO j=1,8

E_QX(i,k,j)=-E_DNDX(i,k,j)
E_QY(i,k,j)=-E_DNDY(i,k,j)
E_QZ(i,k,j)=-E_DNDZ(i,k,j)

ENDDO
ENDDO
ENDDO

DO i=1,E_NumberofGaussPoints
DO k=1,NumberofElements
DO j=1,8

E_BB(i,k,1,j)=E_QX(i,k,j)
E_BB(i,k,2,j)=E_QY(i,k,j)
E_BB(i,k,3,j)=E_QZ(i,k,j)

ENDDO
ENDDO
ENDDO

DO i=1,E_NumberofGaussPoints
DO k=1,NumberofElements
DO i1=1,3
DO i2=1,8
E_BBT(i,k,i2,i1)=E_BB(i,k,i1,i2)
ENDDO
ENDDO
ENDDO
ENDDO

!DO i=1,E_NumberofGaussPoints
!DO k1=1,NumberofElements
!DO i1=1,3
!DO i2=1,3
!IF ( E_sub_sigmasigma_updated(i,k1,i1,i1) <= 1.0D-1 ) THEN
!E_sub_sigmasigma_updated(i,k1,i1,i1)=1.0D-1
!ENDIF
!ENDDO
!ENDDO
!ENDDO
!ENDDO

DO i=1,E_NumberofGaussPoints
DO k=1,NumberofElements
DO i1=1,8
DO i2=1,3
DO i3=1,3
DO i4=1,8

E_KComponent(k,i1,i4)=E_KComponent(k,i1,i4)+&
E_BBT(i,k,i1,i2)*E_sub_sigmasigma(i,k,i2,i3)*E_BB(i,k,i3,i4)*E_Vol(i,k)*E_W(i)

ENDDO ! For i1
ENDDO ! For i2
ENDDO ! For i3 
ENDDO ! For i4
ENDDO
ENDDO

!DO i1=1,8
!DO i4=1,8
!IF(E_KComponent(1,i1,i4) /=0.0D0) THEN
!WRITE(*,*) i1,i4,E_KComponent(1,i1,i4)
!ENDIF
!ENDDO
!ENDDO

END SUBROUTINE E_Components_Hexahedron

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_ReadBoundaryConditions(E_BoundaryIndex,E_Boundaries)
IMPLICIT NONE
INTEGER :: E_IN2=15
INTEGER :: i,j

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

DOUBLE PRECISION, INTENT(OUT), DIMENSION(E_NumberofBoundaryNodes,TotalTimeSteps+1) :: E_Boundaries
INTEGER, INTENT(OUT), DIMENSION(E_NumberofBoundaryNodes,1) :: E_BoundaryIndex

OPEN(E_IN2,file = "3D_Output_Micro_F/3D_Input/Boundaries_E.txt",status ="unknown")
!READ (E_IN2,*) ((E_Boundaries(i,j),j=1,TotalTimeSteps+1),i=1,E_NumberofBoundaryNodes)

DO i=1,E_NumberofBoundaryNodes
READ (E_IN2,*) (E_Boundaries(i,j),j=1,TotalTimeSteps+1)
ENDDO

DO i=1,E_NumberofBoundaryNodes
E_BoundaryIndex(i,1)=E_Boundaries(i,1)
ENDDO
CLOSE(E_IN2)
RETURN
END SUBROUTINE E_ReadBoundaryConditions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE E_ReadControlFlags()
IMPLICIT NONE
INTEGER :: E_IN1=14
INTEGER :: i
DOUBLE PRECISION, DIMENSION (6) :: E_ControlFlags

!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

OPEN(E_IN1,file = '3D_Output_Micro_F/3D_Input/Control_Flags_E.txt',status ='unknown')
DO i=1,6
READ (E_IN1,*) E_ControlFlags(i)
ENDDO
E_NumberofBoundaryNodes=E_ControlFlags(1)
E_DisplacementEquationFlag=E_ControlFlags(2)
E_NumberofGaussPoints=E_ControlFlags(3)
E_Mag=E_ControlFlags(4)
E_MaxDimenofInverseConnec=E_ControlFlags(5)
E_NumberofMaterialModes=E_ControlFlags(6)
ClOSE(E_IN1)

RETURN
END SUBROUTINE E_ReadControlFlags

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Macro_Electrostatic(E_ZeroTimeStep,k,ii,MacroEx_Knot,MacroEy_Knot,MacroEz_Knot,E_sub_sigmasigma_in,E_Boundaries,E_BoundaryIndex,Nodes,&
NodeIndex,CurrentNodes,Connectivities,E_Connectivities,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes,E_Highest_GRCD,&
E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ2_Eff)
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes

INTEGER :: ii1,ii2
INTEGER, INTENT(IN) :: E_ZeroTimeStep,k,ii
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofNodes,3) :: Nodes
INTEGER, INTENT(IN), DIMENSION(NumberofNodes) :: NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofNodes,4) :: CurrentNodes
INTEGER, INTENT(IN), DIMENSION(NumberofElements,17) :: Connectivities
INTEGER, INTENT(IN), DIMENSION(NumberofElements,17) :: E_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION(TotalTimeSteps+1,E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_sigmasigma_in
DOUBLE PRECISION, DIMENSION(E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_sigmasigma
DOUBLE PRECISION, DIMENSION(E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_rhorho
DOUBLE PRECISION, INTENT(IN) :: MacroEx_Knot,MacroEy_Knot,MacroEz_Knot
DOUBLE PRECISION, INTENT(IN), DIMENSION(E_NumberofBoundaryNodes,TotalTimeSteps+1) :: E_Boundaries
INTEGER, INTENT(IN), DIMENSION(E_NumberofBoundaryNodes,1) :: E_BoundaryIndex
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveXNodes) :: PositiveXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveYNodes) :: PositiveYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveZNodes) :: PositiveZNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeXNodes) :: NegativeXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeYNodes) :: NegativeYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeZNodes) :: NegativeZNodes
INTEGER, INTENT(OUT) :: E_Highest_GRCD
DOUBLE PRECISION, INTENT(INOUT) :: E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ2_Eff
DOUBLE PRECISION, DIMENSION(NumberofPositiveXNodes+NumberofPositiveYNodes+NumberofPositiveZNodes,6) :: E_One_Step_DEs
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: E_DNDX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: E_DNDY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: E_DNDZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: E_Vol
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: E_KComponent
INTEGER :: E_NumberofDE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: E_K1GRR
INTEGER, ALLOCATABLE, DIMENSION(:) :: E_K1GRCDD
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: E_K1GRCII
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E_F1Globall
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E_Dis1
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E_EXX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E_EYY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E_EZZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E_JXX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E_JYY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: E_JZZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: E_Sub_EXX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: E_Sub_EYY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: E_Sub_EZZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: E_Sub_JXX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: E_Sub_JYY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: E_Sub_JZZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: E_Vol_Ele

IF (E_ZeroTimeStep==1 .AND. k==1) THEN
E_sub_sigmasigma=E_sub_sigmasigma_in(k,:,:,:,:)
DO ii1=1,E_NumberofGaussPoints
DO ii2=1,NumberofElements
CALL INV(ii1,3,E_sub_sigmasigma(ii1,ii2,:,:),E_sub_rhorho(ii1,ii2,:,:))
ENDDO
ENDDO
ELSE
E_sub_sigmasigma=E_sub_sigmasigma_in(k+1,:,:,:,:)
DO ii1=1,E_NumberofGaussPoints
DO ii2=1,NumberofElements
CALL INV(ii1,3,E_sub_sigmasigma(ii1,ii2,:,:),E_sub_rhorho(ii1,ii2,:,:))
ENDDO
ENDDO
ENDIF

ALLOCATE(E_DNDX(E_NumberofGaussPoints,NumberofElements,8))
ALLOCATE(E_DNDY(E_NumberofGaussPoints,NumberofElements,8))
ALLOCATE(E_DNDZ(E_NumberofGaussPoints,NumberofElements,8))
ALLOCATE(E_KComponent(NumberofElements,8,8))
ALLOCATE(E_Vol(E_NumberofGaussPoints,NumberofElements))
CALL E_Components_Hexahedron(Nodes,Connectivities,E_Connectivities,E_sub_sigmasigma,E_Vol,E_KComponent,E_DNDX,E_DNDY,E_DNDZ)

E_NumberofDE=1*(NumberofPositiveXNodes+NumberofPositiveYNodes+NumberofPositiveZNodes)
ALLOCATE(E_K1GRR(NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec))
ALLOCATE(E_K1GRCDD(NumberofNodes+E_NumberofDE))
ALLOCATE(E_K1GRCII(NumberofNodes+E_NumberofDE,E_MaxDimenofInverseConnec))
ALLOCATE(E_F1Globall(NumberofNodes+E_NumberofDE))
CALL E_KPreprations(Connectivities,E_KComponent,E_NumberofDE,E_K1GRR,E_K1GRCDD,E_K1GRCII,E_F1Globall)

CALL E_OneStepDEs(Nodes,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes,&
MacroEx_Knot,MacroEy_Knot,MacroEz_Knot,E_NumberofDE,E_One_Step_DEs)

ALLOCATE(E_Dis1(NumberofNodes+E_NumberofDE))
IF (E_DisplacementEquationFlag==1) THEN
CALL E_Matrix_DE_Reduced(k,E_NumberofDE,E_One_Step_DEs,E_BoundaryIndex,E_Boundaries,E_K1GRR,E_K1GRCDD,E_K1GRCII,E_F1Globall,Connectivities,&
E_Dis1,E_Highest_GRCD)
ELSEIF (E_DisplacementEquationFlag==0) THEN
CALL E_Matrix_Reduced(k,E_NumberofDE,E_One_Step_DEs,E_BoundaryIndex,E_Boundaries,E_K1GRR,E_K1GRCDD,E_K1GRCII,E_F1Globall,Connectivities,&
E_Dis1,E_Highest_GRCD)
ENDIF

ALLOCATE(E_EXX(NumberofElements))
ALLOCATE(E_EYY(NumberofElements))
ALLOCATE(E_EZZ(NumberofElements))
ALLOCATE(E_JXX(NumberofElements))
ALLOCATE(E_JYY(NumberofElements))
ALLOCATE(E_JZZ(NumberofElements))
ALLOCATE(E_Sub_EXX(NumberofElements,E_NumberofGaussPoints))
ALLOCATE(E_Sub_EYY(NumberofElements,E_NumberofGaussPoints))
ALLOCATE(E_Sub_EZZ(NumberofElements,E_NumberofGaussPoints))
ALLOCATE(E_Sub_JXX(NumberofElements,E_NumberofGaussPoints))
ALLOCATE(E_Sub_JYY(NumberofElements,E_NumberofGaussPoints))
ALLOCATE(E_Sub_JZZ(NumberofElements,E_NumberofGaussPoints))
ALLOCATE(E_Vol_Ele(NumberofElements))
CALL E_EandJ (E_Vol,E_Vol_Ele,E_sub_sigmasigma,Connectivities,Nodes,NodeIndex,E_Dis1,E_EXX,E_EYY,E_EZZ,E_JXX,E_JYY,E_JZZ,E_DNDX,E_DNDY,&
E_DNDZ,E_Sub_EXX,E_Sub_EYY,E_Sub_EZZ,E_Sub_JXX,E_Sub_JYY,E_Sub_JZZ)

CALL E_EffectiveProperties(E_ZeroTimeStep,k,ii,MacroEx_Knot,MacroEy_Knot,MacroEz_Knot,Nodes,E_Vol,E_Sub_EXX,E_Sub_EYY,E_Sub_EZZ,E_Sub_JXX,E_Sub_JYY,&
E_Sub_JZZ,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes,E_sub_sigmasigma,&
E_Sigma_XX_Eff,E_Sigma_YY_Eff,E_Sigma_ZZ_Eff,E_Sigma_XY_Eff,E_Sigma_XZ_Eff,E_Sigma_YZ_Eff,E_Sigma_ZZ2_Eff)

CALL E_FeedForTecplot (E_ZeroTimeStep,k,ii,Connectivities,Nodes,NodeIndex,CurrentNodes,E_Dis1,E_EXX,E_EYY,E_EZZ,E_JXX,E_JYY,E_JZZ,E_sub_sigmasigma,&
E_sub_rhorho)

END SUBROUTINE Macro_Electrostatic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RotatedStrains(Connectivities,BeginofModeIndex,EndofModeIndex,RotationMatrix,Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,&
Sub_StrainsXZ1,Sub_StrainsYZ1,Sub_StrainsXY2,Sub_StrainsXZ2,Sub_StrainsYZ2,Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,&
Sub_StrainsXZ1_R,Sub_StrainsYZ1_R,Sub_StrainsXY2_R,Sub_StrainsXZ2_R,Sub_StrainsYZ2_R,Vol)
IMPLICIT NONE

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER :: i,i1,j1,j2,j3,j4,j5,j6,j7,j8
INTEGER, INTENT(IN), DIMENSION(NumberofElements,17) :: Connectivities
INTEGER, INTENT(IN) :: BeginofModeIndex,EndofModeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofMaterialModes,3,3) :: RotationMatrix 
DOUBLE PRECISION, DIMENSION (NumberofElements,8) :: Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,Sub_StrainsXZ1,Sub_StrainsYZ1,&
Sub_StrainsXY2,Sub_StrainsXZ2,Sub_StrainsYZ2,Sub_StressesXX,Sub_StressesYY,Sub_StressesZZ,Sub_StressesXY,Sub_StressesXZ,Sub_StressesYZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,8) :: Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,&
Sub_StrainsXZ1_R,Sub_StrainsYZ1_R,Sub_StrainsXY2_R,Sub_StrainsXZ2_R,Sub_StrainsYZ2_R

DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,NumberofElements,3,3) :: StrainMatrix
DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,NumberofElements,3,3) :: StrainMatrix_Temp
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofGaussPoints,NumberofElements) :: Vol
DOUBLE PRECISION :: Vol_Avg_Sub_StrainsXX_R,Temp_Vol
INTEGER :: OUT2

OUT2=21

StrainMatrix_Temp=0.0D0

DO i1=1,NumberofGaussPoints
DO i=1,NumberofElements
StrainMatrix(i1,i,1,1)=Sub_StrainsXX(i,i1)
StrainMatrix(i1,i,1,2)=Sub_StrainsXY1(i,i1)
StrainMatrix(i1,i,1,3)=Sub_StrainsXZ1(i,i1)
StrainMatrix(i1,i,2,1)=Sub_StrainsXY2(i,i1)
StrainMatrix(i1,i,2,2)=Sub_StrainsYY(i,i1)
StrainMatrix(i1,i,2,3)=Sub_StrainsYZ1(i,i1)
StrainMatrix(i1,i,3,1)=Sub_StrainsXZ2(i,i1)
StrainMatrix(i1,i,3,2)=Sub_StrainsYZ2(i,i1)
StrainMatrix(i1,i,3,3)=Sub_StrainsZZ(i,i1)
ENDDO
ENDDO

DO i=1,NumberofElements

DO i1=1,NumberofGaussPoints

IF (Connectivities(i,i1+9)>=BeginofModeIndex .AND. Connectivities(i,i1+9)<=EndofModeIndex) THEN

DO j1=1,3
DO j2=1,3
DO j3=1,3
DO j4=1,3
StrainMatrix_Temp(i1,i,j1,j2)=StrainMatrix_Temp(i1,i,j1,j2)+RotationMatrix(Connectivities(i,i1+9)-BeginofModeIndex+1,j1,j3)*&
RotationMatrix(Connectivities(i,i1+9)-BeginofModeIndex+1,j2,j4)*StrainMatrix(i1,i,j3,j4)
ENDDO
ENDDO
ENDDO
ENDDO

ELSE

DO j1=1,3
DO j2=1,3
StrainMatrix_Temp(i1,i,j1,j2)=StrainMatrix(i1,i,j1,j2)
ENDDO
ENDDO

ENDIF

ENDDO

ENDDO !For i

DO i1=1,NumberofGaussPoints
DO i=1,NumberofElements
Sub_StrainsXX_R(i,i1)=StrainMatrix_Temp(i1,i,1,1)
!WRITE(*,*) Sub_StrainsXX_R(i,i1)
Sub_StrainsXY1_R(i,i1)=StrainMatrix_Temp(i1,i,1,2) !StrainMatrix_Temp(i1,i,1,2)
Sub_StrainsXZ1_R(i,i1)=StrainMatrix_Temp(i1,i,1,3) !StrainMatrix_Temp(i1,i,1,3)
Sub_StrainsXY2_R(i,i1)=StrainMatrix_Temp(i1,i,2,1) !StrainMatrix_Temp(i1,i,2,1)
Sub_StrainsYY_R(i,i1)=StrainMatrix_Temp(i1,i,2,2)
Sub_StrainsYZ1_R(i,i1)=StrainMatrix_Temp(i1,i,2,3)
Sub_StrainsXZ2_R(i,i1)=StrainMatrix_Temp(i1,i,3,1) !StrainMatrix_Temp(i1,i,3,1)
Sub_StrainsYZ2_R(i,i1)=StrainMatrix_Temp(i1,i,3,2)
Sub_StrainsZZ_R(i,i1)=StrainMatrix_Temp(i1,i,3,3)
ENDDO
ENDDO


Vol_Avg_Sub_StrainsXX_R=0.0D0
Temp_Vol=0.0D0
DO i1=1,NumberofGaussPoints
DO i=1,NumberofElements
Temp_Vol=Temp_Vol+Vol(i1,1)
Vol_Avg_Sub_StrainsXX_R=Vol_Avg_Sub_StrainsXX_R+Sub_StrainsXX_R(i,i1)*Vol(i1,i)
ENDDO
ENDDO
Vol_Avg_Sub_StrainsXX_R=Vol_Avg_Sub_StrainsXX_R/Temp_Vol

OPEN(OUT2,FILE="3D_Output_Micro_F/3D_Output/Vol_Avg_Strain.txt",STATUS="UNKNOWN")
WRITE(OUT2,*) Vol_Avg_Sub_StrainsXX_R
CLOSE(OUT2)

END SUBROUTINE RotatedStrains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE RotationMatrixSUB(Anisotropic,Nodes,Connectivities,BeginofModeIndex,EndofModeIndex,XYofAxis,MaximumRadius,NN,RotationMatrix,RotationMatrix2)
!IMPLICIT NONE
!INTEGER :: i,j,i1,j1,j2,j3,j4,j5,j6,j7,j8
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!The macroscale common variables
!INTEGER :: TestType
!INTEGER :: DisplacementEquationFlag
!INTEGER :: NumberofNodes
!INTEGER :: NumberofElements
!INTEGER :: NumberofBoundaryNodes
!INTEGER :: NumberofMaterialModes
!INTEGER :: NumberofGaussPoints
!INTEGER :: TotalTimeSteps
!INTEGER :: NumberofPositiveXNodes
!INTEGER :: NumberofPositiveYNodes
!INTEGER :: NumberofPositiveZNodes
!INTEGER :: NumberofNegativeXNodes
!INTEGER :: NumberofNegativeYNodes
!INTEGER :: NumberofNegativeZNodes
!INTEGER :: MaxDimenofInverseConnec
!INTEGER :: NumberofRVETypes
!INTEGER :: ElementShown
!DOUBLE PRECISION :: Mag
!COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
!NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
!NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown
!
!INTEGER, INTENT(IN) :: Anisotropic
!DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: Nodes
!INTEGER, INTENT(IN), DIMENSION (NumberofElements,10) :: Connectivities
!INTEGER, INTENT(IN) :: BeginofModeIndex,EndofModeIndex
!DOUBLE PRECISION, INTENT(IN), DIMENSION(11,2) :: XYofAxis
!DOUBLE PRECISION, INTENT(IN) :: MaximumRadius
!
!DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,EndofModeIndex-BeginofModeIndex+1,3) :: Centroid
!
!DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofGaussPoints,8) :: NN
!
!DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,EndofModeIndex-BeginofModeIndex+1) :: COSTHETA
!DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,EndofModeIndex-BeginofModeIndex+1) :: SINTHETA
!
!DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,7) :: r
!
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofGaussPoints,NumberofMaterialModes,3,3) :: RotationMatrix
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofGaussPoints,NumberofMaterialModes,3,3) :: RotationMatrix2
!
!RotationMatrix=0.0D0
!
!DO i=1,NumberofElements
!
!IF ( Connectivities(i,10) >= BeginofModeIndex .AND. Connectivities(i,10) <= EndofModeIndex) THEN
!
!DO j=1,NumberofGaussPoints
!
!Centroid(j,Connectivities(i,10)-BeginofModeIndex+1,1)=Nodes(Connectivities(i,2),2)*NN(j,1)+Nodes(Connectivities(i,3),2)*NN(j,2)+&
!Nodes(Connectivities(i,4),2)*NN(j,3)+Nodes(Connectivities(i,5),2)*NN(j,4)+Nodes(Connectivities(i,6),2)*NN(j,5)&
!+Nodes(Connectivities(i,7),2)*NN(j,6)+Nodes(Connectivities(i,8),2)*NN(j,7)+Nodes(Connectivities(i,9),2)*NN(j,8)
!
!Centroid(j,Connectivities(i,10)-BeginofModeIndex+1,2)=Nodes(Connectivities(i,2),3)*NN(j,1)+Nodes(Connectivities(i,3),3)*NN(j,2)+&
!Nodes(Connectivities(i,4),3)*NN(j,3)+Nodes(Connectivities(i,5),3)*NN(j,4)+Nodes(Connectivities(i,6),3)*NN(j,5)&
!+Nodes(Connectivities(i,7),3)*NN(j,6)+Nodes(Connectivities(i,8),3)*NN(j,7)+Nodes(Connectivities(i,9),3)*NN(j,8)
!
!Centroid(j,Connectivities(i,10)-BeginofModeIndex+1,3)=Nodes(Connectivities(i,2),4)*NN(j,1)+Nodes(Connectivities(i,3),4)*NN(j,2)+&
!Nodes(Connectivities(i,4),4)*NN(j,3)+Nodes(Connectivities(i,5),4)*NN(j,4)+Nodes(Connectivities(i,6),4)*NN(j,5)&
!+Nodes(Connectivities(i,7),4)*NN(j,6)+Nodes(Connectivities(i,8),4)*NN(j,7)+Nodes(Connectivities(i,9),4)*NN(j,8)
!
!ENDDO
!
!ENDIF
!
!ENDDO
!
!IF (Anisotropic==1) THEN
!DO i=1,EndofModeIndex-BeginofModeIndex+1
!r=0.0D0
!DO i1=1,NumberofGaussPoints
!DO j=1,7
!r(i1,j)=SQRT((Centroid(i1,i,1)-XYofAxis(j,1))**2+(Centroid(i1,i,2)-XYofAxis(j,2))**2)
!IF (r(i1,j) < MaximumRadius) THEN
!COSTHETA(i1,i)= (Centroid(i1,i,1)-XYofAxis(j,1)) /SQRT( (Centroid(i1,i,1)-XYofAxis(j,1))**2+ (Centroid(i1,i,2)-XYofAxis(j,2))**2)
!SINTHETA(i1,i)= (Centroid(i1,i,2)-XYofAxis(j,2)) /SQRT( (Centroid(i1,i,1)-XYofAxis(j,1))**2+ (Centroid(i1,i,2)-XYofAxis(j,2))**2)
!GOTO 1
!ELSE
!ENDIF
!1 ENDDO
!ENDDO
!ENDDO
!ENDIF
!
!!Ori_CC_Temp=0.0D0
!!DO i=1,EndofModeIndex-BeginofModeIndex+1
!!
!!DO j1=1,6
!!DO j2=1,6
!!Ori_CC(i,j1,j2)=CC(i+BeginofModeIndex-1,j1,j2)
!!ENDDO
!!ENDDO
!!
!!Ori_CC_Temp(i,1,1,1,1)=Ori_CC(i,1,1)
!!Ori_CC_Temp(i,1,1,2,2)=Ori_CC(i,1,2)
!!Ori_CC_Temp(i,1,1,3,3)=Ori_CC(i,1,3)
!!Ori_CC_Temp(i,2,2,1,1)=Ori_CC(i,2,1)
!!Ori_CC_Temp(i,2,2,2,2)=Ori_CC(i,2,2)
!!Ori_CC_Temp(i,2,2,3,3)=Ori_CC(i,2,3)
!!Ori_CC_Temp(i,3,3,1,1)=Ori_CC(i,3,1)
!!Ori_CC_Temp(i,3,3,2,2)=Ori_CC(i,3,2)
!!Ori_CC_Temp(i,3,3,3,3)=Ori_CC(i,3,3)
!!Ori_CC_Temp(i,2,3,2,3)=Ori_CC(i,4,4)
!!Ori_CC_Temp(i,1,3,1,3)=Ori_CC(i,5,5)
!!Ori_CC_Temp(i,1,2,1,2)=Ori_CC(i,6,6)
!!
!!ENDDO
!
!!CC_Temp=0.0D0
!DO i1=1,NumberofGaussPoints
!DO i=1,EndofModeIndex-BeginofModeIndex+1
!
!IF (Anisotropic==1) THEN
!!From Cartesian to Cylindrical
!RotationMatrix(i1,i,1,1)=COSTHETA(i1,i)
!RotationMatrix(i1,i,1,2)=SINTHETA(i1,i)
!RotationMatrix(i1,i,2,1)=-SINTHETA(i1,i)
!RotationMatrix(i1,i,2,2)=COSTHETA(i1,i)
!RotationMatrix(i1,i,3,3)=1.0D0
!
!RotationMatrix2(i1,i,1,1)=COSTHETA(i1,i)
!RotationMatrix2(i1,i,1,2)=-SINTHETA(i1,i)
!RotationMatrix2(i1,i,2,1)=SINTHETA(i1,i)
!RotationMatrix2(i1,i,2,2)=COSTHETA(i1,i)
!RotationMatrix2(i1,i,3,3)=1.0D0
!ELSEIF (Anisotropic==2) THEN
!
!ENDIF
!
!!DO j1=1,3
!!DO j2=1,3
!!DO j3=1,3
!!DO j4=1,3
!!DO j5=1,3
!!DO j6=1,3
!!DO j7=1,3
!!DO j8=1,3
!!CC_Temp(i1,i,j1,j2,j3,j4)=CC_Temp(i1,i,j1,j2,j3,j4)+RotationMatrix(i1,i,j1,j5)*RotationMatrix(i1,i,j2,j6)*&
!!RotationMatrix(i1,i,j3,j7)*RotationMatrix(i1,i,j4,j8)*Ori_CC_Temp(i,j5,j6,j7,j8)
!!ENDDO
!!ENDDO
!!ENDDO
!!ENDDO
!!ENDDO
!!ENDDO
!!ENDDO
!!ENDDO
!
!ENDDO
!ENDDO
!
!END SUBROUTINE RotationMatrixSUB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE StrainsandStresses (Vol,Vol_Ele,CC,Sub_CC,Connectivities,Nodes,NodeIndex,DNDX,DNDY,DNDZ,KIC,KII,&
Dis1,StrainsXX,StrainsYY,StrainsZZ,StrainsXY,StrainsXZ,StrainsYZ,StrainsVonMises,StressesXX,StressesYY,StressesZZ,StressesXY,StressesXZ,StressesYZ,&
StressesVonMises,Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,Sub_StrainsXZ1,Sub_StrainsYZ1,&
Sub_StrainsXY2,Sub_StrainsXZ2,Sub_StrainsYZ2,Sub_StressesXX,Sub_StressesYY,Sub_StressesZZ,Sub_StressesXY,Sub_StressesXZ,Sub_StressesYZ,CurrentNodes)
IMPLICIT NONE
INTEGER :: i,j,k,i1,i2,i3,i4

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: Nodes
INTEGER, INTENT(IN), DIMENSION(NumberofNodes) :: NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*NumberofNodes) :: Dis1
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofGaussPoints,NumberofElements) :: Vol
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: Vol_Ele
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofGaussPoints,NumberofElements,8+3) :: DNDX,DNDY,DNDZ
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements,9,24)  :: KIC
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements,9,9)   :: KII
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofMaterialModes,6,6) :: CC
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofGaussPoints,NumberofElements,6,6) :: Sub_CC
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StrainsXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StrainsYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StrainsZZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StrainsXY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StrainsXZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StrainsYZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StrainsVonMises
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StressesXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StressesYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StressesZZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StressesXY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StressesXZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StressesYZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements) :: StressesVonMises
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofNodes,4)  :: CurrentNodes 
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,8) :: Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,Sub_StrainsXZ1,Sub_StrainsYZ1,&
Sub_StrainsXY2,Sub_StrainsXZ2,Sub_StrainsYZ2,Sub_StressesXX,Sub_StressesYY,Sub_StressesZZ,Sub_StressesXY,Sub_StressesXZ,Sub_StressesYZ
DOUBLE PRECISION, DIMENSION(NumberofElements,24) :: Element_U
DOUBLE PRECISION, DIMENSION(NumberofElements,9) :: Element_Alpha
DOUBLE PRECISION, DIMENSION(NumberofElements,9) :: KICU
DOUBLE PRECISION, DIMENSION(NumberofElements,9,9) :: KII_Inverse
!INTEGER :: OUT1,OUT2,OUT3,OUT4

Element_U=0.0D0
Element_Alpha=0.0D0
KICU=0.0D0
KII_Inverse=0.0D0

DO i=1,NumberofElements

DO j=1,8
!Element_U(i,3*(j-1)+1)=Dis1(Connectivities(i,j+1))
!Element_U(i,3*(j-1)+2)=Dis1(NumberofNodes+Connectivities(i,j+1))
!Element_U(i,3*(j-1)+3)=Dis1(2*NumberofNodes+Connectivities(i,j+1))
Element_U(i,j)=Dis1(Connectivities(i,j+1))
Element_U(i,8+j)=Dis1(NumberofNodes+Connectivities(i,j+1))
Element_U(i,16+j)=Dis1(2*NumberofNodes+Connectivities(i,j+1))

ENDDO

DO i1=1,9
DO i2=1,24
KICU(i,i1)=KICU(i,i1)+KIC(i,i1,i2)*Element_U(i,i2)
ENDDO
ENDDO

CALL INV(i,9,KII(i,:,:),KII_Inverse(i,:,:))

DO i1=1,9
DO i2=1,9
Element_Alpha(i,i1)=Element_Alpha(i,i1)-KII_Inverse(i,i1,i2)*KICU(i,i2)
ENDDO
ENDDO

ENDDO


Sub_StrainsXX=0.0D0
Sub_StrainsYY=0.0D0
Sub_StrainsZZ=0.0D0
Sub_StrainsXY1=0.0D0
Sub_StrainsXZ1=0.0D0
Sub_StrainsYZ1=0.0D0
Sub_StrainsXY2=0.0D0
Sub_StrainsXZ2=0.0D0
Sub_StrainsYZ2=0.0D0
StrainsXX=0.0D0
StrainsYY=0.0D0
StrainsZZ=0.0D0
StrainsXY=0.0D0
StrainsXZ=0.0D0
StrainsYZ=0.0D0
StrainsVonMises=0.0D0

Sub_StressesXX=0.0D0
Sub_StressesYY=0.0D0
Sub_StressesZZ=0.0D0
Sub_StressesXY=0.0D0
Sub_StressesXZ=0.0D0
Sub_StressesYZ=0.0D0
StressesXX=0.0D0
StressesYY=0.0D0
StressesZZ=0.0D0
StressesXY=0.0D0
StressesXZ=0.0D0
StressesYZ=0.0D0
StressesVonMises=0.0D0

Vol_Ele=0.0D0


DO k=1,NumberofElements
DO j=1,NumberofGaussPoints

DO i=1,8
Sub_StrainsXX(k,j)=Sub_StrainsXX(k,j)+Dis1(Connectivities(k,i+1))*DNDX(j,k,i)
Sub_StrainsYY(k,j)=Sub_StrainsYY(k,j)+Dis1(NumberofNodes+Connectivities(k,i+1))*DNDY(j,k,i)
Sub_StrainsZZ(k,j)=Sub_StrainsZZ(k,j)+Dis1(2*NumberofNodes+Connectivities(k,i+1))*DNDZ(j,k,i)
Sub_StrainsXY1(k,j)=Sub_StrainsXY1(k,j)+0.5D0*Dis1(Connectivities(k,i+1))*DNDY(j,k,i)
Sub_StrainsXZ1(k,j)=Sub_StrainsXZ1(k,j)+0.5D0*Dis1(Connectivities(k,i+1))*DNDZ(j,k,i)
Sub_StrainsYZ1(k,j)=Sub_StrainsYZ1(k,j)+0.5D0*Dis1(NumberofNodes+Connectivities(k,i+1))*DNDZ(j,k,i)
Sub_StrainsXY2(k,j)=Sub_StrainsXY2(k,j)+0.5D0*Dis1(NumberofNodes+Connectivities(k,i+1))*DNDX(j,k,i)
Sub_StrainsXZ2(k,j)=Sub_StrainsXZ2(k,j)+0.5D0*Dis1(2*NumberofNodes+Connectivities(k,i+1))*DNDX(j,k,i)
Sub_StrainsYZ2(k,j)=Sub_StrainsYZ2(k,j)+0.5D0*Dis1(2*NumberofNodes+Connectivities(k,i+1))*DNDY(j,k,i)
ENDDO

DO i=1,3
!Sub_StrainsXX(k,j)=Sub_StrainsXX(k,j)+Element_Alpha(k,3*(i-1)+1)*DNDX(j,k,i+8)
!Sub_StrainsYY(k,j)=Sub_StrainsYY(k,j)+Element_Alpha(k,3*(i-1)+2)*DNDY(j,k,i+8)
!Sub_StrainsZZ(k,j)=Sub_StrainsZZ(k,j)+Element_Alpha(k,3*(i-1)+3)*DNDZ(j,k,i+8)
!Sub_StrainsXY1(k,j)=Sub_StrainsXY1(k,j)+0.5D0*Element_Alpha(k,3*(i-1)+1)*DNDY(j,k,i+8)
!Sub_StrainsXZ1(k,j)=Sub_StrainsXZ1(k,j)+0.5D0*Element_Alpha(k,3*(i-1)+1)*DNDZ(j,k,i+8)
!Sub_StrainsYZ1(k,j)=Sub_StrainsYZ1(k,j)+0.5D0*Element_Alpha(k,3*(i-1)+2)*DNDZ(j,k,i+8)
!Sub_StrainsXY2(k,j)=Sub_StrainsXY2(k,j)+0.5D0*Element_Alpha(k,3*(i-1)+2)*DNDX(j,k,i+8)
!Sub_StrainsXZ2(k,j)=Sub_StrainsXZ2(k,j)+0.5D0*Element_Alpha(k,3*(i-1)+3)*DNDX(j,k,i+8)
!Sub_StrainsYZ2(k,j)=Sub_StrainsYZ2(k,j)+0.5D0*Element_Alpha(k,3*(i-1)+3)*DNDY(j,k,i+8)

!Sub_StrainsXX(k,j)=Sub_StrainsXX(k,j)+Element_Alpha(k,i)*DNDX(j,k,i+8)
!Sub_StrainsYY(k,j)=Sub_StrainsYY(k,j)+Element_Alpha(k,3+i)*DNDY(j,k,i+8)
!Sub_StrainsZZ(k,j)=Sub_StrainsZZ(k,j)+Element_Alpha(k,6+i)*DNDZ(j,k,i+8)
!Sub_StrainsXY1(k,j)=Sub_StrainsXY1(k,j)+0.5D0*Element_Alpha(k,i)*DNDY(j,k,i+8)
!Sub_StrainsXZ1(k,j)=Sub_StrainsXZ1(k,j)+0.5D0*Element_Alpha(k,i)*DNDZ(j,k,i+8)
!Sub_StrainsYZ1(k,j)=Sub_StrainsYZ1(k,j)+0.5D0*Element_Alpha(k,3+i)*DNDZ(j,k,i+8)
!Sub_StrainsXY2(k,j)=Sub_StrainsXY2(k,j)+0.5D0*Element_Alpha(k,3+i)*DNDX(j,k,i+8)
!Sub_StrainsXZ2(k,j)=Sub_StrainsXZ2(k,j)+0.5D0*Element_Alpha(k,6+i)*DNDX(j,k,i+8)
!Sub_StrainsYZ2(k,j)=Sub_StrainsYZ2(k,j)+0.5D0*Element_Alpha(k,6+i)*DNDY(j,k,i+8)

ENDDO

Sub_StressesXX(k,j)=Sub_CC(j,k,1,1)*Sub_StrainsXX(k,j)+Sub_CC(j,k,1,2)*Sub_StrainsYY(k,j)+&
Sub_CC(j,k,1,3)*Sub_StrainsZZ(k,j)+2.0D0*Sub_CC(j,k,1,4)*(Sub_StrainsYZ1(k,j)+Sub_StrainsYZ2(k,j))+&
2.0D0*Sub_CC(j,k,1,5)*(Sub_StrainsXZ1(k,j)+Sub_StrainsXZ2(k,j))+&
2.0D0*Sub_CC(j,k,1,6)*(Sub_StrainsXY1(k,j)+Sub_StrainsXY2(k,j))
Sub_StressesYY(k,j)=Sub_CC(j,k,1,2)*Sub_StrainsXX(k,j)+Sub_CC(j,k,2,2)*Sub_StrainsYY(k,j)+&
Sub_CC(j,k,2,3)*Sub_StrainsZZ(k,j)+2.0D0*Sub_CC(j,k,2,4)*(Sub_StrainsYZ1(k,j)+Sub_StrainsYZ2(k,j))+&
2.0D0*Sub_CC(j,k,2,5)*(Sub_StrainsXZ1(k,j)+Sub_StrainsXZ2(k,j))+&
2.0D0*Sub_CC(j,k,2,6)*(Sub_StrainsXY1(k,j)+Sub_StrainsXY2(k,j))
Sub_StressesZZ(k,j)=Sub_CC(j,k,1,3)*Sub_StrainsXX(k,j)+Sub_CC(j,k,2,3)*Sub_StrainsYY(k,j)+&
Sub_CC(j,k,3,3)*Sub_StrainsZZ(k,j)+2.0D0*Sub_CC(j,k,3,4)*(Sub_StrainsYZ1(k,j)+Sub_StrainsYZ2(k,j))+&
2.0D0*Sub_CC(j,k,3,5)*(Sub_StrainsXZ1(k,j)+Sub_StrainsXZ2(k,j))+&
2.0D0*Sub_CC(j,k,3,6)*(Sub_StrainsXY1(k,j)+Sub_StrainsXY2(k,j))
Sub_StressesXY(k,j)=Sub_CC(j,k,6,1)*Sub_StrainsXX(k,j)+Sub_CC(j,k,6,2)*Sub_StrainsYY(k,j)+&
Sub_CC(j,k,6,3)*Sub_StrainsZZ(k,j)+2.0D0*Sub_CC(j,k,6,4)*(Sub_StrainsYZ1(k,j)+Sub_StrainsYZ2(k,j))+&
2.0D0*Sub_CC(j,k,6,5)*(Sub_StrainsXZ1(k,j)+Sub_StrainsXZ2(k,j))+&
2.0D0*Sub_CC(j,k,6,6)*(Sub_StrainsXY1(k,j)+Sub_StrainsXY2(k,j))
Sub_StressesXZ(k,j)=Sub_CC(j,k,5,1)*Sub_StrainsXX(k,j)+Sub_CC(j,k,5,2)*Sub_StrainsYY(k,j)+&
Sub_CC(j,k,5,3)*Sub_StrainsZZ(k,j)+2.0D0*Sub_CC(j,k,5,4)*(Sub_StrainsYZ1(k,j)+Sub_StrainsYZ2(k,j))+&
2.0D0*Sub_CC(j,k,5,5)*(Sub_StrainsXZ1(k,j)+Sub_StrainsXZ2(k,j))+&
2.0D0*Sub_CC(j,k,5,6)*(Sub_StrainsXY1(k,j)+Sub_StrainsXY2(k,j))
Sub_StressesYZ(k,j)=Sub_CC(j,k,4,1)*Sub_StrainsXX(k,j)+Sub_CC(j,k,4,2)*Sub_StrainsYY(k,j)+&
Sub_CC(j,k,4,3)*Sub_StrainsZZ(k,j)+2.0D0*Sub_CC(j,k,4,4)*(Sub_StrainsYZ1(k,j)+Sub_StrainsYZ2(k,j))+&
2.0D0*Sub_CC(j,k,4,5)*(Sub_StrainsXZ1(k,j)+Sub_StrainsXZ2(k,j))+&
2.0D0*Sub_CC(j,k,4,6)*(Sub_StrainsXY1(k,j)+Sub_StrainsXY2(k,j))

ENDDO
ENDDO

DO k=1,NumberofElements
DO j=1,NumberofGaussPoints

StrainsXX(k)=StrainsXX(k)+Sub_StrainsXX(k,j)*Vol(j,k)
StrainsYY(k)=StrainsYY(k)+Sub_StrainsYY(k,j)*Vol(j,k)
StrainsZZ(k)=StrainsZZ(k)+Sub_StrainsZZ(k,j)*Vol(j,k)
StrainsXY(k)=StrainsXY(k)+Sub_StrainsXY1(k,j)*Vol(j,k)+Sub_StrainsXY2(k,j)*Vol(j,k)
StrainsXZ(k)=StrainsXZ(k)+Sub_StrainsXZ1(k,j)*Vol(j,k)+Sub_StrainsXZ2(k,j)*Vol(j,k)
StrainsYZ(k)=StrainsYZ(k)+Sub_StrainsYZ1(k,j)*Vol(j,k)+Sub_StrainsYZ2(k,j)*Vol(j,k)

StrainsVonMises(k)=StrainsVonMises(k)+SQRT(0.5D0*( (Sub_StrainsXX(k,j)-Sub_StrainsYY(k,j))**2+&
(Sub_StrainsYY(k,j)-Sub_StrainsZZ(k,j))**2+(Sub_StrainsXX(k,j)-Sub_StrainsZZ(k,j))**2+&
6*( (Sub_StrainsYZ1(k,j)+Sub_StrainsYZ2(k,j))**2+(Sub_StrainsXZ1(k,j)+Sub_StrainsXZ2(k,j))**2+&
(Sub_StrainsXY1(k,j)+Sub_StrainsXY2(k,j))**2) ) )*Vol(j,k)

StressesXX(k)=StressesXX(k)+Sub_StressesXX(k,j)*Vol(j,k)
StressesYY(k)=StressesYY(k)+Sub_StressesYY(k,j)*Vol(j,k)
StressesZZ(k)=StressesZZ(k)+Sub_StressesZZ(k,j)*Vol(j,k)
StressesXY(k)=StressesXY(k)+Sub_StressesXY(k,j)*Vol(j,k)
StressesXZ(k)=StressesXZ(k)+Sub_StressesXZ(k,j)*Vol(j,k)
StressesYZ(k)=StressesYZ(k)+Sub_StressesYZ(k,j)*Vol(j,k)
Vol_Ele(k)=Vol_Ele(k)+Vol(j,k)

StressesVonMises(k)=StressesVonMises(k)+SQRT(0.5D0*( (Sub_StressesXX(k,j)-Sub_StressesYY(k,j))**2+&
(Sub_StressesYY(k,j)-Sub_StressesZZ(k,j))**2+(Sub_StressesXX(k,j)-Sub_StressesZZ(k,j))**2+&
6*(Sub_StressesYZ(k,j)**2+Sub_StressesXZ(k,j)**2+Sub_StressesXY(k,j)**2) ) )*Vol(j,k)

ENDDO
StrainsXX(k)=StrainsXX(k)/Vol_Ele(k)
StrainsYY(k)=StrainsYY(k)/Vol_Ele(k)
StrainsZZ(k)=StrainsZZ(k)/Vol_Ele(k)
StrainsXY(k)=StrainsXY(k)/Vol_Ele(k)
StrainsXZ(k)=StrainsXZ(k)/Vol_Ele(k)
StrainsYZ(k)=StrainsYZ(k)/Vol_Ele(k)

StrainsVonMises(k)=StrainsVonMises(k)/Vol_Ele(k)

StressesXX(k)=StressesXX(k)/Vol_Ele(k)
StressesYY(k)=StressesYY(k)/Vol_Ele(k)
StressesZZ(k)=StressesZZ(k)/Vol_Ele(k)
StressesXY(k)=StressesXY(k)/Vol_Ele(k)
StressesXZ(k)=StressesXZ(k)/Vol_Ele(k)
StressesYZ(k)=StressesYZ(k)/Vol_Ele(k)

StressesVonMises(k)=StressesVonMises(k)/Vol_Ele(k)

ENDDO


DO i=1,NumberofNodes

CurrentNodes(i,1)=Nodes(i,1)
CurrentNodes(i,2)=Nodes(i,2)+Mag*Dis1(i)
CurrentNodes(i,3)=Nodes(i,3)+Mag*Dis1(NumberofNodes+i)
CurrentNodes(i,4)=Nodes(i,4)+Mag*Dis1(2*NumberofNodes+i)

ENDDO

!DO i=1,NumberofElements
!
!WRITE(*,*) StrainsXX(i),StrainsYY(i),StrainsZZ(i),StrainsXY(i),StrainsXZ(i),StrainsYZ(i)
!
!ENDDO

!OUT1=5
!OUT2=6
!OUT3=7
!OUT4=8
!
!OPEN (OUT1,file = '3D_Output_Micro_F/Sub_Strains_11_Microscale.txt',status = 'unknown')
!OPEN (OUT2,file = '3D_Output_Micro_F/Sub_Strains_22_Microscale.txt',status = 'unknown')
!OPEN (OUT3,file = '3D_Output_Micro_F/Sub_Strains_12_Microscale.txt',status = 'unknown')
!OPEN (OUT4,file = '3D_Output_Micro_F/Sub_Strains_33_Microscale.txt',status = 'unknown')
!DO k=1,NumberofElements
!WRITE(OUT1,10) Sub_StrainsXX(k,1),Sub_StrainsXX(k,2),Sub_StrainsXX(k,3),Sub_StrainsXX(k,4),&
!Sub_StrainsXX(k,5),Sub_StrainsXX(k,6),Sub_StrainsXX(k,7),Sub_StrainsXX(k,8)
!WRITE(OUT2,10) Sub_StrainsYY(k,1),Sub_StrainsYY(k,2),Sub_StrainsYY(k,3),Sub_StrainsYY(k,4),&
!Sub_StrainsYY(k,5),Sub_StrainsYY(k,6),Sub_StrainsYY(k,7),Sub_StrainsYY(k,8)
!WRITE(OUT3,10) Sub_StrainsXY1(k,1)+Sub_StrainsXY2(k,1),Sub_StrainsXY1(k,2)+Sub_StrainsXY2(k,2),&
!Sub_StrainsXY1(k,3)+Sub_StrainsXY2(k,3),Sub_StrainsXY1(k,4)+Sub_StrainsXY2(k,4),&
!Sub_StrainsXY1(k,5)+Sub_StrainsXY2(k,5),Sub_StrainsXY1(k,6)+Sub_StrainsXY2(k,6),&
!Sub_StrainsXY1(k,7)+Sub_StrainsXY2(k,7),Sub_StrainsXY1(k,8)+Sub_StrainsXY2(k,8)
!WRITE(OUT4,10) Sub_StrainsZZ(k,1),Sub_StrainsZZ(k,2),Sub_StrainsZZ(k,3),Sub_StrainsZZ(k,4),&
!Sub_StrainsZZ(k,5),Sub_StrainsZZ(k,6),Sub_StrainsZZ(k,7),Sub_StrainsZZ(k,8)
!ENDDO
!CLOSE(OUT1)
!CLOSE(OUT2)
!CLOSE(OUT3)
!CLOSE(OUT4)
!
!10 FORMAT(8(E30.15,1X),/)

END SUBROUTINE StrainsandStresses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Matrix_DE_Reduced(k,NumberofDE,One_Step_DEs,BoundaryIndex,Boundaries,K1GRR,K1GRCDD,K1GRCII,F1Globall,Connectivities,&
Dis1,Highest_GRCD)
IMPLICIT NONE
INTEGER :: i,j,m,n,i1,j1,k1,k2,k3
INTEGER, INTENT(IN) :: k

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(IN) :: NumberofDE
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec) :: K1GRR
INTEGER, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE) :: K1GRCDD
INTEGER, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec) :: K1GRCII
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE) :: F1Globall
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: K1GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: K1GRCD
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: K1GRCI
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: F1Global

INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: K2GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: K2GRRD

DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofDE,8) :: One_Step_DEs
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofBoundaryNodes,TotalTimeSteps+2) :: Boundaries
INTEGER, INTENT(IN), DIMENSION(NumberofBoundaryNodes,2) :: BoundaryIndex
DOUBLE PRECISION, DIMENSION (3*NumberofNodes+NumberofDE) :: Dis
DOUBLE PRECISION, INTENT(OUT),DIMENSION (3*NumberofNodes+NumberofDE) :: Dis1
INTEGER :: Dimen1
DOUBLE PRECISION, DIMENSION(3*NumberofNodes+NumberofDE) :: x
INTEGER :: OUT1=40,OUT2=41
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
!DOUBLE PRECISION :: ElapsedTime,Temp
INTEGER :: TEMP_GRCD
INTEGER, INTENT(OUT) :: Highest_GRCD
INTEGER :: Temp_INT
DOUBLE PRECISION :: Temp
INTEGER :: DUMMY,DUMMY2
INTEGER, ALLOCATABLE, DIMENSION(:) :: IA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: IA
INTEGER, ALLOCATABLE, DIMENSION(:) :: JA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: JA2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: A1

ALLOCATE(K1GR(3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec))
ALLOCATE(K1GRCD(3*NumberofNodes+NumberofDE))
ALLOCATE(K1GRCI(3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec))
ALLOCATE(F1Global(3*NumberofNodes+NumberofDE))

K1GR=K1GRR
K1GRCD=K1GRCDD
K1GRCI=K1GRCII
F1Global=F1Globall

!!!Apply P.B.C.s (Top Right and Bottom Right of F)

DO i=1,NumberofDE

DO j=1,int(One_Step_DEs(i,1))

IF(int(One_Step_DEs(i,3*j))==1) THEN

K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+1 ) = K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+1 )+1

K1GR( 3*(int(One_Step_DEs(i,3*j-1))-1)+1,K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+1 ) ) = One_Step_DEs(i,3*j+1) 

K1GRCI( 3*(int(One_Step_DEs(i,3*j-1))-1)+1, K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+1 ) ) = 3*NumberofNodes+i

ELSEIF(int(One_Step_DEs(i,3*j))==2) THEN

K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+2 ) = K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+2 )+1

K1GR( 3*(int(One_Step_DEs(i,3*j-1))-1)+2,K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+2 ) ) = One_Step_DEs(i,3*j+1) 

K1GRCI( 3*(int(One_Step_DEs(i,3*j-1))-1)+2, K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+2 ) ) = 3*NumberofNodes+i


ELSEIF(int(One_Step_DEs(i,3*j))==3) THEN

K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+3 ) = K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+3 )+1

K1GR( 3*(int(One_Step_DEs(i,3*j-1))-1)+3,K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+3 ) ) = One_Step_DEs(i,3*j+1) 

K1GRCI( 3*(int(One_Step_DEs(i,3*j-1))-1)+3, K1GRCD( 3*(int(One_Step_DEs(i,3*j-1))-1)+3 ) ) = 3*NumberofNodes+i

ENDIF

F1Global(3*NumberofNodes+i)=One_Step_DEs(i,3*int(One_Step_DEs(i,1))+1+1)

ENDDO
ENDDO

!!!Apply P.B.C.s (Bottom left)
DO i=1,NumberofDE
DO j=1,int(One_Step_DEs(i,1))
K1GRCD( 3*NumberofNodes+i ) = K1GRCD( 3*NumberofNodes+i )+1
IF(int(One_Step_DEs(i,3*j))==1) THEN
K1GR( 3*NumberofNodes+i,K1GRCD( 3*NumberofNodes+i ) ) = One_Step_DEs(i,3*j+1) 
K1GRCI( 3*NumberofNodes+i, K1GRCD( 3*NumberofNodes+i ) ) = 3*(int(One_Step_DEs(i,3*j-1))-1)+1
ELSEIF(int(One_Step_DEs(i,3*j))==2) THEN
K1GR( 3*NumberofNodes+i,K1GRCD( 3*NumberofNodes+i ) ) = One_Step_DEs(i,3*j+1) 
K1GRCI( 3*NumberofNodes+i, K1GRCD( 3*NumberofNodes+i ) ) = 3*(int(One_Step_DEs(i,3*j-1))-1)+2
ELSEIF(int(One_Step_DEs(i,3*j))==3) THEN
K1GR( 3*NumberofNodes+i,K1GRCD( 3*NumberofNodes+i ) ) = One_Step_DEs(i,3*j+1) 
K1GRCI( 3*NumberofNodes+i, K1GRCD( 3*NumberofNodes+i ) ) = 3*(int(One_Step_DEs(i,3*j-1))-1)+3
ENDIF
ENDDO
ENDDO

ALLOCATE(K2GR(3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec,2))
ALLOCATE(K2GRRD(3*NumberofNodes+NumberofDE))
K2GR=0
K2GRRD=0

DO i=1,3*NumberofNodes+NumberofDE
DO j=1,K1GRCD(i)
K2GRRD(K1GRCI(i,j))=K2GRRD(K1GRCI(i,j))+1
K2GR(K1GRCI(i,j),K2GRRD(K1GRCI(i,j)),1)=i
K2GR(K1GRCI(i,j),K2GRRD(K1GRCI(i,j)),2)=j
ENDDO
ENDDO

Dis=0.0D0

DO i=1,NumberofBoundaryNodes

IF (BoundaryIndex(i,2)==1) THEN

Dis(3*(BoundaryIndex(i,1)-1)+1)=Boundaries(i,k+2)

ELSEIF (BoundaryIndex(i,2)==2) THEN

Dis(3*(BoundaryIndex(i,1)-1)+2)=Boundaries(i,k+2)

ELSEIF (BoundaryIndex(i,2)==3) THEN

Dis(3*(BoundaryIndex(i,1)-1)+3)=Boundaries(i,k+2)

ENDIF

ENDDO


DO i=1,NumberofBoundaryNodes

IF (BoundaryIndex(i,2)==1) THEN

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+1)
F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+1,m,1) )=F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+1,m,1) )&
-K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+1,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+1,m,2) )*Dis(3*(BoundaryIndex(i,1)-1)+1)

ENDDO

F1Global(3*(BoundaryIndex(i,1)-1)+1)=Dis(3*(BoundaryIndex(i,1)-1)+1)

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+1)
K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+1,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+1,m,2) )=0.0D0
ENDDO

DO m=1,K1GRCD(3*(BoundaryIndex(i,1)-1)+1)
IF ( 3*(BoundaryIndex(i,1)-1)+1 /= K1GRCI(3*(BoundaryIndex(i,1)-1)+1,m) ) THEN
K1GR( 3*(BoundaryIndex(i,1)-1)+1, m )=0.0D0
ELSE
K1GR( 3*(BoundaryIndex(i,1)-1)+1, m )=1.0D0
ENDIF
ENDDO


ELSEIF (BoundaryIndex(i,2)==2) THEN

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+2)
F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+2,m,1) )=F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+2,m,1) )&
-K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+2,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+2,m,2) )*Dis(3*(BoundaryIndex(i,1)-1)+2)

ENDDO

F1Global(3*(BoundaryIndex(i,1)-1)+2)=Dis(3*(BoundaryIndex(i,1)-1)+2)

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+2)
K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+2,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+2,m,2) )=0.0D0
ENDDO

DO m=1,K1GRCD(3*(BoundaryIndex(i,1)-1)+2)
IF ( 3*(BoundaryIndex(i,1)-1)+2 /= K1GRCI(3*(BoundaryIndex(i,1)-1)+2,m) ) THEN
K1GR( 3*(BoundaryIndex(i,1)-1)+2, m )=0.0D0
ELSE
K1GR( 3*(BoundaryIndex(i,1)-1)+2, m )=1.0D0
ENDIF
ENDDO

ELSEIF (BoundaryIndex(i,2)==3) THEN

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+3)
F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+3,m,1) )=F1Global( K2GR(3*(BoundaryIndex(i,1)-1)+3,m,1) )&
-K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+3,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+3,m,2) )*Dis(3*(BoundaryIndex(i,1)-1)+3)

ENDDO

F1Global(3*(BoundaryIndex(i,1)-1)+3)=Dis(3*(BoundaryIndex(i,1)-1)+3)

DO m=1,K2GRRD(3*(BoundaryIndex(i,1)-1)+3)
K1GR( K2GR(3*(BoundaryIndex(i,1)-1)+3,m,1), K2GR(3*(BoundaryIndex(i,1)-1)+3,m,2) )=0.0D0
ENDDO

DO m=1,K1GRCD(3*(BoundaryIndex(i,1)-1)+3)
IF ( 3*(BoundaryIndex(i,1)-1)+3 /= K1GRCI(3*(BoundaryIndex(i,1)-1)+3,m) ) THEN
K1GR( 3*(BoundaryIndex(i,1)-1)+3, m )=0.0D0
ELSE
K1GR( 3*(BoundaryIndex(i,1)-1)+3, m )=1.0D0
ENDIF
ENDDO

ENDIF

ENDDO

DEALLOCATE(K2GR)
DEALLOCATE(K2GRRD)

TEMP_GRCD=K1GRCD(1)
Highest_GRCD=TEMP_GRCD
DO i=2,3*NumberofNodes+NumberofDE
IF ( K1GRCD(i)>TEMP_GRCD ) THEN
TEMP_GRCD=K1GRCD(i)
Highest_GRCD=TEMP_GRCD
ENDIF
ENDDO

!OPEN(33,FILE="3D_Output_Micro_F/Check3.txt",STATUS="UNKNOWN")
!DO i=1,1000  !,MaxDimenofInverseConnec
!WRITE(33,34) (K1GR(i,k1),k1=1,80)
!ENDDO
!34 FORMAT(<80>E15.8,/)
!CLOSE(33)

!CALL CPU_TIME (ElapsedTime)
!WRITE(*,*) "Highest_GRCD  ",Highest_GRCD

DUMMY=0

DO i=1,3*NumberofNodes+NumberofDE
DO j=1,K1GRCD(i)
IF (K1GR(i,j) /= 0.0D0) THEN
DUMMY=DUMMY+1
ELSEIF (K1GR(i,j) == 0.0D0) THEN
GOTO 3
ENDIF
3 ENDDO
ENDDO

ALLOCATE(IA1(DUMMY))
ALLOCATE(JA1(DUMMY))
ALLOCATE(JA2(DUMMY))
ALLOCATE(A1(DUMMY))

DUMMY=0

DO i=1,3*NumberofNodes+NumberofDE
DO j=1,K1GRCD(i)
IF (K1GR(i,j) /= 0.0D0) THEN
DUMMY=DUMMY+1
A1(DUMMY)=K1GR(i,j)
JA1(DUMMY)=K1GRCI(i,j)
JA2(DUMMY)=i
ELSEIF (K1GR(i,j) == 0.0D0) THEN
GOTO 4
ENDIF
4 ENDDO
ENDDO

DEALLOCATE(K1GR)
DEALLOCATE(K1GRCD)
DEALLOCATE(K1GRCI)

!CALL CPU_TIME (ElapsedTime)
!WRITE(20,*) "DUMMY  ",DUMMY,"  Elapsed Time=",ElapsedTime/60.0D0

IA1(1)=1

DUMMY2=1

DO i=2,DUMMY

IF (JA2(i)>JA2(i-1) .AND. i>1) THEN

DUMMY2=DUMMY2+1

IA1(DUMMY2)=i

ENDIF

ENDDO

DUMMY2=DUMMY2+1

IA1(DUMMY2)=DUMMY+1

ALLOCATE(IA(DUMMY2))

DO i=1,DUMMY2
IA(i)=IA1(i)
ENDDO

DEALLOCATE(IA1)

Dimen1=3*NumberofNodes+NumberofDE

CALL HarwellBoeing_Reduced(k,Dimen1,DUMMY,DUMMY2,A1,JA1,JA2,IA,F1Global,x)

!OPEN(27,FILE="3D_Output_Micro_F/A1JA1JA2.txt",STATUS="UNKNOWN")
!DO i=1,6300
!WRITE(27,29) A1(i),JA1(i),JA2(i)
!ENDDO
!CLOSE(27)
!
!29 FORMAT(E15.8,2X,I,2X,I,/)
!
!OPEN(39,FILE="3D_Output_Micro_F/IA.txt",STATUS="UNKNOWN")
!DO i=1,6300
!WRITE(39,33) IA(i)
!ENDDO
!CLOSE(39)
!
!33 FORMAT(I,/)

DEALLOCATE(F1Global)
DEALLOCATE(IA)
DEALLOCATE(JA1)
DEALLOCATE(JA2)
DEALLOCATE(A1)

DO j=1,NumberofNodes
Dis1(j)=x(3*(j-1)+1)
Dis1(NumberofNodes+j)=x(3*(j-1)+2)
Dis1(2*NumberofNodes+j)=x(3*(j-1)+3)
ENDDO

!DO j=1,NumberofNodes
!Dis1(j)=x(3*(j-1)+1)-x(3*(BoundaryIndex(1,1)-1)+1)
!Dis1(NumberofNodes+j)=x(3*(j-1)+2)-x(3*(BoundaryIndex(2,1)-1)+2)
!Dis1(2*NumberofNodes+j)=x(3*(j-1)+3)-x(3*(BoundaryIndex(3,1)-1)+3)
!ENDDO
!
!WRITE(*,*) BoundaryIndex(1,1),BoundaryIndex(2,1),BoundaryIndex(3,1),x(3*(BoundaryIndex(1,1)-1)+1),x(3*(BoundaryIndex(2,1)-1)+2),x(3*(BoundaryIndex(3,1)-1)+3)
!WRITE(*,*) "Dis60",Dis1(60),Dis1(NumberofNodes+60),Dis1(2*NumberofNodes+60)

RETURN
ENDSUBROUTINE Matrix_DE_Reduced

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OneStepDEs(Nodes,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes,MacroStrainsXX,&
MacroStrainsYY,MacroStrainsZZ,MacroStrainsXY1,MacroStrainsXY2,MacroStrainsXZ1,MacroStrainsXZ2,MacroStrainsYZ1,MacroStrainsYZ2,One_Step_DEs)
IMPLICIT NONE
INTEGER :: i,j
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofNodes,3) :: Nodes

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(IN), DIMENSION(NumberofPositiveXNodes) :: PositiveXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveYNodes) :: PositiveYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveZNodes) :: PositiveZNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeXNodes) :: NegativeXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeYNodes) :: NegativeYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeZNodes) :: NegativeZNodes

DOUBLE PRECISION, INTENT(INOUT), DIMENSION(3*(NumberofPositiveXNodes+NumberofPositiveYNodes+NumberofPositiveZNodes),8) :: One_Step_DEs

DOUBLE PRECISION, INTENT(IN) :: MacroStrainsXX,MacroStrainsYY,MacroStrainsZZ,MacroStrainsXY1,MacroStrainsXY2,MacroStrainsXZ1,MacroStrainsXZ2,&
MacroStrainsYZ1,MacroStrainsYZ2

INTEGER :: PosiTemp,NegaTemp

DO j=1,NumberofPositiveXNodes

PosiTemp=PositiveXNodes(j)
NegaTemp=NegativeXNodes(j)

IF ( 3*(PosiTemp-1)+1 < 3*(NegaTemp-1)+1 ) THEN

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=PositiveXNodes(j)
One_Step_DEs(j,3)=1
One_Step_DEs(j,4)=1.0D0
One_Step_DEs(j,5)=NegativeXNodes(j)
One_Step_DEs(j,6)=1
One_Step_DEs(j,7)=-1.0D0
One_Step_DEs(j,8)=MacroStrainsXX*(Nodes(PositiveXNodes(j),2)-Nodes(NegativeXNodes(j),2))

ELSE

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=NegativeXNodes(j)
One_Step_DEs(j,3)=1
One_Step_DEs(j,4)=-1.0D0
One_Step_DEs(j,5)=PositiveXNodes(j)
One_Step_DEs(j,6)=1
One_Step_DEs(j,7)=1.0D0
One_Step_DEs(j,8)=MacroStrainsXX*(Nodes(PositiveXNodes(j),2)-Nodes(NegativeXNodes(j),2))

ENDIF

ENDDO

DO j=NumberofPositiveXNodes+1,NumberofPositiveXNodes+NumberofPositiveXNodes

PosiTemp=PositiveXNodes(j-NumberofPositiveXNodes)
NegaTemp=NegativeXNodes(j-NumberofPositiveXNodes)

IF ( 3*(PosiTemp-1)+2 < 3*(NegaTemp-1)+2 ) THEN


One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=PositiveXNodes(j-NumberofPositiveXNodes)
One_Step_DEs(j,3)=2
One_Step_DEs(j,4)=1.0D0
One_Step_DEs(j,5)=NegativeXNodes(j-NumberofPositiveXNodes)
One_Step_DEs(j,6)=2
One_Step_DEs(j,7)=-1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsXY2*(Nodes(PositiveXNodes(j-NumberofPositiveXNodes),2)-Nodes(NegativeXNodes(j-NumberofPositiveXNodes),2))

ELSE

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=NegativeXNodes(j-NumberofPositiveXNodes)
One_Step_DEs(j,3)=2
One_Step_DEs(j,4)=-1.0D0
One_Step_DEs(j,5)=PositiveXNodes(j-NumberofPositiveXNodes)
One_Step_DEs(j,6)=2
One_Step_DEs(j,7)=1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsXY2*(Nodes(PositiveXNodes(j-NumberofPositiveXNodes),2)-Nodes(NegativeXNodes(j-NumberofPositiveXNodes),2))

ENDIF

ENDDO

DO j=2*NumberofPositiveXNodes+1,3*NumberofPositiveXNodes

PosiTemp=PositiveXNodes(j-2*NumberofPositiveXNodes)
NegaTemp=NegativeXNodes(j-2*NumberofPositiveXNodes)

IF ( 3*(PosiTemp-1)+3 < 3*(NegaTemp-1)+3 ) THEN

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=PositiveXNodes(j-2*NumberofPositiveXNodes)
One_Step_DEs(j,3)=3
One_Step_DEs(j,4)=1.0D0
One_Step_DEs(j,5)=NegativeXNodes(j-2*NumberofPositiveXNodes)
One_Step_DEs(j,6)=3
One_Step_DEs(j,7)=-1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsXZ2*(Nodes(PositiveXNodes(j-2*NumberofPositiveXNodes),2)-Nodes(NegativeXNodes(j-2*NumberofPositiveXNodes),2))

ELSE

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=NegativeXNodes(j-2*NumberofPositiveXNodes) 
One_Step_DEs(j,3)=3
One_Step_DEs(j,4)=-1.0D0
One_Step_DEs(j,5)=PositiveXNodes(j-2*NumberofPositiveXNodes)
One_Step_DEs(j,6)=3
One_Step_DEs(j,7)=+1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsXZ2*(Nodes(PositiveXNodes(j-2*NumberofPositiveXNodes),2)-Nodes(NegativeXNodes(j-2*NumberofPositiveXNodes),2))

ENDIF

ENDDO

DO j=3*NumberofPositiveXNodes+1,3*NumberofPositiveXNodes+NumberofPositiveYNodes

PosiTemp=PositiveYNodes(j-3*NumberofPositiveXNodes)
NegaTemp=NegativeYNodes(j-3*NumberofPositiveXNodes)

IF ( 3*(PosiTemp-1)+1 < 3*(NegaTemp-1)+1 ) THEN

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=PositiveYNodes(j-3*NumberofPositiveXNodes)
One_Step_DEs(j,3)=1
One_Step_DEs(j,4)=1.0D0
One_Step_DEs(j,5)=NegativeYNodes(j-3*NumberofPositiveXNodes)
One_Step_DEs(j,6)=1
One_Step_DEs(j,7)=-1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsXY1*(Nodes(PositiveYNodes(j-3*NumberofPositiveXNodes),3)-Nodes(NegativeYNodes(j-3*NumberofPositiveXNodes),3))

ELSE

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=NegativeYNodes(j-3*NumberofPositiveXNodes)
One_Step_DEs(j,3)=1
One_Step_DEs(j,4)=-1.0D0
One_Step_DEs(j,5)=PositiveYNodes(j-3*NumberofPositiveXNodes)
One_Step_DEs(j,6)=1
One_Step_DEs(j,7)=1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsXY1*(Nodes(PositiveYNodes(j-3*NumberofPositiveXNodes),3)-Nodes(NegativeYNodes(j-3*NumberofPositiveXNodes),3))

ENDIF

ENDDO

DO j=3*NumberofPositiveXNodes+NumberofPositiveYNodes+1,3*NumberofPositiveXNodes+2*NumberofPositiveYNodes

PosiTemp=PositiveYNodes(j-3*NumberofPositiveXNodes-NumberofPositiveYNodes)
NegaTemp=NegativeYNodes(j-3*NumberofPositiveXNodes-NumberofPositiveYNodes)

IF ( 3*(PosiTemp-1)+2 < 3*(NegaTemp-1)+2 ) THEN

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=PositiveYNodes(j-3*NumberofPositiveXNodes-NumberofPositiveYNodes)
One_Step_DEs(j,3)=2
One_Step_DEs(j,4)=1.0D0
One_Step_DEs(j,5)=NegativeYNodes(j-3*NumberofPositiveXNodes-NumberofPositiveYNodes)
One_Step_DEs(j,6)=2
One_Step_DEs(j,7)=-1.0D0
One_Step_DEs(j,8)=MacroStrainsYY*(Nodes(PositiveYNodes(j-3*NumberofPositiveXNodes-NumberofPositiveYNodes),3)&
-Nodes(NegativeYNodes(j-3*NumberofPositiveXNodes-NumberofPositiveYNodes),3))

ELSE

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=NegativeYNodes(j-3*NumberofPositiveXNodes-NumberofPositiveYNodes)
One_Step_DEs(j,3)=2
One_Step_DEs(j,4)=-1.0D0
One_Step_DEs(j,5)=PositiveYNodes(j-3*NumberofPositiveXNodes-NumberofPositiveYNodes)
One_Step_DEs(j,6)=2
One_Step_DEs(j,7)=1.0D0
One_Step_DEs(j,8)=MacroStrainsYY*(Nodes(PositiveYNodes(j-3*NumberofPositiveXNodes-NumberofPositiveYNodes),3)&
-Nodes(NegativeYNodes(j-3*NumberofPositiveXNodes-NumberofPositiveYNodes),3))

ENDIF

ENDDO

DO j=3*NumberofPositiveXNodes+2*NumberofPositiveYNodes+1,3*NumberofPositiveXNodes+3*NumberofPositiveYNodes

PosiTemp=PositiveYNodes(j-3*NumberofPositiveXNodes-2*NumberofPositiveYNodes)
NegaTemp=NegativeYNodes(j-3*NumberofPositiveXNodes-2*NumberofPositiveYNodes)

IF ( 3*(PosiTemp-1)+3 < 3*(NegaTemp-1)+3 ) THEN

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=PositiveYNodes(j-3*NumberofPositiveXNodes-2*NumberofPositiveYNodes)
One_Step_DEs(j,3)=3
One_Step_DEs(j,4)=1.0D0
One_Step_DEs(j,5)=NegativeYNodes(j-3*NumberofPositiveXNodes-2*NumberofPositiveYNodes)
One_Step_DEs(j,6)=3
One_Step_DEs(j,7)=-1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsYZ2*(Nodes(PositiveYNodes(j-3*NumberofPositiveXNodes-2*NumberofPositiveYNodes),3)&
-Nodes(NegativeYNodes(j-3*NumberofPositiveXNodes-2*NumberofPositiveYNodes),3))

ELSE

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=NegativeYNodes(j-3*NumberofPositiveXNodes-2*NumberofPositiveYNodes)
One_Step_DEs(j,3)=3
One_Step_DEs(j,4)=-1.0D0
One_Step_DEs(j,5)=PositiveYNodes(j-3*NumberofPositiveXNodes-2*NumberofPositiveYNodes)
One_Step_DEs(j,6)=3
One_Step_DEs(j,7)=1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsYZ2*(Nodes(PositiveYNodes(j-3*NumberofPositiveXNodes-2*NumberofPositiveYNodes),3)&
-Nodes(NegativeYNodes(j-3*NumberofPositiveXNodes-2*NumberofPositiveYNodes),3))

ENDIF

ENDDO

DO j=3*NumberofPositiveXNodes+3*NumberofPositiveYNodes+1,3*NumberofPositiveXNodes+3*NumberofPositiveYNodes+NumberofPositiveZNodes

PosiTemp=PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes)
NegaTemp=NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes)

IF ( 3*(PosiTemp-1)+1 < 3*(NegaTemp-1)+1 ) THEN

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes)
One_Step_DEs(j,3)=1
One_Step_DEs(j,4)=1.0D0
One_Step_DEs(j,5)=NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes)
One_Step_DEs(j,6)=1
One_Step_DEs(j,7)=-1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsXZ1*(Nodes(PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes),4)&
-Nodes(NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes),4))

ELSE

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes)
One_Step_DEs(j,3)=1
One_Step_DEs(j,4)=-1.0D0
One_Step_DEs(j,5)=PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes)
One_Step_DEs(j,6)=1
One_Step_DEs(j,7)=1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsXZ1*(Nodes(PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes),4)&
-Nodes(NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes),4))

ENDIF

ENDDO

DO j=3*NumberofPositiveXNodes+3*NumberofPositiveYNodes+NumberofPositiveZNodes+1,&
3*NumberofPositiveXNodes+3*NumberofPositiveYNodes+2*NumberofPositiveZNodes

PosiTemp=PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-NumberofPositiveZNodes)
NegaTemp=NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-NumberofPositiveZNodes)

IF ( 3*(PosiTemp-1)+2 < 3*(NegaTemp-1)+2 ) THEN

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-NumberofPositiveZNodes)
One_Step_DEs(j,3)=2
One_Step_DEs(j,4)=1.0D0
One_Step_DEs(j,5)=NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-NumberofPositiveZNodes)
One_Step_DEs(j,6)=2
One_Step_DEs(j,7)=-1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsYZ1*(Nodes(PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-NumberofPositiveZNodes),4)&
-Nodes(NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-NumberofPositiveZNodes),4))

ELSE

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-NumberofPositiveZNodes)
One_Step_DEs(j,3)=2
One_Step_DEs(j,4)=-1.0D0
One_Step_DEs(j,5)=PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-NumberofPositiveZNodes)
One_Step_DEs(j,6)=2
One_Step_DEs(j,7)=1.0D0
One_Step_DEs(j,8)=2.0D0*MacroStrainsYZ1*(Nodes(PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-NumberofPositiveZNodes),4)&
-Nodes(NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-NumberofPositiveZNodes),4))

ENDIF

ENDDO

DO j=3*NumberofPositiveXNodes+3*NumberofPositiveYNodes+2*NumberofPositiveZNodes+1,&
3*NumberofPositiveXNodes+3*NumberofPositiveYNodes+3*NumberofPositiveZNodes

PosiTemp=PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-2*NumberofPositiveZNodes)
NegaTemp=NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-2*NumberofPositiveZNodes)

IF ( 3*(PosiTemp-1)+3 < 3*(NegaTemp-1)+3 ) THEN

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-2*NumberofPositiveZNodes)
One_Step_DEs(j,3)=3
One_Step_DEs(j,4)=1.0D0
One_Step_DEs(j,5)=NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-2*NumberofPositiveZNodes)
One_Step_DEs(j,6)=3
One_Step_DEs(j,7)=-1.0D0
One_Step_DEs(j,8)=MacroStrainsZZ*(Nodes(PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-2*NumberofPositiveZNodes),4)&
-Nodes(NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-2*NumberofPositiveZNodes),4))

ELSE

One_Step_DEs(j,1)=2
One_Step_DEs(j,2)=NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-2*NumberofPositiveZNodes)
One_Step_DEs(j,3)=3
One_Step_DEs(j,4)=-1.0D0
One_Step_DEs(j,5)=PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-2*NumberofPositiveZNodes)
One_Step_DEs(j,6)=3
One_Step_DEs(j,7)=1.0D0
One_Step_DEs(j,8)=MacroStrainsZZ*(Nodes(PositiveZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-2*NumberofPositiveZNodes),4)&
-Nodes(NegativeZNodes(j-3*NumberofPositiveXNodes-3*NumberofPositiveYNodes-2*NumberofPositiveZNodes),4))

ENDIF

ENDDO

ENDSUBROUTINE OneStepDEs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Input_Strains(k,MacroStrainsXX,MacroStrainsYY,MacroStrainsZZ,MacroStrainsXY1,MacroStrainsXY2,MacroStrainsXZ1,MacroStrainsXZ2,&
MacroStrainsYZ1,MacroStrainsYZ2)
IMPLICIT NONE
INTEGER :: IN8=12
INTEGER, INTENT(IN) :: k
DOUBLE PRECISION, DIMENSION(9) :: MacroStrains
DOUBLE PRECISION, INTENT(OUT) :: MacroStrainsXX,MacroStrainsYY,MacroStrainsZZ,MacroStrainsXY1,MacroStrainsXY2,MacroStrainsXZ1,MacroStrainsXZ2,&
MacroStrainsYZ1,MacroStrainsYZ2

OPEN(IN8,file = '3D_Output_Micro_F/3D_Input/Input_Strains.txt',status ='unknown')
READ (IN8,*) MacroStrains(1)
READ (IN8,*) MacroStrains(2)
READ (IN8,*) MacroStrains(3)
READ (IN8,*) MacroStrains(4)
READ (IN8,*) MacroStrains(5)
READ (IN8,*) MacroStrains(6)
READ (IN8,*) MacroStrains(7)
READ (IN8,*) MacroStrains(8)
READ (IN8,*) MacroStrains(9)
ClOSE(IN8)

IF ( k>= 1 .AND. k<=10) THEN
MacroStrainsXX=MacroStrains(1)*k
MacroStrainsYY=MacroStrains(2)*k
MacroStrainsZZ=MacroStrains(3)*k
MacroStrainsXY1=MacroStrains(4)*k
MacroStrainsXZ1=MacroStrains(5)*k
MacroStrainsYZ1=MacroStrains(6)*k
MacroStrainsXY2=MacroStrains(7)*k
MacroStrainsXZ2=MacroStrains(8)*k
MacroStrainsYZ2=MacroStrains(9)*k
ELSEIF (k>=11) THEN
MacroStrainsXX=-MacroStrains(1)*(k-10)
MacroStrainsYY=-MacroStrains(2)*(k-10)
MacroStrainsZZ=-MacroStrains(3)*(k-10)
MacroStrainsXY1=-MacroStrains(4)*(k-10)
MacroStrainsXZ1=-MacroStrains(5)*(k-10)
MacroStrainsYZ1=-MacroStrains(6)*(k-10)
MacroStrainsXY2=-MacroStrains(7)*(k-10)
MacroStrainsXZ2=-MacroStrains(8)*(k-10)
MacroStrainsYZ2=-MacroStrains(9)*(k-10)
ENDIF

END SUBROUTINE Input_Strains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE KPreprations(Connectivities,KComponent,NumberofDE,K1GRR,K1GRCDD,K1GRCII,F1Globall)
IMPLICIT NONE
INTEGER :: i,j,m,n,i1,j1,k1,k2,k3

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements,24,24) :: KComponent
INTEGER, INTENT(IN) :: NumberofDE
DOUBLE PRECISION, INTENT(OUT), DIMENSION (3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec) :: K1GRR
INTEGER, INTENT(OUT), DIMENSION (3*NumberofNodes+NumberofDE) :: K1GRCDD
INTEGER, INTENT(OUT), DIMENSION (3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec) :: K1GRCII
DOUBLE PRECISION, INTENT(OUT), DIMENSION (3*NumberofNodes+NumberofDE) :: F1Globall

K1GRR=0.0D0
K1GRCDD=0
K1GRCII=0

F1Globall=0.0D0

DO j1=1,3
DO i1=1,3
DO j=1,8
DO i=1,8
DO k1=1,NumberofElements

IF ( KComponent(k1,8*(i1-1)+i,8*(j1-1)+j) /= 0.0D0 ) THEN

IF ( K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)==0 ) THEN

K1GRR(3*(Connectivities(k1,i+1)-1)+i1,1)=KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)

K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)=K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)+1

K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,1)=3*(Connectivities(k1,j+1)-1)+j1

GOTO 1

ELSEIF ( K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1) == 1 ) THEN

IF ( 3*(Connectivities(k1,j+1)-1)+j1 == K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,1) ) THEN
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,1)=K1GRR(3*(Connectivities(k1,i+1)-1)+i1,1)+KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ELSEIF ( 3*(Connectivities(k1,j+1)-1)+j1 > K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,1) ) THEN
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,2)=KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)=K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)+1
K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,2)=3*(Connectivities(k1,j+1)-1)+j1
GOTO 1
ELSEIF ( 3*(Connectivities(k1,j+1)-1)+j1 < K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,1) ) THEN
K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)=K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)+1
DO k3=2,K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k3)=K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k3-1)
K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k3)=K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,1)=KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,1)=3*(Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF

ELSEIF ( K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1) > 1 ) THEN

IF ( 3*(Connectivities(k1,j+1)-1)+j1 < K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,1) ) THEN
K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)=K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)+1
DO k3=K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1),2,-1
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k3)=K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k3-1)
K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k3)=K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,1)=KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,1)=3*(Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF

DO k2=1,K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)-1

IF ( 3*(Connectivities(k1,j+1)-1)+j1 == K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k2) ) THEN
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k2)=K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k2)+KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ENDIF

IF (3*(Connectivities(k1,j+1)-1)+j1 > K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k2) .AND. &
3*(Connectivities(k1,j+1)-1)+j1 < K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k2+1) ) THEN
K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)=K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)+1
DO k3=K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1),k2+2,-1
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k3)=K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k3-1)
K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k3)=K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k2+1)=KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k2+1)=3*(Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF  

ENDDO  !k2

DO k2=K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1),K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)

IF ( 3*(Connectivities(k1,j+1)-1)+j1 == K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k2) ) THEN
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k2)=K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k2)+KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ENDIF

IF ( 3*(Connectivities(k1,j+1)-1)+j1 > K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k2)  ) THEN
K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)=K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)+1
DO k3=k2+1,K1GRCDD(3*(Connectivities(k1,i+1)-1)+i1)
K1GRR(3*(Connectivities(k1,i+1)-1)+i1,k3)=KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
K1GRCII(3*(Connectivities(k1,i+1)-1)+i1,k3)=3*(Connectivities(k1,j+1)-1)+j1
ENDDO
GOTO 1
ENDIF  

ENDDO !For k2

1 ENDIF  ! For ( K1GRCD(3*(Connectivities(k1,i+1)-1)+i1)==0 )
 
ENDIF  ! For ( KComponent(k1,8*(i1-1)+i,8*(j1-1)+j) /= 0.0D0 )

ENDDO
ENDDO
ENDDO
ENDDO
ENDDO

ENDSUBROUTINE KPreprations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ConcentricAnisotropic(Anisotropic,Nodes,Connectivities,BeginofModeIndex,EndofModeIndex,XYofAxis,MaximumRadius,CC,Sub_CC,NN,RotationMatrix,&
RotationMatrix2)
IMPLICIT NONE
INTEGER :: i,j,i1,j1,j2,j3,j4,j5,j6,j7,j8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(IN) :: Anisotropic
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: Nodes
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
INTEGER, INTENT(IN) :: BeginofModeIndex,EndofModeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION(11,2) :: XYofAxis
DOUBLE PRECISION, INTENT(IN) :: MaximumRadius
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NumberofMaterialModes,6,6) :: CC
DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NumberofGaussPoints,NumberofElements,6,6) :: Sub_CC
DOUBLE PRECISION, DIMENSION(EndofModeIndex-BeginofModeIndex+1,6,6) :: Ori_CC
DOUBLE PRECISION, DIMENSION(EndofModeIndex-BeginofModeIndex+1,3,3,3,3) :: Ori_CC_Temp
DOUBLE PRECISION, DIMENSION(EndofModeIndex-BeginofModeIndex+1,3,3,3,3) :: CC_Temp

DOUBLE PRECISION, DIMENSION(EndofModeIndex-BeginofModeIndex+1,3) :: Centroid

DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofGaussPoints,8) :: NN

DOUBLE PRECISION, DIMENSION(EndofModeIndex-BeginofModeIndex+1) :: COSTHETA
DOUBLE PRECISION, DIMENSION(EndofModeIndex-BeginofModeIndex+1) :: SINTHETA

DOUBLE PRECISION, DIMENSION(7) :: r

!DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,EndofModeIndex-BeginofModeIndex+1,3,3) :: RotationMatrix
DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofMaterialModes,3,3) :: RotationMatrix
DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofMaterialModes,3,3) :: RotationMatrix2

DOUBLE PRECISION, DIMENSION(EndofModeIndex-BeginofModeIndex+1) :: thetaa,phii

INTEGER :: values(1:8), k
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
DOUBLE PRECISION :: theta,phi,sign

INTEGER :: k1,OUT3

k1=0
OUT3=23

CALL date_and_time(values=values)
CALL random_seed(size=k)
ALLOCATE(seed(1:k))
seed(:) = values(8)
CALL random_seed(put=seed)
CALL random_number(theta)
CALL random_number(phi)
CALL random_number(sign)

RotationMatrix=0.0D0
RotationMatrix2=0.0D0

DO i=1,NumberofElements

DO j=1,NumberofGaussPoints

IF ( Connectivities(i,j+9) >= BeginofModeIndex .AND. Connectivities(i,j+9) <= EndofModeIndex) THEN

Centroid(Connectivities(i,j+9)-BeginofModeIndex+1,1)=Nodes(Connectivities(i,2),2)*NN(j,1)+Nodes(Connectivities(i,3),2)*NN(j,2)+&
Nodes(Connectivities(i,4),2)*NN(j,3)+Nodes(Connectivities(i,5),2)*NN(j,4)+Nodes(Connectivities(i,6),2)*NN(j,5)&
+Nodes(Connectivities(i,7),2)*NN(j,6)+Nodes(Connectivities(i,8),2)*NN(j,7)+Nodes(Connectivities(i,9),2)*NN(j,8)

Centroid(Connectivities(i,j+9)-BeginofModeIndex+1,2)=Nodes(Connectivities(i,2),3)*NN(j,1)+Nodes(Connectivities(i,3),3)*NN(j,2)+&
Nodes(Connectivities(i,4),3)*NN(j,3)+Nodes(Connectivities(i,5),3)*NN(j,4)+Nodes(Connectivities(i,6),3)*NN(j,5)&
+Nodes(Connectivities(i,7),3)*NN(j,6)+Nodes(Connectivities(i,8),3)*NN(j,7)+Nodes(Connectivities(i,9),3)*NN(j,8)

Centroid(Connectivities(i,j+9)-BeginofModeIndex+1,3)=Nodes(Connectivities(i,2),4)*NN(j,1)+Nodes(Connectivities(i,3),4)*NN(j,2)+&
Nodes(Connectivities(i,4),4)*NN(j,3)+Nodes(Connectivities(i,5),4)*NN(j,4)+Nodes(Connectivities(i,6),4)*NN(j,5)&
+Nodes(Connectivities(i,7),4)*NN(j,6)+Nodes(Connectivities(i,8),4)*NN(j,7)+Nodes(Connectivities(i,9),4)*NN(j,8)

ENDIF

ENDDO

ENDDO

IF (Anisotropic==1) THEN
DO i=1,EndofModeIndex-BeginofModeIndex+1
r=0.0D0
!DO i1=1,NumberofGaussPoints
DO j=1,7
r(j)=SQRT((Centroid(i,1)-XYofAxis(j,1))**2+(Centroid(i,2)-XYofAxis(j,2))**2)
IF (r(j) < MaximumRadius) THEN
COSTHETA(i)= (Centroid(i,1)-XYofAxis(j,1)) /SQRT( (Centroid(i,1)-XYofAxis(j,1))**2+ (Centroid(i,2)-XYofAxis(j,2))**2)
SINTHETA(i)= (Centroid(i,2)-XYofAxis(j,2)) /SQRT( (Centroid(i,1)-XYofAxis(j,1))**2+ (Centroid(i,2)-XYofAxis(j,2))**2)
GOTO 1
ELSE
ENDIF
1 ENDDO
!ENDDO
ENDDO
ENDIF

Ori_CC_Temp=0.0D0
DO i=1,EndofModeIndex-BeginofModeIndex+1

DO j1=1,6
DO j2=1,6
Ori_CC(i,j1,j2)=CC(i+BeginofModeIndex-1,j1,j2)
ENDDO
ENDDO

Ori_CC_Temp(i,1,1,1,1)=Ori_CC(i,1,1)
Ori_CC_Temp(i,1,1,2,2)=Ori_CC(i,1,2)
Ori_CC_Temp(i,1,1,3,3)=Ori_CC(i,1,3)
Ori_CC_Temp(i,2,2,1,1)=Ori_CC(i,2,1)
Ori_CC_Temp(i,2,2,2,2)=Ori_CC(i,2,2)
Ori_CC_Temp(i,2,2,3,3)=Ori_CC(i,2,3)
Ori_CC_Temp(i,3,3,1,1)=Ori_CC(i,3,1)
Ori_CC_Temp(i,3,3,2,2)=Ori_CC(i,3,2)
Ori_CC_Temp(i,3,3,3,3)=Ori_CC(i,3,3)

Ori_CC_Temp(i,2,3,2,3)=Ori_CC(i,4,4)
Ori_CC_Temp(i,2,3,3,2)=Ori_CC(i,4,4)
Ori_CC_Temp(i,3,2,2,3)=Ori_CC(i,4,4)
Ori_CC_Temp(i,3,2,3,2)=Ori_CC(i,4,4)

Ori_CC_Temp(i,1,3,1,3)=Ori_CC(i,5,5)
Ori_CC_Temp(i,1,3,3,1)=Ori_CC(i,5,5)
Ori_CC_Temp(i,3,1,1,3)=Ori_CC(i,5,5)
Ori_CC_Temp(i,3,1,3,1)=Ori_CC(i,5,5)

Ori_CC_Temp(i,1,2,1,2)=Ori_CC(i,6,6)
Ori_CC_Temp(i,1,2,2,1)=Ori_CC(i,6,6)
Ori_CC_Temp(i,2,1,1,2)=Ori_CC(i,6,6)
Ori_CC_Temp(i,2,1,2,1)=Ori_CC(i,6,6)

ENDDO


!OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output/theta_phi.txt',status ='unknown')
!DO j1=1,EndofModeIndex-BeginofModeIndex+1
!READ(OUT3,*) i,i1,thetaa(j1),phii(j1)
!ENDDO
!CLOSE(OUT3)


CC_Temp=0.0D0
DO i1=1,NumberofGaussPoints
DO i=1,NumberofElements  !EndofModeIndex-BeginofModeIndex+1

IF ( Connectivities(i,i1+9) >= BeginofModeIndex .AND. Connectivities(i,i1+9) <= EndofModeIndex) THEN

IF (Anisotropic==1) THEN
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,1,1)=COSTHETA(Connectivities(i,9+i1)-BeginofModeIndex+1)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,1,2)=SINTHETA(Connectivities(i,9+i1)-BeginofModeIndex+1)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,2,1)=-SINTHETA(Connectivities(i,9+i1)-BeginofModeIndex+1)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,2,2)=COSTHETA(Connectivities(i,9+i1)-BeginofModeIndex+1)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,3,3)=1.0D0
RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,1,1)=COSTHETA(Connectivities(i,9+i1)-BeginofModeIndex+1)
RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,1,2)=-SINTHETA(Connectivities(i,9+i1)-BeginofModeIndex+1)
RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,2,1)=SINTHETA(Connectivities(i,9+i1)-BeginofModeIndex+1)
RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,2,2)=COSTHETA(Connectivities(i,9+i1)-BeginofModeIndex+1)
RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,3,3)=1.0D0
ELSEIF (Anisotropic==2) THEN
CALL random_number(theta)
CALL random_number(phi)
CALL random_number(sign)
!theta=0.0D0
!phi=0.5D0
!CALL random_seed() 
!CALL random_number(theta)
!CALL random_number(phi)
IF (sign >=0.5D0) THEN
sign=1.0D0
ELSE
sign=-1.0D0
ENDIF
!WRITE(*,*) theta,sign*phi

k1=k1+1

!IF (k1<=EndofModeIndex-BeginofModeIndex+1) THEN
!theta=thetaa(k1) !2*theta*3.141592653589793 !2*acos(theta)
!phi=phii(k1) !0.0D0 !asin(sign*phi)  !asin(sign*phi) !-0.5D0*3.141592653589793 !
!ELSE
theta=2*theta*3.141592653589793
phi=asin(sign*phi) !-0.5D0*3.141592653589793 !asin(sign*phi) !-asin(( 0.99+(1.0D0-0.99)*phi )) !-0.5D0*3.141592653589793 !
!ENDIF

IF (k1==1 .AND. k1<=EndofModeIndex-BeginofModeIndex+1) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output/theta_phi.txt',status ='unknown')
WRITE(OUT3,20) i,i1,theta,phi
ELSEIF ( k1>1 .AND. k1<=EndofModeIndex-BeginofModeIndex+1) THEN
OPEN(OUT3,file = '3D_Output_Micro_F/3D_Output/theta_phi.txt',status ='unknown',position='append')
WRITE(OUT3,20) i,i1,theta,phi
ENDIF
CLOSE(OUT3)
20 FORMAT(I,2X,I,2X,E15.8,2X,E15.8)

RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,1,1)=cos(theta)*cos(phi)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,1,2)=cos(phi)*sin(theta)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,1,3)=sin(phi)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,2,1)=-sin(theta)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,2,2)=cos(theta)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,2,3)=0.0D0
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,3,1)=-cos(theta)*sin(phi)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,3,2)=-sin(theta)*sin(phi)
RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,3,3)=cos(phi)
!phi=2*phi*3.141592653589793 !2*acos(theta)
!theta=acos(sign*theta)
!RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,1,1)=sin(theta)*cos(phi)
!RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,1,2)=-cos(phi)*cos(theta)
!RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,1,3)=sin(phi)
!RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,2,1)=sin(phi)*sin(theta)
!RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,2,2)=-cos(phi)*cos(theta)
!RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,2,3)=-cos(phi)
!RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,3,1)=cos(theta)
!RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,3,2)=sin(theta)
!RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,3,3)=0.0D0

DO j1=1,3
DO j2=1,3
RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,j1,j2)=RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,j2,j1)
ENDDO
ENDDO

!CALL INV(i,3,RotationMatrix(Connectivities(i,9+i1)-BeginofModeIndex+1,:,:),&
!RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,:,:))
ENDIF

DO j1=1,3
DO j2=1,3
DO j3=1,3
DO j4=1,3
DO j5=1,3
DO j6=1,3
DO j7=1,3
DO j8=1,3
CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,j1,j2,j3,j4)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,j1,j2,j3,j4)+&
RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,j1,j5)*RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,j2,j6)*&
RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,j3,j7)*RotationMatrix2(Connectivities(i,9+i1)-BeginofModeIndex+1,j4,j8)*&
Ori_CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,j5,j6,j7,j8)
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO
ENDDO

ENDIF

ENDDO
ENDDO

DO i1=1,NumberofGaussPoints
DO i=1,NumberofElements    !EndofModeIndex-BeginofModeIndex+1

IF ( Connectivities(i,i1+9) >= BeginofModeIndex .AND. Connectivities(i,i1+9) <= EndofModeIndex) THEN
Sub_CC(i1,i,1,1)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,1,1,1)
Sub_CC(i1,i,1,2)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,1,2,2)
Sub_CC(i1,i,1,3)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,1,3,3)
Sub_CC(i1,i,1,4)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,1,2,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,1,3,2)
Sub_CC(i1,i,1,5)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,1,1,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,1,3,1)
Sub_CC(i1,i,1,6)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,1,1,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,1,2,1)

Sub_CC(i1,i,2,1)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,2,1,1)
Sub_CC(i1,i,2,2)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,2,2,2)
Sub_CC(i1,i,2,3)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,2,3,3)
Sub_CC(i1,i,2,4)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,2,2,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,2,3,2)
Sub_CC(i1,i,2,5)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,2,1,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,2,3,1)
Sub_CC(i1,i,2,6)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,2,1,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,2,2,1)

Sub_CC(i1,i,3,1)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,3,1,1)
Sub_CC(i1,i,3,2)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,3,2,2)
Sub_CC(i1,i,3,3)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,3,3,3)
Sub_CC(i1,i,3,4)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,3,2,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,3,3,2)
Sub_CC(i1,i,3,5)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,3,1,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,3,3,1)
Sub_CC(i1,i,3,6)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,3,1,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,3,2,1)

Sub_CC(i1,i,4,1)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,3,1,1) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,2,1,1)
Sub_CC(i1,i,4,2)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,3,2,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,2,2,2)
Sub_CC(i1,i,4,3)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,3,3,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,2,3,3)
Sub_CC(i1,i,4,4)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,3,2,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,3,3,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,2,2,3)+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,2,3,2)
Sub_CC(i1,i,4,5)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,3,1,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,3,3,1) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,2,1,3)+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,2,3,1)
Sub_CC(i1,i,4,6)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,3,1,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,3,2,1) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,2,1,2)+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,2,2,1)

Sub_CC(i1,i,5,1)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,3,1,1) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,1,1,1)
Sub_CC(i1,i,5,2)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,3,2,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,1,2,2)
Sub_CC(i1,i,5,3)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,3,3,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,1,3,3)
Sub_CC(i1,i,5,4)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,3,2,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,3,3,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,1,2,3)+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,1,3,2)
Sub_CC(i1,i,5,5)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,3,1,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,3,3,1) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,1,1,3)+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,1,3,1)
Sub_CC(i1,i,5,6)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,3,1,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,3,2,1) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,1,1,2)+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,3,1,2,1)

Sub_CC(i1,i,6,1)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,2,1,1) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,1,1,1)
Sub_CC(i1,i,6,2)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,2,2,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,1,2,2)
Sub_CC(i1,i,6,3)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,2,3,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,1,3,3)
Sub_CC(i1,i,6,4)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,2,2,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,2,3,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,1,2,3)+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,1,3,2)
Sub_CC(i1,i,6,5)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,2,1,3) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,2,3,1) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,1,1,3)+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,1,3,1)
Sub_CC(i1,i,6,6)=CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,2,1,2) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,1,2,2,1) !+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,1,1,2)+CC_Temp(Connectivities(i,i1+9)-BeginofModeIndex+1,2,1,2,1)
!WRITE(*,*) Sub_CC(i1,i,6,1),Sub_CC(i1,i,6,2),Sub_CC(i1,i,6,3),Sub_CC(i1,i,6,4),Sub_CC(i1,i,6,5),Sub_CC(i1,i,6,6)
ENDIF

ENDDO
ENDDO

ENDSUBROUTINE ConcentricAnisotropic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Components(Nodes,Connectivities,CC,Sub_CC,Vol,KComponent,DNDX,DNDY,DNDZ,TestIndex1,TestIndex2,BeginofModeIndex,EndofModeIndex,&
Anisotropic,RotationMatrix,RotationMatrix2,KIC,KII)
IMPLICIT NONE
INTEGER :: i,j,k,i1,i2,i3,i4,j1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofNodes,4) :: Nodes
INTEGER, INTENT(IN), DIMENSION (NumberofElements,17) :: Connectivities
DOUBLE PRECISION, INTENT(IN),  DIMENSION (NumberofMaterialModes,6,6) :: CC
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofGaussPoints,NumberofElements,6,6) :: Sub_CC
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofGaussPoints,NumberofElements) :: Vol
DOUBLE PRECISION, DIMENSION (NumberofGaussPoints,NumberofElements,6,24) :: BB
DOUBLE PRECISION, DIMENSION (NumberofGaussPoints,NumberofElements,24,6) :: BBT

DOUBLE PRECISION, DIMENSION (NumberofGaussPoints,NumberofElements,6,9) :: BBI
DOUBLE PRECISION, DIMENSION (NumberofGaussPoints,NumberofElements,9,6) :: BBIT

DOUBLE PRECISION, DIMENSION (NumberofGaussPoints,NumberofElements,6,9) :: BBIC

DOUBLE PRECISION, DIMENSION (NumberofGaussPoints,NumberofElements,6,9) :: BBIBar
DOUBLE PRECISION, DIMENSION (NumberofGaussPoints,NumberofElements,9,6) :: BBIBarT

DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,24,24) :: KComponent

DOUBLE PRECISION, DIMENSION (NumberofElements,24,24) :: KCC
DOUBLE PRECISION, DIMENSION (NumberofElements,24,9)  :: KCI
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,9,24)  :: KIC
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,9,9)   :: KII
DOUBLE PRECISION, DIMENSION (NumberofElements,9,9)   :: KII_Inverse
DOUBLE PRECISION, DIMENSION (NumberofElements,9,9)   :: KII2

DOUBLE PRECISION, DIMENSION(NumberofElements,8):: x,y,z

DOUBLE PRECISION, DIMENSION(8) :: Xi,Eta,Mu

DOUBLE PRECISION, DIMENSION(NumberofGaussPoints) :: Xi1,Eta1,Mu1,W

DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,NumberofElements,3,3) :: JJ,InverseJJ

DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,8+3) :: DNDXi,DNDEta,DNDMu
DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,8) :: NN
DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofGaussPoints,NumberofElements,8+3) :: DNDX,DNDY,DNDZ
DOUBLE PRECISION, DIMENSION(NumberofGaussPoints,NumberofElements,8+3) :: QX,QY,QZ

DOUBLE PRECISION :: Delta_A

DOUBLE PRECISION :: a11,a12,a13,a21,a22,a23,a31,a32,a33

DOUBLE PRECISION :: Vol_Tep

DOUBLE PRECISION :: Vf,R,MaximumRadius

INTEGER, INTENT(OUT) :: BeginofModeIndex,EndofModeIndex,Anisotropic

INTEGER, INTENT(OUT) :: TestIndex1,TestIndex2

DOUBLE PRECISION, DIMENSION(11,2) :: XYofAxis

DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofMaterialModes,3,3) :: RotationMatrix
DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofMaterialModes,3,3) :: RotationMatrix2

INTEGER :: IN9=13

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OPEN (IN9,file = '3D_Output_Micro_F/3D_Input/More_Inputs.txt',status ='unknown')
READ (IN9,*) Vf
READ (IN9,*) R
READ (IN9,*) BeginofModeIndex
READ (IN9,*) EndofModeIndex
READ (IN9,*) MaximumRadius !For the subroutine concentric anisotropic
READ (IN9,*) Anisotropic
READ (IN9,*) TestIndex1
READ (IN9,*) TestIndex2
ClOSE(IN9)

!DO i=1,NumberofGaussPoints
!DO j=1,NumberofMaterialModes
!DO i2=1,6
!DO i3=1,6
!Sub_CC(i,j,i2,i3)=CC(j,i2,i3)
!ENDDO
!ENDDO
!ENDDO
!ENDDO

DO i=1,NumberofGaussPoints
DO j=1,NumberofElements
DO i2=1,6
DO i3=1,6
Sub_CC(i,j,i2,i3)=CC(Connectivities(j,i+9),i2,i3)
ENDDO
ENDDO
ENDDO
ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


KII2=0.0D0

JJ=0.0D0
InverseJJ=0.0D0

BB=0.0D0
BBT=0.0D0
KComponent=0.0D0
KCI=0.0D0
KIC=0.0D0
KCC=0.0D0
KII=0.0D0
KII_Inverse=0.0D0
BBI=0.0D0
BBIC=0.0D0
BBIBar=0.0D0
BBIBarT=0.0D0

Xi=(/1.0D0,+1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,-1.0D0,-1.0D0/)
Eta=(/-1.0D0,1.0D0,1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,-1.0D0/)
Mu=(/-1.0D0,-1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,1.0D0,1.0D0/)

DO i=1,NumberofGaussPoints
Xi1(i)=Xi(i)/SQRT(3.0D0)
Eta1(i)=Eta(i)/SQRT(3.0D0)
Mu1(i)=Mu(i)/SQRT(3.0D0)
W(i)=1.0D0
ENDDO

DO k=1,NumberofElements
DO j=1,8
x(k,j)=Nodes(Connectivities(k,j+1),2)
y(k,j)=Nodes(Connectivities(k,j+1),3)
z(k,j)=Nodes(Connectivities(k,j+1),4)
ENDDO
ENDDO

DO i=1,NumberofGaussPoints

DNDXi(i,4)=-1.0D0/8.0D0*(1.0D0-Eta1(i))*(1.0D0-Mu1(i))
DNDEta(i,4)=-1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0-Mu1(i))
DNDMu(i,4)=-1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0-Eta1(i))

DNDXi(i,1)=1.0D0/8.0D0*(1.0D0-Eta1(i))*(1.0D0-Mu1(i))
DNDEta(i,1)=-1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0-Mu1(i))
DNDMu(i,1)=-1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0-Eta1(i))

DNDXi(i,2)=1.0D0/8.0D0*(1.0D0+Eta1(i))*(1.0D0-Mu1(i))
DNDEta(i,2)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0-Mu1(i))
DNDMu(i,2)=-1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0+Eta1(i))

DNDXi(i,3)=-1.0D0/8.0D0*(1.0D0+Eta1(i))*(1.0D0-Mu1(i))
DNDEta(i,3)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0-Mu1(i))
DNDMu(i,3)=-1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0+Eta1(i))

DNDXi(i,8)=-1.0D0/8.0D0*(1.0D0-Eta1(i))*(1.0D0+Mu1(i))
DNDEta(i,8)=-1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0+Mu1(i))
DNDMu(i,8)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0-Eta1(i))

DNDXi(i,5)=1.0D0/8.0D0*(1.0D0-Eta1(i))*(1.0D0+Mu1(i))
DNDEta(i,5)=-1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0+Mu1(i))
DNDMu(i,5)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0-Eta1(i))

DNDXi(i,6)=1.0D0/8.0D0*(1.0D0+Eta1(i))*(1.0D0+Mu1(i))
DNDEta(i,6)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0+Mu1(i))
DNDMu(i,6)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0+Eta1(i))

DNDXi(i,7)=-1.0D0/8.0D0*(1.0D0+Eta1(i))*(1.0D0+Mu1(i))
DNDEta(i,7)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0+Mu1(i))
DNDMu(i,7)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0+Eta1(i))

DNDXi(i,9)=-2.0D0*Xi1(i)
DNDEta(i,9)=0.0D0
DNDMu(i,9)=0.0D0

DNDXi(i,10)=0.0D0
DNDEta(i,10)=-2.0D0*Eta1(i)
DNDMu(i,10)=0.0D0

DNDXi(i,11)=0.0D0
DNDEta(i,11)=0.0D0
DNDMu(i,11)=-2.0D0*Mu1(i)

W(i)=1.0D0

ENDDO

DO i=1,NumberofGaussPoints
DO k=1,NumberofElements

DO j=1,8
JJ(i,k,1,1)=JJ(i,k,1,1)+x(k,j)*DNDXi(i,j)
JJ(i,k,1,2)=JJ(i,k,1,2)+y(k,j)*DNDXi(i,j)
JJ(i,k,1,3)=JJ(i,k,1,3)+z(k,j)*DNDXi(i,j)

JJ(i,k,2,1)=JJ(i,k,2,1)+x(k,j)*DNDEta(i,j)
JJ(i,k,2,2)=JJ(i,k,2,2)+y(k,j)*DNDEta(i,j)
JJ(i,k,2,3)=JJ(i,k,2,3)+z(k,j)*DNDEta(i,j)

JJ(i,k,3,1)=JJ(i,k,3,1)+x(k,j)*DNDMu(i,j)
JJ(i,k,3,2)=JJ(i,k,3,2)+y(k,j)*DNDMu(i,j)
JJ(i,k,3,3)=JJ(i,k,3,3)+z(k,j)*DNDMu(i,j)

!JJ(i,k,1,1)=JJ(i,k,1,1)+x(k,j)*DNDXi(i,j)
!JJ(i,k,1,2)=JJ(i,k,1,2)+x(k,j)*DNDEta(i,j)
!JJ(i,k,1,3)=JJ(i,k,1,3)+x(k,j)*DNDMu(i,j)
!
!JJ(i,k,2,1)=JJ(i,k,2,1)+y(k,j)*DNDXi(i,j)
!JJ(i,k,2,2)=JJ(i,k,2,2)+y(k,j)*DNDEta(i,j)
!JJ(i,k,2,3)=JJ(i,k,2,3)+y(k,j)*DNDMu(i,j)
!
!JJ(i,k,3,1)=JJ(i,k,3,1)+z(k,j)*DNDXi(i,j)
!JJ(i,k,3,2)=JJ(i,k,3,2)+z(k,j)*DNDEta(i,j)
!JJ(i,k,3,3)=JJ(i,k,3,3)+z(k,j)*DNDMu(i,j)
ENDDO

ENDDO
ENDDO

DO i=1,NumberofGaussPoints
DO k=1,NumberofElements

!Delta_A=JJ(i,k,1,1)*(JJ(i,k,2,2)*JJ(i,k,3,3)-JJ(i,k,2,3)*JJ(i,k,3,2))+&
!JJ(i,k,1,2)*(JJ(i,k,3,2)*JJ(i,k,1,3)-JJ(i,k,1,2)*JJ(i,k,3,3))+&
!JJ(i,k,1,3)*(JJ(i,k,1,2)*JJ(i,k,2,3)-JJ(i,k,1,3)*JJ(i,k,2,2))
!
!Vol(i,k)=Delta_A
!
!InverseJJ(i,k,1,1)=(JJ(i,k,2,2)*JJ(i,k,3,3)-JJ(i,k,2,3)*JJ(i,k,3,2))/Delta_A
!InverseJJ(i,k,1,2)=(JJ(i,k,2,3)*JJ(i,k,3,1)-JJ(i,k,2,1)*JJ(i,k,3,3))/Delta_A
!InverseJJ(i,k,1,3)=(JJ(i,k,2,1)*JJ(i,k,2,2)-JJ(i,k,3,1)*JJ(i,k,2,2))/Delta_A
!InverseJJ(i,k,2,1)=(JJ(i,k,3,2)*JJ(i,k,1,3)-JJ(i,k,1,2)*JJ(i,k,3,3))/Delta_A
!InverseJJ(i,k,2,2)=(JJ(i,k,3,3)*JJ(i,k,1,1)-JJ(i,k,3,1)*JJ(i,k,1,3))/Delta_A
!InverseJJ(i,k,2,3)=(JJ(i,k,3,1)*JJ(i,k,1,2)-JJ(i,k,3,2)*JJ(i,k,1,1))/Delta_A
!InverseJJ(i,k,3,1)=(JJ(i,k,1,2)*JJ(i,k,2,3)-JJ(i,k,1,3)*JJ(i,k,2,2))/Delta_A
!InverseJJ(i,k,3,2)=(JJ(i,k,1,3)*JJ(i,k,2,1)-JJ(i,k,2,3)*JJ(i,k,1,1))/Delta_A
!InverseJJ(i,k,3,3)=(JJ(i,k,1,1)*JJ(i,k,2,2)-JJ(i,k,1,2)*JJ(i,k,2,1))/Delta_A

a11=JJ(i,k,1,1)
a12=JJ(i,k,1,2)
a13=JJ(i,k,1,3)
a21=JJ(i,k,2,1)
a22=JJ(i,k,2,2)
a23=JJ(i,k,2,3)
a31=JJ(i,k,3,1)
a32=JJ(i,k,3,2)
a33=JJ(i,k,3,3)

Vol(i,k)=-a13*a22*a31 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 +a11*a22*a33

InverseJJ(i,k,1,1)=(-a23*a32 + a22*a33)/Vol(i,k)
InverseJJ(i,k,1,2)=(a13*a32 - a12*a33)/Vol(i,k)
InverseJJ(i,k,1,3)=(-a13*a22 + a12*a23)/Vol(i,k)
InverseJJ(i,k,2,1)=(a23*a31 - a21*a33)/Vol(i,k)
InverseJJ(i,k,2,2)=(-a13*a31 + a11*a33)/Vol(i,k)
InverseJJ(i,k,2,3)=(a13*a21 - a11*a23)/Vol(i,k)
InverseJJ(i,k,3,1)=(-a22*a31 + a21*a32)/Vol(i,k)
InverseJJ(i,k,3,2)=(a12*a31 - a11*a32)/Vol(i,k)
InverseJJ(i,k,3,3)=(-a12*a21 + a11*a22)/Vol(i,k)

ENDDO
ENDDO

DO i=1,NumberofGaussPoints
DO k=1,NumberofElements
DO j=1,8+3

DNDX(i,k,j)=InverseJJ(i,k,1,1)*DNDXi(i,j)+InverseJJ(i,k,1,2)*DNDEta(i,j)+InverseJJ(i,k,1,3)*DNDMu(i,j)
DNDY(i,k,j)=InverseJJ(i,k,2,1)*DNDXi(i,j)+InverseJJ(i,k,2,2)*DNDEta(i,j)+InverseJJ(i,k,2,3)*DNDMu(i,j)
DNDZ(i,k,j)=InverseJJ(i,k,3,1)*DNDXi(i,j)+InverseJJ(i,k,3,2)*DNDEta(i,j)+InverseJJ(i,k,3,3)*DNDMu(i,j)

ENDDO
ENDDO
ENDDO

DO i=1,NumberofGaussPoints
DO k=1,NumberofElements
DO j=1,8+3

QX(i,k,j)=DNDX(i,k,j)
QY(i,k,j)=DNDY(i,k,j)
QZ(i,k,j)=DNDZ(i,k,j)

ENDDO
ENDDO
ENDDO

DO i=1,NumberofGaussPoints
DO k=1,NumberofElements
DO j=1,8

BB(i,k,1,j)=QX(i,k,j)
BB(i,k,5,j)=QZ(i,k,j)
BB(i,k,6,j)=QY(i,k,j)
BB(i,k,2,8+j)=QY(i,k,j)
BB(i,k,4,8+j)=QZ(i,k,j)
BB(i,k,6,8+j)=QX(i,k,j)
BB(i,k,3,16+j)=QZ(i,k,j)
BB(i,k,4,16+j)=QY(i,k,j)
BB(i,k,5,16+j)=QX(i,k,j)

ENDDO
ENDDO
ENDDO

DO i=1,NumberofGaussPoints
DO k=1,NumberofElements
DO i1=1,6
DO i2=1,24
BBT(i,k,i2,i1)=BB(i,k,i1,i2)
ENDDO
ENDDO
ENDDO
ENDDO


DO i=1,NumberofGaussPoints
DO k=1,NumberofElements
DO j=1,3

BBI(i,k,1,j)=QX(i,k,j+8)
BBI(i,k,5,j)=QZ(i,k,j+8)
BBI(i,k,6,j)=QY(i,k,j+8)
BBI(i,k,2,3+j)=QY(i,k,j+8)
BBI(i,k,4,3+j)=QZ(i,k,j+8)
BBI(i,k,6,3+j)=QX(i,k,j+8)
BBI(i,k,3,6+j)=QZ(i,k,j+8)
BBI(i,k,4,6+j)=QY(i,k,j+8)
BBI(i,k,5,6+j)=QX(i,k,j+8)

ENDDO
ENDDO
ENDDO

DO i=1,NumberofGaussPoints
DO k=1,NumberofElements
DO i1=1,6
DO i2=1,9
BBIT(i,k,i2,i1)=BBI(i,k,i1,i2)
ENDDO
ENDDO
ENDDO
ENDDO

BBIC=0.0D0
DO i=1,NumberofGaussPoints
DO k=1,NumberofElements
DO i1=1,6
DO j1=1,9

Vol_Tep=0.0D0

DO i2=1,NumberofGaussPoints
Vol_Tep=Vol_Tep+Vol(i2,k)
BBIC(i,k,i1,j1)=BBIC(i,k,i1,j1)+BBI(i2,k,i1,j1)*Vol(i2,k)*W(i2)
ENDDO

BBIC(i,k,i1,j1)=-BBIC(i,k,i1,j1)/Vol_Tep/8.0D0
BBIBar(i,k,i1,j1)=BBI(i,k,i1,j1)+BBIC(i,k,i1,j1)

ENDDO
ENDDO
ENDDO
ENDDO

DO i=1,NumberofGaussPoints
DO k=1,NumberofElements
DO i1=1,6
DO i2=1,9
BBIBarT(i,k,i2,i1)=BBIBar(i,k,i1,i2)
ENDDO
ENDDO
ENDDO
ENDDO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (Anisotropic ==1 .OR. Anisotropic==2) THEN !Concentric Transversely Orthotropic or Concentric Random

XYofAxis(1,1)=0.0D0
XYofAxis(1,2)=-SQRT(4.0D0*3.141592653589793*R*R/2.0D0/SQRT(3.0D0)/Vf)
XYofAxis(2,1)=SQRT(4.0D0*3.141592653589793*R*R/2.0D0/SQRT(3.0D0)/Vf)*SQRT(3.0D0)/2.0D0
XYofAxis(2,2)=-0.5D0*SQRT(4.0D0*3.141592653589793*R*R/2.0D0/SQRT(3.0D0)/Vf)
XYofAxis(3,1)=XYofAxis(2,1)
XYofAxis(3,2)=-XYofAxis(2,2)
XYofAxis(4,1)=XYofAxis(1,1)
XYofAxis(4,2)=-XYofAxis(1,2)
XYofAxis(5,1)=-XYofAxis(3,1)
XYofAxis(5,2)=XYofAxis(3,2)
XYofAxis(6,1)=-XYofAxis(2,1)
XYofAxis(6,2)=XYofAxis(2,2)
XYofAxis(7,1)=0.0D0
XYofAxis(7,2)=0.0D0
XYofAxis(8,1)=XYofAxis(2,1)
XYofAxis(8,2)=-1.5D0*SQRT(4.0D0*3.141592653589793*R*R/2.0D0/SQRT(3.0D0)/Vf)
XYofAxis(9,1)=XYofAxis(8,1)
XYofAxis(9,2)=-XYofAxis(8,2)
XYofAxis(10,1)=-XYofAxis(9,1)
XYofAxis(10,2)=XYofAxis(9,2)
XYofAxis(11,1)=-XYofAxis(8,1)
XYofAxis(11,2)=XYofAxis(8,2)

DO i=1,NumberofGaussPoints
NN(i,1)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0-Eta1(i))*(1.0D0-Mu1(i))
NN(i,2)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0+Eta1(i))*(1.0D0-Mu1(i))
NN(i,3)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0+Eta1(i))*(1.0D0-Mu1(i))
NN(i,4)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0-Eta1(i))*(1.0D0-Mu1(i))
NN(i,5)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0-Eta1(i))*(1.0D0+Mu1(i))
NN(i,6)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0+Eta1(i))*(1.0D0+Mu1(i))
NN(i,7)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0+Eta1(i))*(1.0D0+Mu1(i))
NN(i,8)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0-Eta1(i))*(1.0D0+Mu1(i))
ENDDO

CALL ConcentricAnisotropic(Anisotropic,Nodes,Connectivities,BeginofModeIndex,EndofModeIndex,XYofAxis,MaximumRadius,CC,Sub_CC,NN,RotationMatrix,&
RotationMatrix2)
!RotationMatrix: From Cartesian to Cylindrical
!RotationMatrix2: From Cylindrical to Cartesian

!CALL RotationMatrixSUB(Anisotropic,Nodes,Connectivities,BeginofModeIndex,EndofModeIndex,XYofAxis,MaximumRadius,NN,RotationMatrix,RotationMatrix2) 
!RotationMatrix: From Cartesian to Cylindrical
!RotationMatrix2: From Cylindrical to Cartesian

ELSE

RotationMatrix=0.0D0

!DO i1=1,NumberofGaussPoints
DO i=1,EndofModeIndex-BeginofModeIndex+1
!From Cartesian to Cylindrical (From Global to Local)
RotationMatrix(i,1,1)=1.0D0
RotationMatrix(i,2,2)=1.0D0
RotationMatrix(i,3,3)=1.0D0
!From Cylindrical to Cartesian (From Local to Global)
RotationMatrix2(i,1,1)=1.0D0
RotationMatrix2(i,2,2)=1.0D0
RotationMatrix2(i,3,3)=1.0D0
ENDDO
!ENDDO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




DO i=1,NumberofGaussPoints
DO k=1,NumberofElements

DO i1=1,24
DO i2=1,6
DO i3=1,6
DO i4=1,24
KCC(k,i1,i4)=KCC(k,i1,i4)+&
BBT(i,k,i1,i2)*Sub_CC(i,k,i2,i3)*BB(i,k,i3,i4)*Vol(i,k)*W(i)
ENDDO ! For i1
ENDDO ! For i2
ENDDO ! For i3 
ENDDO ! For i4

DO i1=1,24
DO i2=1,6
DO i3=1,6
DO i4=1,9
KCI(k,i1,i4)=KCI(k,i1,i4)+&
BBT(i,k,i1,i2)*Sub_CC(i,k,i2,i3)*BBIBar(i,k,i3,i4)*Vol(i,k)*W(i)
ENDDO ! For i1
ENDDO ! For i2
ENDDO ! For i3 
ENDDO ! For i4

DO i1=1,9
DO i2=1,6
DO i3=1,6
DO i4=1,24
KIC(k,i1,i4)=KIC(k,i1,i4)+&
BBIBarT(i,k,i1,i2)*Sub_CC(i,k,i2,i3)*BB(i,k,i3,i4)*Vol(i,k)*W(i)
ENDDO ! For i1
ENDDO ! For i2
ENDDO ! For i3 
ENDDO ! For i4

DO i1=1,9
DO i2=1,6
DO i3=1,6
DO i4=1,9
KII(k,i1,i4)=KII(k,i1,i4)+&
BBIBarT(i,k,i1,i2)*Sub_CC(i,k,i2,i3)*BBIBar(i,k,i3,i4)*Vol(i,k)*W(i)
ENDDO ! For i1
ENDDO ! For i2
ENDDO ! For i3 
ENDDO ! For i4

ENDDO
ENDDO

DO k=1,NumberofElements

CALL INV(k,9,KII(k,:,:),KII_Inverse(k,:,:))

!DO i1=1,9
!DO i2=1,9
!DO i3=1,9
!
!KII2(k,i1,i3)=KII2(k,i1,i3)+KII(k,i1,i2)*KII_Inverse(k,i2,i3)
!
!ENDDO
!ENDDO
!ENDDO
!
!DO i1=1,9
!DO i3=1,9
!IF(KII2(k,i1,i3) /= 0.0D0) THEN
!WRITE(*,30) k,i1,i3,KII2(k,i1,i3)
!ENDIF
!ENDDO
!ENDDO
!
!30 FORMAT(I,1X,I,1X,I,1X,E30.15,/)

DO i1=1,24
DO i2=1,9
DO i3=1,9
DO i4=1,24
KComponent(k,i1,i4)=KComponent(k,i1,i4)-KCI(k,i1,i2)*KII_Inverse(k,i2,i3)*KIC(k,i3,i4)
ENDDO
ENDDO
ENDDO
ENDDO

DO i1=1,24
DO i2=1,24
KComponent(k,i1,i2)=KCC(k,i1,i2)  !+KComponent(k,i1,i2)
ENDDO
ENDDO

ENDDO


RETURN

END SUBROUTINE Components

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE SixFaces(Nodes,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes)
IMPLICIT NONE
INTEGER :: IN6=10,IN7=11
INTEGER :: i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(OUT), DIMENSION(NumberofPositiveXNodes) :: PositiveXNodes
INTEGER, INTENT(OUT), DIMENSION(NumberofPositiveYNodes) :: PositiveYNodes
INTEGER, INTENT(OUT), DIMENSION(NumberofPositiveZNodes) :: PositiveZNodes
INTEGER, INTENT(OUT), DIMENSION(NumberofNegativeXNodes) :: NegativeXNodes
INTEGER, INTENT(OUT), DIMENSION(NumberofNegativeYNodes) :: NegativeYNodes
INTEGER, INTENT(OUT), DIMENSION(NumberofNegativeZNodes) :: NegativeZNodes


DOUBLE PRECISION, INTENT(INOUT), DIMENSION(NumberofNodes,3) :: Nodes
INTEGER :: Intermediate

OPEN (IN6, file = '3D_Output_Micro_F/3D_Input/Six_faces.txt',status ='unknown')

READ (IN6,*) (PositiveXNodes(i),i=1,NumberofPositiveXNodes)
READ (IN6,*) (PositiveYNodes(i),i=1,NumberofPositiveYNodes)
READ (IN6,*) (PositiveZNodes(i),i=1,NumberofPositiveZNodes)
READ (IN6,*) (NegativeXNodes(i),i=1,NumberofNegativeXNodes)
READ (IN6,*) (NegativeYNodes(i),i=1,NumberofNegativeYNodes)
READ (IN6,*) (NegativeZNodes(i),i=1,NumberofNegativeZNodes)

DO i=1,NumberofPositiveXNodes-1
DO j=i+1,NumberofPositiveXNodes

IF ( Nodes(PositiveXNodes(i),3) > Nodes(PositiveXNodes(j),3)+0.0D0 ) THEN
Intermediate=PositiveXNodes(i)
PositiveXNodes(i)=PositiveXNodes(j)
PositiveXNodes(j)=Intermediate

ELSEIF ( ABS( Nodes(PositiveXNodes(i),3) - Nodes(PositiveXNodes(j),3) ) ==0.0D0 ) THEN

IF ( Nodes(PositiveXNodes(i),4) > Nodes(PositiveXNodes(j),4)+0.0D0 ) THEN
Intermediate=PositiveXNodes(i)
PositiveXNodes(i)=PositiveXNodes(j)
PositiveXNodes(j)=Intermediate
ENDIF

ENDIF

ENDDO
ENDDO


DO i=1,NumberofPositiveYNodes-1
DO j=i+1,NumberofPositiveYNodes

IF ( Nodes(PositiveYNodes(i),2) > Nodes(PositiveYNodes(j),2)+0.0D0 ) THEN
Intermediate=PositiveYNodes(i)
PositiveYNodes(i)=PositiveYNodes(j)
PositiveYNodes(j)=Intermediate

ELSEIF ( ABS( Nodes(PositiveYNodes(i),2) - Nodes(PositiveYNodes(j),2) ) == 0.0D0 ) THEN

IF ( Nodes(PositiveYNodes(i),4) > Nodes(PositiveYNodes(j),4)+0.0D0 ) THEN
Intermediate=PositiveYNodes(i)
PositiveYNodes(i)=PositiveYNodes(j)
PositiveYNodes(j)=Intermediate
ENDIF

ENDIF

ENDDO
ENDDO


DO i=1,NumberofPositiveZNodes-1
DO j=i+1,NumberofPositiveZNodes

IF ( Nodes(PositiveZNodes(i),2) > Nodes(PositiveZNodes(j),2) +0.0D0 ) THEN
Intermediate=PositiveZNodes(i)
PositiveZNodes(i)=PositiveZNodes(j)
PositiveZNodes(j)=Intermediate

ELSEIF ( ABS( Nodes(PositiveZNodes(i),2) - Nodes(PositiveZNodes(j),2) ) ==0.0D0 ) THEN

IF ( Nodes(PositiveZNodes(i),3) > Nodes(PositiveZNodes(j),3)+0.0D0 ) THEN
Intermediate=PositiveZNodes(i)
PositiveZNodes(i)=PositiveZNodes(j)
PositiveZNodes(j)=Intermediate
ENDIF

ENDIF

ENDDO
ENDDO




DO i=1,NumberofNegativeXNodes-1
DO j=i+1,NumberofNegativeXNodes

IF ( Nodes(NegativeXNodes(i),3) > Nodes(NegativeXNodes(j),3)+0.0D0 ) THEN
Intermediate=NegativeXNodes(i)
NegativeXNodes(i)=NegativeXNodes(j)
NegativeXNodes(j)=Intermediate

ELSEIF ( ABS( Nodes(NegativeXNodes(i),3) - Nodes(NegativeXNodes(j),3) ) ==0.0D0 ) THEN

IF ( Nodes(NegativeXNodes(i),4) > Nodes(NegativeXNodes(j),4)+0.0D0 ) THEN
Intermediate=NegativeXNodes(i)
NegativeXNodes(i)=NegativeXNodes(j)
NegativeXNodes(j)=Intermediate
ENDIF

ENDIF

ENDDO
ENDDO


DO i=1,NumberofNegativeYNodes-1
DO j=i+1,NumberofNegativeYNodes

IF ( Nodes(NegativeYNodes(i),2) > Nodes(NegativeYNodes(j),2)+0.0D0 ) THEN
Intermediate=NegativeYNodes(i)
NegativeYNodes(i)=NegativeYNodes(j)
NegativeYNodes(j)=Intermediate

ELSEIF ( ABS( Nodes(NegativeYNodes(i),2) - Nodes(NegativeYNodes(j),2) ) ==0.0D0 ) THEN

IF ( Nodes(NegativeYNodes(i),4) > Nodes(NegativeYNodes(j),4)+0.0D0 ) THEN
Intermediate=NegativeYNodes(i)
NegativeYNodes(i)=NegativeYNodes(j)
NegativeYNodes(j)=Intermediate
ENDIF

ENDIF

ENDDO
ENDDO


DO i=1,NumberofNegativeZNodes-1
DO j=i+1,NumberofNegativeZNodes

IF ( Nodes(NegativeZNodes(i),2) > Nodes(NegativeZNodes(j),2)+0.0D0 ) THEN
Intermediate=NegativeZNodes(i)
NegativeZNodes(i)=NegativeZNodes(j)
NegativeZNodes(j)=Intermediate

ELSEIF ( ABS( Nodes(NegativeZNodes(i),2) - Nodes(NegativeZNodes(j),2) ) ==0.0D0 ) THEN

IF ( Nodes(NegativeZNodes(i),3) > Nodes(NegativeZNodes(j),3)+0.0D0 ) THEN
Intermediate=NegativeZNodes(i)
NegativeZNodes(i)=NegativeZNodes(j)
NegativeZNodes(j)=Intermediate
ENDIF

ENDIF

ENDDO
ENDDO

CLOSE(IN6)

OPEN(IN7,FILE="3D_Output_Micro_F/3D_Output/CheckSixFaces.txt",STATUS="UNKNOWN")

WRITE(IN7,*) "===================X===================="

DO i=1,NumberofPositiveXNodes

IF ( Nodes(PositiveXNodes(i),3)-Nodes(NegativeXNodes(i),3) /=0.0D0 .OR. &
Nodes(PositiveXNodes(i),4)-Nodes(NegativeXNodes(i),4) /=0.0D0) THEN

WRITE(IN7,21) PositiveXNodes(i),Nodes(PositiveXNodes(i),3),Nodes(PositiveXNodes(i),4),&
NegativeXNodes(i),Nodes(NegativeXNodes(i),3),Nodes(NegativeXNodes(i),4),&
Nodes(PositiveXNodes(i),3)-Nodes(NegativeXNodes(i),3),&
Nodes(PositiveXNodes(i),4)-Nodes(NegativeXNodes(i),4)
ENDIF

ENDDO

WRITE(IN7,*) "===================Y===================="

DO i=1,NumberofPositiveYNodes

IF ( Nodes(PositiveYNodes(i),2)-Nodes(NegativeYNodes(i),2) /=0.0D0 .OR. &
Nodes(PositiveYNodes(i),4)-Nodes(NegativeYNodes(i),4) /=0.0D0) THEN

WRITE(IN7,21) PositiveYNodes(i),Nodes(PositiveYNodes(i),2),Nodes(PositiveYNodes(i),4),&
NegativeYNodes(i),Nodes(NegativeYNodes(i),2),Nodes(NegativeYNodes(i),4),&
Nodes(PositiveYNodes(i),2)-Nodes(NegativeYNodes(i),2),&
Nodes(PositiveYNodes(i),4)-Nodes(NegativeYNodes(i),4)
ENDIF

ENDDO

WRITE(IN7,*) "===================Z===================="

DO i=1,NumberofPositiveZNodes

IF ( Nodes(PositiveZNodes(i),2)-Nodes(NegativeZNodes(i),2) /=0.0D0 .OR. &
Nodes(PositiveZNodes(i),3)-Nodes(NegativeZNodes(i),3) /=0.0D0) THEN

WRITE(IN7,21) PositiveZNodes(i),Nodes(PositiveZNodes(i),2),Nodes(PositiveZNodes(i),3),&
NegativeZNodes(i),Nodes(NegativeZNodes(i),2),Nodes(NegativeZNodes(i),3),&
Nodes(PositiveZNodes(i),2)-Nodes(NegativeZNodes(i),2),&
Nodes(PositiveZNodes(i),3)-Nodes(NegativeZNodes(i),3)
ENDIF

ENDDO

21 FORMAT(I10,1X,E15.7,1X,E15.7,1X,I10,1X,E15.7,1X,E15.8,1X,E15.7,1X,E15.8,/)


CLOSE(IN7)

DO i=1,NumberofPositiveXNodes

Nodes(PositiveXNodes(i),3)=Nodes(NegativeXNodes(i),3)
Nodes(PositiveXNodes(i),4)=Nodes(NegativeXNodes(i),4)

ENDDO

DO i=1,NumberofPositiveYNodes

Nodes(PositiveYNodes(i),2)=Nodes(NegativeYNodes(i),2)
Nodes(PositiveYNodes(i),4)=Nodes(NegativeYNodes(i),4)

ENDDO

DO i=1,NumberofPositiveZNodes

Nodes(PositiveZNodes(i),2)=Nodes(NegativeZNodes(i),2)
Nodes(PositiveZNodes(i),3)=Nodes(NegativeZNodes(i),3)

ENDDO

RETURN

END SUBROUTINE SixFaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadMechanicalProperties (CC)
IMPLICIT NONE
INTEGER :: IN5=9
INTEGER :: i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

DOUBLE PRECISION, DIMENSION(NumberofMaterialModes,10) :: MechanicalMaterialProperties

DOUBLE PRECISION :: EE, nu, kappaa12,muu12,EE33,muu23,nuu23,kappatz,mutz,Err,murt,nurt,E11,E22,&
E33,nu12,nu13,nu23,G23,G13,G12,nu21,nu31,nu32

DOUBLE PRECISION :: C11,C22,C12,C66

DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofMaterialModes,6,6) :: CC

MechanicalMaterialProperties=0.0D0
CC=0.0D0

OPEN(IN5,file = '3D_Output_Micro_F/3D_Input/Mechanical_material_properties.txt',status ='unknown')

DO i=1,NumberofMaterialModes

READ(IN5,*) MechanicalMaterialProperties(i,1)

IF (MechanicalMaterialProperties(i,1) == 1) THEN
BACKSPACE IN5
READ(IN5,*) MechanicalMaterialProperties(i,1),MechanicalMaterialProperties(i,2),MechanicalMaterialProperties(i,3)

ELSEIF (MechanicalMaterialProperties(i,1) == 2 .OR. MechanicalMaterialProperties(i,1) == 3) THEN
BACKSPACE IN5
READ(IN5,*) MechanicalMaterialProperties(i,1),MechanicalMaterialProperties(i,2),MechanicalMaterialProperties(i,3),&
MechanicalMaterialProperties(i,4),MechanicalMaterialProperties(i,5),MechanicalMaterialProperties(i,6)

ELSEIF (MechanicalMaterialProperties(i,1) == 4 ) THEN
BACKSPACE IN5
READ(IN5,*) MechanicalMaterialProperties(i,1),MechanicalMaterialProperties(i,2),MechanicalMaterialProperties(i,3),&
MechanicalMaterialProperties(i,4),MechanicalMaterialProperties(i,5),MechanicalMaterialProperties(i,6),&
MechanicalMaterialProperties(i,7),MechanicalMaterialProperties(i,8),MechanicalMaterialProperties(i,9),&
MechanicalMaterialProperties(i,10)

ELSEIF (MechanicalMaterialProperties(i,1) == 5 ) THEN
BACKSPACE IN5
READ(IN5,*) MechanicalMaterialProperties(i,1),MechanicalMaterialProperties(i,2),MechanicalMaterialProperties(i,3),&
MechanicalMaterialProperties(i,4),MechanicalMaterialProperties(i,5)

ENDIF

ENDDO


DO i=1,NumberofMaterialModes

IF (MechanicalMaterialProperties(i,1) == 1) THEN

EE=MechanicalMaterialProperties(i,2)
nu=MechanicalMaterialProperties(i,3)

CC(i,1,1)=EE/(1.0D0+nu)*(1.0D0-nu)/(1.0D0-2.0D0*nu)
CC(i,1,2)=EE/(1.0D0+nu)*nu/(1.0D0-2.0D0*nu)
CC(i,1,3)=EE/(1.0D0+nu)*nu/(1.0D0-2.0D0*nu)
CC(i,2,1)=CC(i,1,2)
CC(i,2,2)=CC(i,1,1)
CC(i,2,3)=CC(i,1,2)
CC(i,3,1)=CC(i,1,3)
CC(i,3,2)=CC(i,1,2)
CC(i,3,3)=CC(i,1,1)
CC(i,4,4)=EE/2.0D0/(1.0D0+nu)
CC(i,5,5)=EE/2.0D0/(1.0D0+nu)
CC(i,6,6)=EE/2.0D0/(1.0D0+nu)


ELSEIF (MechanicalMaterialProperties(i,1) == 2) THEN

kappaa12=MechanicalMaterialProperties(i,2)
muu12=MechanicalMaterialProperties(i,3)
EE33=MechanicalMaterialProperties(i,4)
muu23=MechanicalMaterialProperties(i,5)
nuu23=MechanicalMaterialProperties(i,6)

CC(i,1,1)=muu12+kappaa12
CC(i,1,2)=-muu12+kappaa12
CC(i,1,3)=2.0D0*kappaa12*nuu23
CC(i,2,1)=CC(i,1,2)
CC(i,2,2)=CC(i,1,1)
CC(i,2,3)=CC(i,1,3)
CC(i,3,1)=CC(i,1,3)
CC(i,3,2)=CC(i,2,3)
CC(i,3,3)=EE33+4.0D0*nuu23*nuu23*kappaa12
CC(i,4,4)=muu23
CC(i,5,5)=muu23
CC(i,6,6)=muu12

ELSEIF (MechanicalMaterialProperties(i,1) == 3) THEN

kappatz=MechanicalMaterialProperties(i,2)
mutz=MechanicalMaterialProperties(i,3)
Err=MechanicalMaterialProperties(i,4)
murt=MechanicalMaterialProperties(i,5)
nurt=MechanicalMaterialProperties(i,6)

CC(i,1,1)=Err+4.0D0*nurt**2*kappatz
CC(i,1,2)=2.0D0*kappatz*nurt
CC(i,1,3)=2.0D0*kappatz*nurt
CC(i,2,1)=CC(i,1,2)
CC(i,2,2)=mutz+kappatz
CC(i,2,3)=-mutz+kappatz
CC(i,3,1)=CC(i,1,3)
CC(i,3,2)=CC(i,2,3)
CC(i,3,3)=mutz+kappatz
CC(i,4,4)=mutz
CC(i,5,5)=murt
CC(i,6,6)=murt

ELSEIF (MechanicalMaterialProperties(i,1) == 4 ) THEN

E11=MechanicalMaterialProperties(i,2)
E22=MechanicalMaterialProperties(i,3)
E33=MechanicalMaterialProperties(i,4)
nu12=MechanicalMaterialProperties(i,5)
nu13=MechanicalMaterialProperties(i,6)
nu23=MechanicalMaterialProperties(i,7)
G23=MechanicalMaterialProperties(i,8)
G13=MechanicalMaterialProperties(i,9)
G12=MechanicalMaterialProperties(i,10)

nu21=nu12/E11*E22
nu31=nu13/E11*E33
nu32=nu23/E22*E33

CC(i,1,1)=(E11*(-1.0D0 + nu23*nu32))/(-1.0D0 + nu12*(nu21 + nu23*nu31) + nu23*nu32 + nu13*(nu31 + nu21*nu32))
CC(i,1,2)=-((E22*(nu12 + nu13*nu32))/(-1.0D0 + nu12*(nu21 + nu23*nu31) + nu23*nu32 + nu13*(nu31 + nu21*nu32)))
CC(i,1,3)=-((E33*(nu13 + nu12*nu23))/(-1.0D0 + nu12*(nu21 + nu23*nu31) + nu23*nu32 + nu13*(nu31 + nu21*nu32)))
CC(i,2,2)=(E22*(-1.0D0 + nu13*nu31))/(-1.0D0 + nu12*(nu21 + nu23*nu31) + nu23*nu32 + nu13*(nu31 + nu21*nu32))
CC(i,2,3)=-((E33*(nu13*nu21 + nu23))/(-1.0D0 + nu12*(nu21 + nu23*nu31) + nu23*nu32 + nu13*(nu31 + nu21*nu32)))
CC(i,3,3)=(E33*(-1.0D0 + nu12*nu21))/(-1.0D0 + nu12*(nu21 + nu23*nu31) + nu23*nu32 + nu13*(nu31 + nu21*nu32))
CC(i,2,1)=CC(i,1,2)
CC(i,3,1)=CC(i,1,3)
CC(i,3,2)=CC(i,2,3)
CC(i,4,4)=G23
CC(i,5,5)=G13
CC(i,6,6)=G12

ELSEIF (MechanicalMaterialProperties(i,1) == 5 ) THEN

C11=MechanicalMaterialProperties(i,2)
C22=MechanicalMaterialProperties(i,3)
C12=MechanicalMaterialProperties(i,4)
C66=MechanicalMaterialProperties(i,5)

CC(i,1,1)=C11
CC(i,1,2)=C12
CC(i,1,3)=0.0D0
CC(i,2,2)=C22
CC(i,2,3)=0.0D0
CC(i,3,3)=0.0D0
CC(i,2,1)=C12
CC(i,3,1)=0.0D0
CC(i,3,2)=0.0D0
CC(i,4,4)=0.0D0
CC(i,5,5)=0.0D0
CC(i,6,6)=C66

ENDIF

ENDDO

CLOSE(IN5)


RETURN

END SUBROUTINE ReadMechanicalProperties

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadBoundaryConditions (BoundaryIndex,Boundaries)
IMPLICIT NONE
INTEGER :: IN4=8
INTEGER :: i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofBoundaryNodes,TotalTimeSteps+2) :: Boundaries
INTEGER, INTENT(OUT), DIMENSION(NumberofBoundaryNodes,2) :: BoundaryIndex

OPEN(IN4,file = '3D_Output_Micro_F/3D_Input/Boundaries.txt',status ='unknown')

!READ (IN4,*) ((Boundaries(i,j),j=1,TotalTimeSteps+2),i=1,NumberofBoundaryNodes)

DO i=1,NumberofBoundaryNodes
READ (IN4,*) (Boundaries(i,j),j=1,TotalTimeSteps+2)
ENDDO

DO i=1,NumberofBoundaryNodes

BoundaryIndex(i,1)=Boundaries(i,1)
BoundaryIndex(i,2)=Boundaries(i,2)

ENDDO

CLOSE(IN4)

RETURN
END SUBROUTINE ReadBoundaryConditions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadConnectivities (Connectivities)
IMPLICIT NONE
INTEGER :: IN3=7
INTEGER :: i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(OUT), DIMENSION (NumberofElements,17) :: Connectivities

OPEN(IN3,file = '3D_Output_Micro_F/3D_Input/Connectivities_New.txt',status ='unknown')

READ (IN3,*) ((Connectivities(i,j),j=1,17),i=1,NumberofElements)

CLOSE(IN3)

RETURN

END SUBROUTINE ReadConnectivities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadNodes (NodeIndex,Nodes)
IMPLICIT NONE
INTEGER :: IN2=6
INTEGER :: i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofNodes,4) :: Nodes
INTEGER, INTENT(OUT), DIMENSION(NumberofNodes) :: NodeIndex

OPEN(IN2,file = '3D_Output_Micro_F/3D_Input/The_Nodes.txt',status ='unknown')
READ (IN2,*) ((Nodes(i,j),j=1,4),i=1,NumberofNodes)
DO i=1,NumberofNodes
NodeIndex(i)=Nodes(i,1)
ENDDO
CLOSE(IN2)
RETURN
END SUBROUTINE ReadNodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Macro_Mechanical(k,Nodes,NodeIndex,Connectivities,DNDX,DNDY,DNDZ,Sub_CC,CC,Vol,Boundaries,BoundaryIndex,NumberofDE,CurrentNodes,&
KIC,KII,K1GRR,K1GRCDD,K1GRCII,F1Globall,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,&
NegativeYNodes,NegativeZNodes,TestIndex1,TestIndex2,BeginofModeIndex,EndofModeIndex,Anisotropic,RotationMatrix,Sub_StrainsXX_R,&
Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXZ1_R,Sub_StrainsYZ1_R,Sub_StrainsXY2_R,Sub_StrainsXZ2_R,Sub_StrainsYZ2_R,Highest_GRCD)
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER :: i,j,m
INTEGER, INTENT(IN) :: k
INTEGER, INTENT(IN) :: NumberofDE
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofNodes,3) :: Nodes
INTEGER, INTENT(IN), DIMENSION(NumberofElements,17) :: Connectivities
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveXNodes) :: PositiveXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveYNodes) :: PositiveYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofPositiveZNodes) :: PositiveZNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeXNodes) :: NegativeXNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeYNodes) :: NegativeYNodes
INTEGER, INTENT(IN), DIMENSION(NumberofNegativeZNodes) :: NegativeZNodes

DOUBLE PRECISION :: MacroStrainsXX
DOUBLE PRECISION :: MacroStrainsYY
DOUBLE PRECISION :: MacroStrainsZZ
DOUBLE PRECISION :: MacroStrainsXY1
DOUBLE PRECISION :: MacroStrainsXY2
DOUBLE PRECISION :: MacroStrainsXZ1
DOUBLE PRECISION :: MacroStrainsXZ2
DOUBLE PRECISION :: MacroStrainsYZ1
DOUBLE PRECISION :: MacroStrainsYZ2
DOUBLE PRECISION, DIMENSION(3*(NumberofPositiveXNodes+NumberofPositiveYNodes+NumberofPositiveZNodes),8) :: One_Step_DEs
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofBoundaryNodes,TotalTimeSteps+2) :: Boundaries
INTEGER, INTENT(IN), DIMENSION(NumberofBoundaryNodes,2) :: BoundaryIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements,9,24)  :: KIC
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofElements,9,9)   :: KII
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec) :: K1GRR
INTEGER, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE) :: K1GRCDD
INTEGER, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE,MaxDimenofInverseConnec) :: K1GRCII
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*NumberofNodes+NumberofDE) :: F1Globall
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofGaussPoints,NumberofElements) :: Vol
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofMaterialModes,6,6) :: CC
DOUBLE PRECISION, INTENT(IN), DIMENSION (NumberofGaussPoints,NumberofElements,6,6) :: Sub_CC
INTEGER, INTENT(IN), DIMENSION(NumberofNodes) :: NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofGaussPoints,NumberofElements,8+3) ::  DNDX,DNDY,DNDZ
INTEGER, INTENT(IN) :: TestIndex1,TestIndex2,BeginofModeIndex,EndofModeIndex,Anisotropic
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofMaterialModes,3,3) :: RotationMatrix
INTEGER, INTENT(OUT) :: Highest_GRCD
DOUBLE PRECISION, DIMENSION (NumberofElements) :: Vol_Ele
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StrainsXX
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StrainsYY
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StrainsZZ
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StrainsXY
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StrainsXZ
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StrainsYZ
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StrainsVonMises
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StressesXX
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StressesYY
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StressesZZ
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StressesXY
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StressesXZ
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StressesYZ
DOUBLE PRECISION, DIMENSION (NumberofElements) :: StressesVonMises
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofNodes,4)  :: CurrentNodes 
DOUBLE PRECISION, DIMENSION (NumberofElements,8) :: Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,Sub_StrainsXZ1,Sub_StrainsYZ1,&
Sub_StrainsXY2,Sub_StrainsXZ2,Sub_StrainsYZ2,Sub_StressesXX,Sub_StressesYY,Sub_StressesZZ,Sub_StressesXY,Sub_StressesXZ,Sub_StressesYZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (NumberofElements,8) :: Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,&
Sub_StrainsXZ1_R,Sub_StrainsYZ1_R,Sub_StrainsXY2_R,Sub_StrainsXZ2_R,Sub_StrainsYZ2_R
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: Dis1
DOUBLE PRECISION :: Kappa12,Mu12,Mu23,E33,Nu32,E33_Real,C11,C12,C33

DO m=TestIndex1,TestIndex2

IF (m==1) THEN ! For Mu12

MacroStrainsXX=0.0D0
MacroStrainsYY=0.0D0
MacroStrainsZZ=0.0D0
MacroStrainsXY1=0.005D0*k
MacroStrainsXZ1=0.0D0
MacroStrainsYZ1=0.0D0
MacroStrainsXY2=0.005D0*k
MacroStrainsXZ2=0.0D0
MacroStrainsYZ2=0.0D0

ELSEIF (m==2) THEN ! For Mu23

MacroStrainsXX=0.0D0
MacroStrainsYY=0.0D0
MacroStrainsZZ=0.0D0
MacroStrainsXY1=0.0D0
MacroStrainsXZ1=0.0D0
MacroStrainsYZ1=0.005D0*k
MacroStrainsXY2=0.0D0
MacroStrainsXZ2=0.0D0
MacroStrainsYZ2=0.005D0*k

ELSEIF (m==3) THEN ! For C11 and further for Kappa12

MacroStrainsXX=0.005D0*k
MacroStrainsYY=0.0D0
MacroStrainsZZ=0.0D0
MacroStrainsXY1=0.0D0
MacroStrainsXZ1=0.0D0
MacroStrainsYZ1=0.0D0
MacroStrainsXY2=0.0D0
MacroStrainsXZ2=0.0D0
MacroStrainsYZ2=0.0D0

ELSEIF (m==4) THEN ! For C33

MacroStrainsXX=0.0D0
MacroStrainsYY=0.0D0
MacroStrainsZZ=0.005D0*k
MacroStrainsXY1=0.0D0
MacroStrainsXZ1=0.0D0
MacroStrainsYZ1=0.0D0
MacroStrainsXY2=0.0D0
MacroStrainsXZ2=0.0D0
MacroStrainsYZ2=0.0D0

ELSEIF (m==5) THEN ! For C12 and further for E33 and Nu32

MacroStrainsXX=0.005D0*k
MacroStrainsYY=0.0D0
MacroStrainsZZ=0.005D0*k
MacroStrainsXY1=0.0D0
MacroStrainsXZ1=0.0D0
MacroStrainsYZ1=0.0D0
MacroStrainsXY2=0.0D0
MacroStrainsXZ2=0.0D0
MacroStrainsYZ2=0.0D0


ELSEIF (m==6) THEN

CALL Input_Strains(k,MacroStrainsXX,MacroStrainsYY,MacroStrainsZZ,MacroStrainsXY1,MacroStrainsXY2,MacroStrainsXZ1,MacroStrainsXZ2,&
MacroStrainsYZ1,MacroStrainsYZ2)

ENDIF

CALL OneStepDEs(Nodes,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,NegativeYNodes,NegativeZNodes,&
MacroStrainsXX,MacroStrainsYY,MacroStrainsZZ,MacroStrainsXY1,MacroStrainsXY2,MacroStrainsXZ1,MacroStrainsXZ2,MacroStrainsYZ1,MacroStrainsYZ2,One_Step_DEs)

ALLOCATE(Dis1(3*NumberofNodes+NumberofDE))
IF (DisplacementEquationFlag==1) THEN
CALL Matrix_DE_Reduced(k,NumberofDE,One_Step_DEs,BoundaryIndex,Boundaries,K1GRR,K1GRCDD,K1GRCII,F1Globall,Connectivities,&
Dis1,Highest_GRCD)
ELSEIF(DisplacementEquationFlag==0) THEN
CALL Matrix_Reduced(k,NumberofDE,One_Step_DEs,BoundaryIndex,Boundaries,K1GRR,K1GRCDD,K1GRCII,F1Globall,Connectivities,&
Dis1,Highest_GRCD)
ENDIF

CALL StrainsandStresses (Vol,Vol_Ele,CC,Sub_CC,Connectivities,Nodes,NodeIndex,DNDX,DNDY,DNDZ,KIC,KII,&
Dis1,StrainsXX,StrainsYY,StrainsZZ,StrainsXY,StrainsXZ,StrainsYZ,StrainsVonMises,StressesXX,StressesYY,StressesZZ,StressesXY,StressesXZ,StressesYZ,&
StressesVonMises,Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,Sub_StrainsXZ1,Sub_StrainsYZ1,&
Sub_StrainsXY2,Sub_StrainsXZ2,Sub_StrainsYZ2,Sub_StressesXX,Sub_StressesYY,Sub_StressesZZ,Sub_StressesXY,Sub_StressesXZ,Sub_StressesYZ,CurrentNodes)

IF (k==1) THEN
CALL EffectiveProperties(k,m,MacroStrainsXX,MacroStrainsYY,MacroStrainsZZ,MacroStrainsXY1,MacroStrainsXZ1,MacroStrainsYZ1,MacroStrainsXY2,&
MacroStrainsXZ2,MacroStrainsYZ2,Vol,Vol_Ele,Nodes,Connectivities,CC,PositiveXNodes,PositiveYNodes,PositiveZNodes,NegativeXNodes,&
NegativeYNodes,NegativeZNodes,StrainsXX,StrainsYY,StrainsZZ,StrainsXY,StrainsXZ,StrainsYZ,StressesXX,StressesYY,StressesZZ,&
StressesXY,StressesXZ,StressesYZ,Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,Sub_StrainsXZ1,Sub_StrainsYZ1,&
Sub_StrainsXY2,Sub_StrainsXZ2,Sub_StrainsYZ2,Sub_StressesXX,Sub_StressesYY,Sub_StressesZZ,Sub_StressesXY,Sub_StressesXZ,Sub_StressesYZ,&
Kappa12,Mu12,Mu23,E33,Nu32,E33_Real,C11,C12,C33)
ENDIF

CALL FeedForTecplot (k,m,Connectivities,Nodes,NodeIndex,Dis1,StrainsXX,StrainsYY,StrainsZZ,StrainsXY,StrainsXZ,StrainsYZ,StrainsVonMises,StressesXX,&
StressesYY,StressesZZ,StressesXY,StressesXZ,StressesYZ,StressesVonMises,CurrentNodes)

DEALLOCATE(Dis1)

ENDDO !For m

IF (Anisotropic == 1 .OR. Anisotropic == 2) THEN
CALL RotatedStrains(Connectivities,BeginofModeIndex,EndofModeIndex,RotationMatrix,Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,&
Sub_StrainsXZ1,Sub_StrainsYZ1,Sub_StrainsXY2,Sub_StrainsXZ2,Sub_StrainsYZ2,Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,&
Sub_StrainsXZ1_R,Sub_StrainsYZ1_R,Sub_StrainsXY2_R,Sub_StrainsXZ2_R,Sub_StrainsYZ2_R,Vol)
ELSE
Sub_StrainsXX_R=Sub_StrainsXX
Sub_StrainsYY_R=Sub_StrainsYY
Sub_StrainsZZ_R=Sub_StrainsZZ
Sub_StrainsXY1_R=Sub_StrainsXY1
Sub_StrainsXY2_R=Sub_StrainsXY2
Sub_StrainsXZ1_R=Sub_StrainsXZ1
Sub_StrainsXZ2_R=Sub_StrainsXZ2
Sub_StrainsYZ1_R=Sub_StrainsYZ1
Sub_StrainsYZ2_R=Sub_StrainsYZ2
ENDIF

END SUBROUTINE Macro_Mechanical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_Matrix_DE_Reduced(k,ME_NumberofDE,ME_One_Step_DEs,ME_BoundaryIndex,ME_Boundaries,ME_K1GRR,ME_K1GRCDD,ME_K1GRCII,ME_F1Globall,ME_Connectivities,&
ME_Dis1,ME_Highest_GRCD)
IMPLICIT NONE
INTEGER :: i,j,m,n,i1,j1,k1,k2,k3
INTEGER, INTENT(IN) :: k

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

INTEGER, INTENT(IN) :: ME_NumberofDE
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec) :: ME_K1GRR
INTEGER, INTENT(IN), DIMENSION (M_NumberofNodes+ME_NumberofDE) :: ME_K1GRCDD
INTEGER, INTENT(IN), DIMENSION (M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec) :: ME_K1GRCII
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes+ME_NumberofDE) :: ME_F1Globall
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: ME_K1GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: ME_K1GRCD
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ME_K1GRCI
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: ME_F1Global

INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: ME_K2GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: ME_K2GRRD

DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofDE,6) :: ME_One_Step_DEs
DOUBLE PRECISION, INTENT(IN), DIMENSION (ME_NumberofBoundaryNodes,TotalTimeSteps+1) :: ME_Boundaries
INTEGER, INTENT(IN), DIMENSION(ME_NumberofBoundaryNodes,1) :: ME_BoundaryIndex
DOUBLE PRECISION, DIMENSION (M_NumberofNodes+ME_NumberofDE) :: ME_Dis
DOUBLE PRECISION, INTENT(OUT),DIMENSION (M_NumberofNodes+ME_NumberofDE) :: ME_Dis1
INTEGER :: ME_Dimen1
DOUBLE PRECISION, DIMENSION(M_NumberofNodes+ME_NumberofDE) :: ME_x
INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,10) :: ME_Connectivities
DOUBLE PRECISION :: ElapsedTime,Temp
INTEGER :: ME_TEMP_GRCD
INTEGER, INTENT(OUT) :: ME_Highest_GRCD
INTEGER :: ME_Temp_INT
DOUBLE PRECISION :: ME_Temp
INTEGER :: ME_DUMMY,ME_DUMMY2
INTEGER, ALLOCATABLE, DIMENSION(:) :: ME_IA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: ME_IA
INTEGER, ALLOCATABLE, DIMENSION(:) :: ME_JA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: ME_JA2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ME_A1

ALLOCATE(ME_K1GR(M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec))
ALLOCATE(ME_K1GRCD(M_NumberofNodes+ME_NumberofDE))
ALLOCATE(ME_K1GRCI(M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec))
ALLOCATE(ME_F1Global(M_NumberofNodes+ME_NumberofDE))

ME_K1GR=ME_K1GRR
ME_K1GRCD=ME_K1GRCDD
ME_K1GRCI=ME_K1GRCII
ME_F1Global=ME_F1Globall

!!!Apply P.B.C.s (Top Right and Bottom Right of F)

DO i=1,ME_NumberofDE
DO j=1,int(ME_One_Step_DEs(i,1))

ME_K1GRCD( 1*(int(ME_One_Step_DEs(i,2*j))-1)+1 ) = ME_K1GRCD( 1*(int(ME_One_Step_DEs(i,2*j))-1)+1 )+1

ME_K1GR( 1*(int(ME_One_Step_DEs(i,2*j))-1)+1,ME_K1GRCD( 1*(int(ME_One_Step_DEs(i,2*j))-1)+1 ) ) = ME_One_Step_DEs(i,2*j+1) 

ME_K1GRCI( 1*(int(ME_One_Step_DEs(i,2*j))-1)+1, ME_K1GRCD( 1*(int(ME_One_Step_DEs(i,2*j))-1)+1 ) ) = M_NumberofNodes+i

ME_F1Global(1*M_NumberofNodes+i)=ME_One_Step_DEs(i,2*int(ME_One_Step_DEs(i,1))+1+1)

ENDDO
ENDDO

!!!Apply P.B.C.s (Bottom left)
DO i=1,ME_NumberofDE
DO j=1,int(ME_One_Step_DEs(i,1))
ME_K1GRCD( M_NumberofNodes+i ) = ME_K1GRCD( M_NumberofNodes+i )+1
ME_K1GR( M_NumberofNodes+i,ME_K1GRCD( M_NumberofNodes+i ) ) = ME_One_Step_DEs(i,2*j+1) 
ME_K1GRCI( M_NumberofNodes+i, ME_K1GRCD( M_NumberofNodes+i ) ) = 1*(int(ME_One_Step_DEs(i,2*j))-1)+1
ENDDO
ENDDO

ALLOCATE(ME_K2GR(M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec,2))
ALLOCATE(ME_K2GRRD(M_NumberofNodes+ME_NumberofDE))
ME_K2GR=0
ME_K2GRRD=0

DO i=1,M_NumberofNodes+ME_NumberofDE
DO j=1,ME_K1GRCD(i)
ME_K2GRRD(ME_K1GRCI(i,j))=ME_K2GRRD(ME_K1GRCI(i,j))+1
ME_K2GR(ME_K1GRCI(i,j),ME_K2GRRD(ME_K1GRCI(i,j)),1)=i
ME_K2GR(ME_K1GRCI(i,j),ME_K2GRRD(ME_K1GRCI(i,j)),2)=j
ENDDO
ENDDO

ME_Dis=0.0D0

DO i=1,ME_NumberofBoundaryNodes
ME_Dis(1*(ME_BoundaryIndex(i,1)-1)+1)=ME_Boundaries(i,k+1)
ENDDO


DO i=1,ME_NumberofBoundaryNodes

DO m=1,ME_K2GRRD(1*(ME_BoundaryIndex(i,1)-1)+1)
ME_F1Global( ME_K2GR(1*(ME_BoundaryIndex(i,1)-1)+1,m,1) )=ME_F1Global( ME_K2GR(1*(ME_BoundaryIndex(i,1)-1)+1,m,1) )&
-ME_K1GR( ME_K2GR(1*(ME_BoundaryIndex(i,1)-1)+1,m,1), ME_K2GR(1*(ME_BoundaryIndex(i,1)-1)+1,m,2) )*ME_Dis(1*(ME_BoundaryIndex(i,1)-1)+1)
ENDDO

ME_F1Global(1*(ME_BoundaryIndex(i,1)-1)+1)=ME_Dis(1*(ME_BoundaryIndex(i,1)-1)+1)

DO m=1,ME_K2GRRD(1*(ME_BoundaryIndex(i,1)-1)+1)
ME_K1GR( ME_K2GR(1*(ME_BoundaryIndex(i,1)-1)+1,m,1), ME_K2GR(1*(ME_BoundaryIndex(i,1)-1)+1,m,2) )=0.0D0
ENDDO

DO m=1,ME_K1GRCD(1*(ME_BoundaryIndex(i,1)-1)+1)
IF ( 1*(ME_BoundaryIndex(i,1)-1)+1 /= ME_K1GRCI(1*(ME_BoundaryIndex(i,1)-1)+1,m) ) THEN
ME_K1GR( 1*(ME_BoundaryIndex(i,1)-1)+1, m )=0.0D0
ELSE
ME_K1GR( 1*(ME_BoundaryIndex(i,1)-1)+1, m )=1.0D0
ENDIF
ENDDO

ENDDO

DEALLOCATE(ME_K2GR)
DEALLOCATE(ME_K2GRRD)

ME_TEMP_GRCD=ME_K1GRCD(1)
ME_Highest_GRCD=ME_TEMP_GRCD
DO i=2,1*M_NumberofNodes+ME_NumberofDE
IF ( ME_K1GRCD(i)>ME_TEMP_GRCD ) THEN
ME_TEMP_GRCD=ME_K1GRCD(i)
ME_Highest_GRCD=ME_TEMP_GRCD
ENDIF
ENDDO

ME_DUMMY=0

DO i=1,1*M_NumberofNodes+ME_NumberofDE
DO j=1,ME_K1GRCD(i)
IF (ME_K1GR(i,j) /= 0.0D0) THEN
ME_DUMMY=ME_DUMMY+1
ELSEIF (ME_K1GR(i,j) == 0.0D0) THEN
GOTO 3
ENDIF
3 ENDDO
ENDDO

ALLOCATE(ME_IA1(ME_DUMMY))
ALLOCATE(ME_JA1(ME_DUMMY))
ALLOCATE(ME_JA2(ME_DUMMY))
ALLOCATE(ME_A1(ME_DUMMY))

ME_DUMMY=0

DO i=1,1*M_NumberofNodes+ME_NumberofDE
DO j=1,ME_K1GRCD(i)
IF (ME_K1GR(i,j) /= 0.0D0) THEN
ME_DUMMY=ME_DUMMY+1
ME_A1(ME_DUMMY)=ME_K1GR(i,j)
ME_JA1(ME_DUMMY)=ME_K1GRCI(i,j)
ME_JA2(ME_DUMMY)=i
ELSEIF (ME_K1GR(i,j) == 0.0D0) THEN
GOTO 4
ENDIF
4 ENDDO
ENDDO

!OPEN(33,FILE="3D_Output_Micro_F/3D_Output_Microscale_1/Check2D_KGR_ME.txt",STATUS="UNKNOWN")
!DO i=1,M_NumberofNodes+ME_NumberofDE
!WRITE(33,34) (ME_K1GR(i,k1),k1=1,ME_MaxDimenofInverseConnec)
!ENDDO
!34 FORMAT(<40>E15.8,/)
!WRITE(33,*) "ME_Highest_GRCD  ",ME_Highest_GRCD
!CLOSE(33)

DEALLOCATE(ME_K1GR)
DEALLOCATE(ME_K1GRCD)
DEALLOCATE(ME_K1GRCI)

ME_IA1(1)=1

ME_DUMMY2=1

DO i=2,ME_DUMMY

IF (ME_JA2(i)>ME_JA2(i-1) .AND. i>1) THEN

ME_DUMMY2=ME_DUMMY2+1

ME_IA1(ME_DUMMY2)=i

ENDIF

ENDDO

ME_DUMMY2=ME_DUMMY2+1

ME_IA1(ME_DUMMY2)=ME_DUMMY+1

ALLOCATE(ME_IA(ME_DUMMY2))

DO i=1,ME_DUMMY2
ME_IA(i)=ME_IA1(i)
ENDDO

DEALLOCATE(ME_IA1)

!OPEN(27,FILE="3D_Output_Micro_F/3D_Output_Microscale_1/A1JA1JA2_ME.txt",STATUS="UNKNOWN")
!DO i=1,ME_DUMMY
!WRITE(27,29) ME_A1(i),ME_JA1(i),ME_JA2(i)
!ENDDO
!CLOSE(27)
!
!29 FORMAT(E15.8,2X,I,2X,I,/)
!
!OPEN(39,FILE="3D_Output_Micro_F/3D_Output_Microscale_1/IA_ME.txt",STATUS="UNKNOWN")
!DO i=1,ME_DUMMY2
!WRITE(39,33) ME_IA(i)
!ENDDO
!CLOSE(39)
!
!33 FORMAT(I,/)

ME_Dimen1=1*M_NumberofNodes+ME_NumberofDE

CALL HarwellBoeing_Reduced(k,ME_Dimen1,ME_DUMMY,ME_DUMMY2,ME_A1,ME_JA1,ME_JA2,ME_IA,ME_F1Global,ME_x)

DEALLOCATE(ME_F1Global)
DEALLOCATE(ME_IA)
DEALLOCATE(ME_JA1)
DEALLOCATE(ME_JA2)
DEALLOCATE(ME_A1)

DO j=1,M_NumberofNodes
ME_Dis1(j)=ME_x(1*(j-1)+1)
!WRITE(*,*) ME_Dis1(j)
ENDDO

RETURN
END SUBROUTINE ME_Matrix_DE_Reduced

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_OneStepDEs(M_Nodes,M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes,&
Ex_Knot,Ey_Knot,Ez_Knot,ME_NumberofDE,ME_One_Step_DEs)
IMPLICIT NONE
INTEGER :: i,j
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofNodes,4) :: M_Nodes
INTEGER, INTENT(IN) :: ME_NumberofDE

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveXNodes) :: M_PositiveXNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveYNodes) :: M_PositiveYNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveZNodes) :: M_PositiveZNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeXNodes) :: M_NegativeXNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeYNodes) :: M_NegativeYNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeZNodes) :: M_NegativeZNodes

DOUBLE PRECISION, INTENT(INOUT), DIMENSION(M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+M_NumberofPositiveZNodes,6) :: ME_One_Step_DEs

DOUBLE PRECISION, INTENT(IN) :: Ex_Knot,Ey_Knot,Ez_Knot

INTEGER :: M_PosiTemp,M_NegaTemp

DO j=1,M_NumberofPositiveXNodes

M_PosiTemp=M_PositiveXNodes(j)
M_NegaTemp=M_NegativeXNodes(j)

IF ( 1*(M_PosiTemp-1)+1 < 1*(M_NegaTemp-1)+1 ) THEN
ME_One_Step_DEs(j,1)=2
ME_One_Step_DEs(j,2)=M_PositiveXNodes(j)
ME_One_Step_DEs(j,3)=1.0D0
ME_One_Step_DEs(j,4)=M_NegativeXNodes(j)
ME_One_Step_DEs(j,5)=-1.0D0
ME_One_Step_DEs(j,6)=Ex_Knot*(M_Nodes(M_PositiveXNodes(j),2)-M_Nodes(M_NegativeXNodes(j),2))
ELSE
ME_One_Step_DEs(j,1)=2
ME_One_Step_DEs(j,2)=M_NegativeXNodes(j)
ME_One_Step_DEs(j,3)=-1.0D0
ME_One_Step_DEs(j,4)=M_PositiveXNodes(j)
ME_One_Step_DEs(j,5)=1.0D0
ME_One_Step_DEs(j,6)=Ex_Knot*(M_Nodes(M_PositiveXNodes(j),2)-M_Nodes(M_NegativeXNodes(j),2))
ENDIF

ENDDO

DO j=M_NumberofPositiveXNodes+1,M_NumberofPositiveXNodes+M_NumberofPositiveYNodes

M_PosiTemp=M_PositiveYNodes(j-M_NumberofPositiveXNodes)
M_NegaTemp=M_NegativeYNodes(j-M_NumberofPositiveXNodes)

IF ( 1*(M_PosiTemp-1)+1 < 1*(M_NegaTemp-1)+1 ) THEN
ME_One_Step_DEs(j,1)=2
ME_One_Step_DEs(j,2)=M_PositiveYNodes(j-M_NumberofPositiveXNodes)
ME_One_Step_DEs(j,3)=1.0D0
ME_One_Step_DEs(j,4)=M_NegativeYNodes(j-M_NumberofPositiveXNodes)
ME_One_Step_DEs(j,5)=-1.0D0
ME_One_Step_DEs(j,6)=Ey_Knot*(M_Nodes(M_PositiveYNodes(j-M_NumberofPositiveXNodes),3)-&
M_Nodes(M_NegativeYNodes(j-M_NumberofPositiveXNodes),3))
ELSE
ME_One_Step_DEs(j,1)=2
ME_One_Step_DEs(j,2)=M_NegativeYNodes(j-M_NumberofPositiveXNodes)
ME_One_Step_DEs(j,3)=-1.0D0
ME_One_Step_DEs(j,4)=M_PositiveYNodes(j-M_NumberofPositiveXNodes)
ME_One_Step_DEs(j,5)=1.0D0
ME_One_Step_DEs(j,6)=Ey_Knot*(M_Nodes(M_PositiveYNodes(j-M_NumberofPositiveXNodes),3)-&
M_Nodes(M_NegativeYNodes(j-M_NumberofPositiveXNodes),3))
ENDIF

ENDDO

DO j=M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+1,M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+M_NumberofPositiveZNodes
M_PosiTemp=M_PositiveZNodes(j-M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)
M_NegaTemp=M_NegativeZNodes(j-M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)

IF ( 1*(M_PosiTemp-1)+1 < 1*(M_NegaTemp-1)+1 ) THEN
ME_One_Step_DEs(j,1)=2
ME_One_Step_DEs(j,2)=M_PositiveZNodes(j-M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)
ME_One_Step_DEs(j,3)=1.0D0
ME_One_Step_DEs(j,4)=M_NegativeZNodes(j-M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)
ME_One_Step_DEs(j,5)=-1.0D0
ME_One_Step_DEs(j,6)=Ez_Knot*(M_Nodes(M_PositiveZNodes(j-M_NumberofPositiveXNodes-M_NumberofPositiveYNodes),4)-&
M_Nodes(M_NegativeZNodes(j-M_NumberofPositiveXNodes-M_NumberofPositiveYNodes),4))
ELSE
ME_One_Step_DEs(j,1)=2
ME_One_Step_DEs(j,2)=M_NegativeZNodes(j-M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)
ME_One_Step_DEs(j,3)=-1.0D0
ME_One_Step_DEs(j,4)=M_PositiveZNodes(j-M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)
ME_One_Step_DEs(j,5)=1.0D0
ME_One_Step_DEs(j,6)=Ez_Knot*(M_Nodes(M_PositiveZNodes(j-M_NumberofPositiveXNodes-M_NumberofPositiveYNodes),4)-&
M_Nodes(M_NegativeZNodes(j-M_NumberofPositiveXNodes-M_NumberofPositiveYNodes),4))
ENDIF

ENDDO

RETURN
ENDSUBROUTINE ME_OneStepDEs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_KPreprations(ME_Connectivities,ME_KComponent,ME_NumberofDE,ME_K1GRR,ME_K1GRCDD,ME_K1GRCII,ME_F1Globall)
IMPLICIT NONE
INTEGER :: i,j,m,n,i1,j1,k1,k2,k3

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,10) :: ME_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements,8,8) :: ME_KComponent
INTEGER, INTENT(IN) :: ME_NumberofDE
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec) :: ME_K1GRR
INTEGER, INTENT(OUT), DIMENSION (M_NumberofNodes+ME_NumberofDE) :: ME_K1GRCDD
INTEGER, INTENT(OUT), DIMENSION (M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec) :: ME_K1GRCII
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofNodes+ME_NumberofDE) :: ME_F1Globall

ME_K1GRR=0.0D0
ME_K1GRCDD=0
ME_K1GRCII=0

ME_F1Globall=0.0D0

DO j1=1,1
DO i1=1,1
DO j=1,8
DO i=1,8
DO k1=1,M_NumberofElements

IF ( ME_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j) /= 0.0D0 ) THEN

IF ( ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)==0 ) THEN

ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,1)=ME_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)

ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)=ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)+1

ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,1)=1*(ME_Connectivities(k1,j+1)-1)+j1

GOTO 1

ELSEIF ( ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1) == 1 ) THEN

IF ( 1*(ME_Connectivities(k1,j+1)-1)+j1 == ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,1) ) THEN
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,1)=ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,1)+ME_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ELSEIF ( 1*(ME_Connectivities(k1,j+1)-1)+j1 > ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,1) ) THEN
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,2)=ME_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)=ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)+1
ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,2)=1*(ME_Connectivities(k1,j+1)-1)+j1
GOTO 1
ELSEIF ( 1*(ME_Connectivities(k1,j+1)-1)+j1 < ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,1) ) THEN
ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)=ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)+1
DO k3=2,ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k3)=ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k3-1)
ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k3)=ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,1)=ME_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,1)=1*(ME_Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF

ELSEIF ( ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1) > 1 ) THEN

IF ( 1*(ME_Connectivities(k1,j+1)-1)+j1 < ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,1) ) THEN
ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)=ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)+1
DO k3=ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1),2,-1
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k3)=ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k3-1)
ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k3)=ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,1)=ME_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,1)=1*(ME_Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF

DO k2=1,ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)-1

IF ( 1*(ME_Connectivities(k1,j+1)-1)+j1 == ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k2) ) THEN
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k2)=ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k2)+ME_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ENDIF

IF (1*(ME_Connectivities(k1,j+1)-1)+j1 > ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k2) .AND. &
1*(ME_Connectivities(k1,j+1)-1)+j1 < ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k2+1) ) THEN
ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)=ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)+1
DO k3=ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1),k2+2,-1
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k3)=ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k3-1)
ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k3)=ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k2+1)=ME_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k2+1)=1*(ME_Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF  

ENDDO  !k2

DO k2=ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1),ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)

IF ( 1*(ME_Connectivities(k1,j+1)-1)+j1 == ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k2) ) THEN
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k2)=ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k2)+ME_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ENDIF

IF ( 1*(ME_Connectivities(k1,j+1)-1)+j1 > ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k2)  ) THEN
ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)=ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)+1
DO k3=k2+1,ME_K1GRCDD(1*(ME_Connectivities(k1,i+1)-1)+i1)
ME_K1GRR(1*(ME_Connectivities(k1,i+1)-1)+i1,k3)=ME_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
ME_K1GRCII(1*(ME_Connectivities(k1,i+1)-1)+i1,k3)=1*(ME_Connectivities(k1,j+1)-1)+j1
ENDDO
GOTO 1
ENDIF  

ENDDO !For k2

1 ENDIF  ! For ( M_K1GRCD(3*(M_Connectivities(k1,i+1)-1)+i1)==0 )
 
ENDIF  ! For ( M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j) /= 0.0D0 )

ENDDO
ENDDO
ENDDO
ENDDO
ENDDO

!DO i=1,20
!DO j=1,20
!IF (ME_K1GRR(i,j) /=0.0D0) THEN
!WRITE(*,*) i,j,ME_K1GRR(i,j)
!ENDIF
!ENDDO
!ENDDO

END SUBROUTINE ME_KPreprations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_Components_Hexahedron(M_Nodes,ME_Connectivities,ME_sigmasigma,ME_gg,ME_sub_sigmasigma_updated,ME_Vol,&
ME_KComponent,ME_DNDX,ME_DNDY,ME_DNDZ,ZeroTimeStep,k,M_NumberofDE,Sub_StrainsXX_R,Sub_StrainsYY_R,&
Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R,Sub_StrainsXZ1_R,Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,&
Sub_StrainsYZ2_R,ME_TunnelingIndex,ME_TunnelingIndexDimension,ME_CNTIndex,ME_Critical_Distance,ME_Coordinates2,&
ME_Width1,ME_Width2,M_NumberofElementsPerCNT,ME_ElementsPerCNT,M_Vol,M_NN,M_Dis1,M_CurrentNodes,&
M_Sub_StrainsXX,M_Sub_StrainsYY,M_Sub_StrainsZZ,M_Sub_StrainsXY1,M_Sub_StrainsXZ1,M_Sub_StrainsYZ1,&
M_Sub_StrainsXY2,M_Sub_StrainsXZ2,M_Sub_StrainsYZ2)
IMPLICIT NONE
INTEGER :: i,j,k1,i1,i2,i3,i4

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

!The electrical tunneling control flags in the microscale
INTEGER :: ME_Tunneling
DOUBLE PRECISION :: ME_Vf
DOUBLE PRECISION :: ME_RofCNT
DOUBLE PRECISION :: ME_Lambda
DOUBLE PRECISION :: ME_rho
INTEGER :: ME_BeginningPolymer
INTEGER :: ME_EndingPolymer
INTEGER :: ME_NumberofCNTs
INTEGER :: ME_NumberofCNTs2
INTEGER :: ME_MaxNumberofElementsPerCNT
COMMON/Micro_Electrical_More_Inputs/ME_Tunneling,ME_BeginningPolymer,ME_EndingPolymer,ME_NumberofCNTs,ME_NumberofCNTs2,&
ME_MaxNumberofElementsPerCNT,ME_Vf,ME_RofCNT,ME_Lambda,ME_rho

DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements,8) :: M_Sub_StrainsXX,M_Sub_StrainsYY,M_Sub_StrainsZZ,&
M_Sub_StrainsXY1,M_Sub_StrainsXZ1,M_Sub_StrainsYZ1,M_Sub_StrainsXY2,M_Sub_StrainsXZ2,M_Sub_StrainsYZ2

DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes,4) :: M_Nodes
INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,10) :: ME_Connectivities
DOUBLE PRECISION, INTENT(INOUT),  DIMENSION (ME_NumberofMaterialModes,3,3) :: ME_sigmasigma
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofMaterialModes,6,6) :: ME_gg
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,ME_NumberofMaterialModes,6,6) :: ME_sub_gg
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,ME_NumberofMaterialModes,3,3) :: ME_sub_sigmasigma
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_rhorho
DOUBLE PRECISION, INTENT(OUT), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_sigmasigma_updated
DOUBLE PRECISION, INTENT(OUT), DIMENSION (ME_NumberofGaussPoints,M_NumberofElements) :: ME_Vol
DOUBLE PRECISION, DIMENSION (ME_NumberofGaussPoints,M_NumberofElements,3,8) :: ME_BB
DOUBLE PRECISION, DIMENSION (ME_NumberofGaussPoints,M_NumberofElements,8,3) :: ME_BBT
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements,8,8) :: ME_KComponent

DOUBLE PRECISION, DIMENSION(M_NumberofElements,8) :: ME_x,ME_y,ME_z

DOUBLE PRECISION, DIMENSION(8) :: ME_Xi,ME_Eta,ME_Mu

DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints) :: ME_Xi1,ME_Eta1,ME_Mu1,ME_W

DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_JJ,ME_InverseJJ

DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,8) :: ME_DNDXi,ME_DNDEta,ME_DNDMu
DOUBLE PRECISION, DIMENSION(8,ME_NumberofGaussPoints) :: ME_NN
DOUBLE PRECISION, INTENT(OUT), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,8) :: ME_DNDX,ME_DNDY,ME_DNDZ
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,8) :: ME_QX,ME_QY,ME_QZ

DOUBLE PRECISION :: ME_Delta_A

DOUBLE PRECISION :: ME_a11,ME_a12,ME_a13,ME_a21,ME_a22,ME_a23,ME_a31,ME_a32,ME_a33

DOUBLE PRECISION, DIMENSION(11,2) :: ME_XYofAxis

DOUBLE PRECISION :: ME_MaximumRadius

INTEGER :: ME_BeginofModeIndex,ME_EndofModeIndex

!DOUBLE PRECISION :: ME_Vf,ME_R

INTEGER :: ME_AnisotropicorElectric
!INTEGER, INTENT(OUT) :: ME_TestIndex1,ME_TestIndex2

INTEGER, INTENT(IN) :: ZeroTimeStep,k,M_NumberofDE
DOUBLE PRECISION, INTENT(IN) :: Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R,Sub_StrainsXZ1_R,&
Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,Sub_StrainsYZ2_R    !Changed
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,4,3) :: ME_TunnelingIndex !Changed
INTEGER, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements) :: ME_TunnelingIndexDimension
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs2) :: ME_CNTIndex
DOUBLE PRECISION, INTENT(IN) :: ME_Critical_Distance
DOUBLE PRECISION, INTENT(IN), DIMENSION(9,ME_NumberofCNTs2,2) :: ME_Coordinates2
DOUBLE PRECISION, INTENT(IN) :: ME_Width1,ME_Width2
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs) :: M_NumberofElementsPerCNT
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs,ME_MaxNumberofElementsPerCNT) :: ME_ElementsPerCNT
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofGaussPoints,M_NumberofElements) :: M_Vol
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofGaussPoints,8) :: M_NN
DOUBLE PRECISION, INTENT(IN), DIMENSION(3*M_NumberofNodes+M_NumberofDE) :: M_Dis1
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofNodes,4) :: M_CurrentNodes
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_Tunneling_rhorho

ME_sub_rhorho=0.0D0
ME_sub_sigmasigma_updated=0.0D0
ME_sub_sigmasigma=0.0D0

ME_JJ=0.0D0
ME_InverseJJ=0.0D0

ME_BB=0.0D0
ME_BBT=0.0D0
ME_KComponent=0.0D0

ME_Xi=(/1.0D0,+1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,-1.0D0,-1.0D0/)
ME_Eta=(/-1.0D0,1.0D0,1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,-1.0D0/)
ME_Mu=(/-1.0D0,-1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,1.0D0,1.0D0/)

DO i=1,ME_NumberofGaussPoints
ME_Xi1(i)=ME_Xi(i)/SQRT(3.0D0)
ME_Eta1(i)=ME_Eta(i)/SQRT(3.0D0)
ME_Mu1(i)=ME_Mu(i)/SQRT(3.0D0)
ME_W(i)=1.0D0
ENDDO

DO k1=1,M_NumberofElements
DO j=1,8
ME_x(k1,j)=M_Nodes(ME_Connectivities(k1,j+1),2)
ME_y(k1,j)=M_Nodes(ME_Connectivities(k1,j+1),3)
ME_z(k1,j)=M_Nodes(ME_Connectivities(k1,j+1),4)
ENDDO
ENDDO

DO i=1,ME_NumberofGaussPoints

ME_DNDXi(i,4)=-1.0D0/8.0D0*(1.0D0-ME_Eta1(i))*(1.0D0-ME_Mu1(i))
ME_DNDEta(i,4)=-1.0D0/8.0D0*(1.0D0-ME_Xi1(i))*(1.0D0-ME_Mu1(i))
ME_DNDMu(i,4)=-1.0D0/8.0D0*(1.0D0-ME_Xi1(i))*(1.0D0-ME_Eta1(i))

ME_DNDXi(i,1)=1.0D0/8.0D0*(1.0D0-ME_Eta1(i))*(1.0D0-ME_Mu1(i))
ME_DNDEta(i,1)=-1.0D0/8.0D0*(1.0D0+ME_Xi1(i))*(1.0D0-ME_Mu1(i))
ME_DNDMu(i,1)=-1.0D0/8.0D0*(1.0D0+ME_Xi1(i))*(1.0D0-ME_Eta1(i))

ME_DNDXi(i,2)=1.0D0/8.0D0*(1.0D0+ME_Eta1(i))*(1.0D0-ME_Mu1(i))
ME_DNDEta(i,2)=1.0D0/8.0D0*(1.0D0+ME_Xi1(i))*(1.0D0-ME_Mu1(i))
ME_DNDMu(i,2)=-1.0D0/8.0D0*(1.0D0+ME_Xi1(i))*(1.0D0+ME_Eta1(i))

ME_DNDXi(i,3)=-1.0D0/8.0D0*(1.0D0+ME_Eta1(i))*(1.0D0-ME_Mu1(i))
ME_DNDEta(i,3)=1.0D0/8.0D0*(1.0D0-ME_Xi1(i))*(1.0D0-ME_Mu1(i))
ME_DNDMu(i,3)=-1.0D0/8.0D0*(1.0D0-ME_Xi1(i))*(1.0D0+ME_Eta1(i))

ME_DNDXi(i,8)=-1.0D0/8.0D0*(1.0D0-ME_Eta1(i))*(1.0D0+ME_Mu1(i))
ME_DNDEta(i,8)=-1.0D0/8.0D0*(1.0D0-ME_Xi1(i))*(1.0D0+ME_Mu1(i))
ME_DNDMu(i,8)=1.0D0/8.0D0*(1.0D0-ME_Xi1(i))*(1.0D0-ME_Eta1(i))

ME_DNDXi(i,5)=1.0D0/8.0D0*(1.0D0-ME_Eta1(i))*(1.0D0+ME_Mu1(i))
ME_DNDEta(i,5)=-1.0D0/8.0D0*(1.0D0+ME_Xi1(i))*(1.0D0+ME_Mu1(i))
ME_DNDMu(i,5)=1.0D0/8.0D0*(1.0D0+ME_Xi1(i))*(1.0D0-ME_Eta1(i))

ME_DNDXi(i,6)=1.0D0/8.0D0*(1.0D0+ME_Eta1(i))*(1.0D0+ME_Mu1(i))
ME_DNDEta(i,6)=1.0D0/8.0D0*(1.0D0+ME_Xi1(i))*(1.0D0+ME_Mu1(i))
ME_DNDMu(i,6)=1.0D0/8.0D0*(1.0D0+ME_Xi1(i))*(1.0D0+ME_Eta1(i))

ME_DNDXi(i,7)=-1.0D0/8.0D0*(1.0D0+ME_Eta1(i))*(1.0D0+ME_Mu1(i))
ME_DNDEta(i,7)=1.0D0/8.0D0*(1.0D0-ME_Xi1(i))*(1.0D0+ME_Mu1(i))
ME_DNDMu(i,7)=1.0D0/8.0D0*(1.0D0-ME_Xi1(i))*(1.0D0+ME_Eta1(i))

ME_W(i)=1.0D0

ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!OPEN(101,file = '3D_Electrical_Input/More_Inputs.txt',status ='unknown')
!READ (101,*) Vf
!READ (101,*) R
!READ (101,*) BeginofModeIndex
!READ (101,*) EndofModeIndex
!READ (101,*) MaximumRadius !For the subroutine anisotropic
!READ (101,*) AnisotropicorElectric
!READ (101,*) TestIndex1
!READ (101,*) TestIndex2
!ClOSE(101)
!
!XYofAxis(1,1)=0.0D0
!XYofAxis(1,2)=-SQRT(4.0D0*3.141592653589793*R*R/2.0D0/SQRT(3.0D0)/Vf)
!XYofAxis(2,1)=SQRT(4.0D0*3.141592653589793*R*R/2.0D0/SQRT(3.0D0)/Vf)*SQRT(3.0D0)/2.0D0
!XYofAxis(2,2)=-0.5D0*SQRT(4.0D0*3.141592653589793*R*R/2.0D0/SQRT(3.0D0)/Vf)
!XYofAxis(3,1)=XYofAxis(2,1)
!XYofAxis(3,2)=-XYofAxis(2,2)
!XYofAxis(4,1)=XYofAxis(1,1)
!XYofAxis(4,2)=-XYofAxis(1,2)
!XYofAxis(5,1)=-XYofAxis(3,1)
!XYofAxis(5,2)=XYofAxis(3,2)
!XYofAxis(6,1)=-XYofAxis(2,1)
!XYofAxis(6,2)=XYofAxis(2,2)
!XYofAxis(7,1)=0.0D0
!XYofAxis(7,2)=0.0D0
!XYofAxis(8,1)=XYofAxis(2,1)
!XYofAxis(8,2)=-1.5D0*SQRT(4.0D0*3.141592653589793*R*R/2.0D0/SQRT(3.0D0)/Vf)
!XYofAxis(9,1)=XYofAxis(8,1)
!XYofAxis(9,2)=-XYofAxis(8,2)
!XYofAxis(10,1)=-XYofAxis(9,1)
!XYofAxis(10,2)=XYofAxis(9,2)
!XYofAxis(11,1)=-XYofAxis(8,1)
!XYofAxis(11,2)=XYofAxis(8,2)
!
!DO i=1,NumberofGaussPoints
!NN(1,i)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0-Eta1(i))*(1.0D0-Mu1(i))
!NN(2,i)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0+Eta1(i))*(1.0D0-Mu1(i))
!NN(3,i)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0+Eta1(i))*(1.0D0-Mu1(i))
!NN(4,i)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0-Eta1(i))*(1.0D0-Mu1(i))
!NN(5,i)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0-Eta1(i))*(1.0D0+Mu1(i))
!NN(6,i)=1.0D0/8.0D0*(1.0D0+Xi1(i))*(1.0D0+Eta1(i))*(1.0D0+Mu1(i))
!NN(7,i)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0+Eta1(i))*(1.0D0+Mu1(i))
!NN(8,i)=1.0D0/8.0D0*(1.0D0-Xi1(i))*(1.0D0-Eta1(i))*(1.0D0+Mu1(i))
!ENDDO
!
DO i1=1,ME_NumberofGaussPoints
DO k1=1,ME_NumberofMaterialModes
DO i=1,3
DO j=1,3
ME_sub_sigmasigma(i1,k1,i,j)=ME_sigmasigma(k1,i,j)
ENDDO
ENDDO
ENDDO
ENDDO

DO i1=1,ME_NumberofGaussPoints
DO k1=1,ME_NumberofMaterialModes
DO i=1,6
DO j=1,6
ME_sub_gg(i1,k1,i,j)=ME_gg(k1,i,j)
ENDDO
ENDDO
ENDDO
ENDDO

DO i1=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements
CALL INV(k1,3,ME_sub_sigmasigma(i1,ME_Connectivities(k1,10),:,:),ME_sub_rhorho(i1,k1,:,:))
ENDDO
ENDDO

IF ( k /= 1 .OR. ZeroTimeStep /= 1) THEN
DO i=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements
IF ( k1<ME_BeginningPolymer .OR. k1>ME_EndingPolymer) THEN
ME_sub_rhorho(i,k1,1,1)=ME_sub_rhorho(i,k1,1,1)+ME_sub_gg(i,ME_Connectivities(k1,10),1,1)*M_Sub_StrainsXX(k1,i)+&
ME_sub_gg(i,ME_Connectivities(k1,10),1,2)*M_Sub_StrainsYY(k1,i)+ME_sub_gg(i,ME_Connectivities(k1,10),1,3)*M_Sub_StrainsZZ(k1,i)+&
+ME_sub_gg(i,ME_Connectivities(k1,10),1,4)*(M_Sub_StrainsYZ1(k1,i)+M_Sub_StrainsYZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),1,5)*&
(M_Sub_StrainsXZ1(k1,i)+M_Sub_StrainsXZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),1,6)*(M_Sub_StrainsXY1(k1,i)+M_Sub_StrainsXY2(k1,i))

ME_sub_rhorho(i,k1,2,2)=ME_sub_rhorho(i,k1,2,2)+ME_sub_gg(i,ME_Connectivities(k1,10),2,1)*M_Sub_StrainsXX(k1,i)+&
ME_sub_gg(i,ME_Connectivities(k1,10),2,2)*M_Sub_StrainsYY(k1,i)+ME_sub_gg(i,ME_Connectivities(k1,10),2,3)*M_Sub_StrainsZZ(k1,i)+&
+ME_sub_gg(i,ME_Connectivities(k1,10),2,4)*(M_Sub_StrainsYZ1(k1,i)+M_Sub_StrainsYZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),2,5)*&
(M_Sub_StrainsXZ1(k1,i)+M_Sub_StrainsXZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),2,6)*(M_Sub_StrainsXY1(k1,i)+M_Sub_StrainsXY2(k1,i))

ME_sub_rhorho(i,k1,3,3)=ME_sub_rhorho(i,k1,3,3)+ME_sub_gg(i,ME_Connectivities(k1,10),3,1)*M_Sub_StrainsXX(k1,i)+&
ME_sub_gg(i,ME_Connectivities(k1,10),3,2)*M_Sub_StrainsYY(k1,i)+ME_sub_gg(i,ME_Connectivities(k1,10),3,3)*M_Sub_StrainsZZ(k1,i)+&
+ME_sub_gg(i,ME_Connectivities(k1,10),3,4)*(M_Sub_StrainsYZ1(k1,i)+M_Sub_StrainsYZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),3,5)*&
(M_Sub_StrainsXZ1(k1,i)+M_Sub_StrainsXZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),3,6)*(M_Sub_StrainsXY1(k1,i)+M_Sub_StrainsXY2(k1,i))

ME_sub_rhorho(i,k1,2,3)=ME_sub_rhorho(i,k1,2,3)+ME_sub_gg(i,ME_Connectivities(k1,10),4,1)*M_Sub_StrainsXX(k1,i)+&
ME_sub_gg(i,ME_Connectivities(k1,10),4,2)*M_Sub_StrainsYY(k1,i)+ME_sub_gg(i,ME_Connectivities(k1,10),4,3)*M_Sub_StrainsZZ(k1,i)+&
+ME_sub_gg(i,ME_Connectivities(k1,10),4,4)*(M_Sub_StrainsYZ1(k1,i)+M_Sub_StrainsYZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),4,5)*&
(M_Sub_StrainsXZ1(k1,i)+M_Sub_StrainsXZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),4,6)*(M_Sub_StrainsXY1(k1,i)+M_Sub_StrainsXY2(k1,i))

ME_sub_rhorho(i,k1,1,3)=ME_sub_rhorho(i,k1,1,3)+ME_sub_gg(i,ME_Connectivities(k1,10),5,1)*M_Sub_StrainsXX(k1,i)+&
ME_sub_gg(i,ME_Connectivities(k1,10),5,2)*M_Sub_StrainsYY(k1,i)+ME_sub_gg(i,ME_Connectivities(k1,10),5,3)*M_Sub_StrainsZZ(k1,i)+&
+ME_sub_gg(i,ME_Connectivities(k1,10),5,4)*(M_Sub_StrainsYZ1(k1,i)+M_Sub_StrainsYZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),5,5)*&
(M_Sub_StrainsXZ1(k1,i)+M_Sub_StrainsXZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),5,6)*(M_Sub_StrainsXY1(k1,i)+M_Sub_StrainsXY2(k1,i))

ME_sub_rhorho(i,k1,1,2)=ME_sub_rhorho(i,k1,1,2)+ME_sub_gg(i,ME_Connectivities(k1,10),6,1)*M_Sub_StrainsXX(k1,i)+&
ME_sub_gg(i,ME_Connectivities(k1,10),6,2)*M_Sub_StrainsYY(k1,i)+ME_sub_gg(i,ME_Connectivities(k1,10),6,3)*M_Sub_StrainsZZ(k1,i)+&
+ME_sub_gg(i,ME_Connectivities(k1,10),6,4)*(M_Sub_StrainsYZ1(k1,i)+M_Sub_StrainsYZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),6,5)*&
(M_Sub_StrainsXZ1(k1,i)+M_Sub_StrainsXZ2(k1,i))+ME_sub_gg(i,ME_Connectivities(k1,10),6,6)*(M_Sub_StrainsXY1(k1,i)+M_Sub_StrainsXY2(k1,i))

ME_sub_rhorho(i,k1,3,2)=ME_sub_rhorho(i,k1,2,3)
ME_sub_rhorho(i,k1,3,1)=ME_sub_rhorho(i,k1,1,3)
ME_sub_rhorho(i,k1,2,1)=ME_sub_rhorho(i,k1,1,2)

ENDIF
ENDDO
ENDDO
ENDIF

!
!IF (AnisotropicorElectric == 1) THEN
!CALL Anisotropic(Nodes,Connectivities,BeginofModeIndex,EndofModeIndex,XYofAxis,MaximumRadius,sigmasigma,sub_sigmasigma,NN)
!ELSEIF (AnisotropicorElectric == 2) THEN
!CALL ElectricTunneling(R,Nodes,Connectivities,BeginofModeIndex,EndofModeIndex,XYofAxis,sigmasigma,sub_sigmasigma,NN)
!ELSE
!ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


DO i=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements

DO j=1,8
ME_JJ(i,k1,1,1)=ME_JJ(i,k1,1,1)+ME_x(k1,j)*ME_DNDXi(i,j)
ME_JJ(i,k1,1,2)=ME_JJ(i,k1,1,2)+ME_y(k1,j)*ME_DNDXi(i,j)
ME_JJ(i,k1,1,3)=ME_JJ(i,k1,1,3)+ME_z(k1,j)*ME_DNDXi(i,j)

ME_JJ(i,k1,2,1)=ME_JJ(i,k1,2,1)+ME_x(k1,j)*ME_DNDEta(i,j)
ME_JJ(i,k1,2,2)=ME_JJ(i,k1,2,2)+ME_y(k1,j)*ME_DNDEta(i,j)
ME_JJ(i,k1,2,3)=ME_JJ(i,k1,2,3)+ME_z(k1,j)*ME_DNDEta(i,j)

ME_JJ(i,k1,3,1)=ME_JJ(i,k1,3,1)+ME_x(k1,j)*ME_DNDMu(i,j)
ME_JJ(i,k1,3,2)=ME_JJ(i,k1,3,2)+ME_y(k1,j)*ME_DNDMu(i,j)
ME_JJ(i,k1,3,3)=ME_JJ(i,k1,3,3)+ME_z(k1,j)*ME_DNDMu(i,j)

ENDDO

ENDDO
ENDDO

DO i=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements

ME_a11=ME_JJ(i,k1,1,1)
ME_a12=ME_JJ(i,k1,1,2)
ME_a13=ME_JJ(i,k1,1,3)
ME_a21=ME_JJ(i,k1,2,1)
ME_a22=ME_JJ(i,k1,2,2)
ME_a23=ME_JJ(i,k1,2,3)
ME_a31=ME_JJ(i,k1,3,1)
ME_a32=ME_JJ(i,k1,3,2)
ME_a33=ME_JJ(i,k1,3,3)

ME_Vol(i,k1)=-ME_a13*ME_a22*ME_a31 + ME_a12*ME_a23*ME_a31 + ME_a13*ME_a21*ME_a32&
- ME_a11*ME_a23*ME_a32 - ME_a12*ME_a21*ME_a33 +ME_a11*ME_a22*ME_a33

ME_InverseJJ(i,k1,1,1)=(-ME_a23*ME_a32 + ME_a22*ME_a33)/ME_Vol(i,k1)
ME_InverseJJ(i,k1,1,2)=(ME_a13*ME_a32 - ME_a12*ME_a33)/ME_Vol(i,k1)
ME_InverseJJ(i,k1,1,3)=(-ME_a13*ME_a22 + ME_a12*ME_a23)/ME_Vol(i,k1)
ME_InverseJJ(i,k1,2,1)=(ME_a23*ME_a31 - ME_a21*ME_a33)/ME_Vol(i,k1)
ME_InverseJJ(i,k1,2,2)=(-ME_a13*ME_a31 + ME_a11*ME_a33)/ME_Vol(i,k1)
ME_InverseJJ(i,k1,2,3)=(ME_a13*ME_a21 - ME_a11*ME_a23)/ME_Vol(i,k1)
ME_InverseJJ(i,k1,3,1)=(-ME_a22*ME_a31 + ME_a21*ME_a32)/ME_Vol(i,k1)
ME_InverseJJ(i,k1,3,2)=(ME_a12*ME_a31 - ME_a11*ME_a32)/ME_Vol(i,k1)
ME_InverseJJ(i,k1,3,3)=(-ME_a12*ME_a21 + ME_a11*ME_a22)/ME_Vol(i,k1)

ENDDO
ENDDO

DO i=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements
DO j=1,8

ME_DNDX(i,k1,j)=ME_InverseJJ(i,k1,1,1)*ME_DNDXi(i,j)+ME_InverseJJ(i,k1,1,2)*ME_DNDEta(i,j)+ME_InverseJJ(i,k1,1,3)*ME_DNDMu(i,j)
ME_DNDY(i,k1,j)=ME_InverseJJ(i,k1,2,1)*ME_DNDXi(i,j)+ME_InverseJJ(i,k1,2,2)*ME_DNDEta(i,j)+ME_InverseJJ(i,k1,2,3)*ME_DNDMu(i,j)
ME_DNDZ(i,k1,j)=ME_InverseJJ(i,k1,3,1)*ME_DNDXi(i,j)+ME_InverseJJ(i,k1,3,2)*ME_DNDEta(i,j)+ME_InverseJJ(i,k1,3,3)*ME_DNDMu(i,j)

ENDDO
ENDDO
ENDDO

DO i=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements
DO j=1,8

ME_QX(i,k1,j)=-ME_DNDX(i,k1,j)
ME_QY(i,k1,j)=-ME_DNDY(i,k1,j)
ME_QZ(i,k1,j)=-ME_DNDZ(i,k1,j)

ENDDO
ENDDO
ENDDO

DO i=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements
DO j=1,8

ME_BB(i,k1,1,j)=ME_QX(i,k1,j)
ME_BB(i,k1,2,j)=ME_QY(i,k1,j)
ME_BB(i,k1,3,j)=ME_QZ(i,k1,j)

ENDDO
ENDDO
ENDDO

DO i=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements
DO i1=1,3
DO i2=1,8
ME_BBT(i,k1,i2,i1)=ME_BB(i,k1,i1,i2)
ENDDO
ENDDO
ENDDO
ENDDO

CALL ME_Tunneling_Effect(ZeroTimeStep,k,Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R,Sub_StrainsXZ1_R,&
Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,Sub_StrainsYZ2_R,M_NumberofDE,M_Nodes,ME_Connectivities,ME_sub_rhorho,M_NumberofElementsPerCNT,ME_ElementsPerCNT,&
ME_TunnelingIndex,ME_TunnelingIndexDimension,ME_Width1,ME_Width2,ME_CNTIndex,ME_Critical_Distance,ME_Coordinates2,M_Dis1,M_CurrentNodes,M_NN,M_Vol,&
ME_sub_Tunneling_rhorho)

IF (ME_Tunneling==2 .OR. ME_Tunneling==3) THEN
DO i=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements
IF ( k1>=ME_BeginningPolymer .AND. k1<=ME_EndingPolymer) THEN
ME_sub_rhorho(i,k1,1,1)=ME_sub_Tunneling_rhorho(i,k1,1,1)
ME_sub_rhorho(i,k1,2,2)=ME_sub_Tunneling_rhorho(i,k1,2,2)
ME_sub_rhorho(i,k1,3,3)=ME_sub_Tunneling_rhorho(i,k1,3,3)
ENDIF
ENDDO
ENDDO
ENDIF

DO i1=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements
CALL INV(k1,3,ME_sub_rhorho(i1,k1,:,:),ME_sub_sigmasigma_updated(i1,k1,:,:))
ENDDO
ENDDO

DO i=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements
DO i1=1,3
!IF ( ME_sub_sigmasigma_updated(i,k1,i1,i1) <= 1.0D-13 ) THEN
!ME_sub_sigmasigma_updated(i,k1,i1,i1)=1.0D-13
!ENDIF
IF ( ME_sub_sigmasigma_updated(i,k1,i1,i1) >= 1.0D4 ) THEN
ME_sub_sigmasigma_updated(i,k1,i1,i1)=1.0D-7
ENDIF
ENDDO
ENDDO
ENDDO

DO i=1,ME_NumberofGaussPoints
DO k1=1,M_NumberofElements
DO i1=1,8
DO i2=1,3
DO i3=1,3
DO i4=1,8

ME_KComponent(k1,i1,i4)=ME_KComponent(k1,i1,i4)+&
ME_BBT(i,k1,i1,i2)*ME_sub_sigmasigma_updated(i,k1,i2,i3)*ME_BB(i,k1,i3,i4)*ME_Vol(i,k1)*ME_W(i)

ENDDO ! For i1
ENDDO ! For i2
ENDDO ! For i3 
ENDDO ! For i4
ENDDO
ENDDO

!DO i1=1,8
!DO i4=1,8
!IF(ME_KComponent(1,i1,i4) /=0.0D0) THEN
!WRITE(*,*) i1,i4,ME_KComponent(1,i1,i4)
!ENDIF
!ENDDO
!ENDDO

!DO i=1,ME_NumberofGaussPoints
!DO k1=1,M_NumberofElements
!DO i1=1,3
!IF ( ME_sub_sigmasigma_updated(i,k1,i1,i1) == 1.0D-5) THEN
!ME_sub_sigmasigma_updated(i,k1,i1,i1)=1.0D5
!ENDIF
!ENDDO
!ENDDO
!ENDDO

END SUBROUTINE ME_Components_Hexahedron

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_ReadElectricalProperties (ii,Electrical_material_properties_Microscale_E,Piezo_Coefficients_Microscale_E,ME_sigmasigma,ME_gg)
IMPLICIT NONE
INTEGER :: ME_IN4,ME_IN7
INTEGER :: i,j
INTEGER, INTENT(IN) :: ii
CHARACTER*90, INTENT(IN) :: Electrical_material_properties_Microscale_E
CHARACTER*90, INTENT(IN) :: Piezo_Coefficients_Microscale_E

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

DOUBLE PRECISION, DIMENSION(ME_NumberofMaterialModes) :: ME_sigma11,ME_sigma12,ME_sigma13,ME_sigma22,ME_sigma23,ME_sigma33
DOUBLE PRECISION, INTENT(OUT), DIMENSION(ME_NumberofMaterialModes,3,3) :: ME_sigmasigma

DOUBLE PRECISION, DIMENSION(ME_NumberofMaterialModes) :: ME_rho11,ME_rho12,ME_rho13,ME_rho22,ME_rho23,ME_rho33
DOUBLE PRECISION, DIMENSION(ME_NumberofMaterialModes,3,3) :: ME_rho,ME_sigma
DOUBLE PRECISION, INTENT(OUT), DIMENSION(ME_NumberofMaterialModes,6,6) :: ME_gg

INTEGER :: ME_ModeIndex

ME_rho=0.0D0

ME_sigma=0.0D0

ME_sigmasigma=0.0D0

ME_gg=0.0D0

ME_IN4=14+100*(ii-1)
ME_IN7=17+100*(ii-1)

OPEN(ME_IN4,file = Electrical_material_properties_Microscale_E,status ='unknown')

DO i=1,ME_NumberofMaterialModes

READ(ME_IN4,*) ME_ModeIndex,ME_rho11(i),ME_rho12(i),ME_rho13(i),ME_rho22(i),ME_rho23(i),ME_rho33(i)

ME_rho(i,1,1)=ME_rho11(i)
ME_rho(i,1,2)=ME_rho12(i)
ME_rho(i,1,3)=ME_rho13(i)
ME_rho(i,2,2)=ME_rho22(i)
ME_rho(i,2,3)=ME_rho23(i)
ME_rho(i,3,3)=ME_rho33(i)
ME_rho(i,2,1)=ME_rho(i,1,2)
ME_rho(i,3,1)=ME_rho(i,1,3)
ME_rho(i,3,2)=ME_rho(i,2,3)
CALL INV(i,3,ME_rho(i,:,:),ME_sigma(i,:,:))
ME_sigma11(i)=ME_sigma(i,1,1)
ME_sigma12(i)=ME_sigma(i,1,2)
ME_sigma13(i)=ME_sigma(i,1,3)
ME_sigma22(i)=ME_sigma(i,2,2)
ME_sigma23(i)=ME_sigma(i,2,3)
ME_sigma33(i)=ME_sigma(i,3,3)

ENDDO

DO i=1,ME_NumberofMaterialModes

ME_sigmasigma(i,1,1)=ME_sigma11(i)
ME_sigmasigma(i,1,2)=ME_sigma12(i)
ME_sigmasigma(i,2,1)=ME_sigma12(i)
ME_sigmasigma(i,1,3)=ME_sigma13(i)
ME_sigmasigma(i,3,1)=ME_sigma13(i)
ME_sigmasigma(i,2,2)=ME_sigma22(i)
ME_sigmasigma(i,2,3)=ME_sigma23(i)
ME_sigmasigma(i,3,2)=ME_sigma23(i)
ME_sigmasigma(i,3,3)=ME_sigma33(i)

ENDDO

CLOSE(ME_IN4)

OPEN(ME_IN7,file = Piezo_Coefficients_Microscale_E,status ='unknown')
DO i=1,ME_NumberofMaterialModes
READ(ME_IN7,*) ME_ModeIndex,ME_gg(i,1,1),ME_gg(i,1,2),ME_gg(i,1,3),ME_gg(i,1,4),ME_gg(i,1,5),ME_gg(i,1,6),&
ME_gg(i,2,1),ME_gg(i,2,2),ME_gg(i,2,3),ME_gg(i,2,4),ME_gg(i,2,5),ME_gg(i,2,6),&
ME_gg(i,3,1),ME_gg(i,3,2),ME_gg(i,3,3),ME_gg(i,3,4),ME_gg(i,3,5),ME_gg(i,3,6),&
ME_gg(i,4,1),ME_gg(i,4,2),ME_gg(i,4,3),ME_gg(i,4,4),ME_gg(i,4,5),ME_gg(i,4,6),&
ME_gg(i,5,1),ME_gg(i,5,2),ME_gg(i,5,3),ME_gg(i,5,4),ME_gg(i,5,5),ME_gg(i,5,6),&
ME_gg(i,6,1),ME_gg(i,6,2),ME_gg(i,6,3),ME_gg(i,6,4),ME_gg(i,6,5),ME_gg(i,6,6)
ENDDO
CLOSE(ME_IN7)


RETURN

END SUBROUTINE ME_ReadElectricalProperties

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_ReadBoundaryConditions(ii,Boundaries_Microscale_E,ME_BoundaryIndex,ME_Boundaries)
IMPLICIT NONE
INTEGER :: ME_IN3=13
INTEGER :: i,j
INTEGER, INTENT(IN) :: ii
CHARACTER*90, INTENT(IN) :: Boundaries_Microscale_E

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

DOUBLE PRECISION, INTENT(OUT), DIMENSION(ME_NumberofBoundaryNodes,TotalTimeSteps+2) :: ME_Boundaries
INTEGER, INTENT(OUT), DIMENSION(ME_NumberofBoundaryNodes,2) :: ME_BoundaryIndex

ME_IN3=13+100*(ii-1)
OPEN(ME_IN3,file = Boundaries_Microscale_E,status ='unknown')
!READ (ME_IN3,*) ((ME_Boundaries(i,j),j=1,TotalTimeSteps+1),i=1,ME_NumberofBoundaryNodes)

DO i=1,ME_NumberofBoundaryNodes
READ (ME_IN3,*) (ME_Boundaries(i,j),j=1,TotalTimeSteps+1)
ENDDO

DO i=1,ME_NumberofBoundaryNodes
ME_BoundaryIndex(i,1)=ME_Boundaries(i,1)
ME_BoundaryIndex(i,2)=ME_Boundaries(i,2)
ENDDO
CLOSE(ME_IN3)
RETURN
END SUBROUTINE ME_ReadBoundaryConditions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_ReadConnectivities (ii,Connectivities_Microscale_E,ME_Connectivities)
IMPLICIT NONE
INTEGER :: ME_IN2
INTEGER :: ii
INTEGER :: i,j
CHARACTER*90, INTENT(IN) :: Connectivities_Microscale_E

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

INTEGER, INTENT(OUT), DIMENSION (M_NumberofElements,10) :: ME_Connectivities

ME_IN2=12+100*(ii-1)

OPEN(ME_IN2,file = Connectivities_Microscale_E, status ='unknown')

READ (ME_IN2,*) ((ME_Connectivities(i,j),j=1,10),i=1,M_NumberofElements)

CLOSE(ME_IN2)


RETURN
ENDSUBROUTINE ME_ReadConnectivities

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_ReadControlFlags (ii,Control_Flags_Microscale_E)
IMPLICIT NONE
INTEGER :: i
INTEGER, INTENT(IN) :: ii
INTEGER :: ME_IN1
CHARACTER*90, INTENT(IN) :: Control_Flags_Microscale_E
DOUBLE PRECISION, DIMENSION (5) :: ME_ControlFlags

!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec

ME_IN1=11+100*(ii-1)

OPEN(ME_IN1,file = Control_Flags_Microscale_E,status ='unknown')

DO i=1,5
READ (ME_IN1,*) ME_ControlFlags(i)
ENDDO

ME_NumberofBoundaryNodes=ME_ControlFlags(1)
ME_NumberofMaterialModes=ME_ControlFlags(2)
ME_NumberofGaussPoints=ME_ControlFlags(3)
ME_Mag=ME_ControlFlags(4)
ME_MaxDimenofInverseConnec=ME_ControlFlags(5)
ClOSE(ME_IN1)

RETURN
END SUBROUTINE ME_ReadControlFlags

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ME_MicroFolderNames(ii,Boundaries_Microscale_E,Connectivities_Microscale_E,Control_Flags_Microscale_E,&
Electrical_material_properties_Microscale_E,More_Inputs_Microscale_E,More_Inputs_Microscale_2_E,Effec_Electro_Prop_Microscale_E,&
Piezo_Coefficients_Microscale_E)
IMPLICIT NONE

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(IN) :: ii
CHARACTER*90, INTENT(OUT) :: Boundaries_Microscale_E
CHARACTER*90, INTENT(OUT) :: Connectivities_Microscale_E
CHARACTER*90, INTENT(OUT) :: Control_Flags_Microscale_E
CHARACTER*90, INTENT(OUT) :: Electrical_material_properties_Microscale_E
CHARACTER*90, INTENT(OUT) :: Piezo_Coefficients_Microscale_E
CHARACTER*90, INTENT(OUT) :: More_Inputs_Microscale_E
CHARACTER*90, INTENT(OUT) :: More_Inputs_Microscale_2_E
CHARACTER*90, INTENT(OUT) :: Effec_Electro_Prop_Microscale_E

IF (ii==1) THEN
Boundaries_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_1/Boundaries_Microscale_E.txt"
Connectivities_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_1/Connectivities_Microscale_E.txt"
Control_Flags_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_1/Control_Flags_Microscale_E.txt"
Electrical_material_properties_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_1/Electrical_Microscale_E.txt"
More_Inputs_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_1/More_Inputs_Microscale_E.txt"
More_Inputs_Microscale_2_E="3D_Output_Micro_F/3D_Input_Microscale_1/More_Inputs_Microscale_2_E.txt"
Effec_Electro_Prop_Microscale_E="3D_Output_Micro_F/3D_Output_Microscale_1/Effec_Electro_Prop_Microscale_E.txt"
Piezo_Coefficients_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_1/Piezo_Coef_Microscale_E.txt"

ELSEIF (ii==2) THEN
Boundaries_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_2/Boundaries_Microscale_E.txt"
Connectivities_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_2/Connectivities_Microscale_E.txt"
Control_Flags_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_2/Control_Flags_Microscale_E.txt"
Electrical_material_properties_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_2/Electrical_Microscale_E.txt"
More_Inputs_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_2/More_Inputs_Microscale_E.txt"
More_Inputs_Microscale_2_E="3D_Output_Micro_F/3D_Input_Microscale_2/More_Inputs_Microscale_2_E.txt"
Effec_Electro_Prop_Microscale_E="3D_Output_Micro_F/3D_Output_Microscale_2/Effec_Electro_Prop_Microscale_E.txt"
Piezo_Coefficients_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_2/Piezo_Coef_Microscale_E.txt"

ELSEIF (ii==3) THEN
Boundaries_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_3/Boundaries_Microscale_E.txt"
Connectivities_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_3/Connectivities_Microscale_E.txt"
Control_Flags_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_3/Control_Flags_Microscale_E.txt"
Electrical_material_properties_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_3/Electrical_Microscale_E.txt"
More_Inputs_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_3/More_Inputs_Microscale_E.txt"
More_Inputs_Microscale_2_E="3D_Output_Micro_F/3D_Input_Microscale_3/More_Inputs_Microscale_2_E.txt"
Effec_Electro_Prop_Microscale_E="3D_Output_Micro_F/3D_Output_Microscale_3/Effec_Electro_Prop_Microscale_E.txt"
Piezo_Coefficients_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_3/Piezo_Coef_Microscale_E.txt"

ELSEIF (ii==4) THEN
Boundaries_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_4/Boundaries_Microscale_E.txt"
Connectivities_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_4/Connectivities_Microscale_E.txt"
Control_Flags_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_4/Control_Flags_Microscale_E.txt"
Electrical_material_properties_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_4/Electrical_Microscale_E.txt"
More_Inputs_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_4/More_Inputs_Microscale_E.txt"
More_Inputs_Microscale_2_E="3D_Output_Micro_F/3D_Input_Microscale_4/More_Inputs_Microscale_2_E.txt"
Effec_Electro_Prop_Microscale_E="3D_Output_Micro_F/3D_Output_Microscale_4/Effec_Electro_Prop_Microscale_E.txt"
Piezo_Coefficients_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_4/Piezo_Coef_Microscale_E.txt"

ELSEIF (ii==5) THEN
Boundaries_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_5/Boundaries_Microscale_E.txt"
Connectivities_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_5/Connectivities_Microscale_E.txt"
Control_Flags_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_5/Control_Flags_Microscale_E.txt"
Electrical_material_properties_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_5/Electrical_Microscale_E.txt"
More_Inputs_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_5/More_Inputs_Microscale_E.txt"
More_Inputs_Microscale_2_E="3D_Output_Micro_F/3D_Input_Microscale_5/More_Inputs_Microscale_2_E.txt"
Effec_Electro_Prop_Microscale_E="3D_Output_Micro_F/3D_Output_Microscale_5/Effec_Electro_Prop_Microscale_E.txt"
Piezo_Coefficients_Microscale_E="3D_Output_Micro_F/3D_Input_Microscale_5/Piezo_Coef_Microscale_E.txt"

ENDIF

RETURN
ENDSUBROUTINE ME_MicroFolderNames

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_KPreprations(M_Connectivities,M_KComponent,M_NumberofDE,M_K1GRR,M_K1GRCDD,M_K1GRCII,M_F1Globall)
IMPLICIT NONE
INTEGER :: i,j,m,n,i1,j1,k1,k2,k3

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,10) :: M_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofElements,24,24) :: M_KComponent
INTEGER, INTENT(IN) :: M_NumberofDE
DOUBLE PRECISION, INTENT(OUT), DIMENSION (3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec) :: M_K1GRR
INTEGER, INTENT(OUT), DIMENSION (3*M_NumberofNodes+M_NumberofDE) :: M_K1GRCDD
INTEGER, INTENT(OUT), DIMENSION (3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec) :: M_K1GRCII
DOUBLE PRECISION, INTENT(OUT), DIMENSION (3*M_NumberofNodes+M_NumberofDE) :: M_F1Globall

M_K1GRR=0.0D0
M_K1GRCDD=0
M_K1GRCII=0

M_F1Globall=0.0D0

DO j1=1,3
DO i1=1,3
DO j=1,8
DO i=1,8
DO k1=1,M_NumberofElements

IF ( M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j) /= 0.0D0 ) THEN

IF ( M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)==0 ) THEN

M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,1)=M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)

M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)=M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)+1

M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,1)=3*(M_Connectivities(k1,j+1)-1)+j1

GOTO 1

ELSEIF ( M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1) == 1 ) THEN

IF ( 3*(M_Connectivities(k1,j+1)-1)+j1 == M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,1) ) THEN
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,1)=M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,1)+M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ELSEIF ( 3*(M_Connectivities(k1,j+1)-1)+j1 > M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,1) ) THEN
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,2)=M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)=M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)+1
M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,2)=3*(M_Connectivities(k1,j+1)-1)+j1
GOTO 1
ELSEIF ( 3*(M_Connectivities(k1,j+1)-1)+j1 < M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,1) ) THEN
M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)=M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)+1
DO k3=2,M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k3)=M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k3-1)
M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k3)=M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,1)=M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,1)=3*(M_Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF

ELSEIF ( M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1) > 1 ) THEN

IF ( 3*(M_Connectivities(k1,j+1)-1)+j1 < M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,1) ) THEN
M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)=M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)+1
DO k3=M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1),2,-1
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k3)=M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k3-1)
M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k3)=M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,1)=M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,1)=3*(M_Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF

DO k2=1,M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)-1

IF ( 3*(M_Connectivities(k1,j+1)-1)+j1 == M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k2) ) THEN
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k2)=M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k2)+M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ENDIF

IF (3*(M_Connectivities(k1,j+1)-1)+j1 > M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k2) .AND. &
3*(M_Connectivities(k1,j+1)-1)+j1 < M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k2+1) ) THEN
M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)=M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)+1
DO k3=M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1),k2+2,-1
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k3)=M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k3-1)
M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k3)=M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k3-1)
ENDDO
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k2+1)=M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k2+1)=3*(M_Connectivities(k1,j+1)-1)+j1
GOTO 1
ENDIF  

ENDDO  !k2

DO k2=M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1),M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)

IF ( 3*(M_Connectivities(k1,j+1)-1)+j1 == M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k2) ) THEN
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k2)=M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k2)+M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
GOTO 1
ENDIF

IF ( 3*(M_Connectivities(k1,j+1)-1)+j1 > M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k2)  ) THEN
M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)=M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)+1
DO k3=k2+1,M_K1GRCDD(3*(M_Connectivities(k1,i+1)-1)+i1)
M_K1GRR(3*(M_Connectivities(k1,i+1)-1)+i1,k3)=M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j)
M_K1GRCII(3*(M_Connectivities(k1,i+1)-1)+i1,k3)=3*(M_Connectivities(k1,j+1)-1)+j1
ENDDO
GOTO 1
ENDIF  

ENDDO !For k2

1 ENDIF  ! For ( M_K1GRCD(3*(M_Connectivities(k1,i+1)-1)+i1)==0 )
 
ENDIF  ! For ( M_KComponent(k1,8*(i1-1)+i,8*(j1-1)+j) /= 0.0D0 )

ENDDO
ENDDO
ENDDO
ENDDO
ENDDO

ENDSUBROUTINE M_KPreprations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_StrainsandStresses (M_Vol,M_Vol_Ele,M_CC,M_Sub_CC,M_Connectivities,M_Nodes,M_NodeIndex,M_DNDX,M_DNDY,M_DNDZ,&
M_Dis1,M_StrainsXX,M_StrainsYY,M_StrainsZZ,M_StrainsXY,M_StrainsXZ,M_StrainsYZ,M_StrainsVonMises,M_StressesXX,M_StressesYY,&
M_StressesZZ,M_StressesXY,M_StressesXZ,M_StressesYZ,M_StressesVonMises,M_Sub_StrainsXX,M_Sub_StrainsYY,M_Sub_StrainsZZ,&
M_Sub_StrainsXY1,M_Sub_StrainsXZ1,M_Sub_StrainsYZ1,M_Sub_StrainsXY2,M_Sub_StrainsXZ2,M_Sub_StrainsYZ2,M_Sub_StressesXX,&
M_Sub_StressesYY,M_Sub_StressesZZ,M_Sub_StressesXY,M_Sub_StressesXZ,M_Sub_StressesYZ,M_CurrentNodes)
IMPLICIT NONE
INTEGER :: i,j,k

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,6) :: M_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes,4) :: M_Nodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNodes) :: M_NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*M_NumberofNodes) :: M_Dis1
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofGaussPoints,M_NumberofElements) :: M_Vol
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_Vol_Ele
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofGaussPoints,M_NumberofElements,8+3) :: M_DNDX,M_DNDY,M_DNDZ
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofMaterialModes,6,6) :: M_CC
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofGaussPoints,M_NumberofMaterialModes,6,6) :: M_Sub_CC
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StrainsXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StrainsYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StrainsZZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StrainsXY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StrainsXZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StrainsYZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StrainsVonMises
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StressesXX
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StressesYY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StressesZZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StressesXY
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StressesXZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StressesYZ
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements) :: M_StressesVonMises
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofNodes,4)  :: M_CurrentNodes 

DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements,8) :: M_Sub_StrainsXX,M_Sub_StrainsYY,M_Sub_StrainsZZ,M_Sub_StrainsXY1,M_Sub_StrainsXZ1,M_Sub_StrainsYZ1,&
M_Sub_StrainsXY2,M_Sub_StrainsXZ2,M_Sub_StrainsYZ2,M_Sub_StressesXX,M_Sub_StressesYY,M_Sub_StressesZZ,M_Sub_StressesXY,M_Sub_StressesXZ,M_Sub_StressesYZ

!INTEGER :: M_OUT1,M_OUT2,M_OUT3,M_OUT4

M_Sub_StrainsXX=0.0D0
M_Sub_StrainsYY=0.0D0
M_Sub_StrainsZZ=0.0D0
M_Sub_StrainsXY1=0.0D0
M_Sub_StrainsXZ1=0.0D0
M_Sub_StrainsYZ1=0.0D0
M_Sub_StrainsXY2=0.0D0
M_Sub_StrainsXZ2=0.0D0
M_Sub_StrainsYZ2=0.0D0
M_StrainsXX=0.0D0
M_StrainsYY=0.0D0
M_StrainsZZ=0.0D0
M_StrainsXY=0.0D0
M_StrainsXZ=0.0D0
M_StrainsYZ=0.0D0
M_StrainsVonMises=0.0D0

M_Sub_StressesXX=0.0D0
M_Sub_StressesYY=0.0D0
M_Sub_StressesZZ=0.0D0
M_Sub_StressesXY=0.0D0
M_Sub_StressesXZ=0.0D0
M_Sub_StressesYZ=0.0D0
M_StressesXX=0.0D0
M_StressesYY=0.0D0
M_StressesZZ=0.0D0
M_StressesXY=0.0D0
M_StressesXZ=0.0D0
M_StressesYZ=0.0D0
M_StressesVonMises=0.0D0

M_Vol_Ele=0.0D0


DO k=1,M_NumberofElements
DO j=1,M_NumberofGaussPoints
DO i=1,8

M_Sub_StrainsXX(k,j)=M_Sub_StrainsXX(k,j)+M_Dis1(M_Connectivities(k,i+1))*M_DNDX(j,k,i)
M_Sub_StrainsYY(k,j)=M_Sub_StrainsYY(k,j)+M_Dis1(M_NumberofNodes+M_Connectivities(k,i+1))*M_DNDY(j,k,i)
M_Sub_StrainsZZ(k,j)=M_Sub_StrainsZZ(k,j)+M_Dis1(2*M_NumberofNodes+M_Connectivities(k,i+1))*M_DNDZ(j,k,i)
M_Sub_StrainsXY1(k,j)=M_Sub_StrainsXY1(k,j)+0.5D0*M_Dis1(M_Connectivities(k,i+1))*M_DNDY(j,k,i)
M_Sub_StrainsXZ1(k,j)=M_Sub_StrainsXZ1(k,j)+0.5D0*M_Dis1(M_Connectivities(k,i+1))*M_DNDZ(j,k,i)
M_Sub_StrainsYZ1(k,j)=M_Sub_StrainsYZ1(k,j)+0.5D0*M_Dis1(M_NumberofNodes+M_Connectivities(k,i+1))*M_DNDZ(j,k,i)
M_Sub_StrainsXY2(k,j)=M_Sub_StrainsXY2(k,j)+0.5D0*M_Dis1(M_NumberofNodes+M_Connectivities(k,i+1))*M_DNDX(j,k,i)
M_Sub_StrainsXZ2(k,j)=M_Sub_StrainsXZ2(k,j)+0.5D0*M_Dis1(2*M_NumberofNodes+M_Connectivities(k,i+1))*M_DNDX(j,k,i)
M_Sub_StrainsYZ2(k,j)=M_Sub_StrainsYZ2(k,j)+0.5D0*M_Dis1(2*M_NumberofNodes+M_Connectivities(k,i+1))*M_DNDY(j,k,i)

ENDDO

M_Sub_StressesXX(k,j)=M_Sub_CC(j,M_Connectivities(k,10),1,1)*M_Sub_StrainsXX(k,j)+M_Sub_CC(j,M_Connectivities(k,10),1,2)*M_Sub_StrainsYY(k,j)+&
M_Sub_CC(j,M_Connectivities(k,10),1,3)*M_Sub_StrainsZZ(k,j)
M_Sub_StressesYY(k,j)=M_Sub_CC(j,M_Connectivities(k,10),1,2)*M_Sub_StrainsXX(k,j)+M_Sub_CC(j,M_Connectivities(k,10),2,2)*M_Sub_StrainsYY(k,j)+&
M_Sub_CC(j,M_Connectivities(k,10),2,3)*M_Sub_StrainsZZ(k,j)
M_Sub_StressesZZ(k,j)=M_Sub_CC(j,M_Connectivities(k,10),1,3)*M_Sub_StrainsXX(k,j)+M_Sub_CC(j,M_Connectivities(k,10),2,3)*M_Sub_StrainsYY(k,j)+&
M_Sub_CC(j,M_Connectivities(k,10),3,3)*M_Sub_StrainsZZ(k,j)
M_Sub_StressesXY(k,j)=2.0D0*M_Sub_CC(j,M_Connectivities(k,10),6,6)*(M_Sub_StrainsXY1(k,j)+M_Sub_StrainsXY2(k,j))
M_Sub_StressesXZ(k,j)=2.0D0*M_Sub_CC(j,M_Connectivities(k,10),5,5)*(M_Sub_StrainsXZ1(k,j)+M_Sub_StrainsXZ2(k,j))
M_Sub_StressesYZ(k,j)=2.0D0*M_Sub_CC(j,M_Connectivities(k,10),4,4)*(M_Sub_StrainsYZ1(k,j)+M_Sub_StrainsYZ2(k,j))

ENDDO
ENDDO

DO k=1,M_NumberofElements
DO j=1,M_NumberofGaussPoints

M_StrainsXX(k)=M_StrainsXX(k)+M_Sub_StrainsXX(k,j)*M_Vol(j,k)
M_StrainsYY(k)=M_StrainsYY(k)+M_Sub_StrainsYY(k,j)*M_Vol(j,k)
M_StrainsZZ(k)=M_StrainsZZ(k)+M_Sub_StrainsZZ(k,j)*M_Vol(j,k)
M_StrainsXY(k)=M_StrainsXY(k)+M_Sub_StrainsXY1(k,j)*M_Vol(j,k)+M_Sub_StrainsXY2(k,j)*M_Vol(j,k)
M_StrainsXZ(k)=M_StrainsXZ(k)+M_Sub_StrainsXZ1(k,j)*M_Vol(j,k)+M_Sub_StrainsXZ2(k,j)*M_Vol(j,k)
M_StrainsYZ(k)=M_StrainsYZ(k)+M_Sub_StrainsYZ1(k,j)*M_Vol(j,k)+M_Sub_StrainsYZ2(k,j)*M_Vol(j,k)
M_StrainsVonMises(k)=M_StrainsVonMises(k)+SQRT(0.5D0*( (M_Sub_StrainsXX(k,j)-M_Sub_StrainsYY(k,j))**2+&
(M_Sub_StrainsYY(k,j)-M_Sub_StrainsZZ(k,j))**2+(M_Sub_StrainsXX(k,j)-M_Sub_StrainsZZ(k,j))**2+&
6*( (M_Sub_StrainsYZ1(k,j)+M_Sub_StrainsYZ2(k,j))**2+(M_Sub_StrainsXZ1(k,j)+M_Sub_StrainsXZ2(k,j))**2+&
(M_Sub_StrainsXY1(k,j)+M_Sub_StrainsXY2(k,j))**2) ) )*M_Vol(j,k)

M_StressesXX(k)=M_StressesXX(k)+M_Sub_StressesXX(k,j)*M_Vol(j,k)
M_StressesYY(k)=M_StressesYY(k)+M_Sub_StressesYY(k,j)*M_Vol(j,k)
M_StressesZZ(k)=M_StressesZZ(k)+M_Sub_StressesZZ(k,j)*M_Vol(j,k)
M_StressesXY(k)=M_StressesXY(k)+M_Sub_StressesXY(k,j)*M_Vol(j,k)
M_StressesXZ(k)=M_StressesXZ(k)+M_Sub_StressesXZ(k,j)*M_Vol(j,k)
M_StressesYZ(k)=M_StressesYZ(k)+M_Sub_StressesYZ(k,j)*M_Vol(j,k)
M_Vol_Ele(k)=M_Vol_Ele(k)+M_Vol(j,k)

M_StressesVonMises(k)=M_StressesVonMises(k)+SQRT(0.5D0*( (M_Sub_StressesXX(k,j)-M_Sub_StressesYY(k,j))**2+&
(M_Sub_StressesYY(k,j)-M_Sub_StressesZZ(k,j))**2+(M_Sub_StressesXX(k,j)-M_Sub_StressesZZ(k,j))**2+&
6*(M_Sub_StressesYZ(k,j)**2+M_Sub_StressesXZ(k,j)**2+M_Sub_StressesXY(k,j)**2) ) )*M_Vol(j,k)

ENDDO
M_StrainsXX(k)=M_StrainsXX(k)/M_Vol_Ele(k)
M_StrainsYY(k)=M_StrainsYY(k)/M_Vol_Ele(k)
M_StrainsZZ(k)=M_StrainsZZ(k)/M_Vol_Ele(k)
M_StrainsXY(k)=M_StrainsXY(k)/M_Vol_Ele(k)
M_StrainsXZ(k)=M_StrainsXZ(k)/M_Vol_Ele(k)
M_StrainsYZ(k)=M_StrainsYZ(k)/M_Vol_Ele(k)
M_StrainsVonMises(k)=M_StrainsVonMises(k)/M_Vol_Ele(k)

M_StressesXX(k)=M_StressesXX(k)/M_Vol_Ele(k)
M_StressesYY(k)=M_StressesYY(k)/M_Vol_Ele(k)
M_StressesZZ(k)=M_StressesZZ(k)/M_Vol_Ele(k)
M_StressesXY(k)=M_StressesXY(k)/M_Vol_Ele(k)
M_StressesXZ(k)=M_StressesXZ(k)/M_Vol_Ele(k)
M_StressesYZ(k)=M_StressesYZ(k)/M_Vol_Ele(k)
M_StressesVonMises(k)=M_StressesVonMises(k)/M_Vol_Ele(k)

ENDDO


DO i=1,M_NumberofNodes

M_CurrentNodes(i,1)=M_Nodes(i,1)
M_CurrentNodes(i,2)=M_Nodes(i,2)+M_Mag*M_Dis1(i)
M_CurrentNodes(i,3)=M_Nodes(i,3)+M_Mag*M_Dis1(M_NumberofNodes+i)
M_CurrentNodes(i,4)=M_Nodes(i,4)+M_Mag*M_Dis1(2*M_NumberofNodes+i)

ENDDO

!DO i=1,M_NumberofElements
!
!WRITE(*,*) M_StrainsXX(i),M_StrainsYY(i),M_StrainsZZ(i),M_StrainsXY(i),M_StrainsXZ(i),M_StrainsYZ(i)
!
!ENDDO

!M_OUT1=5
!M_OUT2=6
!M_OUT3=7
!M_OUT4=8
!
!OPEN (M_OUT1,file = '3D_Output_Micro_F/Sub_Strains_11_Microscale.txt',status = 'unknown')
!OPEN (M_OUT2,file = '3D_Output_Micro_F/Sub_Strains_22_Microscale.txt',status = 'unknown')
!OPEN (M_OUT3,file = '3D_Output_Micro_F/Sub_Strains_12_Microscale.txt',status = 'unknown')
!OPEN (M_OUT4,file = '3D_Output_Micro_F/Sub_Strains_33_Microscale.txt',status = 'unknown')
!DO k=1,M_NumberofElements
!WRITE(M_OUT1,10) M_Sub_StrainsXX(k,1),M_Sub_StrainsXX(k,2),M_Sub_StrainsXX(k,3),M_Sub_StrainsXX(k,4),&
!M_Sub_StrainsXX(k,5),M_Sub_StrainsXX(k,6),M_Sub_StrainsXX(k,7),M_Sub_StrainsXX(k,8)
!WRITE(M_OUT2,10) M_Sub_StrainsYY(k,1),M_Sub_StrainsYY(k,2),M_Sub_StrainsYY(k,3),M_Sub_StrainsYY(k,4),&
!M_Sub_StrainsYY(k,5),M_Sub_StrainsYY(k,6),M_Sub_StrainsYY(k,7),M_Sub_StrainsYY(k,8)
!WRITE(M_OUT3,10) M_Sub_StrainsXY1(k,1)+M_Sub_StrainsXY2(k,1),M_Sub_StrainsXY1(k,2)+M_Sub_StrainsXY2(k,2),&
!M_Sub_StrainsXY1(k,3)+M_Sub_StrainsXY2(k,3),M_Sub_StrainsXY1(k,4)+M_Sub_StrainsXY2(k,4),&
!M_Sub_StrainsXY1(k,5)+M_Sub_StrainsXY2(k,5),M_Sub_StrainsXY1(k,6)+M_Sub_StrainsXY2(k,6),&
!M_Sub_StrainsXY1(k,7)+M_Sub_StrainsXY2(k,7),M_Sub_StrainsXY1(k,8)+M_Sub_StrainsXY2(k,8)
!WRITE(M_OUT4,10) M_Sub_StrainsZZ(k,1),M_Sub_StrainsZZ(k,2),M_Sub_StrainsZZ(k,3),M_Sub_StrainsZZ(k,4),&
!M_Sub_StrainsZZ(k,5),M_Sub_StrainsZZ(k,6),M_Sub_StrainsZZ(k,7),M_Sub_StrainsZZ(k,8)
!ENDDO
!CLOSE(M_OUT1)
!CLOSE(M_OUT2)
!CLOSE(M_OUT3)
!CLOSE(M_OUT4)
!
!10 FORMAT(8(E30.15,1X),/)

END SUBROUTINE M_StrainsandStresses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INCLUDE 'mkl_pardiso.f90'
SUBROUTINE HarwellBoeing_Reduced(k,Dimen1,DUMMY,DUMMY2,A1,JA1,JA2,IA,F1Global,x)
IMPLICIT NONE
INTEGER, INTENT(IN) :: Dimen1,DUMMY,DUMMY2
INTEGER, INTENT(IN) :: k
DOUBLE PRECISION, INTENT(IN), DIMENSION (Dimen1) :: F1Global

INTEGER :: i,j

INTEGER, INTENT(IN), DIMENSION(DUMMY2) ::IA
INTEGER, INTENT(IN), DIMENSION(DUMMY) :: JA1
INTEGER, INTENT(IN), DIMENSION(DUMMY) :: JA2
DOUBLE PRECISION, INTENT(IN), DIMENSION(DUMMY) :: A1

!For pardiso_unsymm
INTEGER, DIMENSION(64) :: pt
INTEGER maxfct, mnum, mtype, phase, nrhs, error, msglvl
INTEGER, DIMENSION(64) :: iparm
DOUBLE PRECISION,  Dimension(Dimen1) :: b
DOUBLE PRECISION,  INTENT(OUT), Dimension(Dimen1) :: x
INTEGER idum
DOUBLE PRECISION waltime1, waltime2, ddum
DOUBLE PRECISION :: ElapsedTime

nrhs=1
maxfct=1
mnum=1
     
do i = 1, 64
iparm(i) = 0
end do

iparm(1) = 1 ! no solver default
iparm(2) = 2 ! fill-in reordering from METIS
iparm(3) = 1 ! numbers of processors
iparm(4) = 0 ! no iterative-direct algorithm
iparm(5) = 0 ! no user fill-in reducing permutation
iparm(6) = 0 ! =0 solution on the first n compoments of x
iparm(7) = 0 ! not in use
iparm(8) = 15 ! numbers of iterative refinement steps
iparm(9) = 0 ! not in use
iparm(10) = 16 ! perturbe the pivot elements with 1E-13
iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
iparm(12) = 0 ! not in use
iparm(13) = 1 ! maximum weighted matching algorithm is switched-on (default for non-symmetric)
iparm(14) = 0 ! Output: number of perturbed pivots
iparm(15) = 0 ! not in use
iparm(16) = 0 ! not in use
iparm(17) = 0 ! not in use
iparm(18) = -1 ! Output: number of nonzeros in the factor LU
iparm(19) = -1 ! Output: Mflops for LU factorization
iparm(20) = 0 ! Output: Numbers of CG Iterations
!iparm(28) = 1
!iparm(51) = 1
!iparm(52) = 1
error = 1 ! initialize error flag
msglvl = 0 ! print statistical information
mtype = 11 ! real unsymmetric

do i = 1, 64
   pt(i) = 0
end do

phase = 11
CALL pardiso (pt, maxfct, mnum, mtype, phase, Dimen1, A1, IA, JA1,&
idum, nrhs, iparm, msglvl, ddum, ddum, error)

phase = 22 ! only factorization      
CALL pardiso (pt, maxfct, mnum, mtype, phase, Dimen1, A1, IA, JA1,&
idum, nrhs, iparm, msglvl, ddum, ddum, error)

iparm(8) = 12 ! max numbers of iterative refinement steps

phase = 33 ! only factorization
do i = 1, Dimen1
b(i) = F1Global(i)
end do

CALL pardiso (pt, maxfct, mnum, mtype, phase, Dimen1, A1, IA, JA1,&
idum, nrhs, iparm, msglvl, b, x, error)


phase = -1 ! release internal memory
CALL pardiso (pt, maxfct, mnum, mtype, phase, Dimen1, ddum, idum, idum,&
idum, nrhs, iparm, msglvl, ddum, ddum, error)

!CALL CPU_TIME (ElapsedTime)
!WRITE(20,*) "Right after MKL  ","Elapsed Time=",ElapsedTime/60.0D0
!
!WRITE(20,*) "Error output",error

END SUBROUTINE HarwellBoeing_Reduced

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_Matrix_DE_Reduced(k,M_NumberofDE,M_One_Step_DEs,M_BoundaryIndex,M_Boundaries,M_K1GRR,M_K1GRCDD,M_K1GRCII,M_F1Globall,M_Connectivities,&
M_Dis1,M_Highest_GRCD)
IMPLICIT NONE
INTEGER :: i,j,m,n,i1,j1,k1,k2,k3
INTEGER, INTENT(IN) :: k

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

INTEGER, INTENT(IN) :: M_NumberofDE
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec) :: M_K1GRR
INTEGER, INTENT(IN), DIMENSION (3*M_NumberofNodes+M_NumberofDE) :: M_K1GRCDD
INTEGER, INTENT(IN), DIMENSION (3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec) :: M_K1GRCII
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*M_NumberofNodes+M_NumberofDE) :: M_F1Globall
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: M_K1GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: M_K1GRCD
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: M_K1GRCI
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_F1Global

INTEGER, ALLOCATABLE, DIMENSION (:,:,:) :: M_K2GR
INTEGER, ALLOCATABLE, DIMENSION (:) :: M_K2GRRD

DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofDE,8) :: M_One_Step_DEs
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofBoundaryNodes,TotalTimeSteps+2) :: M_Boundaries
INTEGER, INTENT(IN), DIMENSION(M_NumberofBoundaryNodes,2) :: M_BoundaryIndex
DOUBLE PRECISION, DIMENSION (3*M_NumberofNodes+M_NumberofDE) :: M_Dis
DOUBLE PRECISION, INTENT(OUT),DIMENSION (3*M_NumberofNodes+M_NumberofDE) :: M_Dis1
INTEGER :: M_Dimen1
DOUBLE PRECISION, DIMENSION(3*M_NumberofNodes+M_NumberofDE) :: M_x
INTEGER :: M_OUT1=40,M_OUT2=41
INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,10) :: M_Connectivities
DOUBLE PRECISION :: ElapsedTime,Temp
INTEGER :: M_TEMP_GRCD
INTEGER, INTENT(OUT) :: M_Highest_GRCD
INTEGER :: M_Temp_INT
DOUBLE PRECISION :: M_Temp
INTEGER :: M_DUMMY,M_DUMMY2
INTEGER, ALLOCATABLE, DIMENSION(:) :: M_IA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: M_IA
INTEGER, ALLOCATABLE, DIMENSION(:) :: M_JA1
INTEGER, ALLOCATABLE, DIMENSION(:) :: M_JA2
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: M_A1

ALLOCATE(M_K1GR(3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec))
ALLOCATE(M_K1GRCD(3*M_NumberofNodes+M_NumberofDE))
ALLOCATE(M_K1GRCI(3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec))
ALLOCATE(M_F1Global(3*M_NumberofNodes+M_NumberofDE))

M_K1GR=M_K1GRR
M_K1GRCD=M_K1GRCDD
M_K1GRCI=M_K1GRCII
M_F1Global=M_F1Globall

!!!Apply P.B.C.s (Top Right and Bottom Right of F)

DO i=1,M_NumberofDE

DO j=1,int(M_One_Step_DEs(i,1))

IF(int(M_One_Step_DEs(i,3*j))==1) THEN

M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+1 ) = M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+1 )+1

M_K1GR( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+1,M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+1 ) ) = M_One_Step_DEs(i,3*j+1) 

M_K1GRCI( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+1, M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+1 ) ) = 3*M_NumberofNodes+i

ELSEIF(int(M_One_Step_DEs(i,3*j))==2) THEN

M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+2 ) = M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+2 )+1

M_K1GR( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+2,M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+2 ) ) = M_One_Step_DEs(i,3*j+1) 

M_K1GRCI( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+2, M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+2 ) ) = 3*M_NumberofNodes+i


ELSEIF(int(M_One_Step_DEs(i,3*j))==3) THEN

M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+3 ) = M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+3 )+1

M_K1GR( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+3,M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+3 ) ) = M_One_Step_DEs(i,3*j+1) 

M_K1GRCI( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+3, M_K1GRCD( 3*(int(M_One_Step_DEs(i,3*j-1))-1)+3 ) ) = 3*M_NumberofNodes+i

ENDIF

M_F1Global(3*M_NumberofNodes+i)=M_One_Step_DEs(i,3*int(M_One_Step_DEs(i,1))+1+1)

ENDDO
ENDDO

!!!Apply P.B.C.s (Bottom left)
DO i=1,M_NumberofDE
DO j=1,int(M_One_Step_DEs(i,1))
M_K1GRCD( 3*M_NumberofNodes+i ) = M_K1GRCD( 3*M_NumberofNodes+i )+1
IF(int(M_One_Step_DEs(i,3*j))==1) THEN
M_K1GR( 3*M_NumberofNodes+i,M_K1GRCD( 3*M_NumberofNodes+i ) ) = M_One_Step_DEs(i,3*j+1) 
M_K1GRCI( 3*M_NumberofNodes+i, M_K1GRCD( 3*M_NumberofNodes+i ) ) = 3*(int(M_One_Step_DEs(i,3*j-1))-1)+1
ELSEIF(int(M_One_Step_DEs(i,3*j))==2) THEN
M_K1GR( 3*M_NumberofNodes+i,M_K1GRCD( 3*M_NumberofNodes+i ) ) = M_One_Step_DEs(i,3*j+1) 
M_K1GRCI( 3*M_NumberofNodes+i, M_K1GRCD( 3*M_NumberofNodes+i ) ) = 3*(int(M_One_Step_DEs(i,3*j-1))-1)+2
ELSEIF(int(M_One_Step_DEs(i,3*j))==3) THEN
M_K1GR( 3*M_NumberofNodes+i,M_K1GRCD( 3*M_NumberofNodes+i ) ) = M_One_Step_DEs(i,3*j+1) 
M_K1GRCI( 3*M_NumberofNodes+i, M_K1GRCD( 3*M_NumberofNodes+i ) ) = 3*(int(M_One_Step_DEs(i,3*j-1))-1)+3
ENDIF
ENDDO
ENDDO

ALLOCATE(M_K2GR(3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec,2))
ALLOCATE(M_K2GRRD(3*M_NumberofNodes+M_NumberofDE))
M_K2GR=0
M_K2GRRD=0

DO i=1,3*M_NumberofNodes+M_NumberofDE
DO j=1,M_K1GRCD(i)
M_K2GRRD(M_K1GRCI(i,j))=M_K2GRRD(M_K1GRCI(i,j))+1
M_K2GR(M_K1GRCI(i,j),M_K2GRRD(M_K1GRCI(i,j)),1)=i
M_K2GR(M_K1GRCI(i,j),M_K2GRRD(M_K1GRCI(i,j)),2)=j
ENDDO
ENDDO

M_Dis=0.0D0

DO i=1,M_NumberofBoundaryNodes

IF (M_BoundaryIndex(i,2)==1) THEN

M_Dis(3*(M_BoundaryIndex(i,1)-1)+1)=M_Boundaries(i,k+2)

ELSEIF (M_BoundaryIndex(i,2)==2) THEN

M_Dis(3*(M_BoundaryIndex(i,1)-1)+2)=M_Boundaries(i,k+2)

ELSEIF (M_BoundaryIndex(i,2)==3) THEN

M_Dis(3*(M_BoundaryIndex(i,1)-1)+3)=M_Boundaries(i,k+2)

ENDIF

ENDDO


DO i=1,M_NumberofBoundaryNodes

IF (M_BoundaryIndex(i,2)==1) THEN

DO m=1,M_K2GRRD(3*(M_BoundaryIndex(i,1)-1)+1)
M_F1Global( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+1,m,1) )=M_F1Global( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+1,m,1) )&
-M_K1GR( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+1,m,1), M_K2GR(3*(M_BoundaryIndex(i,1)-1)+1,m,2) )*M_Dis(3*(M_BoundaryIndex(i,1)-1)+1)

ENDDO

M_F1Global(3*(M_BoundaryIndex(i,1)-1)+1)=M_Dis(3*(M_BoundaryIndex(i,1)-1)+1)

DO m=1,M_K2GRRD(3*(M_BoundaryIndex(i,1)-1)+1)
M_K1GR( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+1,m,1), M_K2GR(3*(M_BoundaryIndex(i,1)-1)+1,m,2) )=0.0D0
ENDDO

DO m=1,M_K1GRCD(3*(M_BoundaryIndex(i,1)-1)+1)
IF ( 3*(M_BoundaryIndex(i,1)-1)+1 /= M_K1GRCI(3*(M_BoundaryIndex(i,1)-1)+1,m) ) THEN
M_K1GR( 3*(M_BoundaryIndex(i,1)-1)+1, m )=0.0D0
ELSE
M_K1GR( 3*(M_BoundaryIndex(i,1)-1)+1, m )=1.0D0
ENDIF
ENDDO


ELSEIF (M_BoundaryIndex(i,2)==2) THEN

DO m=1,M_K2GRRD(3*(M_BoundaryIndex(i,1)-1)+2)
M_F1Global( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+2,m,1) )=M_F1Global( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+2,m,1) )&
-M_K1GR( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+2,m,1), M_K2GR(3*(M_BoundaryIndex(i,1)-1)+2,m,2) )*M_Dis(3*(M_BoundaryIndex(i,1)-1)+2)

ENDDO

M_F1Global(3*(M_BoundaryIndex(i,1)-1)+2)=M_Dis(3*(M_BoundaryIndex(i,1)-1)+2)

DO m=1,M_K2GRRD(3*(M_BoundaryIndex(i,1)-1)+2)
M_K1GR( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+2,m,1), M_K2GR(3*(M_BoundaryIndex(i,1)-1)+2,m,2) )=0.0D0
ENDDO

DO m=1,M_K1GRCD(3*(M_BoundaryIndex(i,1)-1)+2)
IF ( 3*(M_BoundaryIndex(i,1)-1)+2 /= M_K1GRCI(3*(M_BoundaryIndex(i,1)-1)+2,m) ) THEN
M_K1GR( 3*(M_BoundaryIndex(i,1)-1)+2, m )=0.0D0
ELSE
M_K1GR( 3*(M_BoundaryIndex(i,1)-1)+2, m )=1.0D0
ENDIF
ENDDO

ELSEIF (M_BoundaryIndex(i,2)==3) THEN

DO m=1,M_K2GRRD(3*(M_BoundaryIndex(i,1)-1)+3)
M_F1Global( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+3,m,1) )=M_F1Global( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+3,m,1) )&
-M_K1GR( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+3,m,1), M_K2GR(3*(M_BoundaryIndex(i,1)-1)+3,m,2) )*M_Dis(3*(M_BoundaryIndex(i,1)-1)+3)

ENDDO

M_F1Global(3*(M_BoundaryIndex(i,1)-1)+3)=M_Dis(3*(M_BoundaryIndex(i,1)-1)+3)

DO m=1,M_K2GRRD(3*(M_BoundaryIndex(i,1)-1)+3)
M_K1GR( M_K2GR(3*(M_BoundaryIndex(i,1)-1)+3,m,1), M_K2GR(3*(M_BoundaryIndex(i,1)-1)+3,m,2) )=0.0D0
ENDDO

DO m=1,M_K1GRCD(3*(M_BoundaryIndex(i,1)-1)+3)
IF ( 3*(M_BoundaryIndex(i,1)-1)+3 /= M_K1GRCI(3*(M_BoundaryIndex(i,1)-1)+3,m) ) THEN
M_K1GR( 3*(M_BoundaryIndex(i,1)-1)+3, m )=0.0D0
ELSE
M_K1GR( 3*(M_BoundaryIndex(i,1)-1)+3, m )=1.0D0
ENDIF
ENDDO

ENDIF

ENDDO

DEALLOCATE(M_K2GR)
DEALLOCATE(M_K2GRRD)

M_TEMP_GRCD=M_K1GRCD(1)
M_Highest_GRCD=M_TEMP_GRCD
DO i=2,3*M_NumberofNodes+M_NumberofDE
IF ( M_K1GRCD(i)>M_TEMP_GRCD ) THEN
M_TEMP_GRCD=M_K1GRCD(i)
M_Highest_GRCD=M_TEMP_GRCD
ENDIF
ENDDO

!OPEN(33,FILE="3D_Output_Micro_F/Check3.txt",STATUS="UNKNOWN")
!DO i=1,1000  !,M_MaxDimenofInverseConnec
!WRITE(33,34) (M_K1GR(i,k1),k1=1,80)
!ENDDO
!34 FORMAT(<80>E15.8,/)
!CLOSE(33)

!CALL CPU_TIME (ElapsedTime)
!WRITE(*,*) "M_Highest_GRCD  ",M_Highest_GRCD

M_DUMMY=0

DO i=1,3*M_NumberofNodes+M_NumberofDE
DO j=1,M_K1GRCD(i)
IF (M_K1GR(i,j) /= 0.0D0) THEN
M_DUMMY=M_DUMMY+1
ELSEIF (M_K1GR(i,j) == 0.0D0) THEN
GOTO 3
ENDIF
3 ENDDO
ENDDO

ALLOCATE(M_IA1(M_DUMMY))
ALLOCATE(M_JA1(M_DUMMY))
ALLOCATE(M_JA2(M_DUMMY))
ALLOCATE(M_A1(M_DUMMY))

M_DUMMY=0

DO i=1,3*M_NumberofNodes+M_NumberofDE
DO j=1,M_K1GRCD(i)
IF (M_K1GR(i,j) /= 0.0D0) THEN
M_DUMMY=M_DUMMY+1
M_A1(M_DUMMY)=M_K1GR(i,j)
M_JA1(M_DUMMY)=M_K1GRCI(i,j)
M_JA2(M_DUMMY)=i
ELSEIF (M_K1GR(i,j) == 0.0D0) THEN
GOTO 4
ENDIF
4 ENDDO
ENDDO

DEALLOCATE(M_K1GR)
DEALLOCATE(M_K1GRCD)
DEALLOCATE(M_K1GRCI)

!CALL CPU_TIME (ElapsedTime)
!WRITE(20,*) "M_DUMMY  ",M_DUMMY,"  Elapsed Time=",ElapsedTime/60.0D0

M_IA1(1)=1

M_DUMMY2=1

DO i=2,M_DUMMY

IF (M_JA2(i)>M_JA2(i-1) .AND. i>1) THEN

M_DUMMY2=M_DUMMY2+1

M_IA1(M_DUMMY2)=i

ENDIF

ENDDO

M_DUMMY2=M_DUMMY2+1

M_IA1(M_DUMMY2)=M_DUMMY+1

ALLOCATE(M_IA(M_DUMMY2))

DO i=1,M_DUMMY2
M_IA(i)=M_IA1(i)
ENDDO

DEALLOCATE(M_IA1)

M_Dimen1=3*M_NumberofNodes+M_NumberofDE

CALL HarwellBoeing_Reduced(k,M_Dimen1,M_DUMMY,M_DUMMY2,M_A1,M_JA1,M_JA2,M_IA,M_F1Global,M_x)

!OPEN(27,FILE="3D_Output_Micro_F/A1JA1JA2.txt",STATUS="UNKNOWN")
!DO i=1,6300
!WRITE(27,29) M_A1(i),M_JA1(i),M_JA2(i)
!ENDDO
!CLOSE(27)
!
!29 FORMAT(E15.8,2X,I,2X,I,/)
!
!OPEN(39,FILE="3D_Output_Micro_F/IA.txt",STATUS="UNKNOWN")
!DO i=1,6300
!WRITE(39,33) M_IA(i)
!ENDDO
!CLOSE(39)
!
!33 FORMAT(I,/)

DEALLOCATE(M_F1Global)
DEALLOCATE(M_IA)
DEALLOCATE(M_JA1)
DEALLOCATE(M_JA2)
DEALLOCATE(M_A1)

DO j=1,M_NumberofNodes
M_Dis1(j)=M_x(3*(j-1)+1) !-M_x(3*(M_BoundaryIndex(1,1)-1)+1)
M_Dis1(M_NumberofNodes+j)=M_x(3*(j-1)+2) !-M_x(3*(M_BoundaryIndex(2,1)-1)+2)
M_Dis1(2*M_NumberofNodes+j)=M_x(3*(j-1)+3) !-M_x(3*(M_BoundaryIndex(3,1)-1)+3)
ENDDO

RETURN
ENDSUBROUTINE M_Matrix_DE_Reduced

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_OneStepDEs(M_Nodes,M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes,StrainsXX,&
StrainsYY,StrainsZZ,StrainsXY1,StrainsXY2,StrainsXZ1,StrainsXZ2,StrainsYZ1,StrainsYZ2,M_One_Step_DEs)
IMPLICIT NONE
INTEGER :: i,j
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofNodes,3) :: M_Nodes

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveXNodes) :: M_PositiveXNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveYNodes) :: M_PositiveYNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveZNodes) :: M_PositiveZNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeXNodes) :: M_NegativeXNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeYNodes) :: M_NegativeYNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeZNodes) :: M_NegativeZNodes

DOUBLE PRECISION, INTENT(INOUT), DIMENSION(3*(M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+M_NumberofPositiveZNodes),8) :: M_One_Step_DEs

DOUBLE PRECISION, INTENT(IN) :: StrainsXX,StrainsYY,StrainsZZ,StrainsXY1,StrainsXY2,StrainsXZ1,StrainsXZ2,StrainsYZ1,StrainsYZ2

INTEGER :: M_PosiTemp,M_NegaTemp

DO j=1,M_NumberofPositiveXNodes

M_PosiTemp=M_PositiveXNodes(j)
M_NegaTemp=M_NegativeXNodes(j)

IF ( 3*(M_PosiTemp-1)+1 < 3*(M_NegaTemp-1)+1 ) THEN

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_PositiveXNodes(j)
M_One_Step_DEs(j,3)=1
M_One_Step_DEs(j,4)=1.0D0
M_One_Step_DEs(j,5)=M_NegativeXNodes(j)
M_One_Step_DEs(j,6)=1
M_One_Step_DEs(j,7)=-1.0D0
M_One_Step_DEs(j,8)=StrainsXX*(M_Nodes(M_PositiveXNodes(j),2)-M_Nodes(M_NegativeXNodes(j),2))

ELSE

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_NegativeXNodes(j)
M_One_Step_DEs(j,3)=1
M_One_Step_DEs(j,4)=-1.0D0
M_One_Step_DEs(j,5)=M_PositiveXNodes(j)
M_One_Step_DEs(j,6)=1
M_One_Step_DEs(j,7)=1.0D0
M_One_Step_DEs(j,8)=StrainsXX*(M_Nodes(M_PositiveXNodes(j),2)-M_Nodes(M_NegativeXNodes(j),2))

ENDIF

ENDDO

DO j=M_NumberofPositiveXNodes+1,M_NumberofPositiveXNodes+M_NumberofPositiveXNodes

M_PosiTemp=M_PositiveXNodes(j-M_NumberofPositiveXNodes)
M_NegaTemp=M_NegativeXNodes(j-M_NumberofPositiveXNodes)

IF ( 3*(M_PosiTemp-1)+2 < 3*(M_NegaTemp-1)+2 ) THEN


M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_PositiveXNodes(j-M_NumberofPositiveXNodes)
M_One_Step_DEs(j,3)=2
M_One_Step_DEs(j,4)=1.0D0
M_One_Step_DEs(j,5)=M_NegativeXNodes(j-M_NumberofPositiveXNodes)
M_One_Step_DEs(j,6)=2
M_One_Step_DEs(j,7)=-1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsXY2*(M_Nodes(M_PositiveXNodes(j-M_NumberofPositiveXNodes),2)-M_Nodes(M_NegativeXNodes(j-M_NumberofPositiveXNodes),2))

ELSE

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_NegativeXNodes(j-M_NumberofPositiveXNodes)
M_One_Step_DEs(j,3)=2
M_One_Step_DEs(j,4)=-1.0D0
M_One_Step_DEs(j,5)=M_PositiveXNodes(j-M_NumberofPositiveXNodes)
M_One_Step_DEs(j,6)=2
M_One_Step_DEs(j,7)=1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsXY2*(M_Nodes(M_PositiveXNodes(j-M_NumberofPositiveXNodes),2)-M_Nodes(M_NegativeXNodes(j-M_NumberofPositiveXNodes),2))

ENDIF

ENDDO

DO j=2*M_NumberofPositiveXNodes+1,3*M_NumberofPositiveXNodes

M_PosiTemp=M_PositiveXNodes(j-2*M_NumberofPositiveXNodes)
M_NegaTemp=M_NegativeXNodes(j-2*M_NumberofPositiveXNodes)

IF ( 3*(M_PosiTemp-1)+3 < 3*(M_NegaTemp-1)+3 ) THEN

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_PositiveXNodes(j-2*M_NumberofPositiveXNodes)
M_One_Step_DEs(j,3)=3
M_One_Step_DEs(j,4)=1.0D0
M_One_Step_DEs(j,5)=M_NegativeXNodes(j-2*M_NumberofPositiveXNodes)
M_One_Step_DEs(j,6)=3
M_One_Step_DEs(j,7)=-1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsXZ2*(M_Nodes(M_PositiveXNodes(j-2*M_NumberofPositiveXNodes),2)-M_Nodes(M_NegativeXNodes(j-2*M_NumberofPositiveXNodes),2))

ELSE

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_NegativeXNodes(j-2*M_NumberofPositiveXNodes) 
M_One_Step_DEs(j,3)=3
M_One_Step_DEs(j,4)=-1.0D0
M_One_Step_DEs(j,5)=M_PositiveXNodes(j-2*M_NumberofPositiveXNodes)
M_One_Step_DEs(j,6)=3
M_One_Step_DEs(j,7)=+1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsXZ2*(M_Nodes(M_PositiveXNodes(j-2*M_NumberofPositiveXNodes),2)-M_Nodes(M_NegativeXNodes(j-2*M_NumberofPositiveXNodes),2))

ENDIF

ENDDO

DO j=3*M_NumberofPositiveXNodes+1,3*M_NumberofPositiveXNodes+M_NumberofPositiveYNodes

M_PosiTemp=M_PositiveYNodes(j-3*M_NumberofPositiveXNodes)
M_NegaTemp=M_NegativeYNodes(j-3*M_NumberofPositiveXNodes)

IF ( 3*(M_PosiTemp-1)+1 < 3*(M_NegaTemp-1)+1 ) THEN

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_PositiveYNodes(j-3*M_NumberofPositiveXNodes)
M_One_Step_DEs(j,3)=1
M_One_Step_DEs(j,4)=1.0D0
M_One_Step_DEs(j,5)=M_NegativeYNodes(j-3*M_NumberofPositiveXNodes)
M_One_Step_DEs(j,6)=1
M_One_Step_DEs(j,7)=-1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsXY1*(M_Nodes(M_PositiveYNodes(j-3*M_NumberofPositiveXNodes),3)-M_Nodes(M_NegativeYNodes(j-3*M_NumberofPositiveXNodes),3))

ELSE

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_NegativeYNodes(j-3*M_NumberofPositiveXNodes)
M_One_Step_DEs(j,3)=1
M_One_Step_DEs(j,4)=-1.0D0
M_One_Step_DEs(j,5)=M_PositiveYNodes(j-3*M_NumberofPositiveXNodes)
M_One_Step_DEs(j,6)=1
M_One_Step_DEs(j,7)=1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsXY1*(M_Nodes(M_PositiveYNodes(j-3*M_NumberofPositiveXNodes),3)-M_Nodes(M_NegativeYNodes(j-3*M_NumberofPositiveXNodes),3))

ENDIF

ENDDO

DO j=3*M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+1,3*M_NumberofPositiveXNodes+2*M_NumberofPositiveYNodes

M_PosiTemp=M_PositiveYNodes(j-3*M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)
M_NegaTemp=M_NegativeYNodes(j-3*M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)

IF ( 3*(M_PosiTemp-1)+2 < 3*(M_NegaTemp-1)+2 ) THEN

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_PositiveYNodes(j-3*M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)
M_One_Step_DEs(j,3)=2
M_One_Step_DEs(j,4)=1.0D0
M_One_Step_DEs(j,5)=M_NegativeYNodes(j-3*M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)
M_One_Step_DEs(j,6)=2
M_One_Step_DEs(j,7)=-1.0D0
M_One_Step_DEs(j,8)=StrainsYY*(M_Nodes(M_PositiveYNodes(j-3*M_NumberofPositiveXNodes-M_NumberofPositiveYNodes),3)&
-M_Nodes(M_NegativeYNodes(j-3*M_NumberofPositiveXNodes-M_NumberofPositiveYNodes),3))

ELSE

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_NegativeYNodes(j-3*M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)
M_One_Step_DEs(j,3)=2
M_One_Step_DEs(j,4)=-1.0D0
M_One_Step_DEs(j,5)=M_PositiveYNodes(j-3*M_NumberofPositiveXNodes-M_NumberofPositiveYNodes)
M_One_Step_DEs(j,6)=2
M_One_Step_DEs(j,7)=1.0D0
M_One_Step_DEs(j,8)=StrainsYY*(M_Nodes(M_PositiveYNodes(j-3*M_NumberofPositiveXNodes-M_NumberofPositiveYNodes),3)&
-M_Nodes(M_NegativeYNodes(j-3*M_NumberofPositiveXNodes-M_NumberofPositiveYNodes),3))

ENDIF

ENDDO

DO j=3*M_NumberofPositiveXNodes+2*M_NumberofPositiveYNodes+1,3*M_NumberofPositiveXNodes+3*M_NumberofPositiveYNodes

M_PosiTemp=M_PositiveYNodes(j-3*M_NumberofPositiveXNodes-2*M_NumberofPositiveYNodes)
M_NegaTemp=M_NegativeYNodes(j-3*M_NumberofPositiveXNodes-2*M_NumberofPositiveYNodes)

IF ( 3*(M_PosiTemp-1)+3 < 3*(M_NegaTemp-1)+3 ) THEN

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_PositiveYNodes(j-3*M_NumberofPositiveXNodes-2*M_NumberofPositiveYNodes)
M_One_Step_DEs(j,3)=3
M_One_Step_DEs(j,4)=1.0D0
M_One_Step_DEs(j,5)=M_NegativeYNodes(j-3*M_NumberofPositiveXNodes-2*M_NumberofPositiveYNodes)
M_One_Step_DEs(j,6)=3
M_One_Step_DEs(j,7)=-1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsYZ2*(M_Nodes(M_PositiveYNodes(j-3*M_NumberofPositiveXNodes-2*M_NumberofPositiveYNodes),3)&
-M_Nodes(M_NegativeYNodes(j-3*M_NumberofPositiveXNodes-2*M_NumberofPositiveYNodes),3))

ELSE

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_NegativeYNodes(j-3*M_NumberofPositiveXNodes-2*M_NumberofPositiveYNodes)
M_One_Step_DEs(j,3)=3
M_One_Step_DEs(j,4)=-1.0D0
M_One_Step_DEs(j,5)=M_PositiveYNodes(j-3*M_NumberofPositiveXNodes-2*M_NumberofPositiveYNodes)
M_One_Step_DEs(j,6)=3
M_One_Step_DEs(j,7)=1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsYZ2*(M_Nodes(M_PositiveYNodes(j-3*M_NumberofPositiveXNodes-2*M_NumberofPositiveYNodes),3)&
-M_Nodes(M_NegativeYNodes(j-3*M_NumberofPositiveXNodes-2*M_NumberofPositiveYNodes),3))

ENDIF

ENDDO

DO j=3*M_NumberofPositiveXNodes+3*M_NumberofPositiveYNodes+1,3*M_NumberofPositiveXNodes+3*M_NumberofPositiveYNodes+M_NumberofPositiveZNodes

M_PosiTemp=M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes)
M_NegaTemp=M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes)

IF ( 3*(M_PosiTemp-1)+1 < 3*(M_NegaTemp-1)+1 ) THEN

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes)
M_One_Step_DEs(j,3)=1
M_One_Step_DEs(j,4)=1.0D0
M_One_Step_DEs(j,5)=M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes)
M_One_Step_DEs(j,6)=1
M_One_Step_DEs(j,7)=-1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsXZ1*(M_Nodes(M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes),4)&
-M_Nodes(M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes),4))

ELSE

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes)
M_One_Step_DEs(j,3)=1
M_One_Step_DEs(j,4)=-1.0D0
M_One_Step_DEs(j,5)=M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes)
M_One_Step_DEs(j,6)=1
M_One_Step_DEs(j,7)=1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsXZ1*(M_Nodes(M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes),4)&
-M_Nodes(M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes),4))

ENDIF

ENDDO

DO j=3*M_NumberofPositiveXNodes+3*M_NumberofPositiveYNodes+M_NumberofPositiveZNodes+1,&
3*M_NumberofPositiveXNodes+3*M_NumberofPositiveYNodes+2*M_NumberofPositiveZNodes

M_PosiTemp=M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-M_NumberofPositiveZNodes)
M_NegaTemp=M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-M_NumberofPositiveZNodes)

IF ( 3*(M_PosiTemp-1)+2 < 3*(M_NegaTemp-1)+2 ) THEN

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-M_NumberofPositiveZNodes)
M_One_Step_DEs(j,3)=2
M_One_Step_DEs(j,4)=1.0D0
M_One_Step_DEs(j,5)=M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-M_NumberofPositiveZNodes)
M_One_Step_DEs(j,6)=2
M_One_Step_DEs(j,7)=-1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsYZ1*(M_Nodes(M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-M_NumberofPositiveZNodes),4)&
-M_Nodes(M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-M_NumberofPositiveZNodes),4))

ELSE

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-M_NumberofPositiveZNodes)
M_One_Step_DEs(j,3)=2
M_One_Step_DEs(j,4)=-1.0D0
M_One_Step_DEs(j,5)=M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-M_NumberofPositiveZNodes)
M_One_Step_DEs(j,6)=2
M_One_Step_DEs(j,7)=1.0D0
M_One_Step_DEs(j,8)=2.0D0*StrainsYZ1*(M_Nodes(M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-M_NumberofPositiveZNodes),4)&
-M_Nodes(M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-M_NumberofPositiveZNodes),4))

ENDIF

ENDDO

DO j=3*M_NumberofPositiveXNodes+3*M_NumberofPositiveYNodes+2*M_NumberofPositiveZNodes+1,&
3*M_NumberofPositiveXNodes+3*M_NumberofPositiveYNodes+3*M_NumberofPositiveZNodes

M_PosiTemp=M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-2*M_NumberofPositiveZNodes)
M_NegaTemp=M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-2*M_NumberofPositiveZNodes)

IF ( 3*(M_PosiTemp-1)+3 < 3*(M_NegaTemp-1)+3 ) THEN

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-2*M_NumberofPositiveZNodes)
M_One_Step_DEs(j,3)=3
M_One_Step_DEs(j,4)=1.0D0
M_One_Step_DEs(j,5)=M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-2*M_NumberofPositiveZNodes)
M_One_Step_DEs(j,6)=3
M_One_Step_DEs(j,7)=-1.0D0
M_One_Step_DEs(j,8)=StrainsZZ*(M_Nodes(M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-2*M_NumberofPositiveZNodes),4)&
-M_Nodes(M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-2*M_NumberofPositiveZNodes),4))

ELSE

M_One_Step_DEs(j,1)=2
M_One_Step_DEs(j,2)=M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-2*M_NumberofPositiveZNodes)
M_One_Step_DEs(j,3)=3
M_One_Step_DEs(j,4)=-1.0D0
M_One_Step_DEs(j,5)=M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-2*M_NumberofPositiveZNodes)
M_One_Step_DEs(j,6)=3
M_One_Step_DEs(j,7)=1.0D0
M_One_Step_DEs(j,8)=StrainsZZ*(M_Nodes(M_PositiveZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-2*M_NumberofPositiveZNodes),4)&
-M_Nodes(M_NegativeZNodes(j-3*M_NumberofPositiveXNodes-3*M_NumberofPositiveYNodes-2*M_NumberofPositiveZNodes),4))

ENDIF

ENDDO

ENDSUBROUTINE M_OneStepDEs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE Macro_Sub_Element(k,ii,ii3,ii1,ii2,M_Nodes,M_NodeIndex,M_DNDX,M_DNDY,M_DNDZ,M_Sub_CC,M_CC,M_Vol,M_Boundaries,M_BoundaryIndex,M_NumberofDE,&
M_K1GRR,M_K1GRCDD,M_K1GRCII,M_F1Globall,M_Connectivities,M_NN,M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,&
M_NegativeYNodes,M_NegativeZNodes,Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R,Sub_StrainsXZ1_R,&
Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,Sub_StrainsYZ2_R,M_Highest_GRCD,ME_sigmasigma,ME_gg,ME_Boundaries,ME_BoundaryIndex,ME_Connectivities,ME_Highest_GRCD,&
ME_Coordinates,M_NumberofElementsPerCNT,ME_ElementsPerCNT,ME_TunnelingIndex,ME_TunnelingIndexDimension,ME_Width1,ME_Width2,ME_CNTIndex,ME_Critical_Distance,&
ME_Coordinates2,E_sub_sigmasigma_in,Effec_Electro_Prop_Microscale_E,RotationMatrix2,BeginofModeIndex,EndofModeIndex,Connectivities)
IMPLICIT NONE

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec
!The electrical tunneling control flags in the microscale
INTEGER :: ME_Tunneling
DOUBLE PRECISION :: ME_Vf
DOUBLE PRECISION :: ME_RofCNT
DOUBLE PRECISION :: ME_Lambda
DOUBLE PRECISION :: ME_rho
INTEGER :: ME_BeginningPolymer
INTEGER :: ME_EndingPolymer
INTEGER :: ME_NumberofCNTs
INTEGER :: ME_NumberofCNTs2
INTEGER :: ME_MaxNumberofElementsPerCNT
COMMON/Micro_Electrical_More_Inputs/ME_Tunneling,ME_BeginningPolymer,ME_EndingPolymer,ME_NumberofCNTs,ME_NumberofCNTs2,&
ME_MaxNumberofElementsPerCNT,ME_Vf,ME_RofCNT,ME_Lambda,ME_rho

INTEGER :: i,j
INTEGER, INTENT(IN) :: k,ii,ii3,ii1,ii2
INTEGER, INTENT(IN) :: M_NumberofDE
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofNodes,3) :: M_Nodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveXNodes) :: M_PositiveXNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveYNodes) :: M_PositiveYNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofPositiveZNodes) :: M_PositiveZNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeXNodes) :: M_NegativeXNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeYNodes) :: M_NegativeYNodes
INTEGER, INTENT(IN), DIMENSION(M_NumberofNegativeZNodes) :: M_NegativeZNodes

DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsXX_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsYY_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsZZ_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsXY1_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsXY2_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsXZ1_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsXZ2_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsYZ1_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsYZ2_R
DOUBLE PRECISION :: StrainsXX
DOUBLE PRECISION :: StrainsYY
DOUBLE PRECISION :: StrainsZZ
DOUBLE PRECISION :: StrainsXY1
DOUBLE PRECISION :: StrainsXY2
DOUBLE PRECISION :: StrainsXZ1
DOUBLE PRECISION :: StrainsXZ2
DOUBLE PRECISION :: StrainsYZ1
DOUBLE PRECISION :: StrainsYZ2
DOUBLE PRECISION, DIMENSION(3*(M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+M_NumberofPositiveZNodes),8) :: M_One_Step_DEs
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofBoundaryNodes,TotalTimeSteps+2) :: M_Boundaries
INTEGER, INTENT(IN), DIMENSION(M_NumberofBoundaryNodes,2) :: M_BoundaryIndex
INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,10) :: M_Connectivities
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec) :: M_K1GRR
INTEGER, INTENT(IN), DIMENSION (3*M_NumberofNodes+M_NumberofDE) :: M_K1GRCDD
INTEGER, INTENT(IN), DIMENSION (3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec) :: M_K1GRCII
DOUBLE PRECISION, INTENT(IN), DIMENSION (3*M_NumberofNodes+M_NumberofDE) :: M_F1Globall
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofGaussPoints,M_NumberofElements) :: M_Vol
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofMaterialModes,6,6) :: M_CC
DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofGaussPoints,M_NumberofMaterialModes,6,6) :: M_Sub_CC
INTEGER, INTENT(IN), DIMENSION(M_NumberofNodes) :: M_NodeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofGaussPoints,M_NumberofElements,8+3) ::  M_DNDX,M_DNDY,M_DNDZ
DOUBLE PRECISION, INTENT(IN), DIMENSION(M_NumberofGaussPoints,8) :: M_NN
INTEGER, INTENT(OUT) :: M_Highest_GRCD
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_Vol_Ele
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StrainsXX
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StrainsYY
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StrainsZZ
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StrainsXY
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StrainsXZ
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StrainsYZ
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StrainsVonMises
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StressesXX
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StressesYY
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StressesZZ
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StressesXY
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StressesXZ
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StressesYZ
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: M_StressesVonMises
DOUBLE PRECISION, DIMENSION (M_NumberofNodes,4)  :: M_CurrentNodes 

DOUBLE PRECISION, DIMENSION (M_NumberofElements,8) :: M_Sub_StrainsXX,M_Sub_StrainsYY,M_Sub_StrainsZZ,M_Sub_StrainsXY1,M_Sub_StrainsXZ1,M_Sub_StrainsYZ1,&
M_Sub_StrainsXY2,M_Sub_StrainsXZ2,M_Sub_StrainsYZ2,M_Sub_StressesXX,M_Sub_StressesYY,M_Sub_StressesZZ,M_Sub_StressesXY,M_Sub_StressesXZ,M_Sub_StressesYZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_Dis1
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofMaterialModes,3,3) :: ME_sigmasigma
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofMaterialModes,6,6) :: ME_gg
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofBoundaryNodes,TotalTimeSteps+1) :: ME_Boundaries
INTEGER, INTENT(IN), DIMENSION(ME_NumberofBoundaryNodes,1) :: ME_BoundaryIndex
INTEGER, INTENT(IN), DIMENSION(M_NumberofElements,10) :: ME_Connectivities
INTEGER :: iiii,ZeroTimeStep
DOUBLE PRECISION :: Ex_Knot,Ey_Knot,Ez_Knot
DOUBLE PRECISION, DIMENSION(M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+M_NumberofPositiveZNodes,6) :: ME_One_Step_DEs
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: ME_Dis1
INTEGER, INTENT(OUT) :: ME_Highest_GRCD
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs) :: M_NumberofElementsPerCNT
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofCNTs2,2) :: ME_Coordinates
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs,ME_MaxNumberofElementsPerCNT) :: ME_ElementsPerCNT
DOUBLE PRECISION, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,4,3) :: ME_TunnelingIndex    !Changed
INTEGER, INTENT(IN), DIMENSION(ME_NumberofGaussPoints,M_NumberofElements) :: ME_TunnelingIndexDimension  !Changed
INTEGER, INTENT(IN), DIMENSION(ME_NumberofCNTs2) :: ME_CNTIndex
DOUBLE PRECISION, INTENT(IN) :: ME_Critical_Distance
DOUBLE PRECISION, INTENT(IN), DIMENSION(9,ME_NumberofCNTs2,2) :: ME_Coordinates2
DOUBLE PRECISION, INTENT(IN) :: ME_Width1,ME_Width2
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_rhorho
DOUBLE PRECISION, DIMENSION(ME_NumberofGaussPoints,M_NumberofElements,3,3) :: ME_sub_Tunneling_rhorho
DOUBLE PRECISION, INTENT(OUT), DIMENSION(TotalTimeSteps+1,E_NumberofGaussPoints,NumberofElements,3,3) :: E_sub_sigmasigma_in
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofMaterialModes,3,3) :: RotationMatrix2
INTEGER, INTENT(IN) :: BeginofModeIndex,EndofModeIndex
INTEGER, INTENT(IN), DIMENSION(NumberofElements,17) :: Connectivities

!DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: ME_sub_sigmasigma
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: ME_sub_sigmasigma_updated
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ME_DNDX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ME_DNDY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ME_DNDZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ME_Vol
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ME_KComponent
INTEGER :: ME_NumberofDE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ME_K1GRR
INTEGER, ALLOCATABLE, DIMENSION(:) :: ME_K1GRCDD
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ME_K1GRCII
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ME_F1Globall

DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: ME_Vol_Ele
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: ME_EXX
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: ME_EYY
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: ME_EZZ
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: ME_JXX
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: ME_JYY
DOUBLE PRECISION, DIMENSION (M_NumberofElements) :: ME_JZZ

DOUBLE PRECISION, DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_EXX
DOUBLE PRECISION, DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_EYY
DOUBLE PRECISION, DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_EZZ
DOUBLE PRECISION, DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_JXX
DOUBLE PRECISION, DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_JYY
DOUBLE PRECISION, DIMENSION (M_NumberofElements,ME_NumberofGaussPoints) :: ME_Sub_JZZ

DOUBLE PRECISION :: ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff
CHARACTER*90, INTENT(IN) :: Effec_Electro_Prop_Microscale_E
DOUBLE PRECISION, DIMENSION(3,3) :: ME_Rotated_Sigma_Eff
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!

!StrainsXX=Sub_StrainsXX_R(ii1,ii2)
!StrainsYY=Sub_StrainsYY_R(ii1,ii2)
!StrainsZZ=Sub_StrainsZZ_R(ii1,ii2)
!StrainsXY1=Sub_StrainsXY1_R(ii1,ii2)
!StrainsXY2=Sub_StrainsXY2_R(ii1,ii2)
!StrainsXZ1=Sub_StrainsXZ1_R(ii1,ii2)
!StrainsXZ2=Sub_StrainsXZ2_R(ii1,ii2)
!StrainsYZ1=Sub_StrainsYZ1_R(ii1,ii2)
!StrainsYZ2=Sub_StrainsYZ2_R(ii1,ii2)

!X_r -> X^{Tildea}_3, X_t -> X^{Tildea}_1, X_z -> X^{Tildea}_2
StrainsZZ=Sub_StrainsXX_R(ii1,ii2)
StrainsXX=Sub_StrainsYY_R(ii1,ii2)
StrainsYY=Sub_StrainsZZ_R(ii1,ii2)
StrainsXZ2=Sub_StrainsXY1_R(ii1,ii2)
StrainsXZ1=Sub_StrainsXY2_R(ii1,ii2)
StrainsYZ2=Sub_StrainsXZ1_R(ii1,ii2)
StrainsYZ1=Sub_StrainsXZ2_R(ii1,ii2)
StrainsXY1=Sub_StrainsYZ1_R(ii1,ii2)
StrainsXY2=Sub_StrainsYZ2_R(ii1,ii2)

!StrainsZZ=0.001D0
!StrainsXX=0.001D0
!StrainsYY=0.001D0
!StrainsXZ2=0.0005D0
!StrainsXZ1=0.0005D0
!StrainsYZ2=0.0005D0
!StrainsYZ1=0.0005D0
!StrainsXY1=0.0005D0
!StrainsXY2=0.0005D0

!WRITE(*,*) StrainsZZ,StrainsXX,StrainsYY,StrainsXZ2,StrainsXZ1,StrainsYZ2,StrainsYZ1,StrainsXY1,StrainsXY2

CALL M_OneStepDEs(M_Nodes,M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes,&
StrainsXX,StrainsYY,StrainsZZ,StrainsXY1,StrainsXY2,StrainsXZ1,StrainsXZ2,StrainsYZ1,StrainsYZ2,M_One_Step_DEs)

ALLOCATE(M_Dis1(3*M_NumberofNodes+M_NumberofDE))
CALL M_Matrix_DE_Reduced(k,M_NumberofDE,M_One_Step_DEs,M_BoundaryIndex,M_Boundaries,M_K1GRR,M_K1GRCDD,M_K1GRCII,M_F1Globall,M_Connectivities,&
M_Dis1,M_Highest_GRCD)
CALL M_StrainsandStresses (M_Vol,M_Vol_Ele,M_CC,M_Sub_CC,M_Connectivities,M_Nodes,M_NodeIndex,M_DNDX,M_DNDY,M_DNDZ,&
M_Dis1,M_StrainsXX,M_StrainsYY,M_StrainsZZ,M_StrainsXY,M_StrainsXZ,M_StrainsYZ,M_StrainsVonMises,M_StressesXX,M_StressesYY,&
M_StressesZZ,M_StressesXY,M_StressesXZ,M_StressesYZ,M_StressesVonMises,M_Sub_StrainsXX,M_Sub_StrainsYY,M_Sub_StrainsZZ,&
M_Sub_StrainsXY1,M_Sub_StrainsXZ1,M_Sub_StrainsYZ1,M_Sub_StrainsXY2,M_Sub_StrainsXZ2,M_Sub_StrainsYZ2,M_Sub_StressesXX,&
M_Sub_StressesYY,M_Sub_StressesZZ,M_Sub_StressesXY,M_Sub_StressesXZ,M_Sub_StressesYZ,M_CurrentNodes)

IF( ii1==ElementShown .AND. ii2==SubElementShown) THEN
CALL M_FeedForTecplot (k,ii,M_Connectivities,M_Nodes,M_NodeIndex,M_Dis1,M_StrainsXX,M_StrainsYY,M_StrainsZZ,M_StrainsXY,M_StrainsXZ,M_StrainsYZ,&
M_StrainsVonMises,M_StressesXX,M_StressesYY,M_StressesZZ,M_StressesXY,M_StressesXZ,M_StressesYZ,M_StressesVonMises,M_CurrentNodes)
ENDIF

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
IF (k==1) THEN
ZeroTimeStep=1
ALLOCATE(ME_sub_sigmasigma_updated(ME_NumberofGaussPoints,M_NumberofElements,3,3))
ALLOCATE(ME_DNDX(ME_NumberofGaussPoints,M_NumberofElements,8))
ALLOCATE(ME_DNDY(ME_NumberofGaussPoints,M_NumberofElements,8))
ALLOCATE(ME_DNDZ(ME_NumberofGaussPoints,M_NumberofElements,8))
ALLOCATE(ME_Vol(ME_NumberofGaussPoints,M_NumberofElements))
ALLOCATE(ME_KComponent(M_NumberofElements,8,8))
CALL ME_Components_Hexahedron(M_Nodes,ME_Connectivities,ME_sigmasigma,ME_gg,ME_sub_sigmasigma_updated,ME_Vol,&
ME_KComponent,ME_DNDX,ME_DNDY,ME_DNDZ,ZeroTimeStep,k,M_NumberofDE,Sub_StrainsXX_R,Sub_StrainsYY_R,&
Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R,Sub_StrainsXZ1_R,Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,&
Sub_StrainsYZ2_R,ME_TunnelingIndex,ME_TunnelingIndexDimension,ME_CNTIndex,ME_Critical_Distance,ME_Coordinates2,&
ME_Width1,ME_Width2,M_NumberofElementsPerCNT,ME_ElementsPerCNT,M_Vol,M_NN,M_Dis1,M_CurrentNodes,&
M_Sub_StrainsXX,M_Sub_StrainsYY,M_Sub_StrainsZZ,M_Sub_StrainsXY1,M_Sub_StrainsXZ1,M_Sub_StrainsYZ1,&
M_Sub_StrainsXY2,M_Sub_StrainsXZ2,M_Sub_StrainsYZ2)

ME_NumberofDE=1*(M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+M_NumberofPositiveZNodes)
ALLOCATE(ME_K1GRR(M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec))
ALLOCATE(ME_K1GRCDD(M_NumberofNodes+ME_NumberofDE))
ALLOCATE(ME_K1GRCII(M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec))
ALLOCATE(ME_F1Globall(M_NumberofNodes+ME_NumberofDE))
CALL ME_KPreprations(ME_Connectivities,ME_KComponent,ME_NumberofDE,ME_K1GRR,ME_K1GRCDD,ME_K1GRCII,ME_F1Globall)

DO iiii=1,5

IF (iiii==1) THEN ! For sigma_xx and sigma_zz (the rule of mixtures)
Ex_Knot=1.0D0
Ey_Knot=0.0D0
Ez_Knot=0.0D0
ELSEIF (iiii==2) THEN ! For sigma_yy
Ex_Knot=0.0D0
Ey_Knot=1.0D0
Ez_Knot=0.0D0
ELSEIF (iiii==3) THEN ! For sigma_xy
Ex_Knot=1.0D0
Ey_Knot=1.0D0
Ez_Knot=0.0D0
ELSEIF (iiii==4) THEN ! For sigma_xz
Ex_Knot=1.0D0
Ey_Knot=0.0D0
Ez_Knot=1.0D0
ELSEIF (iiii==5) THEN ! For sigma_yz
Ex_Knot=0.0D0
Ey_Knot=1.0D0
Ez_Knot=1.0D0
ENDIF

CALL ME_OneStepDEs(M_Nodes,M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes,&
Ex_Knot,Ey_Knot,Ez_Knot,ME_NumberofDE,ME_One_Step_DEs)

ALLOCATE(ME_Dis1(M_NumberofNodes+ME_NumberofDE))
CALL ME_Matrix_DE_Reduced(k,ME_NumberofDE,ME_One_Step_DEs,ME_BoundaryIndex,ME_Boundaries,ME_K1GRR,ME_K1GRCDD,ME_K1GRCII,ME_F1Globall,ME_Connectivities,&
ME_Dis1,ME_Highest_GRCD)

CALL ME_EandJ (ME_Vol,ME_Vol_Ele,ME_sigmasigma,ME_sub_sigmasigma_updated,ME_Connectivities,M_Nodes,M_NodeIndex,ME_Dis1,ME_EXX,ME_EYY,ME_EZZ,ME_JXX,&
ME_JYY,ME_JZZ,ME_DNDX,ME_DNDY,ME_DNDZ,ME_Sub_EXX,ME_Sub_EYY,ME_Sub_EZZ,ME_Sub_JXX,ME_Sub_JYY,ME_Sub_JZZ)

CALL ME_EffectiveProperties(ZeroTimeStep,k,ii,ii3,ii1,ii2,iiii,Ex_Knot,Ey_Knot,Ez_Knot,M_Nodes,ME_Vol,ME_Sub_EXX,ME_Sub_EYY,ME_Sub_EZZ,ME_Sub_JXX,&
ME_Sub_JYY,ME_Sub_JZZ,M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes,ME_sub_sigmasigma_updated,&
ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff,Effec_Electro_Prop_Microscale_E,StrainsZZ)

!IF (iiii==5) THEN
!WRITE(*,*) "11,22,33,12,13,23",ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff
!ENDIF

IF( ii1==ElementShown .AND. ii2==SubElementShown) THEN
CALL ME_FeedForTecplot (ZeroTimeStep,k,ii,iiii,ME_Connectivities,M_Nodes,M_NodeIndex,ME_Dis1,ME_EXX,ME_EYY,ME_EZZ,ME_JXX,ME_JYY,ME_JZZ,M_CurrentNodes,&
ME_sigmasigma,ME_sub_sigmasigma_updated)
ENDIF

DEALLOCATE(ME_Dis1)

ENDDO !For iiii

CALL ME_RotatedConductivities(ii1,ii2,Connectivities,BeginofModeIndex,EndofModeIndex,RotationMatrix2,ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,&
ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff,ME_Rotated_Sigma_Eff)

!E_sub_sigmasigma_in(k,ii1,ii2,1,1)=ME_Sigma_XX_Eff
!E_sub_sigmasigma_in(k,ii1,ii2,2,2)=ME_Sigma_YY_Eff
!E_sub_sigmasigma_in(k,ii1,ii2,3,3)=ME_Sigma_ZZ_Eff
!E_sub_sigmasigma_in(k,ii1,ii2,1,2)=ME_Sigma_XY_Eff
!E_sub_sigmasigma_in(k,ii1,ii2,1,3)=ME_Sigma_XZ_Eff
!E_sub_sigmasigma_in(k,ii1,ii2,2,3)=ME_Sigma_YZ_Eff
!E_sub_sigmasigma_in(k,ii1,ii2,2,1)=E_sub_sigmasigma_in(k,ii1,ii2,1,2)
!E_sub_sigmasigma_in(k,ii1,ii2,3,1)=E_sub_sigmasigma_in(k,ii1,ii2,1,3)
!E_sub_sigmasigma_in(k,ii1,ii2,3,2)=E_sub_sigmasigma_in(k,ii1,ii2,2,3)

E_sub_sigmasigma_in(k,ii2,ii1,1,1)=ME_Rotated_Sigma_Eff(1,1)
E_sub_sigmasigma_in(k,ii2,ii1,2,2)=ME_Rotated_Sigma_Eff(2,2)
E_sub_sigmasigma_in(k,ii2,ii1,3,3)=ME_Rotated_Sigma_Eff(3,3)
E_sub_sigmasigma_in(k,ii2,ii1,1,2)=ME_Rotated_Sigma_Eff(1,2)
E_sub_sigmasigma_in(k,ii2,ii1,1,3)=ME_Rotated_Sigma_Eff(1,3)
E_sub_sigmasigma_in(k,ii2,ii1,2,3)=ME_Rotated_Sigma_Eff(2,3)
E_sub_sigmasigma_in(k,ii2,ii1,2,1)=ME_Rotated_Sigma_Eff(2,1)
E_sub_sigmasigma_in(k,ii2,ii1,3,1)=ME_Rotated_Sigma_Eff(3,1)
E_sub_sigmasigma_in(k,ii2,ii1,3,2)=ME_Rotated_Sigma_Eff(3,2)

DEALLOCATE(ME_sub_sigmasigma_updated)
DEALLOCATE(ME_DNDX)
DEALLOCATE(ME_DNDY)
DEALLOCATE(ME_DNDZ)
DEALLOCATE(ME_Vol)
DEALLOCATE(ME_KComponent)
DEALLOCATE(ME_K1GRR)
DEALLOCATE(ME_K1GRCDD)
DEALLOCATE(ME_K1GRCII)
DEALLOCATE(ME_F1Globall)

ENDIF !For k
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!
ZeroTimeStep=0
ALLOCATE(ME_sub_sigmasigma_updated(ME_NumberofGaussPoints,M_NumberofElements,3,3))
ALLOCATE(ME_DNDX(ME_NumberofGaussPoints,M_NumberofElements,8))
ALLOCATE(ME_DNDY(ME_NumberofGaussPoints,M_NumberofElements,8))
ALLOCATE(ME_DNDZ(ME_NumberofGaussPoints,M_NumberofElements,8))
ALLOCATE(ME_Vol(ME_NumberofGaussPoints,M_NumberofElements))
ALLOCATE(ME_KComponent(M_NumberofElements,8,8))
CALL ME_Components_Hexahedron(M_Nodes,ME_Connectivities,ME_sigmasigma,ME_gg,ME_sub_sigmasigma_updated,ME_Vol,&
ME_KComponent,ME_DNDX,ME_DNDY,ME_DNDZ,ZeroTimeStep,k,M_NumberofDE,Sub_StrainsXX_R,Sub_StrainsYY_R,&
Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R,Sub_StrainsXZ1_R,Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,&
Sub_StrainsYZ2_R,ME_TunnelingIndex,ME_TunnelingIndexDimension,ME_CNTIndex,ME_Critical_Distance,ME_Coordinates2,&
ME_Width1,ME_Width2,M_NumberofElementsPerCNT,ME_ElementsPerCNT,M_Vol,M_NN,M_Dis1,M_CurrentNodes,&
M_Sub_StrainsXX,M_Sub_StrainsYY,M_Sub_StrainsZZ,M_Sub_StrainsXY1,M_Sub_StrainsXZ1,M_Sub_StrainsYZ1,&
M_Sub_StrainsXY2,M_Sub_StrainsXZ2,M_Sub_StrainsYZ2)

ME_NumberofDE=1*(M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+M_NumberofPositiveZNodes)
ALLOCATE(ME_K1GRR(M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec))
ALLOCATE(ME_K1GRCDD(M_NumberofNodes+ME_NumberofDE))
ALLOCATE(ME_K1GRCII(M_NumberofNodes+ME_NumberofDE,ME_MaxDimenofInverseConnec))
ALLOCATE(ME_F1Globall(M_NumberofNodes+ME_NumberofDE))
CALL ME_KPreprations(ME_Connectivities,ME_KComponent,ME_NumberofDE,ME_K1GRR,ME_K1GRCDD,ME_K1GRCII,ME_F1Globall)

DO iiii=1,5

IF (iiii==1) THEN ! For sigma_xx and sigma_zz (the rule of mixtures)
Ex_Knot=1.0D0
Ey_Knot=0.0D0
Ez_Knot=0.0D0
ELSEIF (iiii==2) THEN ! For sigma_yy
Ex_Knot=0.0D0
Ey_Knot=1.0D0
Ez_Knot=0.0D0
ELSEIF (iiii==3) THEN ! For sigma_xy
Ex_Knot=1.0D0
Ey_Knot=1.0D0
Ez_Knot=0.0D0
ELSEIF (iiii==4) THEN ! For sigma_xz
Ex_Knot=1.0D0
Ey_Knot=0.0D0
Ez_Knot=1.0D0
ELSEIF (iiii==5) THEN ! For sigma_yz
Ex_Knot=0.0D0
Ey_Knot=1.0D0
Ez_Knot=1.0D0
ENDIF

CALL ME_OneStepDEs(M_Nodes,M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes,&
Ex_Knot,Ey_Knot,Ez_Knot,ME_NumberofDE,ME_One_Step_DEs)

ALLOCATE(ME_Dis1(M_NumberofNodes+ME_NumberofDE))
CALL ME_Matrix_DE_Reduced(k,ME_NumberofDE,ME_One_Step_DEs,ME_BoundaryIndex,ME_Boundaries,ME_K1GRR,ME_K1GRCDD,ME_K1GRCII,ME_F1Globall,ME_Connectivities,&
ME_Dis1,ME_Highest_GRCD)

CALL ME_EandJ (ME_Vol,ME_Vol_Ele,ME_sigmasigma,ME_sub_sigmasigma_updated,ME_Connectivities,M_Nodes,M_NodeIndex,ME_Dis1,ME_EXX,ME_EYY,ME_EZZ,ME_JXX,&
ME_JYY,ME_JZZ,ME_DNDX,ME_DNDY,ME_DNDZ,ME_Sub_EXX,ME_Sub_EYY,ME_Sub_EZZ,ME_Sub_JXX,ME_Sub_JYY,ME_Sub_JZZ)

CALL ME_EffectiveProperties(ZeroTimeStep,k,ii,ii3,ii1,ii2,iiii,Ex_Knot,Ey_Knot,Ez_Knot,M_Nodes,ME_Vol,ME_Sub_EXX,ME_Sub_EYY,ME_Sub_EZZ,ME_Sub_JXX,&
ME_Sub_JYY,ME_Sub_JZZ,M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes,ME_sub_sigmasigma_updated,&
ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff,Effec_Electro_Prop_Microscale_E,StrainsZZ)

!IF (iiii==5) THEN
!WRITE(*,*) "11,22,33,12,13,23",ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff
!ENDIF

IF( ii1==ElementShown .AND. ii2==SubElementShown) THEN
CALL ME_FeedForTecplot (ZeroTimeStep,k,ii,iiii,ME_Connectivities,M_Nodes,M_NodeIndex,ME_Dis1,ME_EXX,ME_EYY,ME_EZZ,ME_JXX,ME_JYY,ME_JZZ,M_CurrentNodes,&
ME_sigmasigma,ME_sub_sigmasigma_updated)
ENDIF

DEALLOCATE(ME_Dis1)

ENDDO !For iiii

CALL ME_RotatedConductivities(ii1,ii2,Connectivities,BeginofModeIndex,EndofModeIndex,RotationMatrix2,ME_Sigma_XX_Eff,ME_Sigma_YY_Eff,ME_Sigma_ZZ_Eff,&
ME_Sigma_XY_Eff,ME_Sigma_XZ_Eff,ME_Sigma_YZ_Eff,ME_Rotated_Sigma_Eff)

E_sub_sigmasigma_in(k+1,ii2,ii1,1,1)=ME_Rotated_Sigma_Eff(1,1)
E_sub_sigmasigma_in(k+1,ii2,ii1,2,2)=ME_Rotated_Sigma_Eff(2,2)
E_sub_sigmasigma_in(k+1,ii2,ii1,3,3)=ME_Rotated_Sigma_Eff(3,3)
E_sub_sigmasigma_in(k+1,ii2,ii1,1,2)=ME_Rotated_Sigma_Eff(1,2)
E_sub_sigmasigma_in(k+1,ii2,ii1,1,3)=ME_Rotated_Sigma_Eff(1,3)
E_sub_sigmasigma_in(k+1,ii2,ii1,2,3)=ME_Rotated_Sigma_Eff(2,3)
E_sub_sigmasigma_in(k+1,ii2,ii1,2,1)=ME_Rotated_Sigma_Eff(2,1)
E_sub_sigmasigma_in(k+1,ii2,ii1,3,1)=ME_Rotated_Sigma_Eff(3,1)
E_sub_sigmasigma_in(k+1,ii2,ii1,3,2)=ME_Rotated_Sigma_Eff(3,2)

DEALLOCATE(ME_sub_sigmasigma_updated)
DEALLOCATE(ME_DNDX)
DEALLOCATE(ME_DNDY)
DEALLOCATE(ME_DNDZ)
DEALLOCATE(ME_Vol)
DEALLOCATE(ME_KComponent)
DEALLOCATE(ME_K1GRR)
DEALLOCATE(ME_K1GRCDD)
DEALLOCATE(ME_K1GRCII)
DEALLOCATE(ME_F1Globall)
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!

ENDSUBROUTINE Macro_Sub_Element

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INV(k,n,A,Ainv)
INTEGER, INTENT(IN) :: n,k
DOUBLE PRECISION,  intent(in), dimension(n,n) :: A
DOUBLE PRECISION,  intent(out), dimension(size(A,1),size(A,2)) :: Ainv

DOUBLE PRECISION,  dimension(size(A,1)) :: work  ! work array for LAPACK
integer, dimension(size(A,1)) :: ipiv   ! pivot indices
integer :: info

  ! External procedures defined in LAPACK
  external DGETRF
  external DGETRI

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A 
  
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(n, n, Ainv, n, ipiv, info)

  if (info /= 0) then
     write(*,*) k
     stop 'Matrix is numerically singular!'
  end if

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(n, Ainv, n, ipiv, work, n, info)

  if (info /= 0) then
     stop 'Matrix inversion failed!'
  end if
ENDSUBROUTINE INV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_Components(M_Nodes,M_Connectivities,M_CC,M_Sub_CC,M_Vol,M_KComponent,M_DNDX,M_DNDY,M_DNDZ,M_NN)
IMPLICIT NONE
INTEGER :: i,j,k,i1,i2,i3,i4,j1

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

DOUBLE PRECISION, INTENT(IN), DIMENSION (M_NumberofNodes,4) :: M_Nodes
INTEGER, INTENT(IN), DIMENSION (M_NumberofElements,10) :: M_Connectivities
DOUBLE PRECISION, INTENT(IN),  DIMENSION (M_NumberofMaterialModes,6,6) :: M_CC
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofGaussPoints,M_NumberofMaterialModes,6,6) :: M_Sub_CC
DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofGaussPoints,M_NumberofElements) :: M_Vol
DOUBLE PRECISION, DIMENSION (M_NumberofGaussPoints,M_NumberofElements,6,24) :: M_BB
DOUBLE PRECISION, DIMENSION (M_NumberofGaussPoints,M_NumberofElements,24,6) :: M_BBT

DOUBLE PRECISION, DIMENSION (M_NumberofGaussPoints,M_NumberofElements,6,9) :: M_BBI
DOUBLE PRECISION, DIMENSION (M_NumberofGaussPoints,M_NumberofElements,9,6) :: M_BBIT

DOUBLE PRECISION, DIMENSION (M_NumberofGaussPoints,M_NumberofElements,6,9) :: M_BBIC

DOUBLE PRECISION, DIMENSION (M_NumberofGaussPoints,M_NumberofElements,6,9) :: M_BBIBar
DOUBLE PRECISION, DIMENSION (M_NumberofGaussPoints,M_NumberofElements,9,6) :: M_BBIBarT

DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofElements,24,24) :: M_KComponent

DOUBLE PRECISION, DIMENSION (M_NumberofElements,24,24) :: M_KCC
DOUBLE PRECISION, DIMENSION (M_NumberofElements,24,9)  :: M_KCI
DOUBLE PRECISION, DIMENSION (M_NumberofElements,9,24)  :: M_KIC
DOUBLE PRECISION, DIMENSION (M_NumberofElements,9,9)   :: M_KII
DOUBLE PRECISION, DIMENSION (M_NumberofElements,9,9)   :: M_KII_Inverse
DOUBLE PRECISION, DIMENSION (M_NumberofElements,9,9)   :: M_KII2

DOUBLE PRECISION, DIMENSION(M_NumberofElements,8):: M_x,M_y,M_z

DOUBLE PRECISION, DIMENSION(8) :: M_Xi,M_Eta,M_Mu

DOUBLE PRECISION, DIMENSION(M_NumberofGaussPoints) :: M_Xi1,M_Eta1,M_Mu1,M_W

DOUBLE PRECISION, DIMENSION(M_NumberofGaussPoints,M_NumberofElements,3,3) :: M_JJ,M_InverseJJ

DOUBLE PRECISION, DIMENSION(M_NumberofGaussPoints,8+3) :: M_DNDXi,M_DNDEta,M_DNDMu
DOUBLE PRECISION, INTENT(OUT), DIMENSION(M_NumberofGaussPoints,8) :: M_NN
DOUBLE PRECISION, INTENT(OUT), DIMENSION(M_NumberofGaussPoints,M_NumberofElements,8+3) :: M_DNDX,M_DNDY,M_DNDZ
DOUBLE PRECISION, DIMENSION(M_NumberofGaussPoints,M_NumberofElements,8+3) :: M_QX,M_QY,M_QZ

DOUBLE PRECISION :: M_Delta_A

DOUBLE PRECISION :: a11,a12,a13,a21,a22,a23,a31,a32,a33

DOUBLE PRECISION :: M_Vol_Tep

!DOUBLE PRECISION :: M_Vf,M_R,M_MaximumRadius

!INTEGER :: M_BeginofModeIndex,M_EndofModeIndex,M_Anisotropic

!INTEGER, INTENT(OUT) :: M_TestIndex1,M_TestIndex2

!DOUBLE PRECISION, DIMENSION(11,2) :: M_XYofAxis

DO i=1,M_NumberofGaussPoints
DO j=1,M_NumberofMaterialModes
DO i2=1,6
DO i3=1,6
M_Sub_CC(i,j,i2,i3)=M_CC(j,i2,i3)
ENDDO
ENDDO
ENDDO
ENDDO


M_KII2=0.0D0

M_JJ=0.0D0
M_InverseJJ=0.0D0

M_BB=0.0D0
M_BBT=0.0D0
M_KComponent=0.0D0
M_KCI=0.0D0
M_KIC=0.0D0
M_KCC=0.0D0
M_KII=0.0D0
M_KII_Inverse=0.0D0
M_BBI=0.0D0
M_BBIC=0.0D0
M_BBIBar=0.0D0
M_BBIBarT=0.0D0

M_Xi=(/1.0D0,+1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,-1.0D0,-1.0D0/)
M_Eta=(/-1.0D0,1.0D0,1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,-1.0D0/)
M_Mu=(/-1.0D0,-1.0D0,-1.0D0,-1.0D0,1.0D0,1.0D0,1.0D0,1.0D0/)

DO i=1,M_NumberofGaussPoints
M_Xi1(i)=M_Xi(i)/SQRT(3.0D0)
M_Eta1(i)=M_Eta(i)/SQRT(3.0D0)
M_Mu1(i)=M_Mu(i)/SQRT(3.0D0)
M_W(i)=1.0D0
ENDDO

DO k=1,M_NumberofElements
DO j=1,8
M_x(k,j)=M_Nodes(M_Connectivities(k,j+1),2)
M_y(k,j)=M_Nodes(M_Connectivities(k,j+1),3)
M_z(k,j)=M_Nodes(M_Connectivities(k,j+1),4)
ENDDO
ENDDO

DO i=1,M_NumberofGaussPoints

M_DNDXi(i,4)=-1.0D0/8.0D0*(1.0D0-M_Eta1(i))*(1.0D0-M_Mu1(i))
M_DNDEta(i,4)=-1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0-M_Mu1(i))
M_DNDMu(i,4)=-1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0-M_Eta1(i))

M_DNDXi(i,1)=1.0D0/8.0D0*(1.0D0-M_Eta1(i))*(1.0D0-M_Mu1(i))
M_DNDEta(i,1)=-1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0-M_Mu1(i))
M_DNDMu(i,1)=-1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0-M_Eta1(i))

M_DNDXi(i,2)=1.0D0/8.0D0*(1.0D0+M_Eta1(i))*(1.0D0-M_Mu1(i))
M_DNDEta(i,2)=1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0-M_Mu1(i))
M_DNDMu(i,2)=-1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0+M_Eta1(i))

M_DNDXi(i,3)=-1.0D0/8.0D0*(1.0D0+M_Eta1(i))*(1.0D0-M_Mu1(i))
M_DNDEta(i,3)=1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0-M_Mu1(i))
M_DNDMu(i,3)=-1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0+M_Eta1(i))

M_DNDXi(i,8)=-1.0D0/8.0D0*(1.0D0-M_Eta1(i))*(1.0D0+M_Mu1(i))
M_DNDEta(i,8)=-1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0+M_Mu1(i))
M_DNDMu(i,8)=1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0-M_Eta1(i))

M_DNDXi(i,5)=1.0D0/8.0D0*(1.0D0-M_Eta1(i))*(1.0D0+M_Mu1(i))
M_DNDEta(i,5)=-1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0+M_Mu1(i))
M_DNDMu(i,5)=1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0-M_Eta1(i))

M_DNDXi(i,6)=1.0D0/8.0D0*(1.0D0+M_Eta1(i))*(1.0D0+M_Mu1(i))
M_DNDEta(i,6)=1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0+M_Mu1(i))
M_DNDMu(i,6)=1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0+M_Eta1(i))

M_DNDXi(i,7)=-1.0D0/8.0D0*(1.0D0+M_Eta1(i))*(1.0D0+M_Mu1(i))
M_DNDEta(i,7)=1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0+M_Mu1(i))
M_DNDMu(i,7)=1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0+M_Eta1(i))

M_DNDXi(i,9)=-2.0D0*M_Xi1(i)
M_DNDEta(i,9)=0.0D0
M_DNDMu(i,9)=0.0D0

M_DNDXi(i,10)=0.0D0
M_DNDEta(i,10)=-2.0D0*M_Eta1(i)
M_DNDMu(i,10)=0.0D0

M_DNDXi(i,11)=0.0D0
M_DNDEta(i,11)=0.0D0
M_DNDMu(i,11)=-2.0D0*M_Mu1(i)

M_W(i)=1.0D0

ENDDO

DO i=1,M_NumberofGaussPoints
M_NN(i,1)=1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0-M_Eta1(i))*(1.0D0-M_Mu1(i))
M_NN(i,2)=1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0+M_Eta1(i))*(1.0D0-M_Mu1(i))
M_NN(i,3)=1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0+M_Eta1(i))*(1.0D0-M_Mu1(i))
M_NN(i,4)=1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0-M_Eta1(i))*(1.0D0-M_Mu1(i))
M_NN(i,5)=1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0-M_Eta1(i))*(1.0D0+M_Mu1(i))
M_NN(i,6)=1.0D0/8.0D0*(1.0D0+M_Xi1(i))*(1.0D0+M_Eta1(i))*(1.0D0+M_Mu1(i))
M_NN(i,7)=1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0+M_Eta1(i))*(1.0D0+M_Mu1(i))
M_NN(i,8)=1.0D0/8.0D0*(1.0D0-M_Xi1(i))*(1.0D0-M_Eta1(i))*(1.0D0+M_Mu1(i))
ENDDO

DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements

DO j=1,8
M_JJ(i,k,1,1)=M_JJ(i,k,1,1)+M_x(k,j)*M_DNDXi(i,j)
M_JJ(i,k,1,2)=M_JJ(i,k,1,2)+M_y(k,j)*M_DNDXi(i,j)
M_JJ(i,k,1,3)=M_JJ(i,k,1,3)+M_z(k,j)*M_DNDXi(i,j)

M_JJ(i,k,2,1)=M_JJ(i,k,2,1)+M_x(k,j)*M_DNDEta(i,j)
M_JJ(i,k,2,2)=M_JJ(i,k,2,2)+M_y(k,j)*M_DNDEta(i,j)
M_JJ(i,k,2,3)=M_JJ(i,k,2,3)+M_z(k,j)*M_DNDEta(i,j)

M_JJ(i,k,3,1)=M_JJ(i,k,3,1)+M_x(k,j)*M_DNDMu(i,j)
M_JJ(i,k,3,2)=M_JJ(i,k,3,2)+M_y(k,j)*M_DNDMu(i,j)
M_JJ(i,k,3,3)=M_JJ(i,k,3,3)+M_z(k,j)*M_DNDMu(i,j)

!M_JJ(i,k,1,1)=M_JJ(i,k,1,1)+M_x(k,j)*M_DNDXi(i,j)
!M_JJ(i,k,1,2)=M_JJ(i,k,1,2)+M_x(k,j)*M_DNDEta(i,j)
!M_JJ(i,k,1,3)=M_JJ(i,k,1,3)+M_x(k,j)*M_DNDMu(i,j)
!
!M_JJ(i,k,2,1)=M_JJ(i,k,2,1)+M_y(k,j)*M_DNDXi(i,j)
!M_JJ(i,k,2,2)=M_JJ(i,k,2,2)+M_y(k,j)*M_DNDEta(i,j)
!M_JJ(i,k,2,3)=M_JJ(i,k,2,3)+M_y(k,j)*M_DNDMu(i,j)
!
!M_JJ(i,k,3,1)=M_JJ(i,k,3,1)+M_z(k,j)*M_DNDXi(i,j)
!M_JJ(i,k,3,2)=M_JJ(i,k,3,2)+M_z(k,j)*M_DNDEta(i,j)
!M_JJ(i,k,3,3)=M_JJ(i,k,3,3)+M_z(k,j)*M_DNDMu(i,j)
ENDDO

ENDDO
ENDDO

DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements

!M_Delta_A=M_JJ(i,k,1,1)*(M_JJ(i,k,2,2)*M_JJ(i,k,3,3)-M_JJ(i,k,2,3)*M_JJ(i,k,3,2))+&
!M_JJ(i,k,1,2)*(M_JJ(i,k,3,2)*M_JJ(i,k,1,3)-M_JJ(i,k,1,2)*M_JJ(i,k,3,3))+&
!M_JJ(i,k,1,3)*(M_JJ(i,k,1,2)*M_JJ(i,k,2,3)-M_JJ(i,k,1,3)*M_JJ(i,k,2,2))
!
!M_Vol(i,k)=M_Delta_A
!
!M_InverseJJ(i,k,1,1)=(M_JJ(i,k,2,2)*M_JJ(i,k,3,3)-M_JJ(i,k,2,3)*M_JJ(i,k,3,2))/M_Delta_A
!M_InverseJJ(i,k,1,2)=(M_JJ(i,k,2,3)*M_JJ(i,k,3,1)-M_JJ(i,k,2,1)*M_JJ(i,k,3,3))/M_Delta_A
!M_InverseJJ(i,k,1,3)=(M_JJ(i,k,2,1)*M_JJ(i,k,2,2)-M_JJ(i,k,3,1)*M_JJ(i,k,2,2))/M_Delta_A
!M_InverseJJ(i,k,2,1)=(M_JJ(i,k,3,2)*M_JJ(i,k,1,3)-M_JJ(i,k,1,2)*M_JJ(i,k,3,3))/M_Delta_A
!M_InverseJJ(i,k,2,2)=(M_JJ(i,k,3,3)*M_JJ(i,k,1,1)-M_JJ(i,k,3,1)*M_JJ(i,k,1,3))/M_Delta_A
!M_InverseJJ(i,k,2,3)=(M_JJ(i,k,3,1)*M_JJ(i,k,1,2)-M_JJ(i,k,3,2)*M_JJ(i,k,1,1))/M_Delta_A
!M_InverseJJ(i,k,3,1)=(M_JJ(i,k,1,2)*M_JJ(i,k,2,3)-M_JJ(i,k,1,3)*M_JJ(i,k,2,2))/M_Delta_A
!M_InverseJJ(i,k,3,2)=(M_JJ(i,k,1,3)*M_JJ(i,k,2,1)-M_JJ(i,k,2,3)*M_JJ(i,k,1,1))/M_Delta_A
!M_InverseJJ(i,k,3,3)=(M_JJ(i,k,1,1)*M_JJ(i,k,2,2)-M_JJ(i,k,1,2)*M_JJ(i,k,2,1))/M_Delta_A

a11=M_JJ(i,k,1,1)
a12=M_JJ(i,k,1,2)
a13=M_JJ(i,k,1,3)
a21=M_JJ(i,k,2,1)
a22=M_JJ(i,k,2,2)
a23=M_JJ(i,k,2,3)
a31=M_JJ(i,k,3,1)
a32=M_JJ(i,k,3,2)
a33=M_JJ(i,k,3,3)

M_Vol(i,k)=-a13*a22*a31 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 +a11*a22*a33

M_InverseJJ(i,k,1,1)=(-a23*a32 + a22*a33)/M_Vol(i,k)
M_InverseJJ(i,k,1,2)=(a13*a32 - a12*a33)/M_Vol(i,k)
M_InverseJJ(i,k,1,3)=(-a13*a22 + a12*a23)/M_Vol(i,k)
M_InverseJJ(i,k,2,1)=(a23*a31 - a21*a33)/M_Vol(i,k)
M_InverseJJ(i,k,2,2)=(-a13*a31 + a11*a33)/M_Vol(i,k)
M_InverseJJ(i,k,2,3)=(a13*a21 - a11*a23)/M_Vol(i,k)
M_InverseJJ(i,k,3,1)=(-a22*a31 + a21*a32)/M_Vol(i,k)
M_InverseJJ(i,k,3,2)=(a12*a31 - a11*a32)/M_Vol(i,k)
M_InverseJJ(i,k,3,3)=(-a12*a21 + a11*a22)/M_Vol(i,k)

ENDDO
ENDDO

DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements
DO j=1,8+3

M_DNDX(i,k,j)=M_InverseJJ(i,k,1,1)*M_DNDXi(i,j)+M_InverseJJ(i,k,1,2)*M_DNDEta(i,j)+M_InverseJJ(i,k,1,3)*M_DNDMu(i,j)
M_DNDY(i,k,j)=M_InverseJJ(i,k,2,1)*M_DNDXi(i,j)+M_InverseJJ(i,k,2,2)*M_DNDEta(i,j)+M_InverseJJ(i,k,2,3)*M_DNDMu(i,j)
M_DNDZ(i,k,j)=M_InverseJJ(i,k,3,1)*M_DNDXi(i,j)+M_InverseJJ(i,k,3,2)*M_DNDEta(i,j)+M_InverseJJ(i,k,3,3)*M_DNDMu(i,j)

ENDDO
ENDDO
ENDDO

DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements
DO j=1,8+3

M_QX(i,k,j)=M_DNDX(i,k,j)
M_QY(i,k,j)=M_DNDY(i,k,j)
M_QZ(i,k,j)=M_DNDZ(i,k,j)

ENDDO
ENDDO
ENDDO

DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements
DO j=1,8

M_BB(i,k,1,j)=M_QX(i,k,j)
M_BB(i,k,5,j)=M_QZ(i,k,j)
M_BB(i,k,6,j)=M_QY(i,k,j)
M_BB(i,k,2,8+j)=M_QY(i,k,j)
M_BB(i,k,4,8+j)=M_QZ(i,k,j)
M_BB(i,k,6,8+j)=M_QX(i,k,j)
M_BB(i,k,3,16+j)=M_QZ(i,k,j)
M_BB(i,k,4,16+j)=M_QY(i,k,j)
M_BB(i,k,5,16+j)=M_QX(i,k,j)

ENDDO
ENDDO
ENDDO

DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements
DO i1=1,6
DO i2=1,24
M_BBT(i,k,i2,i1)=M_BB(i,k,i1,i2)
ENDDO
ENDDO
ENDDO
ENDDO


DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements
DO j=1,3

M_BBI(i,k,1,j)=M_QX(i,k,j+8)
M_BBI(i,k,5,j)=M_QZ(i,k,j+8)
M_BBI(i,k,6,j)=M_QY(i,k,j+8)
M_BBI(i,k,2,3+j)=M_QY(i,k,j+8)
M_BBI(i,k,4,3+j)=M_QZ(i,k,j+8)
M_BBI(i,k,6,3+j)=M_QX(i,k,j+8)
M_BBI(i,k,3,6+j)=M_QZ(i,k,j+8)
M_BBI(i,k,4,6+j)=M_QY(i,k,j+8)
M_BBI(i,k,5,6+j)=M_QX(i,k,j+8)

ENDDO
ENDDO
ENDDO

DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements
DO i1=1,6
DO i2=1,9
M_BBIT(i,k,i2,i1)=M_BBI(i,k,i1,i2)
ENDDO
ENDDO
ENDDO
ENDDO

M_BBIC=0.0D0
DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements
DO i1=1,6
DO j1=1,9

M_Vol_Tep=0.0D0

DO i2=1,M_NumberofGaussPoints
M_Vol_Tep=M_Vol_Tep+M_Vol(i2,k)
M_BBIC(i,k,i1,j1)=M_BBIC(i,k,i1,j1)+M_BBI(i2,k,i1,j1)*M_Vol(i2,k)*M_W(i2)
ENDDO

M_BBIC(i,k,i1,j1)=-M_BBIC(i,k,i1,j1)/M_Vol_Tep/8.0D0
M_BBIBar(i,k,i1,j1)=M_BBI(i,k,i1,j1)+M_BBIC(i,k,i1,j1)

ENDDO
ENDDO
ENDDO
ENDDO

DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements
DO i1=1,6
DO i2=1,9
M_BBIBarT(i,k,i2,i1)=M_BBIBar(i,k,i1,i2)
ENDDO
ENDDO
ENDDO
ENDDO



DO i=1,M_NumberofGaussPoints
DO k=1,M_NumberofElements

DO i1=1,24
DO i2=1,6
DO i3=1,6
DO i4=1,24
M_KCC(k,i1,i4)=M_KCC(k,i1,i4)+&
M_BBT(i,k,i1,i2)*M_Sub_CC(i,M_Connectivities(k,10),i2,i3)*M_BB(i,k,i3,i4)*M_Vol(i,k)*M_W(i)
ENDDO ! For i1
ENDDO ! For i2
ENDDO ! For i3 
ENDDO ! For i4

DO i1=1,24
DO i2=1,6
DO i3=1,6
DO i4=1,9
M_KCI(k,i1,i4)=M_KCI(k,i1,i4)+&
M_BBT(i,k,i1,i2)*M_Sub_CC(i,M_Connectivities(k,10),i2,i3)*M_BBIBar(i,k,i3,i4)*M_Vol(i,k)*M_W(i)
ENDDO ! For i1
ENDDO ! For i2
ENDDO ! For i3 
ENDDO ! For i4

DO i1=1,9
DO i2=1,6
DO i3=1,6
DO i4=1,24
M_KIC(k,i1,i4)=M_KIC(k,i1,i4)+&
M_BBIBarT(i,k,i1,i2)*M_Sub_CC(i,M_Connectivities(k,10),i2,i3)*M_BB(i,k,i3,i4)*M_Vol(i,k)*M_W(i)
ENDDO ! For i1
ENDDO ! For i2
ENDDO ! For i3 
ENDDO ! For i4

DO i1=1,9
DO i2=1,6
DO i3=1,6
DO i4=1,9
M_KII(k,i1,i4)=M_KII(k,i1,i4)+&
M_BBIBarT(i,k,i1,i2)*M_Sub_CC(i,M_Connectivities(k,10),i2,i3)*M_BBIBar(i,k,i3,i4)*M_Vol(i,k)*M_W(i)
ENDDO ! For i1
ENDDO ! For i2
ENDDO ! For i3 
ENDDO ! For i4

ENDDO
ENDDO

DO k=1,M_NumberofElements

CALL INV(k,9,M_KII(k,:,:),M_KII_Inverse(k,:,:))

!DO i1=1,9
!DO i2=1,9
!DO i3=1,9
!
!M_KII2(k,i1,i3)=M_KII2(k,i1,i3)+M_KII(k,i1,i2)*M_KII_Inverse(k,i2,i3)
!
!ENDDO
!ENDDO
!ENDDO
!
!DO i1=1,9
!DO i3=1,9
!IF(M_KII2(k,i1,i3) /= 0.0D0) THEN
!WRITE(*,30) k,i1,i3,M_KII2(k,i1,i3)
!ENDIF
!ENDDO
!ENDDO
!
!30 FORMAT(I,1X,I,1X,I,1X,E30.15,/)

DO i1=1,24
DO i2=1,9
DO i3=1,9
DO i4=1,24
M_KComponent(k,i1,i4)=M_KComponent(k,i1,i4)-M_KCI(k,i1,i2)*M_KII_Inverse(k,i2,i3)*M_KIC(k,i3,i4)
ENDDO
ENDDO
ENDDO
ENDDO

DO i1=1,24
DO i2=1,24
M_KComponent(k,i1,i2)=M_KComponent(k,i1,i2)+M_KCC(k,i1,i2)
ENDDO
ENDDO

ENDDO


RETURN
END SUBROUTINE M_Components

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_SixFaces(ii,Six_faces_Microscale,Six_faces_Microscale_Check,M_Nodes,M_PositiveXNodes,M_PositiveYNodes,&
M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes)
IMPLICIT NONE
INTEGER, INTENT(IN) :: ii
CHARACTER*90, INTENT(IN) :: Six_faces_Microscale
CHARACTER*90, INTENT(IN) :: Six_faces_Microscale_Check
INTEGER :: M_IN6,M_OUT6
INTEGER :: i,j

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

INTEGER, INTENT(OUT), DIMENSION(M_NumberofPositiveXNodes) :: M_PositiveXNodes
INTEGER, INTENT(OUT), DIMENSION(M_NumberofPositiveYNodes) :: M_PositiveYNodes
INTEGER, INTENT(OUT), DIMENSION(M_NumberofPositiveZNodes) :: M_PositiveZNodes
INTEGER, INTENT(OUT), DIMENSION(M_NumberofNegativeXNodes) :: M_NegativeXNodes
INTEGER, INTENT(OUT), DIMENSION(M_NumberofNegativeYNodes) :: M_NegativeYNodes
INTEGER, INTENT(OUT), DIMENSION(M_NumberofNegativeZNodes) :: M_NegativeZNodes


DOUBLE PRECISION, INTENT(INOUT), DIMENSION(M_NumberofNodes,3) :: M_Nodes
INTEGER :: M_Intermediate

M_IN6=10+100*(ii-1)
M_OUT6=50+100*(ii-1)
OPEN (M_IN6, file = Six_faces_Microscale,status ='unknown')

READ (M_IN6,*) (M_PositiveXNodes(i),i=1,M_NumberofPositiveXNodes)
READ (M_IN6,*) (M_PositiveYNodes(i),i=1,M_NumberofPositiveYNodes)
READ (M_IN6,*) (M_PositiveZNodes(i),i=1,M_NumberofPositiveZNodes)
READ (M_IN6,*) (M_NegativeXNodes(i),i=1,M_NumberofNegativeXNodes)
READ (M_IN6,*) (M_NegativeYNodes(i),i=1,M_NumberofNegativeYNodes)
READ (M_IN6,*) (M_NegativeZNodes(i),i=1,M_NumberofNegativeZNodes)

DO i=1,M_NumberofPositiveXNodes-1
DO j=i+1,M_NumberofPositiveXNodes

IF ( M_Nodes(M_PositiveXNodes(i),3) > M_Nodes(M_PositiveXNodes(j),3)+0.0D0 ) THEN
M_Intermediate=M_PositiveXNodes(i)
M_PositiveXNodes(i)=M_PositiveXNodes(j)
M_PositiveXNodes(j)=M_Intermediate

ELSEIF ( ABS( M_Nodes(M_PositiveXNodes(i),3) - M_Nodes(M_PositiveXNodes(j),3) ) ==0.0D0 ) THEN

IF ( M_Nodes(M_PositiveXNodes(i),4) > M_Nodes(M_PositiveXNodes(j),4)+0.0D0 ) THEN
M_Intermediate=M_PositiveXNodes(i)
M_PositiveXNodes(i)=M_PositiveXNodes(j)
M_PositiveXNodes(j)=M_Intermediate
ENDIF

ENDIF

ENDDO
ENDDO


DO i=1,M_NumberofPositiveYNodes-1
DO j=i+1,M_NumberofPositiveYNodes

IF ( M_Nodes(M_PositiveYNodes(i),2) > M_Nodes(M_PositiveYNodes(j),2)+0.0D0 ) THEN
M_Intermediate=M_PositiveYNodes(i)
M_PositiveYNodes(i)=M_PositiveYNodes(j)
M_PositiveYNodes(j)=M_Intermediate

ELSEIF ( ABS( M_Nodes(M_PositiveYNodes(i),2) - M_Nodes(M_PositiveYNodes(j),2) ) == 0.0D0 ) THEN

IF ( M_Nodes(M_PositiveYNodes(i),4) > M_Nodes(M_PositiveYNodes(j),4)+0.0D0 ) THEN
M_Intermediate=M_PositiveYNodes(i)
M_PositiveYNodes(i)=M_PositiveYNodes(j)
M_PositiveYNodes(j)=M_Intermediate
ENDIF

ENDIF

ENDDO
ENDDO


DO i=1,M_NumberofPositiveZNodes-1
DO j=i+1,M_NumberofPositiveZNodes

IF ( M_Nodes(M_PositiveZNodes(i),2) > M_Nodes(M_PositiveZNodes(j),2) +0.0D0 ) THEN
M_Intermediate=M_PositiveZNodes(i)
M_PositiveZNodes(i)=M_PositiveZNodes(j)
M_PositiveZNodes(j)=M_Intermediate

ELSEIF ( ABS( M_Nodes(M_PositiveZNodes(i),2) - M_Nodes(M_PositiveZNodes(j),2) ) ==0.0D0 ) THEN

IF ( M_Nodes(M_PositiveZNodes(i),3) > M_Nodes(M_PositiveZNodes(j),3)+0.0D0 ) THEN
M_Intermediate=M_PositiveZNodes(i)
M_PositiveZNodes(i)=M_PositiveZNodes(j)
M_PositiveZNodes(j)=M_Intermediate
ENDIF

ENDIF

ENDDO
ENDDO




DO i=1,M_NumberofNegativeXNodes-1
DO j=i+1,M_NumberofNegativeXNodes

IF ( M_Nodes(M_NegativeXNodes(i),3) > M_Nodes(M_NegativeXNodes(j),3)+0.0D0 ) THEN
M_Intermediate=M_NegativeXNodes(i)
M_NegativeXNodes(i)=M_NegativeXNodes(j)
M_NegativeXNodes(j)=M_Intermediate

ELSEIF ( ABS( M_Nodes(M_NegativeXNodes(i),3) - M_Nodes(M_NegativeXNodes(j),3) ) ==0.0D0 ) THEN

IF ( M_Nodes(M_NegativeXNodes(i),4) > M_Nodes(M_NegativeXNodes(j),4)+0.0D0 ) THEN
M_Intermediate=M_NegativeXNodes(i)
M_NegativeXNodes(i)=M_NegativeXNodes(j)
M_NegativeXNodes(j)=M_Intermediate
ENDIF

ENDIF

ENDDO
ENDDO


DO i=1,M_NumberofNegativeYNodes-1
DO j=i+1,M_NumberofNegativeYNodes

IF ( M_Nodes(M_NegativeYNodes(i),2) > M_Nodes(M_NegativeYNodes(j),2)+0.0D0 ) THEN
M_Intermediate=M_NegativeYNodes(i)
M_NegativeYNodes(i)=M_NegativeYNodes(j)
M_NegativeYNodes(j)=M_Intermediate

ELSEIF ( ABS( M_Nodes(M_NegativeYNodes(i),2) - M_Nodes(M_NegativeYNodes(j),2) ) ==0.0D0 ) THEN

IF ( M_Nodes(M_NegativeYNodes(i),4) > M_Nodes(M_NegativeYNodes(j),4)+0.0D0 ) THEN
M_Intermediate=M_NegativeYNodes(i)
M_NegativeYNodes(i)=M_NegativeYNodes(j)
M_NegativeYNodes(j)=M_Intermediate
ENDIF

ENDIF

ENDDO
ENDDO


DO i=1,M_NumberofNegativeZNodes-1
DO j=i+1,M_NumberofNegativeZNodes

IF ( M_Nodes(M_NegativeZNodes(i),2) > M_Nodes(M_NegativeZNodes(j),2)+0.0D0 ) THEN
M_Intermediate=M_NegativeZNodes(i)
M_NegativeZNodes(i)=M_NegativeZNodes(j)
M_NegativeZNodes(j)=M_Intermediate

ELSEIF ( ABS( M_Nodes(M_NegativeZNodes(i),2) - M_Nodes(M_NegativeZNodes(j),2) ) ==0.0D0 ) THEN

IF ( M_Nodes(M_NegativeZNodes(i),3) > M_Nodes(M_NegativeZNodes(j),3)+0.0D0 ) THEN
M_Intermediate=M_NegativeZNodes(i)
M_NegativeZNodes(i)=M_NegativeZNodes(j)
M_NegativeZNodes(j)=M_Intermediate
ENDIF

ENDIF

ENDDO
ENDDO

CLOSE(M_IN6)

OPEN(M_OUT6,FILE= Six_faces_Microscale_Check, STATUS="UNKNOWN")

WRITE(M_OUT6,*) "===================X===================="

DO i=1,M_NumberofPositiveXNodes

IF ( M_Nodes(M_PositiveXNodes(i),3)-M_Nodes(M_NegativeXNodes(i),3) /=0.0D0 .OR. &
M_Nodes(M_PositiveXNodes(i),4)-M_Nodes(M_NegativeXNodes(i),4) /=0.0D0) THEN

WRITE(M_OUT6,21) M_PositiveXNodes(i),M_Nodes(M_PositiveXNodes(i),3),M_Nodes(M_PositiveXNodes(i),4),&
M_NegativeXNodes(i),M_Nodes(M_NegativeXNodes(i),3),M_Nodes(M_NegativeXNodes(i),4),&
M_Nodes(M_PositiveXNodes(i),3)-M_Nodes(M_NegativeXNodes(i),3),&
M_Nodes(M_PositiveXNodes(i),4)-M_Nodes(M_NegativeXNodes(i),4)
ENDIF

ENDDO

WRITE(M_OUT6,*) "===================Y===================="

DO i=1,M_NumberofPositiveYNodes

IF ( M_Nodes(M_PositiveYNodes(i),2)-M_Nodes(M_NegativeYNodes(i),2) /=0.0D0 .OR. &
M_Nodes(M_PositiveYNodes(i),4)-M_Nodes(M_NegativeYNodes(i),4) /=0.0D0) THEN

WRITE(M_OUT6,21) M_PositiveYNodes(i),M_Nodes(M_PositiveYNodes(i),2),M_Nodes(M_PositiveYNodes(i),4),&
M_NegativeYNodes(i),M_Nodes(M_NegativeYNodes(i),2),M_Nodes(M_NegativeYNodes(i),4),&
M_Nodes(M_PositiveYNodes(i),2)-M_Nodes(M_NegativeYNodes(i),2),&
M_Nodes(M_PositiveYNodes(i),4)-M_Nodes(M_NegativeYNodes(i),4)
ENDIF

ENDDO

WRITE(M_OUT6,*) "===================Z===================="

DO i=1,M_NumberofPositiveZNodes

IF ( M_Nodes(M_PositiveZNodes(i),2)-M_Nodes(M_NegativeZNodes(i),2) /=0.0D0 .OR. &
M_Nodes(M_PositiveZNodes(i),3)-M_Nodes(M_NegativeZNodes(i),3) /=0.0D0) THEN

WRITE(M_OUT6,21) M_PositiveZNodes(i),M_Nodes(M_PositiveZNodes(i),2),M_Nodes(M_PositiveZNodes(i),3),&
M_NegativeZNodes(i),M_Nodes(M_NegativeZNodes(i),2),M_Nodes(M_NegativeZNodes(i),3),&
M_Nodes(M_PositiveZNodes(i),2)-M_Nodes(M_NegativeZNodes(i),2),&
M_Nodes(M_PositiveZNodes(i),3)-M_Nodes(M_NegativeZNodes(i),3)
ENDIF

ENDDO

21 FORMAT(I10,1X,E15.7,1X,E15.7,1X,I10,1X,E15.7,1X,E15.8,1X,E15.7,1X,E15.8,/)

CLOSE(M_OUT6)

DO i=1,M_NumberofPositiveXNodes

M_Nodes(M_PositiveXNodes(i),3)=M_Nodes(M_NegativeXNodes(i),3)
M_Nodes(M_PositiveXNodes(i),4)=M_Nodes(M_NegativeXNodes(i),4)

ENDDO

DO i=1,M_NumberofPositiveYNodes

M_Nodes(M_PositiveYNodes(i),2)=M_Nodes(M_NegativeYNodes(i),2)
M_Nodes(M_PositiveYNodes(i),4)=M_Nodes(M_NegativeYNodes(i),4)

ENDDO

DO i=1,M_NumberofPositiveZNodes

M_Nodes(M_PositiveZNodes(i),2)=M_Nodes(M_NegativeZNodes(i),2)
M_Nodes(M_PositiveZNodes(i),3)=M_Nodes(M_NegativeZNodes(i),3)

ENDDO

RETURN

END SUBROUTINE M_SixFaces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_ReadMechanicalProperties (ii,Mechanical_material_properties_Microscale,M_CC)
IMPLICIT NONE
INTEGER, INTENT(IN) :: ii
CHARACTER*90, INTENT(IN) :: Mechanical_material_properties_Microscale
INTEGER :: M_IN5
INTEGER :: i,j

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

DOUBLE PRECISION, DIMENSION(M_NumberofMaterialModes,10) :: M_MechanicalMaterialProperties

DOUBLE PRECISION :: M_EE, M_nu, M_kappaa12,M_muu12,M_EE33,M_muu23,M_nuu23,M_kappatz,M_mutz,M_Err,M_murt,M_nurt,M_E11,M_E22,&
M_E33,M_nu12,M_nu13,M_nu23,M_G23,M_G13,M_G12,M_nu21,M_nu31,M_nu32

DOUBLE PRECISION :: M_C11,M_C22,M_C12,M_C66

DOUBLE PRECISION, INTENT(OUT), DIMENSION(M_NumberofMaterialModes,6,6) :: M_CC

M_MechanicalMaterialProperties=0.0D0
M_CC=0.0D0

M_IN5=9+100*(ii-1)
OPEN(M_IN5,file = Mechanical_material_properties_Microscale,status ='unknown')

DO i=1,M_NumberofMaterialModes

READ(M_IN5,*) M_MechanicalMaterialProperties(i,1)

IF (M_MechanicalMaterialProperties(i,1) == 1) THEN
BACKSPACE M_IN5
READ(M_IN5,*) M_MechanicalMaterialProperties(i,1),M_MechanicalMaterialProperties(i,2),M_MechanicalMaterialProperties(i,3)

ELSEIF (M_MechanicalMaterialProperties(i,1) == 2 .OR. M_MechanicalMaterialProperties(i,1) == 3) THEN
BACKSPACE M_IN5
READ(M_IN5,*) M_MechanicalMaterialProperties(i,1),M_MechanicalMaterialProperties(i,2),M_MechanicalMaterialProperties(i,3),&
M_MechanicalMaterialProperties(i,4),M_MechanicalMaterialProperties(i,5),M_MechanicalMaterialProperties(i,6)

ELSEIF (M_MechanicalMaterialProperties(i,1) == 4 ) THEN
BACKSPACE M_IN5
READ(M_IN5,*) M_MechanicalMaterialProperties(i,1),M_MechanicalMaterialProperties(i,2),M_MechanicalMaterialProperties(i,3),&
M_MechanicalMaterialProperties(i,4),M_MechanicalMaterialProperties(i,5),M_MechanicalMaterialProperties(i,6),&
M_MechanicalMaterialProperties(i,7),M_MechanicalMaterialProperties(i,8),M_MechanicalMaterialProperties(i,9),&
M_MechanicalMaterialProperties(i,10)

ELSEIF (M_MechanicalMaterialProperties(i,1) == 5 ) THEN
BACKSPACE M_IN5
READ(M_IN5,*) M_MechanicalMaterialProperties(i,1),M_MechanicalMaterialProperties(i,2),M_MechanicalMaterialProperties(i,3),&
M_MechanicalMaterialProperties(i,4),M_MechanicalMaterialProperties(i,5)

ENDIF

ENDDO


DO i=1,M_NumberofMaterialModes

IF (M_MechanicalMaterialProperties(i,1) == 1) THEN

M_EE=M_MechanicalMaterialProperties(i,2)
M_nu=M_MechanicalMaterialProperties(i,3)

M_CC(i,1,1)=M_EE/(1.0D0+M_nu)*(1.0D0-M_nu)/(1.0D0-2.0D0*M_nu)
M_CC(i,1,2)=M_EE/(1.0D0+M_nu)*M_nu/(1.0D0-2.0D0*M_nu)
M_CC(i,1,3)=M_EE/(1.0D0+M_nu)*M_nu/(1.0D0-2.0D0*M_nu)
M_CC(i,2,1)=M_CC(i,1,2)
M_CC(i,2,2)=M_CC(i,1,1)
M_CC(i,2,3)=M_CC(i,1,2)
M_CC(i,3,1)=M_CC(i,1,3)
M_CC(i,3,2)=M_CC(i,1,2)
M_CC(i,3,3)=M_CC(i,1,1)
M_CC(i,4,4)=M_EE/2.0D0/(1.0D0+M_nu)
M_CC(i,5,5)=M_EE/2.0D0/(1.0D0+M_nu)
M_CC(i,6,6)=M_EE/2.0D0/(1.0D0+M_nu)


ELSEIF (M_MechanicalMaterialProperties(i,1) == 2) THEN

M_kappaa12=M_MechanicalMaterialProperties(i,2)
M_muu12=M_MechanicalMaterialProperties(i,3)
M_EE33=M_MechanicalMaterialProperties(i,4)
M_muu23=M_MechanicalMaterialProperties(i,5)
M_nuu23=M_MechanicalMaterialProperties(i,6)

M_CC(i,1,1)=M_muu12+M_kappaa12
M_CC(i,1,2)=-M_muu12+M_kappaa12
M_CC(i,1,3)=2.0D0*M_kappaa12*M_nuu23
M_CC(i,2,1)=M_CC(i,1,2)
M_CC(i,2,2)=M_CC(i,1,1)
M_CC(i,2,3)=M_CC(i,1,3)
M_CC(i,3,1)=M_CC(i,1,3)
M_CC(i,3,2)=M_CC(i,2,3)
M_CC(i,3,3)=M_EE33+4.0D0*M_nuu23*M_nuu23*M_kappaa12
M_CC(i,4,4)=M_muu23
M_CC(i,5,5)=M_muu23
M_CC(i,6,6)=M_muu12

ELSEIF (M_MechanicalMaterialProperties(i,1) == 3) THEN

M_kappatz=M_MechanicalMaterialProperties(i,2)
M_mutz=M_MechanicalMaterialProperties(i,3)
M_Err=M_MechanicalMaterialProperties(i,4)
M_murt=M_MechanicalMaterialProperties(i,5)
M_nurt=M_MechanicalMaterialProperties(i,6)

M_CC(i,1,1)=M_Err+4.0D0*M_nurt**2*M_kappatz
M_CC(i,1,2)=2.0D0*M_kappatz*M_nurt
M_CC(i,1,3)=2.0D0*M_kappatz*M_nurt
M_CC(i,2,1)=M_CC(i,1,2)
M_CC(i,2,2)=M_mutz+M_kappatz
M_CC(i,2,3)=-M_mutz+M_kappatz
M_CC(i,3,1)=M_CC(i,1,3)
M_CC(i,3,2)=M_CC(i,2,3)
M_CC(i,3,3)=M_mutz+M_kappatz
M_CC(i,4,4)=M_mutz
M_CC(i,5,5)=M_murt
M_CC(i,6,6)=M_murt

ELSEIF (M_MechanicalMaterialProperties(i,1) == 4 ) THEN

M_E11=M_MechanicalMaterialProperties(i,2)
M_E22=M_MechanicalMaterialProperties(i,3)
M_E33=M_MechanicalMaterialProperties(i,4)
M_nu12=M_MechanicalMaterialProperties(i,5)
M_nu13=M_MechanicalMaterialProperties(i,6)
M_nu23=M_MechanicalMaterialProperties(i,7)
M_G23=M_MechanicalMaterialProperties(i,8)
M_G13=M_MechanicalMaterialProperties(i,9)
M_G12=M_MechanicalMaterialProperties(i,10)

M_nu21=M_nu12/M_E11*M_E22
M_nu31=M_nu13/M_E11*M_E33
M_nu32=M_nu23/M_E22*M_E33

M_CC(i,1,1)=(M_E11*(-1.0D0 + M_nu23*M_nu32))/(-1.0D0 + M_nu12*(M_nu21 + M_nu23*M_nu31) + M_nu23*M_nu32 + M_nu13*(M_nu31 + M_nu21*M_nu32))
M_CC(i,1,2)=-((M_E22*(M_nu12 + M_nu13*M_nu32))/(-1.0D0 + M_nu12*(M_nu21 + M_nu23*M_nu31) + M_nu23*M_nu32 + M_nu13*(M_nu31 + M_nu21*M_nu32)))
M_CC(i,1,3)=-((M_E33*(M_nu13 + M_nu12*M_nu23))/(-1.0D0 + M_nu12*(M_nu21 + M_nu23*M_nu31) + M_nu23*M_nu32 + M_nu13*(M_nu31 + M_nu21*M_nu32)))
M_CC(i,2,2)=(M_E22*(-1.0D0 + M_nu13*M_nu31))/(-1.0D0 + M_nu12*(M_nu21 + M_nu23*M_nu31) + M_nu23*M_nu32 + M_nu13*(M_nu31 + M_nu21*M_nu32))
M_CC(i,2,3)=-((M_E33*(M_nu13*M_nu21 + M_nu23))/(-1.0D0 + M_nu12*(M_nu21 + M_nu23*M_nu31) + M_nu23*M_nu32 + M_nu13*(M_nu31 + M_nu21*M_nu32)))
M_CC(i,3,3)=(M_E33*(-1.0D0 + M_nu12*M_nu21))/(-1.0D0 + M_nu12*(M_nu21 + M_nu23*M_nu31) + M_nu23*M_nu32 + M_nu13*(M_nu31 + M_nu21*M_nu32))
M_CC(i,2,1)=M_CC(i,1,2)
M_CC(i,3,1)=M_CC(i,1,3)
M_CC(i,3,2)=M_CC(i,2,3)
M_CC(i,4,4)=M_G23
M_CC(i,5,5)=M_G13
M_CC(i,6,6)=M_G12

ELSEIF (M_MechanicalMaterialProperties(i,1) == 5 ) THEN

M_C11=M_MechanicalMaterialProperties(i,2)
M_C22=M_MechanicalMaterialProperties(i,3)
M_C12=M_MechanicalMaterialProperties(i,4)
M_C66=M_MechanicalMaterialProperties(i,5)

M_CC(i,1,1)=M_C11
M_CC(i,1,2)=M_C12
M_CC(i,1,3)=0.0D0
M_CC(i,2,2)=M_C22
M_CC(i,2,3)=0.0D0
M_CC(i,3,3)=0.0D0
M_CC(i,2,1)=M_C12
M_CC(i,3,1)=0.0D0
M_CC(i,3,2)=0.0D0
M_CC(i,4,4)=0.0D0
M_CC(i,5,5)=0.0D0
M_CC(i,6,6)=M_C66

ENDIF

ENDDO

CLOSE(M_IN5)


RETURN
END SUBROUTINE M_ReadMechanicalProperties

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_ReadBoundaryConditions(ii,Boundaries_Microscale,M_BoundaryIndex,M_Boundaries)
IMPLICIT NONE
INTEGER :: M_IN4
INTEGER :: i,j
INTEGER, INTENT(IN) :: ii
CHARACTER*90, INTENT(IN) :: Boundaries_Microscale

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

DOUBLE PRECISION, INTENT(OUT), DIMENSION(M_NumberofBoundaryNodes,TotalTimeSteps+2) :: M_Boundaries
INTEGER, INTENT(OUT), DIMENSION(M_NumberofBoundaryNodes,2) :: M_BoundaryIndex

M_IN4=8+100*(ii-1)
OPEN(M_IN4,file = Boundaries_Microscale,status ='unknown')
!READ (M_IN4,*) ((M_Boundaries(i,j),j=1,TotalTimeSteps+2),i=1,M_NumberofBoundaryNodes)

DO i=1,M_NumberofBoundaryNodes
READ (M_IN4,*) (M_Boundaries(i,j),j=1,TotalTimeSteps+2)
ENDDO

DO i=1,M_NumberofBoundaryNodes
M_BoundaryIndex(i,1)=M_Boundaries(i,1)
M_BoundaryIndex(i,2)=M_Boundaries(i,2)
ENDDO
CLOSE(M_IN4)
RETURN
END SUBROUTINE M_ReadBoundaryConditions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_ReadConnectivities (ii,Connectivities_Microscale,M_Connectivities)
IMPLICIT NONE
INTEGER :: M_IN3
INTEGER :: ii
INTEGER :: i,j
CHARACTER*90, INTENT(IN) :: Connectivities_Microscale

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

INTEGER, INTENT(OUT), DIMENSION (M_NumberofElements,10) :: M_Connectivities

M_IN3=7+100*(ii-1)

OPEN(M_IN3,file = Connectivities_Microscale, status ='unknown')

READ (M_IN3,*) ((M_Connectivities(i,j),j=1,10),i=1,M_NumberofElements)

CLOSE(M_IN3)


RETURN
ENDSUBROUTINE M_ReadConnectivities


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_ReadNodes (ii,The_Nodes_Microscale,M_NodeIndex,M_Nodes)
IMPLICIT NONE
INTEGER, INTENT(IN) :: ii
INTEGER :: M_IN2
CHARACTER*90, INTENT(IN) :: The_Nodes_Microscale
INTEGER :: i,j

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

DOUBLE PRECISION, INTENT(OUT), DIMENSION (M_NumberofNodes,4) :: M_Nodes
INTEGER, INTENT(OUT), DIMENSION(M_NumberofNodes) :: M_NodeIndex

M_IN2=6+100*(ii-1)

OPEN(M_IN2,file = The_Nodes_Microscale, status ='unknown')
READ (M_IN2,*) ((M_Nodes(i,j),j=1,4),i=1,M_NumberofNodes)
DO i=1,M_NumberofNodes
M_NodeIndex(i)=M_Nodes(i,1)
ENDDO
CLOSE(M_IN2)

RETURN
ENDSUBROUTINE M_ReadNodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE M_ReadControlFlags (ii,Control_Flags_Microscale) 
IMPLICIT NONE
INTEGER :: i
INTEGER, INTENT(IN) :: ii
INTEGER :: M_IN1
CHARACTER*90, INTENT(IN) :: Control_Flags_Microscale
DOUBLE PRECISION, DIMENSION (13) :: M_ControlFlags

!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec

M_IN1=5+100*(ii-1)

OPEN(M_IN1,file = Control_Flags_Microscale,status ='unknown')

DO i=1,13
READ (M_IN1,*) M_ControlFlags(i)
ENDDO

M_NumberofNodes=M_ControlFlags(1)
M_NumberofElements=M_ControlFlags(2)
M_NumberofBoundaryNodes=M_ControlFlags(3)
M_NumberofMaterialModes=M_ControlFlags(4)
M_NumberofGaussPoints=M_ControlFlags(5)
M_Mag=M_ControlFlags(6)
M_NumberofPositiveXNodes=M_ControlFlags(7)
M_NumberofPositiveYNodes=M_ControlFlags(8)
M_NumberofPositiveZNodes=M_ControlFlags(9)
M_NumberofNegativeXNodes=M_ControlFlags(10)
M_NumberofNegativeYNodes=M_ControlFlags(11)
M_NumberofNegativeZNodes=M_ControlFlags(12)
M_MaxDimenofInverseConnec=M_ControlFlags(13)
ClOSE(M_IN1)

RETURN
END SUBROUTINE M_ReadControlFlags


SUBROUTINE M_MicroFolderNames(ii,Boundaries_Microscale,Connectivities_Microscale,Control_Flags_Microscale,Input_Strains_Microscale,&
Mechanical_material_properties_Microscale,Six_faces_Microscale,The_Nodes_Microscale,Six_faces_Microscale_Check)
IMPLICIT NONE

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER, INTENT(IN) :: ii
CHARACTER*90, INTENT(OUT) :: Boundaries_Microscale
CHARACTER*90, INTENT(OUT) :: Connectivities_Microscale
CHARACTER*90, INTENT(OUT) :: Control_Flags_Microscale
CHARACTER*90, INTENT(OUT) :: Input_Strains_Microscale
CHARACTER*90, INTENT(OUT) :: Mechanical_material_properties_Microscale
CHARACTER*90, INTENT(OUT) :: Six_faces_Microscale
CHARACTER*90, INTENT(OUT) :: The_Nodes_Microscale
CHARACTER*90, INTENT(OUT) :: Six_faces_Microscale_Check

IF (ii==1) THEN
Boundaries_Microscale="3D_Output_Micro_F/3D_Input_Microscale_1/Boundaries_Microscale.txt"
Connectivities_Microscale="3D_Output_Micro_F/3D_Input_Microscale_1/Connectivities_Microscale.txt"
Control_Flags_Microscale="3D_Output_Micro_F/3D_Input_Microscale_1/Control_Flags_Microscale.txt"
Input_Strains_Microscale="3D_Output_Micro_F/3D_Input_Microscale_1/Input_Strains_Microscale.txt"
Mechanical_material_properties_Microscale="3D_Output_Micro_F/3D_Input_Microscale_1/Mechanical_Microscale.txt"
Six_faces_Microscale="3D_Output_Micro_F/3D_Input_Microscale_1/Six_faces_Microscale.txt"
The_Nodes_Microscale="3D_Output_Micro_F/3D_Input_Microscale_1/The_Nodes_Microscale.txt"
Six_faces_Microscale_Check="3D_Output_Micro_F/3D_Output_Microscale_1/CheckSixFaces_Microscale.txt"
ELSEIF (ii==2) THEN
Boundaries_Microscale="3D_Output_Micro_F/3D_Input_Microscale_2/Boundaries_Microscale.txt"
Connectivities_Microscale="3D_Output_Micro_F/3D_Input_Microscale_2/Connectivities_Microscale.txt"
Control_Flags_Microscale="3D_Output_Micro_F/3D_Input_Microscale_2/Control_Flags_Microscale.txt"
Input_Strains_Microscale="3D_Output_Micro_F/3D_Input_Microscale_2/Input_Strains_Microscale.txt"
Mechanical_material_properties_Microscale="3D_Output_Micro_F/3D_Input_Microscale_2/Mechanical_Microscale.txt"
Six_faces_Microscale="3D_Output_Micro_F/3D_Input_Microscale_2/Six_faces_Microscale.txt"
The_Nodes_Microscale="3D_Output_Micro_F/3D_Input_Microscale_2/The_Nodes_Microscale.txt"
Six_faces_Microscale_Check="3D_Output_Micro_F/3D_Output_Microscale_2/CheckSixFaces_Microscale.txt"
ELSEIF (ii==3) THEN
Boundaries_Microscale="3D_Output_Micro_F/3D_Input_Microscale_3/Boundaries_Microscale.txt"
Connectivities_Microscale="3D_Output_Micro_F/3D_Input_Microscale_3/Connectivities_Microscale.txt"
Control_Flags_Microscale="3D_Output_Micro_F/3D_Input_Microscale_3/Control_Flags_Microscale.txt"
Input_Strains_Microscale="3D_Output_Micro_F/3D_Input_Microscale_3/Input_Strains_Microscale.txt"
Mechanical_material_properties_Microscale="3D_Output_Micro_F/3D_Input_Microscale_3/Mechanical_Microscale.txt"
Six_faces_Microscale="3D_Output_Micro_F/3D_Input_Microscale_3/Six_faces_Microscale.txt"
The_Nodes_Microscale="3D_Output_Micro_F/3D_Input_Microscale_3/The_Nodes_Microscale.txt"
Six_faces_Microscale_Check="3D_Output_Micro_F/3D_Output_Microscale_3/CheckSixFaces_Microscale.txt"
ELSEIF (ii==4) THEN
Boundaries_Microscale="3D_Output_Micro_F/3D_Input_Microscale_4/Boundaries_Microscale.txt"
Connectivities_Microscale="3D_Output_Micro_F/3D_Input_Microscale_4/Connectivities_Microscale.txt"
Control_Flags_Microscale="3D_Output_Micro_F/3D_Input_Microscale_4/Control_Flags_Microscale.txt"
Input_Strains_Microscale="3D_Output_Micro_F/3D_Input_Microscale_4/Input_Strains_Microscale.txt"
Mechanical_material_properties_Microscale="3D_Output_Micro_F/3D_Input_Microscale_4/Mechanical_Microscale.txt"
Six_faces_Microscale="3D_Output_Micro_F/3D_Input_Microscale_4/Six_faces_Microscale.txt"
The_Nodes_Microscale="3D_Output_Micro_F/3D_Input_Microscale_4/The_Nodes_Microscale.txt"
Six_faces_Microscale_Check="3D_Output_Micro_F/3D_Output_Microscale_4/CheckSixFaces_Microscale.txt"
ELSEIF (ii==5) THEN
Boundaries_Microscale="3D_Output_Micro_F/3D_Input_Microscale_5/Boundaries_Microscale.txt"
Connectivities_Microscale="3D_Output_Micro_F/3D_Input_Microscale_5/Connectivities_Microscale.txt"
Control_Flags_Microscale="3D_Output_Micro_F/3D_Input_Microscale_5/Control_Flags_Microscale.txt"
Input_Strains_Microscale="3D_Output_Micro_F/3D_Input_Microscale_5/Input_Strains_Microscale.txt"
Mechanical_material_properties_Microscale="3D_Output_Micro_F/3D_Input_Microscale_5/Mechanical_Microscale.txt"
Six_faces_Microscale="3D_Output_Micro_F/3D_Input_Microscale_5/Six_faces_Microscale.txt"
The_Nodes_Microscale="3D_Output_Micro_F/3D_Input_Microscale_5/The_Nodes_Microscale.txt"
Six_faces_Microscale_Check="3D_Output_Micro_F/3D_Output_Microscale_5/CheckSixFaces_Microscale.txt"
ENDIF

RETURN
END SUBROUTINE M_MicroFolderNames

SUBROUTINE Macro_Each_RVE_Type(k,ii,ElementType,ElementType2,ElementTypeIndex,Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,&
Sub_StrainsXY1_R,Sub_StrainsXY2_R,Sub_StrainsXZ1_R,Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,Sub_StrainsYZ2_R,M_Highest_GRCD,ME_Highest_GRCD,&
E_sub_sigmasigma_in,RotationMatrix2,BeginofModeIndex,EndofModeIndex,Connectivities)
IMPLICIT NONE
!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown
!The macroscale electrostatic control flags
INTEGER :: E_NumberofBoundaryNodes
INTEGER :: E_DisplacementEquationFlag
INTEGER :: E_NumberofGaussPoints
INTEGER :: E_MaxDimenofInverseConnec
INTEGER :: E_NumberofMaterialModes
DOUBLE PRECISION :: E_Mag
COMMON/E_Control_Flags/ E_Mag,E_NumberofBoundaryNodes,E_DisplacementEquationFlag,E_NumberofGaussPoints,E_MaxDimenofInverseConnec,&
E_NumberofMaterialModes
!The control flags in the microscale
INTEGER :: M_NumberofElements
INTEGER :: M_NumberofNodes
INTEGER :: M_NumberofBoundaryNodes
INTEGER :: M_NumberofMaterialModes
INTEGER :: M_NumberofGaussPoints
INTEGER :: M_NumberofPositiveXNodes
INTEGER :: M_NumberofPositiveYNodes
INTEGER :: M_NumberofPositiveZNodes
INTEGER :: M_NumberofNegativeXNodes
INTEGER :: M_NumberofNegativeYNodes
INTEGER :: M_NumberofNegativeZNodes
INTEGER :: M_MaxDimenofInverseConnec
DOUBLE PRECISION :: M_Mag
COMMON/M_Control_Flags/M_Mag,M_NumberofElements,M_NumberofNodes,M_NumberofBoundaryNodes,M_NumberofMaterialModes,M_NumberofGaussPoints,&
M_NumberofPositiveXNodes,M_NumberofPositiveYNodes,M_NumberofPositiveZNodes,M_NumberofNegativeXNodes,M_NumberofNegativeYNodes,&
M_NumberofNegativeZNodes,M_MaxDimenofInverseConnec
!The electrical control flags in the microscale
INTEGER :: ME_NumberofBoundaryNodes
INTEGER :: ME_NumberofMaterialModes
INTEGER :: ME_NumberofGaussPoints
INTEGER :: ME_MaxDimenofInverseConnec
DOUBLE PRECISION :: ME_Mag
COMMON/ME_Control_Flags/ME_Mag,ME_NumberofBoundaryNodes,ME_NumberofMaterialModes,ME_NumberofGaussPoints,ME_MaxDimenofInverseConnec
!The electrical tunneling control flags in the microscale
INTEGER :: ME_Tunneling
DOUBLE PRECISION :: ME_Vf
DOUBLE PRECISION :: ME_RofCNT
DOUBLE PRECISION :: ME_Lambda
DOUBLE PRECISION :: ME_rho
INTEGER :: ME_BeginningPolymer
INTEGER :: ME_EndingPolymer
INTEGER :: ME_NumberofCNTs
INTEGER :: ME_NumberofCNTs2
INTEGER :: ME_MaxNumberofElementsPerCNT
COMMON/Micro_Electrical_More_Inputs/ME_Tunneling,ME_BeginningPolymer,ME_EndingPolymer,ME_NumberofCNTs,ME_NumberofCNTs2,&
ME_MaxNumberofElementsPerCNT,ME_Vf,ME_RofCNT,ME_Lambda,ME_rho

INTEGER :: ii1,ii2
INTEGER, INTENT(IN) :: k,ii
INTEGER, INTENT(IN), DIMENSION(NumberofElements,9) :: ElementType
INTEGER, INTENT(IN), DIMENSION(NumberofRVETypes,NumberofElements*NumberofGaussPoints,2) :: ElementType2
INTEGER, INTENT(IN), DIMENSION(NumberofRVETypes) :: ElementTypeIndex
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsXX_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsYY_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsZZ_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsXY1_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsXY2_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsXZ1_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsXZ2_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsYZ1_R
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofElements,8) :: Sub_StrainsYZ2_R

! The nodes
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: M_Nodes
INTEGER, ALLOCATABLE, DIMENSION(:) :: M_NodeIndex
! Connectivities
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: M_Connectivities
! Boundaries
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: M_Boundaries
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: M_BoundaryIndex
! MechanicalMaterialProperties
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: M_CC
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) :: M_Sub_CC
! Matrix
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: M_Vol
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:,:) :: M_BB
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_A1,M_B1,M_C1,M_D1,M_A2,M_B2,M_C2,M_D2,M_A3,M_B3,M_C3,M_D3,M_A4,M_B4,M_C4,M_D4
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:,:) :: M_KComponent
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: M_K1GRR
INTEGER, ALLOCATABLE, DIMENSION (:) :: M_K1GRCDD
INTEGER, ALLOCATABLE, DIMENSION (:,:) :: M_K1GRCII
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_F1Globall
!Matrix solver
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_Dis1
!StrainsandStresses
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StrainsXX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StrainsYY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StrainsZZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StrainsXY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StrainsXZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StrainsYZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StressesXX
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StressesYY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StressesZZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StressesXY
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StressesXZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StressesYZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:) :: M_StressesVonMises
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: M_CurrentNodes 
INTEGER, ALLOCATABLE, DIMENSION (:) :: M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes
INTEGER :: M_NumberofDE
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: M_One_Step_DEs
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) ::  M_DNDX,M_DNDY,M_DNDZ
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: M_NN
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: M_Vol_Ele
DOUBLE PRECISION :: M_Kappa12,M_Mu12,M_Mu23,M_E33,M_Nu32
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: M_Sub_StrainsXX,M_Sub_StrainsYY,M_Sub_StrainsZZ,M_Sub_StrainsXY1,M_Sub_StrainsXZ1,M_Sub_StrainsYZ1,&
M_Sub_StrainsXY2,M_Sub_StrainsXZ2,M_Sub_StrainsYZ2,M_Sub_StressesXX,M_Sub_StressesYY,M_Sub_StressesZZ,M_Sub_StressesXY,M_Sub_StressesXZ,M_Sub_StressesYZ
DOUBLE PRECISION :: ElapsedTime
INTEGER :: M_TestIndex1,M_TestIndex2
INTEGER, INTENT(OUT) :: M_Highest_GRCD
DOUBLE PRECISION, INTENT(OUT), DIMENSION(TotalTimeSteps+1,E_NumberofGaussPoints,NumberofElements,3,3) ::E_sub_sigmasigma_in
DOUBLE PRECISION, INTENT(IN), DIMENSION(NumberofMaterialModes,3,3) :: RotationMatrix2
INTEGER, INTENT(IN) :: BeginofModeIndex,EndofModeIndex
INTEGER, INTENT(IN), DIMENSION(NumberofElements,17) :: Connectivities

CHARACTER*90 :: Boundaries_Microscale
CHARACTER*90 :: Connectivities_Microscale
CHARACTER*90 :: Control_Flags_Microscale
CHARACTER*90 :: Input_Strains_Microscale
CHARACTER*90 :: Mechanical_material_properties_Microscale
CHARACTER*90 :: Six_faces_Microscale
CHARACTER*90 :: The_Nodes_Microscale
CHARACTER*90 :: Six_faces_Microscale_Check

CHARACTER*90 :: Boundaries_Microscale_E
CHARACTER*90 :: Connectivities_Microscale_E
CHARACTER*90 :: Control_Flags_Microscale_E
CHARACTER*90 :: Electrical_material_properties_Microscale_E
CHARACTER*90 :: Piezo_Coefficients_Microscale_E
CHARACTER*90 :: More_Inputs_Microscale_E
CHARACTER*90 :: More_Inputs_Microscale_2_E
CHARACTER*90 :: Effec_Electro_Prop_Microscale_E

INTEGER, ALLOCATABLE, DIMENSION (:,:) :: ME_Connectivities
DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: ME_Boundaries
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ME_BoundaryIndex
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ME_sigmasigma
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ME_gg
INTEGER, INTENT(OUT) :: ME_Highest_GRCD
INTEGER, ALLOCATABLE, DIMENSION(:) :: M_NumberofElementsPerCNT
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ME_Coordinates
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ME_ElementsPerCNT
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: ME_TunnelingIndex
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ME_TunnelingIndexDimension
DOUBLE PRECISION :: ME_Width1,ME_Width2
INTEGER, ALLOCATABLE, DIMENSION(:) :: ME_CNTIndex
DOUBLE PRECISION :: ME_Critical_Distance
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: ME_Coordinates2
INTEGER :: ii3,ii4,ii3_Temp

CALL M_MicroFolderNames(ii,Boundaries_Microscale,Connectivities_Microscale,Control_Flags_Microscale,Input_Strains_Microscale,&
Mechanical_material_properties_Microscale,Six_faces_Microscale,The_Nodes_Microscale,Six_faces_Microscale_Check)
CALL M_ReadControlFlags (ii,Control_Flags_Microscale)

CALL ME_MicroFolderNames(ii,Boundaries_Microscale_E,Connectivities_Microscale_E,Control_Flags_Microscale_E,&
Electrical_material_properties_Microscale_E,More_Inputs_Microscale_E,More_Inputs_Microscale_2_E,Effec_Electro_Prop_Microscale_E,&
Piezo_Coefficients_Microscale_E)
CALL ME_ReadControlFlags (ii,Control_Flags_Microscale_E)

ALLOCATE(M_Sub_StrainsXX(M_NumberofElements,8))
ALLOCATE(M_Sub_StrainsYY(M_NumberofElements,8))
ALLOCATE(M_Sub_StrainsZZ(M_NumberofElements,8))
ALLOCATE(M_Sub_StrainsXY1(M_NumberofElements,8))
ALLOCATE(M_Sub_StrainsXZ1(M_NumberofElements,8))
ALLOCATE(M_Sub_StrainsYZ1(M_NumberofElements,8))
ALLOCATE(M_Sub_StrainsXY2(M_NumberofElements,8))
ALLOCATE(M_Sub_StrainsXZ2(M_NumberofElements,8))
ALLOCATE(M_Sub_StrainsYZ2(M_NumberofElements,8))
ALLOCATE(M_Sub_StressesXX(M_NumberofElements,8))
ALLOCATE(M_Sub_StressesYY(M_NumberofElements,8))
ALLOCATE(M_Sub_StressesZZ(M_NumberofElements,8))
ALLOCATE(M_Sub_StressesXY(M_NumberofElements,8))
ALLOCATE(M_Sub_StressesXZ(M_NumberofElements,8))
ALLOCATE(M_Sub_StressesYZ(M_NumberofElements,8))
! The microscale nodes
ALLOCATE(M_Nodes(M_NumberofNodes,4))
ALLOCATE(M_NodeIndex(M_NumberofNodes))
! Microscale Connectivities
ALLOCATE(M_Connectivities(M_NumberofElements,10))
! Microscale Boundaries
ALLOCATE(M_Boundaries(M_NumberofBoundaryNodes,TotalTimeSteps+2))
ALLOCATE(M_BoundaryIndex(M_NumberofBoundaryNodes,2))
! Microscale MechanicalMaterialProperties
ALLOCATE(M_CC(M_NumberofMaterialModes,6,6))
ALLOCATE(M_Sub_CC(M_NumberofGaussPoints,M_NumberofMaterialModes,6,6))
! Microscale Matrix
ALLOCATE(M_Vol(M_NumberofGaussPoints,M_NumberofElements))
ALLOCATE(M_Vol_Ele(M_NumberofElements))
ALLOCATE(M_BB(M_NumberofGaussPoints,M_NumberofElements,6,24))
ALLOCATE(M_A1(M_NumberofElements))
ALLOCATE(M_B1(M_NumberofElements))
ALLOCATE(M_C1(M_NumberofElements))
ALLOCATE(M_D1(M_NumberofElements))
ALLOCATE(M_A2(M_NumberofElements))
ALLOCATE(M_B2(M_NumberofElements))
ALLOCATE(M_C2(M_NumberofElements))
ALLOCATE(M_D2(M_NumberofElements))
ALLOCATE(M_A3(M_NumberofElements))
ALLOCATE(M_B3(M_NumberofElements))
ALLOCATE(M_C3(M_NumberofElements))
ALLOCATE(M_D3(M_NumberofElements))
ALLOCATE(M_A4(M_NumberofElements))
ALLOCATE(M_B4(M_NumberofElements))
ALLOCATE(M_C4(M_NumberofElements))
ALLOCATE(M_D4(M_NumberofElements))
ALLOCATE(M_KComponent(M_NumberofElements,24,24))
ALLOCATE(M_Dis1(3*M_NumberofNodes))
!StrainsandStresses
ALLOCATE(M_StrainsXX(M_NumberofElements))
ALLOCATE(M_StrainsYY(M_NumberofElements))
ALLOCATE(M_StrainsZZ(M_NumberofElements))
ALLOCATE(M_StrainsXY(M_NumberofElements))
ALLOCATE(M_StrainsXZ(M_NumberofElements))
ALLOCATE(M_StrainsYZ(M_NumberofElements))
ALLOCATE(M_StressesXX(M_NumberofElements))
ALLOCATE(M_StressesYY(M_NumberofElements))
ALLOCATE(M_StressesZZ(M_NumberofElements))
ALLOCATE(M_StressesXY(M_NumberofElements))
ALLOCATE(M_StressesXZ(M_NumberofElements))
ALLOCATE(M_StressesYZ(M_NumberofElements))
ALLOCATE(M_StressesVonMises(M_NumberofElements))
ALLOCATE(M_CurrentNodes(M_NumberofNodes,4)) 
ALLOCATE(M_PositiveXNodes(M_NumberofPositiveXNodes))
ALLOCATE(M_PositiveYNodes(M_NumberofPositiveYNodes))
ALLOCATE(M_PositiveZNodes(M_NumberofPositiveZNodes))
ALLOCATE(M_NegativeXNodes(M_NumberofNegativeXNodes))
ALLOCATE(M_NegativeYNodes(M_NumberofNegativeYNodes))
ALLOCATE(M_NegativeZNodes(M_NumberofNegativeZNodes))
ALLOCATE(M_DNDX(M_NumberofGaussPoints,M_NumberofElements,8+3))
ALLOCATE(M_DNDY(M_NumberofGaussPoints,M_NumberofElements,8+3))
ALLOCATE(M_DNDZ(M_NumberofGaussPoints,M_NumberofElements,8+3))
ALLOCATE(M_NN(M_NumberofGaussPoints,8))
ALLOCATE(M_One_Step_DEs(3*(M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+M_NumberofPositiveZNodes),8))

CALL M_ReadNodes (ii,The_Nodes_Microscale,M_NodeIndex,M_Nodes)
CALL M_ReadConnectivities (ii,Connectivities_Microscale,M_Connectivities)
CALL M_ReadBoundaryConditions(ii,Boundaries_Microscale,M_BoundaryIndex,M_Boundaries)
CALL M_ReadMechanicalProperties (ii,Mechanical_material_properties_Microscale,M_CC)
CALL M_SixFaces(ii,Six_faces_Microscale,Six_faces_Microscale_Check,M_Nodes,M_PositiveXNodes,M_PositiveYNodes,&
M_PositiveZNodes,M_NegativeXNodes,M_NegativeYNodes,M_NegativeZNodes)
CALL M_Components(M_Nodes,M_Connectivities,M_CC,M_Sub_CC,M_Vol,M_KComponent,M_DNDX,M_DNDY,M_DNDZ,M_NN)
M_NumberofDE=3*(M_NumberofPositiveXNodes+M_NumberofPositiveYNodes+M_NumberofPositiveZNodes)
ALLOCATE(M_K1GRR(3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec))
ALLOCATE(M_K1GRCDD(3*M_NumberofNodes+M_NumberofDE))
ALLOCATE(M_K1GRCII(3*M_NumberofNodes+M_NumberofDE,M_MaxDimenofInverseConnec))
ALLOCATE(M_F1Globall(3*M_NumberofNodes+M_NumberofDE))
CALL M_KPreprations(M_Connectivities,M_KComponent,M_NumberofDE,M_K1GRR,M_K1GRCDD,M_K1GRCII,M_F1Globall)


CALL ME_MoreInputs (ii,More_Inputs_Microscale_E)
ALLOCATE(M_NumberofElementsPerCNT(ME_NumberofCNTs))
ALLOCATE(ME_Coordinates(ME_NumberofCNTs2,2))
ALLOCATE(ME_ElementsPerCNT(ME_NumberofCNTs,ME_MaxNumberofElementsPerCNT))
CALL ME_MoreInputs2 (ii,More_Inputs_Microscale_2_E,M_NumberofElementsPerCNT,ME_Coordinates,ME_ElementsPerCNT)
ALLOCATE(ME_Connectivities(M_NumberofElements,10))
CALL ME_ReadConnectivities (ii,Connectivities_Microscale_E,ME_Connectivities)
ALLOCATE(ME_Boundaries(ME_NumberofBoundaryNodes,TotalTimeSteps+1))
ALLOCATE(ME_BoundaryIndex(ME_NumberofBoundaryNodes,1))
CALL ME_ReadBoundaryConditions(ii,Boundaries_Microscale_E,ME_BoundaryIndex,ME_Boundaries)
ALLOCATE(ME_sigmasigma(ME_NumberofMaterialModes,3,3))
ALLOCATE(ME_gg(ME_NumberofMaterialModes,6,6))
CALL ME_ReadElectricalProperties (ii,Electrical_material_properties_Microscale_E,Piezo_Coefficients_Microscale_E,ME_sigmasigma,ME_gg)

ALLOCATE(ME_TunnelingIndex(ME_NumberofGaussPoints,M_NumberofElements,4,3))
ALLOCATE(ME_TunnelingIndexDimension(ME_NumberofGaussPoints,M_NumberofElements))
ALLOCATE(ME_CNTIndex(ME_NumberofCNTs2))
ALLOCATE(ME_Coordinates2(9,ME_NumberofCNTs2,2))
CALL ME_Prep_Tunneling(M_Nodes,ME_Connectivities,M_NN,ME_Coordinates,M_NumberofElementsPerCNT,ME_ElementsPerCNT,ME_TunnelingIndex,&
ME_TunnelingIndexDimension,ME_Width1,ME_Width2,ME_CNTIndex,ME_Critical_Distance,ME_Coordinates2)

OPEN(51,FILE="3D_Output_Micro_F/3D_Output/Test.txt",STATUS="UNKNOWN")
WRITE(51,*) "START"
CLOSE(51)
IF (TestType==4) THEN
ii3_Temp=1
ELSE
ii3_Temp=ElementTypeIndex(ii)
ENDIF
!$omp parallel private(ii1,ii2,ii3)
!$omp DO
DO ii3=1,ii3_Temp  !2008,2008
ii1=ElementType2(ii,ii3,1)
ii2=ElementType2(ii,ii3,2)
CALL Macro_Sub_Element(k,ii,ii3,ii1,ii2,M_Nodes,M_NodeIndex,M_DNDX,M_DNDY,M_DNDZ,M_Sub_CC,M_CC,M_Vol,M_Boundaries,M_BoundaryIndex,M_NumberofDE,&
M_K1GRR,M_K1GRCDD,M_K1GRCII,M_F1Globall,M_Connectivities,M_NN,M_PositiveXNodes,M_PositiveYNodes,M_PositiveZNodes,M_NegativeXNodes,&
M_NegativeYNodes,M_NegativeZNodes,Sub_StrainsXX_R,Sub_StrainsYY_R,Sub_StrainsZZ_R,Sub_StrainsXY1_R,Sub_StrainsXY2_R,Sub_StrainsXZ1_R,&
Sub_StrainsXZ2_R,Sub_StrainsYZ1_R,Sub_StrainsYZ2_R,M_Highest_GRCD,ME_sigmasigma,ME_gg,ME_Boundaries,ME_BoundaryIndex,ME_Connectivities,ME_Highest_GRCD,&
ME_Coordinates,M_NumberofElementsPerCNT,ME_ElementsPerCNT,ME_TunnelingIndex,ME_TunnelingIndexDimension,ME_Width1,ME_Width2,ME_CNTIndex,&
ME_Critical_Distance,ME_Coordinates2,E_sub_sigmasigma_in,Effec_Electro_Prop_Microscale_E,RotationMatrix2,BeginofModeIndex,EndofModeIndex,&
Connectivities)
OPEN(51,FILE="3D_Output_Micro_F/3D_Output/Test.txt",STATUS="UNKNOWN",POSITION="APPEND")
WRITE(51,*) "k,ii3,ii1,ii2",k,ii3,ii1,ii2
CLOSE(51)
ENDDO
!$omp END DO
!$omp end parallel
IF (k==TotalTimeSteps .AND. ElementTypeIndex(ii) > 0 ) THEN
CALL M_TecplotGenerator(ii)
DO ii3=1,5
CALL ME_TecplotGenerator(ii,ii3)
ENDDO
ENDIF

END SUBROUTINE Macro_Each_RVE_Type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadElementType(ElementType,ElementType2,ElementTypeIndex)
IMPLICIT NONE
INTEGER :: i,j,i1,IN_ElementType=210

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

INTEGER,INTENT(OUT),DIMENSION(NumberofElements,9) :: ElementType
INTEGER, INTENT(OUT), DIMENSION(NumberofRVETypes,NumberofElements*NumberofGaussPoints,2) :: ElementType2
INTEGER, INTENT(OUT), DIMENSION(NumberofRVETypes) :: ElementTypeIndex

OPEN(IN_ElementType,file = '3D_Output_Micro_F/3D_Input/Element_Type.txt',status ='unknown')

READ(IN_ElementType,*) ((ElementType(i,j),j=1,9),i=1,NumberofElements)

CLOSE(IN_ElementType)

ElementType2=0
ElementTypeIndex=0

DO i=1,NumberofElements
DO j=1,8

DO i1=1,NumberofRVETypes
IF (ElementType(i,j+1) == i1) THEN
ElementTypeIndex(i1)=ElementTypeIndex(i1)+1
ElementType2(i1,ElementTypeIndex(i1),1)=i
ElementType2(i1,ElementTypeIndex(i1),2)=j
GOTO 1
ENDIF
1 ENDDO

ENDDO
ENDDO

RETURN
ENDSUBROUTINE ReadElementType

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE ReadStrains (Sub_StrainsXX,Sub_StrainsYY,Sub_StrainsZZ,Sub_StrainsXY1,Sub_StrainsXY2,Sub_StrainsXZ1,Sub_StrainsXZ2,&
!Sub_StrainsYZ1,Sub_StrainsYZ2)
!IMPLICIT NONE
!INTEGER :: IN11=9,IN12=10,IN13=11,IN14=12,IN15=13,IN16=14,IN17=15,IN18=16,IN19=17
!INTEGER :: i,j
!
!INTEGER :: NumberofNodes
!INTEGER :: NumberofElements
!INTEGER :: NumberofBoundaryNodes
!INTEGER :: NumberofMaterialModes
!INTEGER :: NumberofGaussPoints
!INTEGER :: TotalTimeSteps
!INTEGER :: NumberofPositiveXNodes
!INTEGER :: NumberofPositiveYNodes
!INTEGER :: NumberofPositiveZNodes
!INTEGER :: NumberofNegativeXNodes
!INTEGER :: NumberofNegativeYNodes
!INTEGER :: NumberofNegativeZNodes
!INTEGER :: MaxDimenofInverseConnec
!INTEGER :: NumberofRVETypes
!INTEGER :: ElementShown
!DOUBLE PRECISION :: Mag
!COMMON/Control_Flags/ Mag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,NumberofGaussPoints,TotalTimeSteps,&
!NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,NumberofNegativeYNodes,NumberofNegativeZNodes,&
!MaxDimenofInverseConnec,NumberofRVETypes,ElementShown
!
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofElements,8) :: Sub_StrainsXX
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofElements,8) :: Sub_StrainsYY
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofElements,8) :: Sub_StrainsZZ
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofElements,8) :: Sub_StrainsXY1
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofElements,8) :: Sub_StrainsXY2
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofElements,8) :: Sub_StrainsXZ1
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofElements,8) :: Sub_StrainsXZ2
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofElements,8) :: Sub_StrainsYZ1
!DOUBLE PRECISION, INTENT(OUT), DIMENSION(NumberofElements,8) :: Sub_StrainsYZ2
!
!OPEN(IN11,file = '3D_Output_Micro_F/3D_Output/Sub_Strains_11.txt',status ='unknown')
!OPEN(IN12,file = '3D_Output_Micro_F/3D_Output/Sub_Strains_22.txt',status ='unknown')
!OPEN(IN13,file = '3D_Output_Micro_F/3D_Output/Sub_Strains_33.txt',status ='unknown')
!OPEN(IN14,file = '3D_Output_Micro_F/3D_Output/Sub_Strains_12_1.txt',status ='unknown')
!OPEN(IN15,file = '3D_Output_Micro_F/3D_Output/Sub_Strains_12_2.txt',status ='unknown')
!OPEN(IN16,file = '3D_Output_Micro_F/3D_Output/Sub_Strains_13_1.txt',status ='unknown')
!OPEN(IN17,file = '3D_Output_Micro_F/3D_Output/Sub_Strains_13_2.txt',status ='unknown')
!OPEN(IN18,file = '3D_Output_Micro_F/3D_Output/Sub_Strains_23_1.txt',status ='unknown')
!OPEN(IN19,file = '3D_Output_Micro_F/3D_Output/Sub_Strains_23_2.txt',status ='unknown')
!
!READ(IN11,*) ((Sub_StrainsXX(i,j),j=1,8),i=1,NumberofElements)
!READ(IN12,*) ((Sub_StrainsYY(i,j),j=1,8),i=1,NumberofElements)
!READ(IN13,*) ((Sub_StrainsZZ(i,j),j=1,8),i=1,NumberofElements)
!READ(IN14,*) ((Sub_StrainsXY1(i,j),j=1,8),i=1,NumberofElements)
!READ(IN15,*) ((Sub_StrainsXY2(i,j),j=1,8),i=1,NumberofElements)
!READ(IN16,*) ((Sub_StrainsXZ1(i,j),j=1,8),i=1,NumberofElements)
!READ(IN17,*) ((Sub_StrainsXZ2(i,j),j=1,8),i=1,NumberofElements)
!READ(IN18,*) ((Sub_StrainsYZ1(i,j),j=1,8),i=1,NumberofElements)
!READ(IN19,*) ((Sub_StrainsYZ2(i,j),j=1,8),i=1,NumberofElements)
!
!CLOSE(IN11)
!CLOSE(IN12)
!CLOSE(IN13)
!CLOSE(IN14)
!CLOSE(IN15)
!CLOSE(IN16)
!CLOSE(IN17)
!CLOSE(IN18)
!CLOSE(IN19)
!
!RETURN
!
!END SUBROUTINE ReadStrains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ReadControlFlags()
IMPLICIT NONE
INTEGER :: IN1=5
INTEGER :: i
DOUBLE PRECISION, DIMENSION (18) :: ControlFlags

!The control flags
INTEGER :: TestType
INTEGER :: DisplacementEquationFlag
INTEGER :: NumberofNodes
INTEGER :: NumberofElements
INTEGER :: NumberofBoundaryNodes
INTEGER :: NumberofMaterialModes
INTEGER :: NumberofGaussPoints
INTEGER :: TotalTimeSteps
INTEGER :: NumberofPositiveXNodes
INTEGER :: NumberofPositiveYNodes
INTEGER :: NumberofPositiveZNodes
INTEGER :: NumberofNegativeXNodes
INTEGER :: NumberofNegativeYNodes
INTEGER :: NumberofNegativeZNodes
INTEGER :: MaxDimenofInverseConnec
INTEGER :: NumberofRVETypes
INTEGER :: ElementShown
INTEGER :: SubElementShown
DOUBLE PRECISION :: Mag
COMMON/Control_Flags/ Mag,TestType,DisplacementEquationFlag,NumberofNodes,NumberofElements,NumberofBoundaryNodes,NumberofMaterialModes,&
NumberofGaussPoints,TotalTimeSteps,NumberofPositiveXNodes,NumberofPositiveYNodes,NumberofPositiveZNodes,NumberofNegativeXNodes,&
NumberofNegativeYNodes,NumberofNegativeZNodes,MaxDimenofInverseConnec,NumberofRVETypes,ElementShown,SubElementShown

OPEN(IN1,file = '3D_Output_Micro_F/3D_Input/Control_Flags.txt',status ='unknown')
DO i=1,19
READ (IN1,*) ControlFlags(i)
ENDDO
TestType=ControlFlags(1)
DisplacementEquationFlag=ControlFlags(2)
NumberofNodes=ControlFlags(3)
NumberofElements=ControlFlags(4)
NumberofBoundaryNodes=ControlFlags(5)
NumberofMaterialModes=ControlFlags(6)
NumberofGaussPoints=ControlFlags(7)
TotalTimeSteps=ControlFlags(8)
Mag=ControlFlags(9)
NumberofPositiveXNodes=ControlFlags(10)
NumberofPositiveYNodes=ControlFlags(11)
NumberofPositiveZNodes=ControlFlags(12)
NumberofNegativeXNodes=ControlFlags(13)
NumberofNegativeYNodes=ControlFlags(14)
NumberofNegativeZNodes=ControlFlags(15)
MaxDimenofInverseConnec=ControlFlags(16)
NumberofRVETypes=ControlFlags(17)
ElementShown=ControlFlags(18)
SubElementShown=ControlFlags(19)
ClOSE(IN1)

RETURN
END SUBROUTINE ReadControlFlags