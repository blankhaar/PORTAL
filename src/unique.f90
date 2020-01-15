Module SortUnique
contains
  Recursive Subroutine MergeSort(temp, Begin, Finish, list)
    ! 1st 3 arguments are input, 4th is output sorted list
    implicit none
    integer(kind=4),intent(inout) :: Begin,list(:),temp(:)
    integer(kind=4),intent(in) :: Finish
    integer(kind=4) :: Middle
    if (Finish-Begin<2) then    !if run size =1
       return                   !it is sorted
    else
       ! split longer runs into halves
       Middle = (Finish+Begin)/2
       ! recursively sort both halves from list into temp
       call MergeSort(list, Begin, Middle, temp)
       call MergeSort(list, Middle, Finish, temp)
       ! merge sorted runs from temp into list
       call Merge(temp, Begin, Middle, Finish, list)
     endif
  End Subroutine MergeSort

  Subroutine Merge(list, Begin, Middle, Finish, temp)
    implicit none
    integer(kind=4),intent(inout) :: list(:),temp(:)
    integer(kind=4),intent(in) ::Begin,Middle,Finish
    integer(kind=4)    :: kx,ky,kz
    ky=Begin
    kz=Middle
    !! While there are elements in the left or right runs...
    do kx=Begin,Finish-1
       !! If left run head exists and is <= existing right run head.
       if (ky.lt.Middle.and.(kz.ge.Finish.or.list(ky).le.list(kz))) then
          temp(kx)=list(ky)
          ky=ky+1
       else
          temp(kx)=list(kz)
          kz = kz + 1
       end if
    end do

  End Subroutine Merge

  Function Unique(list)
    !! usage sortedlist=Unique(list)
    implicit none
    integer(kind=4) :: strt,fin,N
    integer(kind=4), intent(inout) :: list(:)
    integer(kind=4), allocatable  :: unique(:),work(:)
    logical,allocatable :: mask(:)
    ! sort
    work=list;strt=1;N=size(list);fin=N+1
    call MergeSort(work,strt,fin,list) 
    ! cull duplicate indices
    allocate(mask(N));
    mask=.false.
    mask(1:N-1)=list(1:N-1)==list(2:N)
    unique=pack(list,.not.mask)

  End Function Unique

End Module SortUnique

