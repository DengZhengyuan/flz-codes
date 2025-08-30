subroutine USR1
  ! 每个时间步调用：把 5 个 y=const 水平线上的整行 cell 逐个写入各自 CSV
  use run,       only: TIME
  use compar,    only: istart3, iend3, jstart3, jend3, kstart3
  use functions, only: FUNIJK
  use geometry,  only: DX, DY
  use cutcell,   only: XG_E, YG_N
  use rxns,      only: SPECIES_ALIAS_g
  use fldvar,    only: U_G, V_G, U_S, V_S, EP_G, EP_S, X_G, RO_G
  implicit none

  integer, parameter :: NPLANES = 5
  double precision, save :: PLANES(NPLANES) = (/0.4d0,0.6d0,1.0d0,1.5d0,2.0d0/)
  integer,          save :: Jline(NPLANES) = 0
  integer,          save :: unitno(NPLANES)
  character(len=128),save :: fname(NPLANES)
  logical,          save :: inited = .false.
  integer,          save :: NG_O3 = -1      ! O3 的气相组分下标

  integer :: p, i, j, k, ijk, ng
  double precision :: xc, yc, dyc, y0c

  if (.not. inited) then
     ! ===== 等距网格：用公式映射每个目标高度到 J 行 =====
     dyc = DY(jstart3)
     y0c = YG_N(jstart3) - 0.5d0*DY(jstart3)
     do p = 1, NPLANES
        Jline(p) = nint((PLANES(p) - y0c)/dyc) + jstart3
        if (Jline(p) < jstart3) Jline(p) = jstart3
        if (Jline(p) > jend3  ) Jline(p) = jend3
     end do

     ! ===== 通过别名查找 O3 的下标 =====
     NG_O3 = -1
     do ng = 1, size(SPECIES_ALIAS_g)
        if (upper(SPECIES_ALIAS_g(ng)) == 'O3') then
           NG_O3 = ng
           exit
        end if
     end do
     if (NG_O3 < 0) then
        write(*,*) 'WARNING: O3 alias not found. Set NG_O3 manually!'
        NG_O3 = 2   ! 按你 mfx 配置兜底为 2
     end if

     ! ===== 打开 5 个输出文件（各写各的）=====
     do p = 1, NPLANES
        write(fname(p),'(A,F0.2,A)') 'line_y', PLANES(p), 'm.csv'
        unitno(p) = 70 + p
        open(unit=unitno(p), file=fname(p), status='unknown', position='append', recl=100000)
        write(unitno(p),'(A)') 'time,x,y,rho_g,u_g,v_g,ep_g,ep_s,u_s,v_s,Xg_O3'
     end do

     inited = .true.
  end if

  k = kstart3   ! 2D：唯一的 k 层

  ! ===== 对每个 y 平面，写整行（该 J 行上所有 I）=====
  do p = 1, NPLANES
     j = Jline(p)
     do i = istart3, iend3
        ijk = FUNIJK(i, j, k)
        xc  = XG_E(i) - 0.5d0*DX(i)
        yc  = YG_N(j) - 0.5d0*DY(j)

        !$omp critical
        write(unitno(p),'(*(G0.8,:,","))') &
               TIME, xc, yc, &
               RO_G(ijk), &
               U_G(ijk), V_G(ijk), &
               EP_G(ijk), EP_S(ijk,1), &
               U_S(ijk,1), V_S(ijk,1), X_G(ijk, NG_O3)
        !$omp end critical
     end do
     flush(unitno(p))
  end do

contains
  pure function upper(s) result(t)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: t
    integer :: n, c
    do n = 1, len(s)
       c = iachar(s(n:n))
       if (c >= 97 .and. c <= 122) then
          t(n:n) = achar(c - 32)
       else
          t(n:n) = s(n:n)
       end if
    end do
  end function upper
end subroutine USR1