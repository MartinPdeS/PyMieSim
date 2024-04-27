subroutine bjndd ( n, x, bj, dj, fj )

  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) bj(n+1)
  real ( kind = 8 ) bs
  real ( kind = 8 ) dj(n+1)
  real ( kind = 8 ) f
  real ( kind = 8 ) f0
  real ( kind = 8 ) f1
  real ( kind = 8 ) fj(n+1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mt
  integer ( kind = 4 ) nt
  real ( kind = 8 ) x

  do nt = 1, 900
    mt = int ( 0.5D+00 * log10 ( 6.28D+00 * nt ) &
      - nt * log10 ( 1.36D+00 * abs ( x ) / nt ) )
    if ( 20 < mt ) then
      exit
    end if
  end do

  m = nt
  bs = 0.0D+00
  f0 = 0.0D+00
  f1 = 1.0D-35
  do k = m, 0, -1
    f = 2.0D+00 * ( k + 1.0D+00 ) * f1 / x - f0
    if ( k <= n ) then
      bj(k+1) = f
    end if
    if ( k == 2 * int ( k / 2 ) ) then
      bs = bs + 2.0D+00 * f
    end if
    f0 = f1
    f1 = f
  end do

  do k = 0, n
    bj(k+1) = bj(k+1) / ( bs - f )
  end do

  dj(1) = -bj(2)
  fj(1) = -1.0D+00 * bj(1) - dj(1) / x
  do k = 1, n
    dj(k+1) = bj(k) - k * bj(k+1) / x
    fj(k+1) = ( k * k / ( x * x ) - 1.0D+00 ) * bj(k+1) - dj(k+1) / x
  end do

  return
end

subroutine jdzo(nt, n, m, zo) bind(C, name="jdzo")
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), intent(in) :: nt
    integer(c_int), dimension(1400), intent(out) :: n, m
    real(c_double), dimension(1400), intent(out) :: zo

  real ( kind = 8 ) bj(101)
  real ( kind = 8 ) dj(101)
  real ( kind = 8 ) fj(101)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) m1(70)
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) n1(70)
  integer ( kind = 4 ) nm
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xm
  real ( kind = 8 ) zoc(70)

  if ( nt < 600 ) then
    xm = -1.0D+00 + 2.248485D+00 * real ( nt, kind = 8 ) ** 0.5D+00 &
      - 0.0159382D+00 * nt + 3.208775D-04 * real ( nt, kind = 8 ) ** 1.5D+00
    nm = int ( 14.5D+00 + 0.05875D+00 * nt )
    mm = int ( 0.02D+00 * nt ) + 6
  else
    xm = 5.0D+00 + 1.445389D+00 * ( real ( nt, kind = 8 ) ) ** 0.5D+00 &
      + 0.01889876D+00 * nt &
      - 2.147763D-04 * ( real ( nt, kind = 8 ) ) ** 1.5D+00
    nm = int ( 27.8D+00 + 0.0327D+00 * nt )
    mm = int ( 0.01088D+00 * nt ) + 10
  end if

  l0 = 0

  do i = 1, 70
    x1 = 0.407658D+00 + 0.4795504D+00 * ( real ( i - 1, kind = 8 ) ) ** 0.5D+00 + 0.983618D+00 * ( i - 1 )
    x2 = 1.99535D+00 + 0.8333883 * ( real ( i - 1, kind = 8 ) ) ** 0.5D+00 + 0.984584D+00 * ( i - 1 )
    l1 = 0

    do j = 1, mm
      if ( i == 1 .and. j == 1 ) then

        l1 = l1 + 1
        n1(l1) = i - 1
        m1(l1) = j
        if ( i == 1 ) then
          m1(l1) = j - 1
        end if
        zoc(l1) = x

        if ( i <= 15 ) then
          x1 = x + 3.057D+00 + 0.0122D+00 * ( i - 1 ) + ( 1.555D+00 + 0.41575D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
        else
          x1 = x + 2.918D+00 + 0.01924D+00 * ( i - 1 ) + ( 6.26D+00 + 0.13205D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
        end if

      else
        x = x1

        do
          call bjndd ( i, x, bj, dj, fj )

          x0 = x
          x = x - dj(i) / fj(i)

          if ( xm < x1 ) then
            exit
          end if

          if ( abs ( x - x0 ) <= 1.0D-10 ) then
            l1 = l1 + 1
            n1(l1) = i - 1
            m1(l1) = j
            if ( i == 1 ) then
              m1(l1) = j - 1
            end if
            zoc(l1) = x

            if ( i <= 15 ) then
              x1 = x + 3.057D+00 + 0.0122D+00 * ( i - 1 ) + ( 1.555D+00 + 0.41575D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
            else
              x1 = x + 2.918D+00 + 0.01924D+00 * ( i - 1 ) + ( 6.26D+00 + 0.13205D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
            end if
            exit
          end if

        end do

      end if

      x = x2

      do

        call bjndd ( i, x, bj, dj, fj )
        x0 = x
        x = x - bj(i) / dj(i)

        if ( xm < x ) then
          exit
        end if

        if ( abs ( x - x0 ) <= 1.0D-10 ) then
          exit
        end if

      end do

      if ( x <= xm ) then

        l1 = l1 + 1
        n1(l1) = i - 1
        m1(l1) = j
        zoc(l1) = x
        if ( i <= 15 ) then
          x2 = x + 3.11D+00 + 0.0138D+00 * ( i - 1 ) + ( 0.04832D+00 + 0.2804D+00 * ( i - 1 ) ) / ( j + 1 ) ** 2
        else
          x2 = x + 3.001D+00 + 0.0105D+00 * ( i - 1 ) &
            + ( 11.52D+00 + 0.48525D+00 * ( i - 1 ) ) / ( j + 3 ) ** 2
        end if

      end if

    end do

    l = l0 + l1
    l2 = l
    do
      if ( l0 == 0 ) then
        do k = 1, l
          zo(k) = zoc(k)
          n(k) = n1(k)
          m(k) = m1(k)
        end do
        l1 = 0
      else if ( l0 /= 0 ) then
        if ( zoc(l1) .le. zo(l0) ) then
          zo(l0+l1) = zo(l0)
          n(l0+l1) = n(l0)
          m(l0+l1) = m(l0)
          l0 = l0 - 1
        else
          zo(l0+l1) = zoc(l1)
          n(l0+l1) = n1(l1)
          m(l0+l1) = m1(l1)
          l1 = l1 - 1
        end if
      end if

      if ( l1 == 0 ) then
        exit
      end if

    end do

    l0 = l2
  PRINT *, i
  end do

  return
end