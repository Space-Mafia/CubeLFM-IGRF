        function cross(x1,x2)
! ... Parameter variables ..............................................        
        dimension x1(2)
!> x1 - first 2d vector to take cross product with.       
        dimension x2(2)
!> x2 - second 2d vector to take cross product with.       
! ... Begin ............................................................
        cross = x1(1)*x2(2) - x2(1)*x1(2)
        return
        end
