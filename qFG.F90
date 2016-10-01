subroutine qv(q,i,j)
use Global
real q(4)
integer i,j

q(1) = rho(i,j)
q(2) = rho(i,j)*u(i,j)
q(3) = rho(i,j)*v(i,j)
q(4) = E(i,j)

end subroutine

subroutine Fv(F,i,j)
implicit none
use Global
    real F(4,0:N,0:M)
    integer i,j

    F(1,i,j) = rho(i,j)*u(i,j)
    F(2,i,j) = rho(i,j)*u(i,j)*u(i,j) + p(i,j) - txx(i,j)
    F(3,i,j) = rho(i,j)*u(i,j)*v(i,j) - txy(i,j)
    F(4,i,j) = (E(i,j)+p(i,j)-txx(i,j))*u(i,j) - txy(i,j)*v(i,j)+Qdx(i,j)

end subroutine


subroutine Gv(G,i,j)
implicit none
use Global
    real G(4,0:N,0:M)
    integer i,j

    G(1,i,j) = rho(i,j)*v(i,j)
    G(2,i,j) = rho(i,j)*u(i,j)*v(i,j) - txy(i,j)
    G(3,i,j) = rho(i,j)*v(i,j)*v(i,j) + p(i,j) - tyy(i,j)
    G(4,i,j) = (E(i,j) + p(i,j) - tyy(i,j))*v(i,j) - txy(i,j)*u(i,j)+Qdy(i,j)

end subroutine