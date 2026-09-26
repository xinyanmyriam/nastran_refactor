// Standard wedge stiffness
// Nodes: (0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,0,1),(0,1,1)
// Natural coords (r,s,t), r,s in [0,1], r+s<=1, t in [0,1]
// Shape functions:
// N1=(1-r-s)(1-t), N2=r(1-t), N3=s(1-t)
// N4=(1-r-s)t, N5=r t, N6=s t
// Jacobian: x=r, y=s, z=t -> J=I
// B matrix from derivatives
// K = integral B^T D B dV