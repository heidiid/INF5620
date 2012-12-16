from dolfin import *
import matplotlib.pylab as plt
import numpy
 
plt.ion()

N = 100 # number of subintervals
p_degree = 1 # polyn. degree of approx.

eps = Constant(1E-3)
alfa = Constant(0.5)
dt = Constant(1E-2)
h = Constant(-2.)
T = 10E-2

mesh = UnitInterval(N)
mesh.coordinates()[:] = mesh.coordinates()*pi

V = FunctionSpace(mesh, "Lagrange", p_degree)

(eta, u) = TrialFunction(W)  # solution

(phi, psi) = TestFunction(W) # test

u_ = Function(V)
eta_ = Function(V)

u_ = interpolate(Constant(0.0), V) # initial conditions
eta_ = interpolate(Expression("sin(x[0])"), V)

a = -dt*eta*Dx(psi,0)*dx + eta*phi*dx - dt*(h+alfa*eta_)*u*Dx(psi,0)*dx + u*psi*dx - dt*alfa/2.*u_*u*Dx(psi,0)*dx \
    + eps/3.*h**2*Dx(u,0)*Dx(psi,0)*dx
l = u_*psi*dx + eps/3.*h**2*Dx(u_,0)*Dx(psi,0)*dx + eta_*phi*dx

bc_eta = DirichletBC(W.sub(0), Constant(0.), DomainBoundary())
bc_u = DirichletBC(W.sub(1), Constant(0.), DomainBoundary())
bcs = [bc_eta, bc_u]

U = Function(W)
u = Function(V)
eta = Function(V)

t = 0
x = mesh.coordinates()
line, = plt.plot(x,x)  
y_eta = eta_.vector().array()
line.set_ydata(y_eta)           

while t < T:
	t += float(dt)
	print t 
	y_eta = eta_.vector().array()
	line.set_ydata(y_eta)
	
	solve(a==l, U, bcs)	
	(eta, u) = U.split(deepcopy=True)
	
	eta_.assign(eta)
	u_.assign(u)

	y_eta = eta_.vector().array()
	line.set_ydata(y_eta)
	#print y_eta
#	y_u = u_.vector().array()
#	line.set_ydata(y_u)
	plt.draw()

raw_input("")

