from dolfin import *
import matplotlib.pylab as plt
import numpy
 
plt.ion()
 
#while(t < 5.0):                 # we'll limit ourselves to 5 seconds.
#                                # set this to while(True) if you want to loop forever
#    t = time.time() - starttime # find out how long the script has been running
#    y = -2*sin(x)*sin(t)        # just a function for a standing wave
#                                # replace this with any function you want to animate
#                                # for instance, y = sin(x-t)
# 
#    line.set_ydata(y)           # update the plot data
#    draw()                      # redraw the canvas

N = 1000 # number of subintervals
p_degree = 1 # polyn. degree of approx.

eps = Constant(1E-6)
alfa = Constant(1E-3)
dt = Constant(1E-2)
h = Constant(-2.)
T = 2

mesh = UnitInterval(N)
mesh.coordinates()[:] = mesh.coordinates()*pi

Veta = FunctionSpace(mesh, "Lagrange", p_degree)
Vu = FunctionSpace(mesh, "Lagrange", p_degree)

eta = TrialFunction(Veta)  # solution
u = TrialFunction(Vu)  # solution

phi = TestFunction(Veta) # test
psi = TestFunction(Vu) # test

u_ = Function(Vu)
eta_ = Function(Veta)

u_ = interpolate(Constant(0.0), Vu) # initial conditions
eta_ = interpolate(Expression("sin(x[0])"), Veta)

a_eta = eta*psi*dx - dt*alfa*eta*u_*Dx(psi, 0)*dx
l_eta = eta_*psi*dx + dt*h*u_*Dx(psi, 0)*dx

a_u = u*phi*dx - dt*alfa/2.*u_*u*Dx(phi, 0)*dx + eps/3.*h**2*Dx(u,0)*Dx(phi,0)*dx
l_u = eta_*phi*dx + u_*phi*dx + eps/3.*h**2*Dx(u_,0)*Dx(phi,0)*dx

bc_eta = DirichletBC(Veta, Constant(0.), DomainBoundary())

u = Function(Vu)
eta = Function(Veta)

t = 0
x = mesh.coordinates()
line, = plt.plot(x,x)  
y_eta = eta_.vector().array()
line.set_ydata(y_eta)           

while t < T:
	t += float(dt)

	solve(a_eta==l_eta, eta, bc_eta)
	eta_.assign(eta)

	solve(a_u==l_u, u)
	u_.assign(u)

	print t 
#	y_eta = eta_.vector().array()
#	line.set_ydata(y_eta)
	
	#print y_eta
	y_u = u_.vector().array()
	line.set_ydata(y_u)
	plt.draw()

raw_input("")

