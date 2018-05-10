import sympy
from sympy import sin, cos, diff
x,t = sympy.symbols('y t')
v=sin(x)
rho=cos(x)+4+.01*x**2*sin(t)
m=rho*v
sympy.diff(m,x)
sympy.diff(-m,x)
rho_s=diff(rho,t)+diff(m,x)
p0=sympy.symbols('p0')
p=rho*p0
v_s=rho*diff(v,t)+v*diff(v,x)+diff(p,x)

m_s=rho*v_s

from boutdata import mms
print('[rho]')
print('solution = ',mms.exprToStr(rho))
print('source = ',mms.exprToStr(rho_s))
#print('[m]')
#print('solution =',mms.exprToStr(m))
#print('source =',mms.exprToStr(m_s))
print('[v_n]')
print('solution =',mms.exprToStr(v))
print('source =',mms.exprToStr(v_s))
