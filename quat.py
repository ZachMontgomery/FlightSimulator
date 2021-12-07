from math import asin, atan2, sin, cos, sqrt

def QuatMult(a,b):
	a0,a1,a2,a3=a
	b0,b1,b2,b3=b
	return a0*b0-a1*b1-a2*b2-a3*b3,a0*b1+a1*b0+a2*b3-a3*b2,a0*b2-a1*b3+a2*b0+a3*b1,a0*b3+a1*b2-a2*b1+a3*b0

def Euler2Quat(a):
	phi2=a[0]/2.
	theta2=a[1]/2.
	psi2=a[2]/2.
	sb=sin(phi2)
	cb=cos(phi2)
	se=sin(theta2)
	ce=cos(theta2)
	sh=sin(psi2)
	ch=cos(psi2)
	cbch=cb*ch
	sbsh=sb*sh
	sbch=sb*ch
	cbsh=cb*sh
	return cbch*ce+sbsh*se,sbch*ce-cbsh*se,cbch*se+sbsh*ce,cbsh*ce-sbch*se

def Quat2Euler(e):
	e0,e1,e2,e3=e
	test=e0*e2-e1*e3
	if test>=.5:
		a=2.*asin(e1/0.70710678118654757),1.5707963267948966,0.
	elif test<=-.5:
		a=2.*asin(e1/0.70710678118654757),-1.5707963267948966,0.
	else:
		c1=e0*e0-e2*e2
		c2=e1*e1-e3*e3
		a=atan2(2.*(e0*e1+e2*e3),c1-c2),asin(2.*test),atan2(2.*(e0*e3+e1*e2),c1+c2)
	return a

def NormalizeQuaternion(e):
	e0,e1,e2,e3=e
	d=1.5-0.5*(e0*e0+e1*e1+e2*e2+e3*e3)
	return e0*d,e1*d,e2*d,e3*d

def Body2Fixed(v,e):
	e0,e1,e2,e3=e
	v0,v1,v2=v
	c0=v0*e1+v1*e2+v2*e3
	c1=v0*e0-v1*e3+v2*e2
	c2=v1*e0-v2*e1+v0*e3
	c3=v1*e1+v2*e0-v0*e2
	return e0*c1+e1*c0+e2*c3-e3*c2,e0*c2-e1*c3+e2*c0+e3*c1,e0*c3+e1*c2-e2*c1+e3*c0

def Fixed2Body(v,e):
	e0,e1,e2,e3=e
	v0,v1,v2=v
	c0=-e1*v0-e2*v1-e3*v2
	c1=e3*v1+e0*v0-e2*v2
	c2=e1*v2+e0*v1-e3*v0
	c3=e2*v0+e0*v2-e1*v1
	return c1*e0-c0*e1+c2*e3-c3*e2,c2*e0-c0*e2-c1*e3+c3*e1,c1*e2-c0*e3-c2*e1+c3*e0

