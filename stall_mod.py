import numpy as np

d2r = np.pi / 180.
as1 = 10. * d2r
as2 = 30. * d2r
Wref = 2800.
rhoref = .0023769
Vref = 180.
Sw = 185.
bw = 33.
c_w = Sw / bw
CLref = 2. * Wref / rhoref / Vref ** 2. / Sw
CLa = 4.4
CLde = .35
deref = 5. * d2r
CL_as1 = CLref + CLa * as1
x = 1. * CL_as1
CDde = 0.01
CDa = 0.35
CDaa = 0.1
CDref = 0.05
CD_2 = CDaa / 2. / CLa ** 2.
CD_1 = (CDa - 2. * CD_2 * CLref * CLa) / CLa
CD_0 = CDref - CD_1 * CLref - CD_2 * CLref ** 2.
zbp = -1.
Cmref = -zbp * CDref / c_w
Cma = -.68
Cmde = -.92
m = -0.4
CYb = -.56
CYdr = .155

def CL(a,de):
	if a < -5. * d2r:
		cl = -CL(-a-10.*d2r,-de+7.*d2r)
	elif a <= as1:
		cl = CLref + CLa * a + CLde * (de - deref)
	elif a < as2:
		mat = np.array([[   as2**3., as2**2., as2, 1.],
						[   as1**3., as1**2., as1, 1.],
						[3.*as2**2.,  2.*as2,  1., 0.],
						[3.*as1**2.,  2.*as1,  1., 0.]])
		b = np.array([x*np.cos(as2), CLref + CLa * as1, -x*np.sin(as2), CLa])
		A, B, C, D = np.linalg.solve(mat,b)
		cl = A * a ** 3. + B * a ** 2. + C * a + D + ((as2 - a) / (as2 - as1)) ** 2. * CLde * (de - deref)
	else:
		cl = x * np.cos(a)
	return cl

def CD(a,de,cl):
	if a < -5. * d2r:
		cd = CD(-a-10.*d2r,de,CL(-a-10.*d2r,de))
	elif a <= as1:
		cd = CD_0 + CD_1 * cl + CD_2 * cl ** 2. + CDde * (de - deref)
	elif a < as2:
		mat = np.array([[   as2**3., as2**2., as2, 1.],
						[   as1**3., as1**2., as1, 1.],
						[3.*as2**2.,  2.*as2,  1., 0.],
						[3.*as1**2.,  2.*as1,  1., 0.]])
		CLk = CL(as1,de)
		b = np.array([2. * np.sin(as2), CD_0 + CD_1 * CLk + CD_2 * CLk ** 2., 2. * np.cos(as2), CD_1 * CLa + 2. * CD_2 * CLk * CLa])
		A, B, C, D = np.linalg.solve(mat,b)
		cd = A * a ** 3. + B * a ** 2. + C * a + D + ((as2 - a) / (as2 - as1)) ** 2. * CDde * (de - deref)
	else:
		cd = 2. * np.sin(a)
	return cd

def CM(a,de,cl,cd):
	if a < -5.*d2r:
		de1 = -de+7.*d2r
		de2 = -de+0.*d2r
		cm = -CM(-a-10*d2r,-de,CL(-a-10.*d2r,de1),-CD(-a-10.*d2r,de,CL(-a-10.*d2r,de))) + .3
	elif a <= as1:
		cm = Cmref + Cma / CLa * (cl * np.cos(a) - CLref + cd * np.sin(a)) + Cmde * (de - deref)
	elif a < as2:
		asmin = 20. * d2r
		mat = np.array([[     as1**5.,      as1**4.,      as1**3.,   as1**2.,   as1, 1.],
						[  5.*as1**4.,   4.*as1**3.,   3.*as1**2.,    2.*as1,    1., 0.],
						[     as2**5.,      as2**4.,      as2**3.,   as2**2.,   as2, 1.],
						[  5.*as2**4.,   4.*as2**3.,   3.*as2**2.,    2.*as2,    1., 0.],
						[   asmin**5.,    asmin**4.,    asmin**3., asmin**2., asmin, 1.],
						[5.*asmin**4., 4.*asmin**3., 3.*asmin**2.,  2.*asmin,    1., 0.]])
		CLk = CL(as1,de)
		CDk = CD(as1,de,CLk)
		b = np.array([Cmref + Cma / CLa * (CLk * np.cos(as1) - CLref + CDk * np.sin(as1)) + Cmde * (de - deref),
					  Cma/CLa*(-CLk*np.sin(as1) + CLa*np.cos(as1) + CDk*np.cos(as1) + (CD_1*CLa + 2.*CD_2*CLk*CLa)*np.sin(as1)),
					  m * np.sin(as2),
					  m * np.cos(as2),
					  -0.5 + .5**2.*Cmde*(de-deref),
					  0.])
		A, B, C, D, E, F = np.linalg.solve(mat,b)
		cm = A*a**5. + B*a**4. + C*a**3. + D*a**2. + E*a + F #+ ((as2 - a) / (as2 - as1)) ** 2. * Cmde * (de - deref)
	else:
		cm = m*np.sin(a)
	return cm

def CS_func(b, dr):
	
	d2r = np.pi / 180.
	bs1 = 10. * d2r
	bs2 = 30. * d2r
	x = CYb * bs1
	
	if b < 0. * d2r:
		cs = -CS_func(-b-0.*d2r,-dr+0.*d2r)
	elif b <= bs1:
		cs = CYb * b + CYdr * dr
	elif b < bs2:
		mat = np.array([[   bs2**3., bs2**2., bs2, 1.],
						[   bs1**3., bs1**2., bs1, 1.],
						[3.*bs2**2.,  2.*bs2,  1., 0.],
						[3.*bs1**2.,  2.*bs1,  1., 0.]])
		bmat = np.array([x*np.cos(bs2), CYb * bs1, -x*np.sin(bs2), CYb])
		A, B, C, D = np.linalg.solve(mat,bmat)
		cs = A * b ** 3. + B * b ** 2. + C * b + D + ((bs2 - b) / (bs2 - bs1)) ** 2. * CYdr * dr
	else:
		cs = x * np.cos(b)
	return cs


res_a = 10000
res_de = 7
CLp = np.zeros((res_a,res_de))
CDp = np.zeros((res_a,res_de))
CMp = np.zeros((res_a,res_de))
alpha = np.zeros(res_a)
CSp = np.zeros((res_a,res_de))

for j in range(res_de):
	de = float(j) / float(res_de-1) * 30. - 15.
	de *= d2r
	for i in range(res_a):
		alpha[i] = float(i) / float(res_a-1) * 180. - 90.
		alpha[i] *= d2r
		CLp[i,j] = CL(alpha[i],de)
		CDp[i,j] = CD(alpha[i],de,CLp[i,j])
		CMp[i,j] = CM(alpha[i],de,CLp[i,j],CDp[i,j])
		CSp[i,j] = CS_func(alpha[i],de)

import matplotlib.pyplot as plt

plt.figure()
plt.plot(alpha/d2r,CLp,'b')

# plt.figure()
plt.plot(alpha/d2r,CDp,'r')

# plt.figure()
plt.plot(alpha/d2r,CMp,'g')
plt.axis([-90.,90.,-1.5,2.1])
plt.grid()

plt.figure()
plt.plot(alpha/d2r,CSp)
plt.grid()

plt.show()

