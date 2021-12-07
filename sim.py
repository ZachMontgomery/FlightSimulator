# Zach Montgomery
# A01621681
# simulator module

import numpy as np
import quat
import json
import stdatmos
import scipy.linalg

class vehicle:
	
	def text(self, word):
		l = int((100 - len(word)) / 2)
		print('\n'+'='*100)
		print('{}{}'.format(l*' ',word))
		print('='*100,'\n')
		
	
	# initialization routine
	def __init__(self, filename):
		
		self.text('Initializing Aircraft')
		
		# import json file as dictionary
		#*******************************************************************
		data = json.load(open(filename))
		
		# output dictionary entries into variables
		#*******************************************************************
		self.Sw       = data['aircraft']['wing_area']
		self.bw       = data['aircraft']['wing_span']
		self.c_w      = self.Sw / self.bw
		self.zbp      = data['aircraft']['thrust']['offset']
		self.xbp, self.ybp = 0., 0.
		self.T0       = data['aircraft']['thrust']['T0']
		self.T1       = data['aircraft']['thrust']['T1']
		self.T2       = data['aircraft']['thrust']['T2']
		self.a        = data['aircraft']['thrust']['a']
		
		self.V0       = data['initial']['airspeed']
		self.Altitude = data['initial']['altitude']
		self.W       = data['initial']['weight']
		self.climb    = data['initial']['climb'] *np.pi/180.
		self.bank     = data['initial']['bank'] *np.pi/180.
		self.psi0     = data['initial']['heading'] *np.pi/180.
		self.Phi1     = data['initial']['latitude'] *np.pi/180.
		self.Psi1     = data['initial']['longitude'] *np.pi/180.
		
		self.V_ref    = data['reference']['airspeed']
		self.rho_ref  = data['reference']['density']
		self.de_ref   = data['reference']['elevator'] * np.pi / 180.
		self.L_ref    = data['reference']['lift']
		self.CD_ref   = data['reference']['CD']
		self.Ixx      = data['reference']['Ixx']
		self.Iyy      = data['reference']["Iyy"]
		self.Izz      = data['reference']["Izz"]
		self.Ixy      = data['reference']["Ixy"]
		self.Ixz      = data['reference']["Ixz"]
		self.Iyz      = data['reference']["Iyz"]
		self.hxb      = data['reference']['hx']
		self.hyb      = data['reference']['hy']
		self.hzb      = data['reference']['hz']
		self.CLa      = data['reference']["CL,a"]
		self.CDa      = data['reference']["CD,a"]
		self.CDaa     = data['reference']['CD,a,a']
		self.Cma      = data['reference']["Cm,a"]
		self.CYb      = data['reference']["CY,b"]
		self.Clb      = data['reference']["Cl,b"]
		self.Cnb      = data['reference']["Cn,b"]
		self.CLq_     = data['reference']["CL,q"]
		self.CDq_     = data['reference']["CD,q"]
		self.Cmq_     = data['reference']["Cm,q"]
		self.CYp_     = data['reference']["CY,p"]
		self.Clp_     = data['reference']["Cl,p"]
		self.Cnp_     = data['reference']["Cn,p"]
		self.CYr_     = data['reference']["CY,r"]
		self.Clr_     = data['reference']["Cl,r"]
		self.Cnr_     = data['reference']["Cn,r"]
		self.CLde     = data['reference']["CL,de"]
		self.CDde     = data['reference']["CD,de"]
		self.Cmde     = data['reference']["Cm,de"]
		self.CYda     = data['reference']["CY,da"]
		self.Clda     = data['reference']["Cl,da"]
		self.Cnda     = data['reference']["Cn,da"]
		self.CYdr     = data['reference']["CY,dr"]
		self.Cldr     = data['reference']["Cl,dr"]
		self.Cndr     = data['reference']["Cn,dr"]
		self.CD3      = data['reference']['CD3']
		
		# set simulator variables
		#*******************************************************************
		self.t = 0.
		self.dt = .0166666666666667
		self.g = 32.2
		self.constant_density = True
		self.auto_pilot = False
		
		# initialize state vector
		#*******************************************************************
		#self.states = np.zeros(13)
		#self.states[0] = self.V0
		#g = 32.2
		#self.states[3] = -g / self.V0 * np.sin(self.theta0) * np.sin(self.phi0) / np.cos(self.phi0)
		#self.states[4] = g / self.V0 * np.sin(self.phi0) ** 2. * np.cos(self.theta0) / np.cos(self.phi0)
		#self.states[5] = g / self.V0 * np.sin(self.phi0) * np.cos(self.theta0)
		#self.states[8] = -self.zf0
		#self.states[9:13] = quat.Euler2Quat([self.phi0,self.theta0,self.psi0])
		
		# inertia tensor
		#*******************************************************************
		self.I = np.zeros((3,3))
		self.I[0,0] = self.Ixx
		self.I[1,1] = self.Iyy
		self.I[2,2] = self.Izz
		self.I[0,1] = -self.Ixy
		self.I[1,0] = -self.Ixy
		self.I[0,2] = -self.Ixz
		self.I[2,0] = -self.Ixz
		self.I[1,2] = -self.Iyz
		self.I[2,1] = -self.Iyz
		self.I_inv = np.linalg.inv(self.I)
		
		# angular momentum tensor
		#*******************************************************************
		self.h = np.zeros((3,3))
		self.h[0,1] = -self.hzb
		self.h[0,2] = self.hyb
		self.h[1,0] = self.hzb
		self.h[1,2] = -self.hxb
		self.h[2,0] = -self.hyb
		self.h[2,1] = self.hxb
	
	# trim 1124
	def trim(self, flag, controller_flag):
		
		
		def calc_R(self,x):
			tau,alpha,beta,de,da,dr=x
			
			const = .5 * self.rho * self.V0 ** 2. * self.Sw
			
			f1 = self.T + const * (self.CL * self.Salpha - self.CS * self.Sbeta - self.CD * self.u / self.V0) - self.W * self.Stheta + self.W / self.g * (self.r * self.v - self.q * self.w)
			
			f2 = const * (self.CS * self.Cbeta - self.CD * self.v / self.V0) + self.W * self.Sphi * self.Ctheta + self.W / self.g * (self.p * self.w - self.r * self.u)
			
			f3 = -const * (self.CL * self.Calpha + self.CD * self.w / self.V0) + self.W * self.Cphi * self.Ctheta + self.W / self.g * (self.q * self.u - self.p * self.w)
			
			f4 = const * self.bw * self.Cl - self.hzb * self.q + self.hyb * self.r + (self.Iyy - self.Izz) * self.q * self.r + self.Iyz * (self.q ** 2. - self.r ** 2.) + self.Ixz * self.p * self.q - self.Ixy * self.p * self.r
			
			f5 = self.zbp * self.T + const * self.c_w * self.Cm + self.hzb * self.p - self.hxb * self.r + (self.Izz - self.Ixx) * self.p * self.r + self.Ixz * (self.r**2. - self.p**2.) + self.Ixy * self.q * self.r - self.Iyz * self.p * self.q
			
			f6 = const * self.bw * self.Cn - self.hyb * self.p + self.hxb * self.q + (self.Ixx - self.Iyy) * self.p * self.q + self.Ixy * (self.p**2. - self.q**2.) + self.Iyz * self.p * self.r - self.Ixz * self.q * self.r
			
			return np.array([f1,f2,f3,f4,f5,f6])
		
		def calc_jacobian(self,x):
			tau,alpha,beta,de,da,dr=x
			
			CLa = self.CLa + self.CLq_ * self.q_a
			CLb = self.CLq_ * self.q_b
			
			CSa = self.CYp_ * self.p_a + self.CYr_ * self.r_a
			CSb = self.CYb + self.CYp_ * self.p_b + self.CYr_ * self.r_b
			
			CDa = self.CD1 * CLa + 2. * self.CD2 * self.CL * CLa + 2. * self.CD3 * self.CS * CSa + self.CDq_ * self.q_a
			CDb = self.CD1 * CLb + 2. * self.CD2 * self.CL * CLb + 2. * self.CD3 * self.CS * CSb + self.CDq_ * self.q_b
			CDde = self.CD1 * self.CLde + 2. * self.CD2 * self.CL * self.CLde + self.CDde
			CDda = 2. * self.CD3 * self.CS * self.CYda
			CDdr = 2. * self.CD3 * self.CS * self.CYdr
			
			Cla = self.Clp_ * self.p_a + self.Clr_ / self.CL_ref * (self.CL * self.r_a + CLa * self.r_)
			Clb = self.Clb + self.Clp_ * self.p_b + self.Clr_ / self.CL_ref * (self.CL * self.r_b + CLb * self.r_)
			
			Cma = self.Cma / self.CLa * (self.CL * self.ua / self.V0 + CLa * self.u / self.V0 + self.CD * self.wa / self.V0 + CDa * self.w / self.V0) + self.Cmq_ * self.q_a
			Cmb = self.Cma / self.CLa * (self.CL * self.ub / self.V0 + CLb * self.u / self.V0 + self.CD * self.wb / self.V0 + CDb * self.w / self.V0) + self.Cmq_ * self.q_b
			
			Cna = self.Cnb / self.CYb / self.V0 * (self.CS * self.ua + CSa * self.u - self.CD * self.va - CDa * self.v) + self.Cnp_ / self.CL_ref * (self.CL * self.p_a + CLa * self.p_) + self.Cnr_ * self.r_a
			Cnb = self.Cnb / self.CYb / self.V0 * (self.CS * self.ub + CSb * self.u - self.CD * self.vb - CDb * self.v) + self.Cnp_ / self.CL_ref * (self.CL * self.p_b + CLb * self.p_) + self.Cnr_ * self.r_b
			
			k = .5 * self.rho * self.V0 ** 2. * self.Sw
			
			J = np.zeros((6,6))
			
			#tau
			J[0,0] = (self.rho / self.rho_sea) ** self.a * (self.T0 + self.T1 * self.V0 + self.T2 * self.V0 ** 2.)
			J[4,0] = (self.rho / self.rho_sea) ** self.a * (self.T0 + self.T1 * self.V0 + self.T2 * self.V0 ** 2.) * self.zbp
			
			#alpha
			J[0,1] = k * (self.CL * self.Calpha + self.Salpha * CLa - CSa * self.Sbeta - self.CD * self.ua / self.V0 - CDa * self.u / self.V0) - self.W * self.Ctheta * self.thetaa + self.W / self.g * (self.r * self.va + self.ra * self.v - self.q * self.wa - self.qa * self.w)
			
			J[1,1] = k * (CSa * self.Cbeta - self.CD * self.va / self.V0 - CDa * self.v / self.V0) - self.W * self.Sphi * self.Stheta * self.thetaa + self.W / self.g * (self.p * self.wa + self.pa * self.w - self.r * self.ua - self.ra * self.u)
			
			J[2,1] = k * (self.CL * self.Salpha - CLa * self.Calpha - self.CD * self.wa / self.V0 - CDa * self.w / self.V0) - self.W * self.Cphi * self.Stheta * self.thetaa + self.W / self.g * (self.q * self.ua + self.qa * self.u - self.p * self.wa - self.pa * self.w)
			
			J[3,1] = k * self.bw * Cla - self.hzb * self.qa + self.hyb * self.ra + (self.Iyy - self.Izz) * (self.q * self.ra + self.qa * self.r) + 2. * self.Iyz * (self.q * self.qa - self.r * self.ra) + self.Ixz * (self.p * self.qa + self.pa * self.q) - self.Ixy * (self.p * self.ra + self.pa * self.r)
			
			J[4,1] = k * self.c_w * Cma + self.hzb * self.pa - self.hxb * self.ra + (self.Izz - self.Ixx) * (self.p * self.ra + self.pa * self.r) + 2. * self.Ixz * (self.r * self.ra - self.p * self.pa) + self.Ixy * (self.q * self.ra + self.qa * self.r) - self.Iyz * (self.p * self.qa + self.pa * self.q)
			
			J[5,1] = k * self.bw * Cna - self.hyb * self.pa + self.hxb * self.qa + (self.Ixx - self.Iyy) * (self.p * self.qa + self.pa * self.q) + 2. * self.Ixy * (self.p * self.pa - self.q * self.qa) + self.Iyz * (self.p * self.ra + self.pa * self.r) - self.Ixz * (self.q * self.ra + self.qa * self.r)
			
			#beta
			J[0,2] = k * (CLb * self.Salpha - self.CS * self.Cbeta - CSa * self.Sbeta - self.CD * self.ub / self.V0 - CDb * self.u / self.V0) - self.W * self.Ctheta * self.thetab + self.W / self.g * (self.r * self.vb + self.rb * self.v - self.q * self.wb - self.qb * self.w)
			
			J[1,2] = k * (-self.CS * self.Sbeta + CSb * self.Cbeta - self.CD * self.vb / self.V0 - CDb * self.v / self.V0) - self.W * self.Sphi * self.Stheta * self.thetab + self.W / self.g * (self.p * self.wb + self.pb * self.w - self.r * self.ub - self.rb * self.u)
			
			J[2,2] = k * (-CLb * self.Calpha - self.CD * self.wb / self.V0 - CDb * self.w / self.V0) - self.W * self.Cphi * self.Stheta * self.thetab + self.W / self.g * (self.q * self.ub + self.qb * self.u - self.p * self.wb - self.pb * self.w)
			
			J[3,2] = k * self.bw * Clb - self.hzb * self.qb + self.hyb * self.rb + (self.Iyy - self.Izz) * (self.q * self.rb + self.qb * self.r) + 2. * self.Iyz * (self.q * self.qb - self.r * self.rb) + self.Ixz * (self.p * self.qb + self.pb * self.q) - self.Ixy * (self.p * self.rb + self.pb * self.r)
			
			J[4,2] = k * self.c_w * Cmb + self.hzb * self.pb - self.hxb * self.rb + (self.Izz - self.Ixx) * (self.p * self.rb + self.pb * self.r) + 2. * self.Ixz * (self.r * self.rb - self.p * self.pb) + self.Ixy * (self.q * self.rb + self.qb * self.r) - self.Iyz * (self.p * self.qb + self.pb * self.q)
			
			J[5,2] = k * self.bw * Cnb - self.hyb * self.pb + self.hxb * self.qb + (self.Ixx - self.Iyy) * (self.p * self.qb + self.pb * self.q) + 2. * self.Ixy * (self.p * self.pb - self.q * self.qb) + self.Iyz * (self.p * self.rb + self.pb * self.r) - self.Ixz * (self.q * self.rb + self.qb * self.r)
			
			#de
			J[0,3] = k * (self.Salpha * self.CLde - self.u / self.V0 * CDde)
			J[1,3] = -k * self.v / self.V0 * CDde
			J[2,3] = k * (-self.Calpha * self.CLde - self.w / self.V0 * CDde)
			J[4,3] = k * self.c_w * (self.Cma / self.CLa * (self.u / self.V0 * self.CLde + self.w / self.V0 * CDde) + self.Cmde)
			
			#da
			J[0,4] = k * (-self.CYda * self.Sbeta - CDda * self.u / self.V0)
			J[1,4] = k * (self.CYda * self.Cbeta - CDda * self.v / self.V0)
			J[2,4] = -k * CDda * self.w / self.V0
			J[3,4] = k * self.bw * self.Clda
			J[4,4] = k * self.c_w * self.Cma / self.CLa * self.w / self.V0 * CDda
			J[5,4] = k * self.bw * (self.Cnb / self.CYb / self.V0 * ( self.u * self.CYda - self.v * CDda) + self.Cnda)
			
			#dr
			J[0,5] = k * (-self.CYdr * self.Sbeta - CDdr * self.u / self.V0)
			J[1,5] = k * (self.CYdr * self.Cbeta - CDdr * self.v / self.V0)
			J[2,5] = -k * CDdr * self.w / self.V0
			J[3,5] = k * self.bw * self.Cldr
			J[4,5] = k * self.c_w * self.Cma / self.CLa * self.w / self.V0 * CDdr
			J[5,5] = k * self.bw * (self.Cnb / self.CYb / self.V0 * ( self.u * self.CYdr - self.v * CDdr) + self.Cndr)
			
			return J
		
		self.text('Begin Trim Solver')
		
		# set initial guess
		tau = 0.0
		alpha = 0.
		beta = 0.
		de = self.de_ref
		da = 0.
		dr = 0.
		xn = np.array([tau,alpha,beta,de,da,dr])
		
		# prepare and start the iterative loop
		tol = 1.e-14
		maxit = 200
		ea = 1.
		count = 0
		print('Iteration		Error')
		while ea >= tol:
			# save the new guess as the old guess
			xo = np.copy(xn)
			
			# calculate theta and other unknown constants
			self.calc_constants(xo)
			
			# calculate residual vector
			R = calc_R(self, xo)
			
			# calculate jacobian
			J = calc_jacobian(self, xo)
			
			# update guess
			dx = np.linalg.solve(J,-R)
			relax = 1.
			xn = xo + relax * dx
			
			# calculate error
			ea = abs( (xn[1] - xo[1]) / xn[1] )
			
			# display iteration
			count += 1
			print('{:>9d}		{:8.2E}'.format(count, ea))
			if count >= maxit: break
		
		# update constants
		self.calc_constants(xn)
		
		if flag:
			# set trim state and controls
			self.states = np.zeros(13)
			
			self.states[0] = self.V0 * self.Calpha * self.Cbeta / np.sqrt(1. - self.Salpha**2. * self.Sbeta**2.)
			self.states[1] = self.V0 * self.Calpha * self.Sbeta / np.sqrt(1. - self.Salpha**2. * self.Sbeta**2.)
			self.states[2] = self.V0 * self.Salpha * self.Cbeta / np.sqrt(1. - self.Salpha**2. * self.Sbeta**2.)
			
			self.states[3] = self.p
			self.states[4] = self.q
			self.states[5] = self.r
			
			self.states[8] = -self.Altitude
			
			self.states[9:13] = quat.Euler2Quat([self.bank,self.theta,self.psi0])
		
		if controller_flag:
			states = np.zeros(12)
			
			states[0] = self.V0 * self.Calpha * self.Cbeta / np.sqrt(1. - self.Salpha**2. * self.Sbeta**2.)
			states[1] = self.V0 * self.Calpha * self.Sbeta / np.sqrt(1. - self.Salpha**2. * self.Sbeta**2.)
			states[2] = self.V0 * self.Salpha * self.Cbeta / np.sqrt(1. - self.Salpha**2. * self.Sbeta**2.)
			
			states[3] = self.p
			states[4] = self.q
			states[5] = self.r
			
			states[8] = -self.Altitude
			
			states[9] = self.bank
			states[10] = self.theta
			states[11] = self.psi0
			
			control = np.array([xn[0],xn[3],xn[4],xn[5]])
			
			return states, control
		else:
			self.control = np.array([xn[0],xn[3],xn[4],xn[5]])
		
		#display trim results
		r2d = 180. / np.pi
		print('\n\nTrim Results (deg except Tau)\n')
		print('Tau   = {:16.12f}'.format(xn[0]))
		print('alpha = {:16.12f}'.format(xn[1]*r2d))
		print('beta  = {:16.12f}'.format(xn[2]*r2d))
		print('de    = {:16.12f}'.format(xn[3]*r2d))
		print('da    = {:16.12f}'.format(xn[4]*r2d))
		print('dr    = {:16.12f}'.format(xn[5]*r2d))
		print('theta = {:16.12f}\n\n'.format(self.theta*r2d))
	
	def calc_theta(self, climb):
		self.Sg = np.sin(climb)
		
		denom = 1. - self.Salpha**2. * self.Sbeta**2.
		knum = self.Sphi * self.Calpha * self.Sbeta + self.Cphi * self.Salpha * self.Cbeta
		
		denoma = -2. * self.Salpha * self.Sbeta**2.
		denomb = -2. * self.Salpha**2. * self.Sbeta
		knuma = self.Cphi * self.Cbeta * self.Calpha - self.Sphi * self.Sbeta * self.Salpha
		knumb = self.Sphi * self.Calpha * self.Cbeta - self.Cphi * self.Salpha * self.Sbeta
		
		k = knum / np.sqrt(denom)
		# print('k is {}'.format(k))
		ka = -.5 * knum * denom ** -1.5 * denoma + knuma / np.sqrt(denom)
		kb = -.5 * knum * denom ** -1.5 * denomb + knumb / np.sqrt(denom)
		
		m = self.Calpha * self.Cbeta / np.sqrt(denom)
		# print('m is {}'.format(m))
		ma = -.5 * self.Calpha * self.Cbeta * denom ** -1.5 * denoma - self.Salpha * self.Cbeta / np.sqrt(denom)
		mb = -.5 * self.Calpha * self.Cbeta * denom ** -1.5 * denomb - self.Calpha * self.Sbeta / np.sqrt(denom)
		
		# test = k**4.+k**2.*m**2.-k**2.*self.Sg**2.
		
		xnum1 = self.Sg * m + k * np.sqrt(k**2. + m**2. - self.Sg**2.)
		# if test > 0.:
		xnum1a = self.Sg * ma + .5 * k * (k**2.+m**2.-self.Sg**2.)**-.5*(2.*k*ka+2.*m*ma-self.Sg**2.) + ka * (k**2.+m**2.-self.Sg**2.)**.5
		xnum1b = self.Sg * mb + .5 * k * (k**2.+m**2.-self.Sg**2.)**-.5*(2.*k*kb+2.*m*mb-self.Sg**2.) + kb * (k**2.+m**2.-self.Sg**2.)**.5
		# else:
			# print(test)
			# input()
		xnum2 = self.Sg * m - k * np.sqrt(k**2. + m**2. - self.Sg**2.)
		# if test > 0.:
		xnum2a = self.Sg * ma - .5 * k * (k**2.+m**2.-self.Sg**2.)**-.5*(2.*k*ka+2.*m*ma-self.Sg**2.) + ka * (k**2.+m**2.-self.Sg**2.)**.5
		xnum2b = self.Sg * mb - .5 * k * (k**2.+m**2.-self.Sg**2.)**-.5*(2.*k*kb+2.*m*mb-self.Sg**2.) + kb * (k**2.+m**2.-self.Sg**2.)**.5
		# else:
			# print(test)
			# input()
		
		xdenom = k**2. + m**2.
		xdenoma = 2.*k*ka + 2.*m*ma
		xdenomb = 2.*k*kb + 2.*m*mb
		
		x1 = xnum1 / xdenom
		x1a = -xnum1*xdenom**-2.*xdenoma + xnum1a / xdenom
		x1b = -xnum1*xdenom**-2.*xdenomb + xnum1b / xdenom
		x2 = xnum2 / xdenom
		x2a = -xnum2*xdenom**-2.*xdenoma + xnum2a / xdenom
		x2b = -xnum2*xdenom**-2.*xdenomb + xnum2b / xdenom
		
		t1 = np.arcsin(x1)
		t1a = (1.-x1**2.)**-.5 * x1a
		t1b = (1.-x1**2.)**-.5 * x1b
		t2 = np.arcsin(x2)
		t2a = (1.-x2**2.)**-.5 * x2a
		t2b = (1.-x2**2.)**-.5 * x2b
		
		e1 = abs(np.sin(t1) * m - np.cos(t1) * k - self.Sg)
		e2 = abs(np.sin(t2) * m - np.cos(t2) * k - self.Sg)
		
		if e1 < e2:
			t = t1
			ta = t1a
			tb = t1b
			if e1 > 1.e-16: print('Not a good theta value. Error is {}'.format(e1))
		else:
			t = t2
			ta = t2a
			tb = t2b
			if e2 > 1.e-16: print('Not a good theta value. Error is {}'.format(e2))
		
		return t, ta, tb
	
	def calc_constants(self, x):
		tau,alpha,beta,de,da,dr=x
		
		# calculate trig values for the angles
		self.Sphi = np.sin(self.bank)
		self.Cphi = np.cos(self.bank)
		self.Salpha = np.sin(alpha)
		self.Calpha = np.cos(alpha)
		self.Sbeta = np.sin(beta)
		self.Cbeta = np.cos(beta)
		
		
		
		# calculate theta and its derivatives and trig values
		self.theta, self.thetaa, self.thetab = self.calc_theta(self.climb)
		self.Stheta = np.sin(self.theta)
		self.Ctheta = np.cos(self.theta)
		
		
		
		# calculate uvw and their derivatives
		denom = 1. - self.Salpha**2. * self.Sbeta**2.
		denoma = -2. * self.Salpha * self.Calpha * self.Sbeta ** 2.
		denomb = -2. * self.Salpha ** 2. * self.Sbeta * self.Cbeta
		
		self.u = self.V0 * self.Calpha * self.Cbeta / np.sqrt(denom)
		self.v = self.V0 * self.Calpha * self.Sbeta / np.sqrt(denom)
		self.w = self.V0 * self.Salpha * self.Cbeta / np.sqrt(denom)
		
		self.ua = -.5 * self.V0 * self.Calpha * self.Cbeta * denom ** -1.5 * denoma - self.V0 * self.Salpha * self.Cbeta / np.sqrt(denom)
		self.ub = -.5 * self.V0 * self.Calpha * self.Cbeta * denom ** -1.5 * denomb - self.V0 * self.Calpha * self.Sbeta / np.sqrt(denom)
		
		self.va = -.5 * self.V0 * self.Calpha * self.Sbeta * denom ** -1.5 * denoma - self.V0 * self.Salpha * self.Sbeta / np.sqrt(denom)
		self.vb = -.5 * self.V0 * self.Calpha * self.Sbeta * denom ** -1.5 * denomb + self.V0 * self.Calpha * self.Cbeta / np.sqrt(denom)
		
		self.wa = -.5 * self.V0 * self.Salpha * self.Cbeta * denom ** -1.5 * denoma + self.V0 * self.Calpha * self.Cbeta / np.sqrt(denom)
		self.wb = -.5 * self.V0 * self.Salpha * self.Cbeta * denom ** -1.5 * denomb - self.V0 * self.Salpha * self.Sbeta / np.sqrt(denom)
		
		
		
		# calculate omega and its derivatives
		omega = self.g * self.Sphi * self.Ctheta * np.sqrt(1.-self.Salpha**2.*self.Sbeta**2.) / self.V0 / (self.Salpha * self.Cbeta * self.Stheta + self.Calpha * self.Cbeta * self.Cphi * self.Ctheta)
		
		denom = self.Stheta * self.w + self.Cphi * self.Ctheta * self.u
		denoma = self.Stheta * self.wa + self.w * self.Ctheta * self.thetaa + self.Cphi * (self.Ctheta * self.ua - self.u * self.Stheta * self.thetaa)
		denomb = self.Stheta * self.wb + self.w * self.Ctheta * self.thetab + self.Cphi * (self.Ctheta * self.ub - self.u * self.Stheta * self.thetab)
		
		omegaa = -1. * self.g * self.Sphi * self.Ctheta * denom ** -2. * denoma - self.g * self.Sphi * self.Stheta * self.thetaa / denom
		omegab = -1. * self.g * self.Sphi * self.Ctheta * denom ** -2. * denomb - self.g * self.Sphi * self.Stheta * self.thetab / denom
		
		
		
		# calculate pqr and their derivatives
		self.p = -self.Stheta * omega
		self.q = self.Sphi * self.Ctheta * omega
		self.r = self.Cphi * self.Ctheta * omega
		
		self.p_ = self.p * self.bw / 2. / self.V0
		self.q_ = self.q * self.c_w / 2. / self.V0
		self.r_ = self.r * self.bw / 2. / self.V0
		
		self.pa = -self.Stheta * omegaa - omega * self.Ctheta * self.thetaa
		self.pb = -self.Stheta * omegab - omega * self.Ctheta * self.thetab
		self.qa = self.Sphi * (self.Ctheta * omegaa - omega * self.Stheta * self.thetaa)
		self.qb = self.Sphi * (self.Ctheta * omegab - omega * self.Stheta * self.thetab)
		self.ra = self.Cphi * (self.Ctheta * omegaa - omega * self.Stheta * self.thetaa)
		self.rb = self.Cphi * (self.Ctheta * omegab - omega * self.Stheta * self.thetab)
		
		self.p_a = self.pa * self.bw / 2. / self.V0
		self.p_b = self.pb * self.bw / 2. / self.V0
		self.q_a = self.qa * self.c_w / 2. / self.V0
		self.q_b = self.qb * self.c_w / 2. / self.V0
		self.r_a = self.ra * self.bw / 2. / self.V0
		self.r_b = self.rb * self.bw / 2. / self.V0
		
		
		
		# calculate unknown constants
		self.CL_ref = self.L_ref * 2. / self.rho_ref / self.V_ref ** 2. / self.Sw
		self.Cm_ref = -self.zbp * self.CD_ref / self.c_w
		self.CD2 = self.CDaa / 2. / self.CLa **2.
		self.CD1 = (self.CDa - 2. * self.CD2 * self.CL_ref * self.CLa) / self.CLa
		self.CD0 = self.CD_ref - self.CD1 * self.CL_ref - self.CD2 * self.CL_ref ** 2.
		
		
		
		# calculate densities
		self.rho = stdatmos.statee(self.Altitude)[3]
		self.rho_sea = stdatmos.statee(0.)[3]
		if self.constant_density:
			self.rho = self.rho_ref
			self.rho_sea = self.rho_ref
		
		
		
		# other values
		self.CL = self.CL_ref + self.CLa * alpha + self.CLq_ * self.q_ + self.CLde * (de - self.de_ref)
		self.CS = self.CYb * beta + self.CYp_ * self.p_ + self.CYr_ * self.r_ + self.CYda * da + self.CYdr * dr
		self.CD = self.CD0 + self.CD1 * self.CL + self.CD2 * self.CL**2. + self.CD3 * self.CS**2. + self.CDq_ * self.q_ + self.CDde * (de - self.de_ref)
		self.Cl = self.Clb * beta + self.Clp_ * self.p_ + self.Clr_ / self.CL_ref * self.CL * self.r_ + self.Clda * da + self.Cldr * dr
		self.Cm = self.Cm_ref + self.Cma / self.CLa * (self.CL * self.u / self.V0 - self.CL_ref + self.CD * self.w / self.V0) + self.Cmq_ * self.q_ + self.Cmde * (de - self.de_ref)
		self.Cn = self.Cnb / self.CYb * (self.CS * self.u / self.V0 - self.CD * self.v / self.V0) + self.Cnp_ / self.CL_ref * self.CL * self.p_ + self.Cnr_ * self.r_ + self.Cnda * da + self.Cndr * dr
		self.T = tau * (self.rho / self.rho_sea) ** self.a * (self.T0 + self.T1 * self.V0 + self.T2 * self.V0**2.)
	
	def trim_hunsaker(self, V0, phi, climb):
		
		CT0 = (self.rho/self.rho_sea)**self.a*(self.T0+self.T1*V0+self.T2*V0**2.)*2./self.rho/V0**2./self.Sw
		
		tau, alpha, beta, d_de, da, dr = 0., 0., 0., 0., 0., 0.
		
		self.Sphi = np.sin(phi)
		self.Cphi = np.cos(phi)
		self.Salpha = np.sin(alpha)
		self.Calpha = np.cos(alpha)
		self.Sbeta = np.sin(beta)
		self.Cbeta = np.cos(beta)
		
		
		
		# calculate theta and its derivatives and trig values
		theta = self.calc_theta(climb)[0]
		
		
		xo = np.zeros(6)
		ea = np.ones(6)
		count = 0
		
		while max(ea) > 1.e-15:
			denom = np.sqrt(1.-np.sin(alpha)**2.*np.sin(beta)**2.)
			u = V0/denom*np.cos(alpha)*np.cos(beta)
			v = V0/denom*np.cos(alpha)*np.sin(beta)
			w = V0/denom*np.sin(alpha)*np.cos(beta)
			
			k = self.g*np.sin(phi)*np.cos(theta)/(np.sin(theta)*w+np.cos(phi)*np.cos(theta)*u)
			p = -k * np.sin(theta)
			q = k * np.sin(phi) * np.cos(theta)
			r = k * np.cos(phi) * np.cos(theta)
			
			p_ = self.bw / 2. / V0 * p
			q_ = self.c_w / 2. / V0 * q
			r_ = self.bw / 2. / V0 * r
			
			CL = self.CL_ref + self.CLa * alpha + self.CLq_ * q_ + self.CLde * d_de
			CS = self.CYb * beta + self.CYp_ * p_ + self.CYr_ * r_ + self.CYda * da + self.CYdr * dr
			CD = self.CD0 + self.CD1 * CL + self.CD2 * CL ** 2. + self.CD3 * CS ** 2. + self.CDq_ * q_ + self.CDde * d_de
			
			A = np.zeros((6,6))
			A[0,0] = CT0
			A[0,3] = -self.CDde * u / V0
			A[1,2] = self.CYb * np.cos(beta)
			A[1,4] = self.CYda * np.cos(beta)
			A[1,5] = self.CYdr * np.cos(beta)
			A[2,1] = -self.CLa * np.cos(alpha)
			A[2,3] = -self.CLde * np.cos(alpha)
			A[3,2] = self.bw * self.Clb
			A[3,4] = self.bw * self.Clda
			A[3,5] = self.bw * self.Cldr
			A[4,1] = self.zbp * CT0
			A[4,3] = self.c_w * self.Cmde
			A[5,4] = self.bw * self.Cnda
			A[5,5] = self.bw * self.Cndr
			
			b1 = np.zeros(6)
			b1[0] = -CL*np.sin(alpha)+CS*np.sin(beta)+(self.CD0+self.CD1*CL+self.CD2*CL**2.+self.CD3*CS**2.+self.CDq_*q_)*u/V0
			b1[1] = (-self.CYp_ * p_ - self.CYr_ * r_) * np.cos(beta) + CD * v / V0
			b1[2] = (self.CL_ref + self.CLq_ * q_) * np.cos(alpha) + CD * w / V0
			b1[3] = -self.bw * (self.Clp_ * p_ + self.Clr_ / self.CL_ref * CL * r_)
			b1[4] = -self.c_w * (self.Cm_ref + self.Cma / self.CLa * (CL * u / V0 - self.CL_ref + CD * w / V0) + self.Cmq_ * q_)
			b1[5] = -self.bw * (self.Cnb / self.CYb * (CS * u / V0 - CD * v / V0) + self.Cnp_ / self.CL_ref * CL * p_ + self.Cnr_ * r_)
			
			b2 = np.zeros(6)
			b2[0] = np.sin(theta) - (r * v - q * w) / self.g
			b2[1] = -np.sin(phi) * np.cos(theta) - (p * w - r * u) / self.g
			b2[2] = -np.cos(phi) * np.cos(theta) - (q * u - p * v) / self.g
			b2[3] = ( self.hzb*q-self.hyb*r-(self.Iyy-self.Izz)*q*r-self.Iyz*(q*q-r*r)-self.Ixz*p*q+self.Ixy*p*r)/self.W
			b2[4] = (-self.hzb*p+self.hxb*r-(self.Izz-self.Ixx)*p*r-self.Ixz*(r*r-p*p)-self.Ixy*q*r+self.Iyz*p*q)/self.W
			b2[5] = ( self.hyb*p-self.hxb*q-(self.Ixx-self.Iyy)*p*q-self.Ixy*(p*p-q*q)-self.Iyz*p*r+self.Ixz*q*r)/self.W
			
			xn = np.linalg.solve(A,b1+b2)
			
			for i in range(6):
				ea[i] = abs((xn[i]-xo[i])/xn[i])
			
			xo = np.copy(xn)
			count += 1
			
			print('Iteration {}   Error {}'.format(count, max(ea)))
		
		states = np.array([u,v,w,p,q,r,0.,0.,-500.,phi,theta,0.])
		control = np.array([xn[0],xn[3]+self.de_ref,xn[4],xn[5]])
		return states, control
	
	# controls
	def Control(self):
		#self.control = np.zeros(4)
		#self.control[2] = 1. * np.pi / 180.
		pass
	
	# lift coefficient
	def CL_func(self, a, de, q_):
		
		d2r = np.pi / 180.
		as1 = 10. * d2r
		as2 = 30. * d2r
		x = self.CL_ref + self.CLa * as1
		
		if a < -5. * d2r:
			cl = -self.CL_func(-a-10.*d2r,-de+7.*d2r,q_)
		elif a <= as1:
			cl = self.CL_ref + self.CLa * a + self.CLq_ * q_ + self.CLde * (de - self.de_ref)
		elif a < as2:
			mat = np.array([[   as2**3., as2**2., as2, 1.],
							[   as1**3., as1**2., as1, 1.],
							[3.*as2**2.,  2.*as2,  1., 0.],
							[3.*as1**2.,  2.*as1,  1., 0.]])
			b = np.array([x*np.cos(as2), self.CL_ref + self.CLa * as1, -x*np.sin(as2), self.CLa])
			A, B, C, D = np.linalg.solve(mat,b)
			cl = A * a ** 3. + B * a ** 2. + C * a + D + ((as2 - a) / (as2 - as1)) ** 2. * self.CLde * (de - self.de_ref)
		else:
			cl = x * np.cos(a)
		return cl
	
	# drag coefficient
	def CD_func(self, a,de,cl,q_,cs):
		
		d2r = np.pi / 180.
		as1 = 10. * d2r
		as2 = 30. * d2r
		
		if a < -5. * d2r:
			cd = self.CD_func(-a-10.*d2r,de,self.CL_func(-a-10.*d2r,de,q_),q_,cs)
		elif a <= as1:
			cd = self.CD0 + self.CD1 * cl + self.CD2 * cl**2. + self.CD3 * cs**2. + self.CDq_ * q_ + self.CDde * (de - self.de_ref)
		elif a < as2:
			mat = np.array([[   as2**3., as2**2., as2, 1.],
							[   as1**3., as1**2., as1, 1.],
							[3.*as2**2.,  2.*as2,  1., 0.],
							[3.*as1**2.,  2.*as1,  1., 0.]])
			CLk = self.CL_func(as1,de,q_)
			b = np.array([2. * np.sin(as2), self.CD0 + self.CD1 * CLk + self.CD2 * CLk ** 2., 2. * np.cos(as2), self.CD1 * self.CLa + 2. * self.CD2 * CLk * self.CLa])
			A, B, C, D = np.linalg.solve(mat,b)
			cd = A * a ** 3. + B * a ** 2. + C * a + D + ((as2 - a) / (as2 - as1)) ** 2. * self.CDde * (de - self.de_ref)
		else:
			cd = 2. * np.sin(a)
		return cd
	
	# moment coefficient
	def CM_func(self, a,de,cl,cd,q_,cs,u,w,V):
		
		d2r = np.pi / 180.
		as1 = 10. * d2r
		as2 = 30. * d2r
		m = -0.4
		
		if a < -5.*d2r:
			de1 = -de+7.*d2r
			de2 = -de+0.*d2r
			cm = -self.CM_func(-a-10*d2r,-de,self.CL_func(-a-10.*d2r,de1,q_),-self.CD_func(-a-10.*d2r,de,self.CL_func(-a-10.*d2r,de,q_),q_,cs),q_,cs,u,w,V) + .3
		elif a <= as1:
			cm = self.Cm_ref + self.Cma / self.CLa * (cl * u / V - self.CL_ref + cd * w / V) + self.Cmq_ * q_ + self.Cmde * (de - self.de_ref)
		elif a < as2:
			asmin = 20. * d2r
			mat = np.array([[     as1**5.,      as1**4.,      as1**3.,   as1**2.,   as1, 1.],
							[  5.*as1**4.,   4.*as1**3.,   3.*as1**2.,    2.*as1,    1., 0.],
							[     as2**5.,      as2**4.,      as2**3.,   as2**2.,   as2, 1.],
							[  5.*as2**4.,   4.*as2**3.,   3.*as2**2.,    2.*as2,    1., 0.],
							[   asmin**5.,    asmin**4.,    asmin**3., asmin**2., asmin, 1.],
							[5.*asmin**4., 4.*asmin**3., 3.*asmin**2.,  2.*asmin,    1., 0.]])
			CLk = self.CL_func(as1,de,q_)
			CDk = self.CD_func(as1,de,CLk,q_,cs)
			b = np.array([self.Cm_ref + self.Cma / self.CLa * (CLk * np.cos(as1) - self.CL_ref + CDk * np.sin(as1)) + self.Cmde * (de - self.de_ref),
						self.Cma/self.CLa*(-CLk*np.sin(as1) + self.CLa*np.cos(as1) + CDk*np.cos(as1) + (self.CD1*self.CLa + 2.*self.CD2*CLk*self.CLa)*np.sin(as1)),
						  m * np.sin(as2),
						  m * np.cos(as2),
						  -0.5 + .5**2.*self.Cmde*(de-self.de_ref),
						  0.])
			A, B, C, D, E, F = np.linalg.solve(mat,b)
			cm = A*a**5. + B*a**4. + C*a**3. + D*a**2. + E*a + F #+ ((as2 - a) / (as2 - as1)) ** 2. * Cmde * (de - deref)
		else:
			cm = m*np.sin(a)
		return cm
	
	# side coefficient
	def CS_func(self, b, da, dr, p_, r_):
		
		d2r = np.pi / 180.
		bs1 = 10. * d2r
		bs2 = 30. * d2r
		x = self.CYb * bs1
		
		if b < 0. * d2r:
			cs = -self.CS_func(-b-0.*d2r,da,-dr+0.*d2r,p_,r_)
		elif b <= bs1:
			cs = self.CYb * b + self.CYp_ * p_ + self.CYr_ * r_ + self.CYda * da + self.CYdr * dr
		elif b < bs2:
			mat = np.array([[   bs2**3., bs2**2., bs2, 1.],
							[   bs1**3., bs1**2., bs1, 1.],
							[3.*bs2**2.,  2.*bs2,  1., 0.],
							[3.*bs1**2.,  2.*bs1,  1., 0.]])
			bmat = np.array([x*np.cos(bs2), self.CYb * bs1, -x*np.sin(bs2), self.CYb])
			A, B, C, D = np.linalg.solve(mat,bmat)
			cs = A * b ** 3. + B * b ** 2. + C * b + D + ((bs2 - b) / (bs2 - bs1)) ** 2. * self.CYdr * dr
		else:
			cs = x * np.cos(b)
		return cs
	
	# aerodynamic forces and moments
	def aero(self,states):
		u,v,w,p,q,r,xf,yf,zf,e0,ex,ey,ez=states
		tau,de,da,dr=self.control
		
		if self.constant_density:
			self.rho = self.rho_ref
			self.rho_sea = self.rho_ref
		else:
			self.rho = stdatmos.statee(-zf)[3]
			self.rho_sea = stdatmos.statee(0.)[3]
		
		V = np.sqrt(u*u+v*v+w*w)
		const = 0.5 * self.rho * V ** 2. * self.Sw
		
		alpha = np.arctan2(w,u)
		beta = np.arctan2(v,u)
		
		p_ = self.bw / 2. / V * p
		q_ = self.c_w / 2. / V * q
		r_ = self.bw / 2. / V * r
		
		# CL = self.CL_ref + self.CLa * alpha + self.CLq_ * q_ + self.CLde * (de - self.de_ref)
		CL = self.CL_func(alpha,de,q_)
		# CS = self.CYb * beta + self.CYp_ * p_ + self.CYr_ * r_ + self.CYda * da + self.CYdr * dr
		CS = self.CS_func(beta,da,dr,p_,r_)
		# CD = self.CD0 + self.CD1 * CL + self.CD2 * CL**2. + self.CD3 * CS**2. + self.CDq_ * q_ + self.CDde * (de - self.de_ref)
		CD = self.CD_func(alpha,de,CL,q_,CS)
		Cl = self.Clb * beta + self.Clp_ * p_ + self.Clr_ / self.CL_ref * CL * r_ + self.Clda * da + self.Cldr * dr
		# Cm = self.Cm_ref + self.Cma / self.CLa * (CL * u / V - self.CL_ref + CD * w / V) + self.Cmq_ * q_ + self.Cmde * (de - self.de_ref)
		Cm = self.CM_func(alpha,de,CL,CD,q_,CS,u,w,V)
		Cn = self.Cnb / self.CYb * (CS * u / V - CD * v / V) + self.Cnp_ / self.CL_ref * CL * p_ + self.Cnr_ * r_ + self.Cnda * da + self.Cndr * dr
		T = tau * (self.rho / self.rho_sea) ** self.a * (self.T0 + self.T1 * V + self.T2 * V**2.)
		
		X = const * (CL * np.sin(alpha) - CS * np.sin(beta) - CD * u / V)
		Y = const * (CS * np.cos(beta) - CD * v / V)
		Z = const * (-CL * np.cos(alpha) - CD * w / V)
		
		l = const * self.bw * Cl
		m = const * self.c_w * Cm
		n = const * self.bw * Cn
		
		Txb = T
		Tyb = 0.
		Tzb = 0.
		
		
		
		# print(X,Y,Z,l,m,n)
		# print(CL,p_,q_,r_)
		# print(alpha,beta)
		# print(self.rho)
		# print(de,da,dr)
		# print(const)
		# input()
		return X,Y,Z,l,m,n,Txb,Tyb,Tzb
	
	# derivitive function used in the integrator function
	def f(self,states):
		u,v,w,p,q,r,xf,yf,zf,e0,ex,ey,ez=states
		
		X,Y,Z,l,m,n,Txb,Tyb,Tzb=self.aero(states)
		
		self.gForces = np.array([(X+Txb)/self.W,(Y+Tyb)/self.W,-(Z+Tzb)/self.W])
		
		u_dot = (2. * (ex * ez - ey * e0) + (X + Txb) / self.W) * self.g + r * v - q * w
		v_dot = (2. * (ey * ez + ex * e0) + (Y + Tyb) / self.W) * self.g + p * w - r * u
		w_dot = (ez * ez + e0 * e0 - ex * ex - ey * ey + (Z + Tzb) / self.W) * self.g + q * u - p * v
		
		b = np.zeros(3)
		b[0] = l + Tzb * self.ybp - Tyb * self.zbp + (self.Iyy - self.Izz) * q * r + self.Iyz * (q*q - r*r) + self.Ixz * p * q - self.Ixy * p * r
		b[1] = m + Txb * self.zbp - Tzb * self.xbp + (self.Izz - self.Ixx) * p * r + self.Ixz * (r*r - p*p) + self.Ixy * q * r - self.Iyz * p * q
		b[2] = n + Tyb * self.xbp - Txb * self.ybp + (self.Ixx - self.Iyy) * p * q + self.Ixy * (p*p - q*q) + self.Iyz * p * r - self.Ixz * q * r
		b += np.dot(self.h, np.array([p,q,r]))
		[p_dot, q_dot, r_dot] = np.dot(self.I_inv, b)
		
		xf_dot, yf_dot, zf_dot = quat.Body2Fixed([u,v,w],[e0,ex,ey,ez])
		
		e0_dot = .5 * (-ex*p-ey*q-ez*r)
		ex_dot = .5 * (e0*p-ez*q+ey*r)
		ey_dot = .5 * (ez*p+e0*q-ex*r)
		ez_dot = .5 * (-ey*p+ex*q+e0*r)
		
		return np.array([u_dot, v_dot, w_dot, p_dot, q_dot, r_dot, xf_dot, yf_dot, zf_dot, e0_dot, ex_dot, ey_dot, ez_dot])
	
	# 4th order Runge-Kutta Integrator
	def rk4(self):
		
		k1 = np.zeros(13)
		k2 = np.zeros(13)
		k3 = np.zeros(13)
		k4 = np.zeros(13)
		
		k1 = self.f(self.states)
		k2 = self.f(self.states + self.dt / 2. * k1)
		k3 = self.f(self.states + self.dt / 2. * k2)
		k4 = self.f(self.states + self.dt * k3)
		
		self.states += self.dt / 6. * (k1 + 2. * k2 + 2. * k3 + k4)
		self.t += self.dt
		self.states[9:13] = quat.NormalizeQuaternion(self.states[9:13])
	
	def calc_A_and_B(self, states, control):
		
		A = np.zeros((12,12))
		B = np.zeros((12,4))
		F = np.zeros((6,12))
		
		u,v,w,p,q,r,xf,yf,zf,phi,theta,psi=states
		tau,de,da,dr=control
		
		V = u #np.sqrt(u*u+v*v+w*w)
		# dVdu = u/V
		# dVdv = v/V
		# dVdw = w/V
		
		p_ = self.bw / 2. / V * p
		q_ = self.c_w / 2. / V * q
		r_ = self.bw / 2. / V * r
		
		alpha = w / u
		beta = v / u
		
		if self.constant_density:
			self.rho = self.rho_ref
			self.rho_sea = self.rho_ref
		else:
			self.rho = stdatmos.statee(-zf)[3]
			self.rho_sea = stdatmos.statee(0.)[3]
		
		
		k = .5 * self.rho * V ** 2. * self.Sw
		
		# lift coefficient
		CL = self.CL_ref + self.CLa * alpha + self.CLq_ * q_ + self.CLde * (de - self.de_ref)
		CLu = -self.CLa * w / u **2. - self.CLq_ * self.c_w * q / 2. / u **2.
		CLv = 0.
		CLw = self.CLa / u
		CLp = 0.
		CLq = self.CLq_ * self.c_w / 2. / V
		CLr = 0.
		CLxf = CLyf = CLzf = 0.
		CLphi = CLtheta = CLpsi = 0.
		
		CLtau = 0.
		CLde = self.CLde
		CLda = 0.
		CLdr = 0.
		
		# side coefficient
		CS = self.CYb * beta + self.CYp_ * p_ + self.CYr_ * r_ + self.CYda * da + self.CYdr * dr
		CSu = -self.CYb * v / u ** 2. - self.CYp_ * self.bw * p / 2. / u **2. - self.CYr_ * self.bw * r / 2. / u **2.
		CSv = self.CYb / u
		CSw = 0.
		CSp = self.CYp_ * self.bw / 2. / V
		CSq = 0.
		CSr = self.CYr_ * self.bw / 2. / V
		CSxf = CSyf = CSzf = 0.
		CSphi = CStheta = CSpsi = 0.
		
		CStau = 0.
		CSde = 0.
		CSda = self.CYda
		CSdr = self.CYdr
		
		# drag coefficient
		CD = self.CD0 + self.CD1 * CL + self.CD2 * CL**2. + self.CD3 * CS**2. + self.CDq_ * q_ + self.CDde * (de - self.de_ref)
		CDu = self.CD1 * CLu + 2. * self.CD2 * CL * CLu + 2. * self.CD3 * CS * CSu - self.CDq_ * self.c_w * q / 2. / u **2.
		CDv = self.CD1 * CLv + 2. * self.CD2 * CL * CLv + 2. * self.CD3 * CS * CSv
		CDw = self.CD1 * CLw + 2. * self.CD2 * CL * CLw + 2. * self.CD3 * CS * CSw
		CDp = self.CD1 * CLp + 2. * self.CD2 * CL * CLp + 2. * self.CD3 * CS * CSp
		CDq = self.CD1 * CLq + 2. * self.CD2 * CL * CLq + 2. * self.CD3 * CS * CSq + self.CDq_ * self.c_w / 2. / V
		CDr = self.CD1 * CLr + 2. * self.CD2 * CL * CLr + 2. * self.CD3 * CS * CSr
		CDxf = CDyf = CDzf = 0.
		CDphi = CDtheta = CDpsi = 0.
		
		CDtau = 0.
		CDde = self.CD1 * CLde + 2. * self.CD2 * CL * CLde + self.CDde
		CDda = 2. * self.CD3 * CS * CSda
		CDdr = 2. * self.CD3 * CS * CSdr
		
		# roll coefficient
		Cl = self.Clb * beta + self.Clp_ * p_ + self.Clr_ / self.CL_ref * CL * r_ + self.Clda * da + self.Cldr * dr
		Clu = -self.Clb * v / u **2. + self.Clr_ / self.CL_ref * (r_ * CLu - CL * self.bw * r / 2. / u **2.) - self.Clp_ * self.bw * p / 2. / u **2.
		Clv = self.Clb / u + self.Clr_ / self.CL_ref * r_ * CLv
		Clw = self.Clr_ / self.CL_ref * r_ * CLw
		Clp = self.Clp_ * self.bw / 2. / V + self.Clr_ / self.CL_ref * r_ * CLp
		Clq = self.Clr_ / self.CL_ref * r_ * CLq
		Clr = self.Clr_ / self.CL_ref * (CL * self.bw / 2. / V + r_ * CLr)
		Clxf = Clyf = Clzf = 0.
		Clphi = Cltheta = Clpsi = 0.
		
		Cltau = 0.
		Clde = self.Clr_ / self.CL_ref * r_ * CLde
		Clda = self.Clda
		Cldr = self.Cldr
		
		# pitch coefficient
		Cm = self.Cm_ref + self.Cma / self.CLa * (CL - self.CL_ref + CD * w / V) + self.Cmq_ * q_ + self.Cmde * (de - self.de_ref)
		Cmu = self.Cma / self.CLa * (CLu + CDu * w / u - CD * w / u **2.) - self.Cmq_ * self.c_w * q / 2. / u **2.
		Cmv = self.Cma / self.CLa * (CLv + CDv * w / u)
		Cmw = self.Cma / self.CLa * (CLw + CDw * w / u + CD / u)
		Cmp = self.Cma / self.CLa * (CLp + CDp * w / u)
		Cmq = self.Cma / self.CLa * (CLq + CDq * w / u) + self.Cmq_ * self.c_w / 2. / V
		Cmr = self.Cma / self.CLa * (CLr + CDr * w / u)
		Cmxf = Cmyf = Cmzf = 0.
		Cmphi = Cmtheta = Cmpsi = 0.
		
		Cmtau = 0.
		Cmde = self.Cma / self.CLa * (CLde + CDde * w / u) + self.Cmde
		Cmda = self.Cma / self.CLa * w / u * CDda
		Cmdr = self.Cma / self.CLa * w / u * CDdr
		
		# yaw coefficient
		Cn = self.Cnb / self.CYb * (CS - CD * v / V) + self.Cnp_ / self.CL_ref * CL * p_ + self.Cnr_ * r_ + self.Cnda * da + self.Cndr * dr
		Cnu = self.Cnb / self.CYb * (CSu - CDu * v / u + CD * v / u **2.) + self.Cnp_ / self.CL_ref * (CLu * p_ - CL * self.bw * p / 2. / u**2.)
		Cnv = self.Cnb / self.CYb * (CSv - CDv * v / u - CD / u) + self.Cnp_ / self.CL_ref * CLv * p_
		Cnw = self.Cnb / self.CYb * (CSw - CDw * v / u) + self.Cnp_ / self.CL_ref * CLw * p_
		Cnp = self.Cnb / self.CYb * (CSp - CDp * v / u) + self.Cnp_ / self.CL_ref * (CLp * p_ + CL * self.bw / 2. / V)
		Cnq = self.Cnb / self.CYb * (CSq - CDq * v / u) + self.Cnp_ / self.CL_ref * CLq * p_
		Cnr = self.Cnb / self.CYb * (CSr - CDr * v / u) + self.Cnp_ / self.CL_ref * CLr * p_ + self.Cnr_ * self.bw / 2. / V
		Cnxf = Cnyf = Cnzf = 0.
		Cnphi = Cntheta = Cnpsi = 0.
		
		Cntau = 0.
		Cnde = -self.Cnb / self.CYb * CDde * v / u + self.Cnp_ / self.CL_ref * p_ * CLde
		Cnda = self.Cnb / self.CYb * (CSda - v / u * CDda) + self.Cnda
		Cndr = self.Cnb / self.CYb * (CSdr - v / u * CDdr) + self.Cndr
		
		# forces and moments vector
		# row 1
		F[0,0] = k * (CLu * w / u - CL * w / u ** 2. - CSu * v / u + CS * v / u ** 2. - CDu) + self.rho * u * self.Sw * (CL * w/u - CS * v/u -CD) + tau * (self.rho / self.rho_sea) ** self.a * (self.T1 + 2. * self.T2 * u)
		F[0,1] = k * (CLv * w / u - CSv * v / u - CS / u - CDv)
		F[0,2] = k * (CLw * w / u - CL / u - CDw)
		F[0,3] = k * (CLp * w / u - CSp * v / u - CDp)
		F[0,4] = k * (CLq * w / u - CSq * v / u - CDq)
		F[0,5] = k * (CLr * w / u - CSr * v / u - CDr)
		# row 2
		F[1,0] = k * (CSu - CDu * v / u + CD * v / u **2.) + self.rho * u * self.Sw * (CS - CD * v / u)
		F[1,1] = k * (CSv - CDv * v / u - CD / u)
		F[1,2] = k * (CSw - CDw * v / u)
		F[1,3] = k * (CSp - CDp * v / u)
		F[1,4] = k * (CSq - CDq * v / u)
		F[1,5] = k * (CSr - CDr * v / u)
		# row 3
		F[2,0] = k * (-CLu - CDu * w / u + CD * w / u **2.) + self.rho * u * self.Sw * (-CL - CD * w / u)
		F[2,1] = k * (-CLv - CDv * w / u)
		F[2,2] = k * (-CLw - CDw * w / u - CD / u)
		F[2,3] = k * (-CLp - CDp * w / u)
		F[2,4] = k * (-CLq - CDq * w / u)
		F[2,5] = k * (-CLr - CDr * w / u)
		# row 4
		F[3,0] = k * self.bw * Clu + self.rho * u * self.Sw * self.bw * Cl
		F[3,1] = k * self.bw * Clv
		F[3,2] = k * self.bw * Clw
		F[3,3] = k * self.bw * Clp
		F[3,4] = k * self.bw * Clq
		F[3,5] = k * self.bw * Clr
		# row 5
		F[4,0] = k * self.c_w * Cmu + self.rho * u * self.Sw * self.c_w * Cm + self.zbp * tau * (self.rho / self.rho_sea) ** self.a * (self.T1 + 2. * self.T2 * u)
		F[4,1] = k * self.c_w * Cmv
		F[4,2] = k * self.c_w * Cmw
		F[4,3] = k * self.c_w * Cmp
		F[4,4] = k * self.c_w * Cmq
		F[4,5] = k * self.c_w * Cmr
		# row 6
		F[5,0] = k * self.bw * Cnu + self.rho * u * self.Sw * self.bw * Cn
		F[5,1] = k * self.bw * Cnv
		F[5,2] = k * self.bw * Cnw
		F[5,3] = k * self.bw * Cnp
		F[5,4] = k * self.bw * Cnq
		F[5,5] = k * self.bw * Cnr
		
		# Calculate the A matrix
		# =========================================================================================
		# row 1
		A[0,0] = self.g / self.W * F[0,0]
		A[0,1] = self.g / self.W * F[0,1] + r
		A[0,2] = self.g / self.W * F[0,2] - q
		A[0,3] = self.g / self.W * F[0,3]
		A[0,4] = self.g / self.W * F[0,4] - w
		A[0,5] = self.g / self.W * F[0,5] + v
		A[0,6:9] = 0.
		A[0,9] = 0.
		A[0,10] = -self.g * np.cos(theta)
		A[0,11] = 0.
		# row 2
		A[1,0] = self.g / self.W * F[1,0] - r
		A[1,1] = self.g / self.W * F[1,1]
		A[1,2] = self.g / self.W * F[1,2] + p
		A[1,3] = self.g / self.W * F[1,3] + w
		A[1,4] = self.g / self.W * F[1,4]
		A[1,5] = self.g / self.W * F[1,5] - u
		A[1,6:9] = 0.
		A[1,9] = self.g * np.cos(phi) * np.cos(theta)
		A[1,10] = -self.g * np.sin(phi) * np.sin(theta)
		A[1,11] = 0.
		# row 3
		A[2,0] = self.g / self.W * F[2,0] + q
		A[2,1] = self.g / self.W * F[2,1] - p
		A[2,2] = self.g / self.W * F[2,2]
		A[2,3] = self.g / self.W * F[2,3] - v
		A[2,4] = self.g / self.W * F[2,4] + u
		A[2,5] = self.g / self.W * F[2,5]
		A[2,6:9] = 0.
		A[2,9] = -self.g * np.sin(phi) * np.cos(theta)
		A[2,10] = -self.g * np.cos(phi) * np.sin(theta)
		A[2,11] = 0.
		# rows 4-6
		A[3:6,0] = np.dot(self.I_inv,F[3:6,0])
		A[3:6,1] = np.dot(self.I_inv,F[3:6,1])
		A[3:6,2] = np.dot(self.I_inv,F[3:6,2])
		x = np.array([self.Ixz*q-self.Ixy*r, (self.Izz-self.Ixx)*r-2.*self.Ixz*p-self.Iyz*q, (self.Ixx-self.Iyy)*q+2.*self.Ixy*p+self.Iyz*r])
		xx = F[3:6,3] + np.dot(self.h,[1.,0.,0.]) + x
		A[3:6,3] = np.dot(self.I_inv,xx)
		x = np.array([(self.Iyy-self.Izz)*r+2.*self.Iyz*q+self.Ixz*p, self.Ixy*r-self.Iyz*p, (self.Ixx-self.Iyy)*p-2.*self.Ixy*q-self.Ixz*r])
		xx = F[3:6,4] + np.dot(self.h,[0.,1.,0.]) + x
		A[3:6,4] = np.dot(self.I_inv,xx)
		x = np.array([(self.Iyy-self.Izz)*q-2.*self.Iyz*r-self.Ixy*p, (self.Izz-self.Ixx)*p+2.*self.Ixz*r+self.Ixy*q, self.Iyz*p-self.Ixz*q])
		xx = F[3:6,5] + np.dot(self.h,[0.,0.,1.]) + x
		A[3:6,5] = np.dot(self.I_inv,xx)
		# row 7
		A[6,0] = np.cos(theta) * np.cos(psi)
		A[6,1] = np.sin(phi) * np.sin(theta) * np.cos(psi) - np.cos(phi) * np.sin(psi)
		A[6,2] = np.cos(phi) * np.sin(theta) * np.cos(psi) + np.sin(phi) * np.sin(psi)
		A[6,3:9] = 0.
		A[6,9]=v*(np.cos(phi)*np.sin(theta)*np.cos(psi)+np.sin(phi)*np.sin(psi))+w*(-np.sin(phi)*np.sin(theta)*np.cos(psi)+np.cos(phi)*np.sin(psi))
		A[6,10]=-u*np.sin(theta)*np.cos(psi)+v*np.sin(phi)*np.cos(theta)*np.cos(psi)+w*np.cos(phi)*np.cos(theta)*np.cos(psi)
		A[6,11]=-u*np.cos(theta)*np.sin(psi)-v*(np.sin(phi)*np.sin(theta)*np.sin(psi)+np.cos(phi)*np.cos(psi))+w*(np.sin(phi)*np.cos(psi)-np.cos(phi)*np.sin(theta)*np.sin(psi))
		# row 8
		A[7,0] = np.cos(theta) * np.sin(psi)
		A[7,1] = np.sin(phi) * np.sin(theta) * np.sin(psi) + np.cos(phi) * np.cos(psi)
		A[7,2] = np.cos(phi) * np.sin(theta) * np.sin(psi) - np.sin(phi) * np.cos(psi)
		A[7,3:9] = 0.
		A[7,9]=v*(np.cos(phi)*np.sin(theta)*np.sin(psi)-np.sin(phi)*np.cos(psi))-w*(np.sin(phi)*np.sin(theta)*np.sin(psi)+np.cos(phi)*np.cos(psi))
		A[7,10]=-u*np.sin(theta)*np.sin(psi)+v*np.sin(phi)*np.cos(theta)*np.sin(psi)+w*np.cos(phi)*np.cos(theta)*np.sin(psi)
		A[7,11]=u*np.cos(theta)*np.cos(psi)+v*(np.sin(phi)*np.sin(theta)*np.cos(psi)-np.cos(phi)*np.sin(psi))-w*(np.cos(phi)*np.sin(theta)*np.cos(psi)+np.sin(phi)*np.sin(psi))
		# row 9
		A[8,0] = -np.sin(theta)
		A[8,1] = np.sin(phi) * np.cos(theta)
		A[8,2] = np.cos(phi) * np.cos(theta)
		A[8,3:9] = 0.
		A[8,9] = v * np.cos(phi) * np.cos(theta) - w * np.sin(phi) * np.cos(theta)
		A[8,10] = -u * np.cos(theta) - v * np.sin(phi) * np.sin(theta) - w * np.cos(phi) * np.sin(theta)
		A[8,11] = 0.
		# row 10
		A[9,0:3] = 0.
		A[9,3] = 1.
		A[9,4] = np.sin(phi) * np.tan(theta)
		A[9,5] = np.cos(phi) * np.tan(theta)
		A[9,6:9] = 0.
		A[9,9] = q * np.cos(phi) * np.tan(theta) - r * np.sin(phi) * np.tan(theta)
		A[9,10] = q * np.sin(phi) / np.cos(theta) ** 2. + r * np.cos(phi) / np.cos(theta) ** 2.
		A[9,11] = 0.
		# row 11
		A[10,0:3] = 0.
		A[10,3] = 0.
		A[10,4] = np.cos(phi)
		A[10,5] = -np.sin(phi)
		A[10,6:9] = 0.
		A[10,9] = -q * np.sin(phi) - r * np.cos(phi)
		A[10,10:12] = 0.
		# row 12
		A[11,0:3] = 0.
		A[11,3] = 0.
		A[11,4] = np.sin(phi) / np.cos(theta)
		A[11,5] = np.cos(phi) / np.cos(theta)
		A[11,6:9] = 0.
		A[11,9] = q * np.cos(phi) / np.cos(theta) - r * np.sin(phi) / np.cos(theta)
		A[11,10] = q * np.sin(phi) * np.tan(theta) / np.cos(theta) + r * np.cos(phi) * np.tan(theta) / np.cos(theta)
		A[11,11] = 0.
		
		
		# Calculate the B matrix
		# =========================================================================================
		# row 1
		B[0,0] = self.g / self.W * (self.rho / self.rho_sea) ** self.a * (self.T0 + self.T1 * u + self.T2 * u ** 2.)
		B[0,1] = self.g / self.W * k * (CLde * alpha - CSde * beta - CDde)
		B[0,2] = self.g / self.W * k * (CLda * alpha - CSda * beta - CDda)
		B[0,3] = self.g / self.W * k * (CLdr * alpha - CSdr * beta - CDdr)
		# row 2
		B[1,0] = 0.
		B[1,1] = self.g / self.W * k * (CSde - CDde * beta)
		B[1,2] = self.g / self.W * k * (CSda - CDda * beta)
		B[1,3] = self.g / self.W * k * (CSdr - CDdr * beta)
		# row 3
		B[2,0] = 0.
		B[2,1] = self.g / self.W * k * (-CLde - CDde * alpha)
		B[2,2] = self.g / self.W * k * (-CLda - CDda * alpha)
		B[2,3] = self.g / self.W * k * (-CLdr - CDdr * alpha)
		# row 4-6
		x = np.array([0.,self.zbp * (self.rho/self.rho_sea)**self.a*(self.T0+self.T1*u+self.T2*u**2.),0.])
		B[3:6,0] = np.dot(self.I_inv,x)
		x = k * np.array([self.bw * Clde, self.c_w * Cmde, self.bw * Cnde])
		B[3:6,1] = np.dot(self.I_inv,x)
		x = k * np.array([self.bw * Clda, self.c_w * Cmda, self.bw * Cnda])
		B[3:6,2] = np.dot(self.I_inv,x)
		x = k * np.array([self.bw * Cldr, self.c_w * Cmdr, self.bw * Cndr])
		B[3:6,3] = np.dot(self.I_inv,x)
		# row 7-12
		B[6:12,:] = 0.
		
		# print('\n\nA matrix \n\n{}\n\nB matrix\n\n{}'.format(A,B))
		
		return A, B
	
	def calc_gain(self):
		# print('\n'*100)
		self.text('Enter in values for Autopilot')
		
		self.V0 = float(input('Enter desired value for airspeed (ft/sec): '))
		self.climb = float(input('Enter desired value for climb (deg): ')) * np.pi / 180.
		if self.climb == 0.:
			self.Altitude = float(input('Enter desired value for altitude (ft): '))
			if self.Altitude >= 0.:
				self.alt_flag = True
			else:
				self.alt_flag = False
		else:
			self.alt_flag = False
		self.temp = float(input('Enter desired value for bank (deg): ')) * np.pi / 180.
		if self.temp == 0.:
			self.psi0 = float(input('Enter desired value for heading (deg): ')) * np.pi / 180.
			if self.psi0 > -180. * np.pi / 180. and self.psi0 <= 180. * np.pi / 180.:
				self.head_flag = True
			else:
				self.head_flag = False
		else:
			self.head_flag = False
		
		self.bank = self.temp
		# self.bank = 0.
		
		# self.xo, self.uo = self.trim_hunsaker(self.V0, self.phi, self.climb)
		self.xo, self.uo = self.trim(False, True)
		
		# print(self.xo,self.uo)
		
		A, B = self.calc_A_and_B(self.xo, self.uo)
		
		Q = np.zeros((12,12))
		
		u_err = .005
		v_err = .005
		w_err = .005
		p_err = 1.*np.pi/180.
		q_err = 1.*np.pi/180.
		r_err = 1.*np.pi/180.
		x_err = 10000.
		y_err = 10000.
		if self.alt_flag:
			z_err = .01
			theta_err = .005*np.pi/180.
		else:
			z_err = 10000.
			theta_err = .005*np.pi/180.
		if self.head_flag:
			phi_err = .01 * np.pi / 180.
			psi_err = .005 * np.pi / 180.
		else:
			phi_err = .005*np.pi/180.
			psi_err = 1000.
		Q[0,0] = 1. / u_err**2.
		Q[1,1] = 1. / v_err**2.
		Q[2,2] = 1. / w_err**2.
		Q[3,3] = 1. / p_err**2.
		Q[4,4] = 1. / q_err**2.
		Q[5,5] = 1. / r_err**2.
		Q[6,6] = 1. / x_err**2.
		Q[7,7] = 1. / y_err**2.
		Q[8,8] = 1. / z_err**2.
		Q[9,9] = 1. / phi_err**2.
		Q[10,10] = 1. / theta_err**2.
		Q[11,11] = 1. / psi_err**2.
		
		max_throttle_setting = .001
		max_aileron_deflection = .08*np.pi/180.
		max_elevator_deflection = .08*np.pi/180.
		max_rudder_deflection = .08*np.pi/180.
		R = np.zeros((4,4))
		R[0,0] = 1. / max_throttle_setting ** 2.
		R[1,1] = 1. / max_elevator_deflection ** 2.
		R[2,2] = 1. / max_aileron_deflection ** 2.
		R[3,3] = 1. / max_rudder_deflection ** 2.
		
		P = scipy.linalg.solve_continuous_are(A,B,Q,R)
		self.K = np.dot(np.linalg.inv(R), np.dot(np.matrix.transpose(B), P))
		
		for i in range(4):
			for j in range(12):
				if abs(self.K[i,j]) < 1.e-15:
					self.K[i,j] = 0.
					# pass
		
		# self.K[:,6:9] = 0.
		
		self.text('Gain Matrix')
		print(self.K)
		
		
		
		
		
		self.auto_pilot = True
		
		self.text('Begin Autopilot')
	


