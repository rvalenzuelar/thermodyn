# -*- coding: utf-8 -*-
"""
	Thermodynamic functions

	Raul Valenzuela
	August, 2015
"""

import numpy as np
import pandas as pd
from numpy.polynomial.polynomial import polyval
import bisect
import sys

class meteo(object):
	def __init__(self,**kwargs):
		for key,value in kwargs.iteritems():
			''' convert input value to numpy array '''
			if isinstance(value,list) or isinstance(value,int) or isinstance(value,float):
				value=np.asarray(value)
			if key == 'C':
				self.C = value # [°C]
			elif key == 'K':
				self.K = value # [K]
			elif key == 'theta':
				self.theta = value # [K]			
			elif key == 'Dewp':
				self.Dewp = value	# [°C]
			elif key == 'relh':
				self.relh = value	# [%]
			elif key == 'mixing_ratio':
				self.mixing_ratio = value # [kg/kg]
			elif key == 'mb' or key == 'hPa':
				self.pressure = value # [mb] or [hPa]
			elif key == 'bar':
				self.pressure = value # [mb]
			elif key == 'Pa':
				self.pressure = value # [hPa]
			elif key == 'agl_m':
				self.height = value # m altitude above ground level			
			elif key == 'm':
				self.elevation = value # m station elevation

		''' constants '''
		self.cp = 1004.6 # [J K-1 kg-1] specific heat at const press dry air
		self.Lv = 2.25E6 # [J kg-1] latent heat of evaporation
		self.Rd = 287.0 # [J K-1 kg-1] gas constant for dry air
		self.Rm = 461.4 # [J K-1 kg-1] gas constant for moist air
		self.p0 = 1000. # [hPa] mean sea level pressure 
		self.g  = 9.8 # [m s-1]

def parse_args(**kwargs):

	return meteo(**kwargs)

def sea_level_press(**kwargs):
	meteo=parse_args(**kwargs)	
	check_K=hasattr(meteo,'K')
	check_p=hasattr(meteo,'pressure')
	check_m=hasattr(meteo,'m')

	return meteo.pressure*np.exp((meteo.g*meteo.elevation)/(meteo.Rd*meteo.K))

def sat_vapor_press_lowe(**kwargs):
	"""  sat_wv_press = f(C[K]) [mb]
		valid for liquid water and -50°C<T<50°C
		error within 1%
		Lowe, 1977, JAM
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_K=hasattr(meteo,'K')
	if check_C:
		coef=[	6.107799961,
				4.436518521E-1,
				1.428945805E-2,
				2.650648471E-4,
				3.031240396E-6,
				2.034080948E-8,
				6.136820929E-11]
		return polyval(meteo.C,coef)				
	elif check_K:
		coef =[	6984.505294,
				-188.9039310,
				2.133357675,
				-1.288580973E-2,
				4.393587233E-5,
				-8.023923082E-8,
				6.136820929E-11]
		return polyval(meteo.K,coef)
	else:
		print "Error: check input arguments\n"

def sat_vapor_press_cc(**kwargs):
	"""  sat_wv_press = f(C[K]) [mb]
		from Clasius-Clapyron equation
		Bhoren and Albrecht (1998)
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_K=hasattr(meteo,'K')
	es0=6.11 # [mb]
	if check_C:
		Tk=meteo.C+273.15
		return es0*np.exp(19.84-5417/Tk)
	elif check_K:
		return es0*np.exp(19.84-5417/meteo.K)
	else:
		print "Error: check input arguments\n"

def sat_mix_ratio(**kwargs):
	"""  sat_mix_ratio = f(C,hPa {mb}) [kg/kg]
		Saucier, 1989, p.11
		Bohren and Albrecht, 1998, p.186
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_p=hasattr(meteo,'pressure')	
	if check_C and check_p:
		es = sat_vapor_press_cc(C=meteo.C)
		p=meteo.pressure 
		return (0.622*es)/(p - es ) # [kg/kg]
	else:
		print "Error: check input arguments\n"

def relative_humidity(**kwargs):
	""" 	relative_humidity = f(C,Dewp) [%]
		Lawrence, 2005, BAMS
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_Dewp=hasattr(meteo,'Dewp')		
	if check_C and check_relh:
		relh = np.asarray(100-5*(meteo.C- meteo.Dewp)) #[%]
		relh[relh>100.0] = 100.0
		return relh	
	else:
		print "\nError: check input arguments\n"

def dew_point(**kwargs):
	""" 	Dewp = f(C,relative_humidity) [C]
		Lawrence, 2005, BAMS
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_relh=hasattr(meteo,'relh')		
	if check_C and check_relh:
		relh =meteo.relh
		relh[relh>100.0] = 100.0		
		dewp = meteo.C - np.asarray((100. - relh)/5.) #[%]
		return dewp
	else:
		print "\nError: check input arguments\n"		

def virtual_temperature(**kwargs):
	""" Compute virtual temperature
		Tv = f(C{K}, mixing_ratio) [C] or [K]
		Thetav = f(theta, mixing_ratio) [K]
		Saucier, 1989, p.12
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_K=hasattr(meteo,'K')
	check_theta=hasattr(meteo,'theta')
	check_mixingr=hasattr(meteo,'mixing_ratio')
	if  check_C and check_mixingr:
		return meteo.C*(1+0.61*meteo.mixing_ratio)
	elif check_K and check_mixingr:
		return meteo.K*(1+0.61*meteo.mixing_ratio)
	elif check_theta and check_mixingr:
		return meteo.theta*(1+0.61*meteo.mixing_ratio)
	else:
		print "\nError: check input arguments\n"


def lcl_temperature(**kwargs):
	""" Lifting condensation level
		temperature
		Bolton, 1980, MWR (Eq22)
	"""
	meteo=parse_args(**kwargs)	
	check_K=hasattr(meteo,'K')	
	check_relh=hasattr(meteo,'relh')	
	if check_K and check_relh:
		Tk=meteo.K 
		relh=meteo.relh # [%]
		a = 1/(Tk-55)
		b = np.log(relh/100.0)/2840.
		return (1/(a-b))+55

def theta1(**kwargs):
	""" Compute potential temperature
		theta = f(C {K}, hPa {mb}) [K]
		Wallace&Hobbs, 2006, p.85		
	"""
	meteo=parse_args(**kwargs)		
	check_C=hasattr(meteo,'C')
	check_K=hasattr(meteo,'K')
	check_p=hasattr(meteo,'pressure')	
	c = meteo.Rd/meteo.cp 	
	if check_K and check_p:
		quotient=np.divide(meteo.p0,meteo.pressure)		
		return meteo.K*np.power(quotient,c) # [K]
	elif check_C and check_p:
		Tk=meteo.C+273.15
		quotient=np.divide(meteo.p0,meteo.pressure)		
		return Tk*np.power(quotient,c) # [K]
	else:
		print "\nError in theta1: check input arguments\n"


def theta2(**kwargs):
	""" Compute potential temperature
		theta = f(C {K}, hPa {mb}, mixing_ratio) [K]
		mixing_ratio in [kg/kg]
		Bolton, 1980, MWR	
	"""
	meteo=parse_args(**kwargs)	
	check_C = hasattr(meteo,'C')
	check_K = hasattr(meteo,'K')
	check_p = hasattr(meteo,'pressure')	
	check_mixingr = hasattr(meteo,'mixing_ratio')	
	if check_K and check_p and check_mixingr:
		Tk=meteo.K
		p=meteo.pressure
		MR=meteo.mixing_ratio 
		quotient=np.divide(meteo.p0, p)
		power=0.2584*(1-0.28*MR)
		return Tk*np.power(quotient, power) # [K]
	elif check_C and check_p and check_mixingr:
		Tk=meteo.C+273.15
		p=meteo.pressure
		MR=meteo.mixing_ratio
		quotient=np.divide(meteo.p0, p)
		power=0.2584*(1-0.28*MR)
		return Tk*np.power(quotient, power) # [K]
	else:
		print "\nError in theta2: check input arguments\n"

def theta_equiv1(**kwargs):
	""" Compute equivalent potential temperature
		theta_equiv = f(C{K}, hPa{mb}) [K]
		Saucier, 1989, p.14
		Wallace&Hobbs, 2006, p.85
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_K=hasattr(meteo,'K')	
	check_p=hasattr(meteo,'pressure')
	if check_K and check_p:
		Tc=meteo.K-273.15
		satmixr=sat_mix_ratio(C=Tc, hPa=meteo.pressure)
		th=theta1(K=meteo.K, hPa=meteo.pressure)
		exp=np.exp( np.divide( meteo.Lv*satmixr, meteo.cp*meteo.K ) )
		return th*exp
	elif check_C and check_p:
		satmixr=sat_mix_ratio(C=meteo.C, hPa=meteo.pressure)
		th=theta1(C=meteo.C, hPa=meteo.pressure)
		Tk=meteo.C+273.15
		exp=np.exp( np.divide( meteo.Lv*satmixr, meteo.cp*Tk ) )
		return th*exp
	else:
		print "\nError in theta_equiv1: check input arguments\n"	

def theta_equiv2(**kwargs):
	""" Compute equivalent potential temperature
		theta_equiv = f(C{K}, hPa{mb},mixing_ratio) [K]
		mixing_ratio in [kg/kg]
		Bolton, 1980, MWR (Eq43)
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_K=hasattr(meteo,'K')	
	check_p=hasattr(meteo,'pressure')
	check_relh=hasattr(meteo,'relh')
	check_mixingr = hasattr(meteo,'mixing_ratio')

	if check_K and check_p and check_mixingr and check_relh:
		Tk=meteo.K
		P=meteo.pressure
		MR=meteo.mixing_ratio
		RH=meteo.relh
		th=theta2(K=Tk,hPa=P,mixing_ratio=MR)
		TL=lcl_temperature(K=Tk,relh=RH)
		a=(3.376/TL)-0.00254
		b=1E3*MR*(1+0.81*MR)
		return th*np.exp(a*b)
	elif check_C and check_p and check_mixingr and check_relh:
		Tk=meteo.C+273.15
		P=meteo.pressure
		MR=meteo.mixing_ratio
		RH=meteo.relh
		th=theta2(K=Tk,hPa=P,mixing_ratio=MR)
		TL=lcl_temperature(K=Tk,relh=RH)
		a=(3.376/TL)-0.00254
		b=1E3*MR*(1+0.81*MR)
		return th*np.exp(a*b)
	else:
		print "\nError in theta_equiv2: check input arguments\n"	

def bv_freq_dry(**kwargs):
	"""	Compute dry Brunt-Vaisala frequency (N^2)
		Bohren and Albrecht (1998)		
	"""
	meteo=parse_args(**kwargs)	
	check_theta=hasattr(meteo,'theta')
	check_height=hasattr(meteo,'height')		

	''' creates layer field'''
	layer = make_layer(meteo.height,depth_m=kwargs['depth_m'],centered=kwargs['centered'])
	
	''' creates dataframe '''	
	d = {'theta':meteo.theta,'layer':layer}
	df=pd.DataFrame(d,index=meteo.height)

	''' creates group '''
	grp=df.groupby('layer')

	''' mean value of layers '''
	grp_mean = grp.mean()
	theta_mean = grp_mean.theta

	''' bottom values of layers '''
	grp_bot = grp.first()
	theta_bot = grp_bot.theta
	z_bot = grp.apply(get_min_hgt)

	''' top values of layers '''
	grp_top = grp.last()
	theta_top = grp_top.theta
	z_top = grp.apply(get_max_hgt)

	''' differentials '''
	dTheta = theta_top - theta_bot
	dZ = z_top - z_bot

	''' Brunt-Vaisala frequency '''
	quotient = dTheta.values/dZ.values
	bvf_raw = (meteo.g/theta_mean)*quotient

	''' get layer center '''
	hgt = grp.apply(get_layer_center)
	
	''' if last value in row is nan then drop it'''
	tail= bvf_raw.tail(1).values[0]
	if np.isnan(tail):
		d={'bvf_dry':bvf_raw[:-1].values}
		bvf=pd.DataFrame(d,index=hgt[:-1].values)
		bvf.index.name=['Height']
	else:
		d={'bvf_dry':bvf_raw.values}
		bvf=pd.DataFrame(d,index=hgt.values)
		bvf.index.name='Height'

	return bvf

def bv_freq_moist(**kwargs):
	"""	Compute moist Brunt-Vaisala frequency (N^2)
		Durran and Klemp (1982)		
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_K=hasattr(meteo,'K')
	check_hgt=hasattr(meteo,'height')		
	check_press=hasattr(meteo,'pressure')		
	check_mixingr=hasattr(meteo,'mixing_ratio')

	''' creates layer field '''
	layer = make_layer(meteo.height,depth_m=kwargs['depth_m'],centered=kwargs['centered'])

	''' add fields  '''
	sat_mixr = sat_mix_ratio(C=meteo.K-273.15,hPa=meteo.pressure)
	th = theta2(K=meteo.K,hPa=meteo.pressure,mixing_ratio=meteo.mixing_ratio)

	''' creates dataframe '''
	if check_C:
		Tk=meteo.C+273.15
		d = {'temp':Tk,'press':meteo.pressure,'sat_mixr':sat_mixr,'layer':layer,'theta':th}
	elif check_K:
		d = {'temp':meteo.K,'press':meteo.pressure,'sat_mixr':sat_mixr,'layer':layer,'theta':th}
	df=pd.DataFrame(d,index=meteo.height)

	''' creates group '''
	grp=df.groupby('layer')

	''' mean value of layers '''
	grp_mean = grp.mean()
	temp_mean = grp_mean.temp
	press_mean = grp_mean.press
	satmixr_mean = grp_mean.sat_mixr

	''' bottom values of layer '''
	grp_bot = grp.first()
	theta_bot = grp_bot.theta
	lnTheta_bot = np.log(theta_bot)
	sat_mixr_bot = grp_bot.sat_mixr
	z_bot = grp.apply(get_min_hgt)

	''' top values of layer '''	
	grp_top = grp.last()
	theta_top = grp_top.theta
	lnTheta_top = np.log(theta_top)
	sat_mixr_top = grp_top.sat_mixr
	z_top = grp.apply(get_max_hgt)	

	''' differentials '''
	dZ = z_top - z_bot
	dlnTheta = lnTheta_top - lnTheta_bot
	dQs = sat_mixr_top - sat_mixr_bot
	dQ = dQs

	''' Brunt-Vaisala frequency '''
	eps=meteo.Rd/meteo.Rm
	Lv=meteo.Lv
	cp=meteo.cp
	Rd=meteo.Rd
	g=meteo.g
	
	F1 = 1+((Lv*satmixr_mean)/(Rd*temp_mean))
	F2 = 1+(eps*np.power(Lv,2)*satmixr_mean)/(cp*Rd*np.power(temp_mean,2))
	F3 = dlnTheta/dZ
	F4 = (Lv/(cp*temp_mean))*dQ/dZ
	F5 = dQ/dZ
	bvf_raw = g*((F1/F2)*(F3+F4)-F5)

	''' get layer center '''
	hgt = grp.apply(get_layer_center)
	
	''' if last value in row is nan then drop it'''
	tail= bvf_raw.tail(1).values[0]
	if np.isnan(tail):
		d={'bvf_moist':bvf_raw[:-1].values}
		bvf=pd.DataFrame(d,index=hgt[:-1].values)
		bvf.index.name=['Height']
	else:
		d={'bvf_moist':bvf_raw.values}
		bvf=pd.DataFrame(d,index=hgt.values)
		bvf.index.name='Height'

	return bvf

""" 
	supporting functions 
------------------------------------
"""
def get_min_hgt(x):
	''' x is a pandas group instance '''
	return min(x.index)

def get_max_hgt(x):
	''' x is a pandas group instance '''
	return max(x.index)

def get_layer_center(x):
	''' x is a pandas group '''
	values = np.squeeze(x.index)
	s = values.shape
	if s:
		target = int(x.name)
		out = values[find_nearest2(values,target)]
		return out
	else:
		return values

def make_layer(height,**kwargs):
	''' makes a new field layer 
		so that values can be grouped 
		later for layer-based calculations
		(i.e. Brunt-Vaisala freq)
	'''
	centered=False
	for key,value in kwargs.iteritems():
		if key == 'depth_m': 
			depth = value
		elif key == 'centered':
			centered = value
	bottom=min(height)
	top=max(height)
	layer_value=range(0,int(top)+depth,depth)
	layers = [layer_value[bisect.bisect_left(layer_value,item)] for item in height]

	if centered:
		layer_value=range(0,int(top)+depth/2,depth/2)
		layers_half = [layer_value[bisect.bisect_left(layer_value,item)] for item in height]
		f = lambda x,y: x if x == y else y - depth/2
		return pd.Series(list(map(f,layers,layers_half)),index=height)
	else:
		return layers

def find_nearest2(array,target):

	""" See stackoverflow answer from Bi Rico """
	''' array must be sorted '''
	idx = array.searchsorted(target)
	idx = np.clip(idx, 1, len(array)-1)
	left = array[idx-1]
	right = array[idx]
	idx -= target - left < right - target
	return idx

