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

class meteo(object):
	def __init__(self,**kwargs):
		for key,value in kwargs.iteritems():
			if isinstance(value,list) or isinstance(value,int):
				value=np.asarray(value)
			if key == 'C':
				self.C = value # [°C]
			elif key == 'K':
				self.K = value # [K]
			elif key == 'theta':
				self.theta = value # [K]			
			elif key == 'Dewp':
				self.Dewp = value	# [°C]
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
		self.cp = 1004.6	 # [J K-1 kg-1] specific heat at const press dry air
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

def sat_wv_press(**kwargs):
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

def sat_mix_ratio(**kwargs):
	"""  sat_mix_ratio = f(C,hPa {mb}) [kg/kg]
		Saucier, 1989, p.11
		Bohren and Albrecht, 1998, p.186
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_p=hasattr(meteo,'pressure')	
	if check_C and check_p:
		es = sat_wv_press(C=meteo.C)
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
	check_dewp=hasattr(meteo,'Dewp')		
	if check_C and check_dewp:
		relh = np.asarray(100-5*(meteo.C- meteo.Dewp)) #[%]
		relh[relh>100.0] = 100.0
		return relh	
	else:
		print "Error: check input arguments\n"

def dew_point(**kwargs):
	""" 	relative_humidity = f(C,Dewp) [%]
		Dewp = f(C,relative_humidity) [C]
		Lawrence, 2005, BAMS
	"""
	meteo=parse_args(**kwargs)	
	check_C=hasattr(meteo,'C')
	check_dewp=hasattr(meteo,'Dewp')		
	if check_C and check_dewp:
		relh = np.asarray(100-5*(meteo.C- meteo.Dewp)) #[%]
		relh[relh>100.0] = 100.0
		return relh	
	else:
		print "Error: check input arguments\n"		

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
		print "Error: check input arguments\n"

def theta(**kwargs):
	""" Compute potential temperature
		theta = f(C {K}, hPa {mb}) [K]
	"""
	meteo=parse_args(**kwargs)	
	c = meteo.Rd/meteo.cp 
	check_C=hasattr(meteo,'C')
	check_K=hasattr(meteo,'K')
	quotient=np.divide(meteo.p0,meteo.pressure)
	if check_K:
		return meteo.K*np.power(quotient,c) # [K]
	elif check_C:
		Tk=meteo.C+273.15
		return Tk*np.power(quotient,c) # [K]
	else:
		print "Error: check input arguments\n"

def theta_equiv(**kwargs):
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
		th=theta(K=meteo.K, hPa=meteo.pressure)
		exp=np.exp( np.divide( meteo.Lv*satmixr, meteo.cp*meteo.K ) )
		return th*exp
	elif check_C and check_p:
		satmixr=sat_mix_ratio(C=meteo.C, hPa=meteo.pressure)
		th=theta(C=meteo.C, hPa=meteo.pressure)
		Tk=meteo.C+273.15
		exp=np.exp( np.divide( meteo.Lv*satmixr, meteo.cp*Tk ) )
		return th*exp		

def bv_freq_dry(**kwargs):
	"""	Compute dry Brunt-Vaisala frequency (N^2)
		Bohren and Albrecht (1998)		
	"""
	meteo=parse_args(**kwargs)	
	check_theta=hasattr(meteo,'theta')
	check_height=hasattr(meteo,'height')		

	''' creates layer field'''
	layer = make_layer(meteo.height,depth_m=kwargs['depth_m'],centered=kwargs['centered'])
	
	''' creates dataframes '''	
	d_theta = {'theta':meteo.theta,'layer':layer}
	d_hgt = {'hgt':meteo.height,'layer':layer}
	grp_theta=pd.DataFrame(d_theta).groupby('layer')
	grp_hgt=pd.DataFrame(d_hgt).groupby('layer')


	''' mean value of the layer '''
	theta_mean = grp_theta.mean()

	''' differentials '''
	theta_first = grp_theta.first()
	theta_last = grp_theta.last()
	dTheta = theta_last - theta_first
	z_first = grp_hgt.first()
	z_last = grp_hgt.last()
	dZ = z_last - z_first

	''' Brunt-Vaisala frequency '''
	quotient = dTheta.values/dZ.values
	bvf_raw = (meteo.g/theta_mean)*quotient
	bvf_raw.columns=['bvf_dry']	
	foo= bvf_raw.tail(1).values[0]

	# print bvf_raw



	''' get layer center '''
	hgt = grp_hgt.apply(get_layer_center)
	hgt.columns=['Height']
	
	''' if last value in row is nan then drop it'''
	if np.isnan(foo[0]):
		bvf=bvf_raw[:-1]
		# bvf.index=hgt[:-1]
	else:
		bvf=bvf_raw
		# bvf.index=hgt

	# print bvf

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

	''' creates layer field '''
	layer = make_layer(meteo.height,depth_m=kwargs['depth_m'],centered=kwargs['centered'])

	''' add fields  '''
	sat_mixr = sat_mix_ratio(C=meteo.K-273.15,hPa=meteo.pressure)
	th = theta(K=meteo.K,hPa=meteo.pressure)

	''' creates dataframe '''
	if check_C:
		d = {'temp':meteo.C,'press':meteo.pressure,'sat_mixr':sat_mixr,'layer':layer,'theta':th}
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

	''' bot  values of layer '''
	grp_bot = grp.first()
	theta_bot = grp_bot.theta
	sat_mixr_bot = grp_bot.sat_mixr
	z_bot = grp.apply(get_min_hgt)

	''' top values of layer '''	
	grp_top = grp.last()
	theta_top = grp_top.theta
	sat_mixr_top = grp_top.sat_mixr
	z_top = grp.apply(get_max_hgt)	

	''' differentials '''
	dZ = z_top - z_bot
	dLogTheta = np.log(theta_top - theta_bot)
	dQs = sat_mixr_top - sat_mixr_bot
	dQw = dQs

	''' Brunt-Vaisala frequency '''
	eps=meteo.Rd/meteo.Rm
	Lv=meteo.Lv
	cp=meteo.cp
	Rd=meteo.Rd

	F1 = 1+(Lv*satmixr_mean*np.power(Rd*temp_mean,-1))
	F2 = 1+(eps*np.power(Lv,2)*satmixr_mean*np.power(cp*Rd*np.power(temp_mean,2),-1))
	F3 = dLogTheta/dZ
	F4 = Lv*np.power(cp*temp_mean,-1)*dQs/dZ
	F5 = dQw/dZ
	bvf = meteo.g*((F1/F2)*(F3+F4)-F5)
	# print bvf
	# exit()

	bvf_raw = (meteo.g/theta_mean)*quotient
	bvf_raw.columns=['bvf_dry']	
	foo= bvf_raw.tail(1).values[0]

	''' get layer center '''
	hgt = grp_hgt.apply(get_layer_center)
	hgt.columns=['Height']
	
	''' if last value in row is nan then drop it'''
	if np.isnan(foo[0]):
		bvf=bvf_raw[:-1]
		bvf.index=hgt[:-1]
	else:
		bvf=bvf_raw
		bvf.index=hgt

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
	''' x is a pandas group instance '''
	values = np.squeeze(x.values)
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

