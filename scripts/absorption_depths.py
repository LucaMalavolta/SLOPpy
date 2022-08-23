import numpy as np
import matplotlib.pyplot as plt
import cPickle as pickle

#w,f1,f2,f3,fall,sigma=np.loadtxt('/home/sicilia/Astro/Wyttenbach/HD189733_589289_1200_transmission_spectrum_TransSpec_atmcorr.rdb',unpack=True,skiprows=2)
#w,t1,sigma=np.loadtxt('/home/sicilia/Astro/Wyttenbach/transmission_spectrum_1.rdb',unpack=True)
#w,t2,sigma=np.loadtxt('/home/sicilia/Astro/Wyttenbach/transmission_spectrum_2.rdb',unpack=True)
#w,t3,sigma=np.loadtxt('/home/sicilia/Astro/Wyttenbach/transmission_spectrum_3.rdb',unpack=True)
#w,t,sigma=np.loadtxt('/home/sicilia/Astro/HD189733/list/average_transmission_spectrum2.txt',unpack=True)
#w,t,sigma=np.loadtxt('/home/sicilia/Astro/Wyttenbach/test.rdb',unpack=True)
#n, w, t, sigma = np.loadtxt('/home/sicilia/Astro/Wyttenbach/Wyttenbach_by_Casasayas.txt',unpack=True)
#our_average =  pickle.load(open('/home/sicilia/Astro/HD189733/list/HD1897333_HARPS_transmission_planetRF_average.p', "rb"))
#data_night1 = pickle.load(open('/home/sicilia/Astro/HD189733/list/HD1897333_HARPS_2006-09-07_transmission_planetRF_second_correction.p',"rb"))
#data_night2 = pickle.load(open('/home/sicilia/Astro/HD189733/list/HD1897333_HARPS_2007-07-19_transmission_planetRF_second_correction.p',"rb"))
data_night3 = pickle.load(open('/home/sicilia/Astro/HD189733/list/HD1897333_HARPS_2007-08-28_transmission_planetRF_second_correction.p',"rb"))

'''
w=our_average['wave']
t=(our_average['average']-1.)*100
t[t>1.1]=0.
sigma=our_average['average_err']
'''
w=data_night3['wave']
t=(data_night3['average']-1.)*100
t[t>1.1]=0.
sigma=data_night3['average_err']


#t = (t + 1)
#sigma = sigma/100

#central band
#C6_D2 = lambda wl: ((wl >= 5888.416) and (wl <= 5891.396)) or ((wl >= 5894.389) and (wl <= 5897.369))

w1=5889.906
w2=5895.879
wc=5892.89

#shifting the spectrum and not the passbands
#w = w - 0.04
#w = w + 0.16
#w = w + 0.075
#w = w + 0.05

#blue band & red band
B = (w >= 5874.89) & (w <= 5886.89)
R = (w >= 5898.89) & (w <= 5910.89)

#red band for Casasayas
#R = (w >= 5898.89) & (w <= 5907.89)

#alla c6 aggiungo e sottraggo 1.48
C12 = (w >= 5886.8925) & (w <= 5898.8925)
C6_D2 = (w >= 5888.426) & (w <= 5891.386)
C6_D1 = (w >= 5894.399) & (w <= 5897.359)
C3_D2 = (w >= 5889.156) & (w <= 5890.656)
C3_D1 = (w >= 5895.129) & (w <= 5896.629)
C1_5_D2 = (w >= 5889.531) & (w <= 5890.281)
C1_5_D1 = (w >= 5895.504) & (w <= 5896.254)
C0_75_D2 = (w >= 5889.718) & (w <= 5890.094)
C0_75_D1 = (w >= 5895.691) & (w <= 5896.067)
C0_375_D2 = (w >= 5889.812) & (w <= 5890.000)
C0_375_D1 = (w >= 5895.785) & (w <= 5895.973)

'''
#including -0.16 A
B = (w >= 5874.73) & (w <= 5886.73)
R = (w >= 5898.73) & (w <= 5910.73)

C12 = (w >= 5886.73) & (w <= 5898.73)
C6_D2 = (w >= 5888.246) & (w <= 5891.246)
C6_D1 = (w >= 5894.219) & (w <= 5897.219)
C3_D2 = (w >= 5888.996) & (w <= 5890.496)
C3_D1 = (w >= 5894.969) & (w <= 5896.469)
C1_5_D2 = (w >= 5889.371) & (w <= 5890.121)
C1_5_D1 = (w >= 5895.344) & (w <= 5896.094)
C0_75_D2 = (w >= 5889.558) & (w <= 5889.934)
C0_75_D1 = (w >= 5895.531) & (w <= 5895.907)
C0_375_D2 = (w >= 5889.652) & (w <= 5889.840)
C0_375_D1 = (w >= 5895.625) & (w <= 5895.813)


#sottraendo le bande come dice W e meno 0.16
C12 = (w >= 5886.73) & (w <= 5898.73)
C6_D2 = (w >= 5886.746) & (w <= 5892.746)
C6_D1 = (w >= 5892.719) & (w <= 5898.719)
C3_D2 = (w >= 5888.246) & (w <= 5891.246)
C3_D1 = (w >= 5894.219) & (w <= 5897.219)
C1_5_D2 = (w >= 5888.996) & (w <= 5890.496)
C1_5_D1 = (w >= 5894.969) & (w <= 5896.469)
C0_75_D2 = (w >= 5889.371) & (w <= 5890.121)
C0_75_D1 = (w >= 5895.344) & (w <= 5896.094)
C0_375_D2 = (w >= 5889.558) & (w <= 5889.934)
C0_375_D1 = (w >= 5895.531) & (w <= 5895.907)

'''

flux_C12, sum_weights_C12 = np.average(t[C12],axis=0,weights=1/sigma[C12]**2,returned=True)
flux_C6_D2, sum_weights_C6_D2 = np.average(t[C6_D2],axis=0,weights=1/sigma[C6_D2]**2,returned=True)
#flux_C6_D2, sum_weights_C6_D2 = np.average((np.array(filter(C6_D2, t))),weights=1/np.array(filter(C6_D2, sigma)**2),returned=True)
flux_C6_D1, sum_weights_C6_D1 = np.average(t[C6_D1],axis=0,weights=1/sigma[C6_D1]**2,returned=True)
flux_C3_D2, sum_weights_C3_D2 = np.average(t[C3_D2],axis=0,weights=1/sigma[C3_D2]**2,returned=True)
flux_C3_D1, sum_weights_C3_D1 = np.average(t[C3_D1],axis=0,weights=1/sigma[C3_D1]**2,returned=True)
flux_C1_5_D2, sum_weights_C1_5_D2 = np.average(t[C1_5_D2],axis=0,weights=1/sigma[C1_5_D2]**2,returned=True)
flux_C1_5_D1, sum_weights_C1_5_D1 = np.average(t[C1_5_D1],axis=0,weights=1/sigma[C1_5_D1]**2,returned=True)
flux_C0_75_D2, sum_weights_C0_75_D2 = np.average(t[C0_75_D2],axis=0,weights=1/sigma[C0_75_D2]**2,returned=True)
flux_C0_75_D1, sum_weights_C0_75_D1 = np.average(t[C0_75_D1],axis=0,weights=1/sigma[C0_75_D1]**2,returned=True)
flux_C0_375_D2, sum_weights_C0_375_D2 = np.average(t[C0_375_D2],axis=0,weights=1/sigma[C0_375_D2]**2,returned=True)
flux_C0_375_D1, sum_weights_C0_375_D1 = np.average(t[C0_375_D1],axis=0,weights=1/sigma[C0_375_D1]**2,returned=True)
flux_B, sum_weights_B = np.average(t[B],axis=0,weights=1/sigma[B]**2,returned=True)
flux_R, sum_weights_R = np.average(t[R],axis=0,weights=1/sigma[R]**2,returned=True)

deltaC12 = flux_C12 - (flux_B + flux_R)/2
deltaC6_D2 = flux_C6_D2 - (flux_B + flux_R)/2
deltaC6_D1 = flux_C6_D1 - (flux_B + flux_R)/2
deltaC3_D2 = flux_C3_D2 - (flux_B + flux_R)/2
deltaC3_D1 = flux_C3_D1 - (flux_B + flux_R)/2
deltaC1_5_D2 = flux_C1_5_D2 - (flux_B + flux_R)/2
deltaC1_5_D1 = flux_C1_5_D1 - (flux_B + flux_R)/2
deltaC0_75_D2 = flux_C0_75_D2 - (flux_B + flux_R)/2
deltaC0_75_D1 = flux_C0_75_D1 - (flux_B + flux_R)/2
deltaC0_375_D2 = flux_C0_375_D2 - (flux_B + flux_R)/2
deltaC0_375_D1 = flux_C0_375_D1 - (flux_B + flux_R)/2

delta_medio_6 = (deltaC6_D2 + deltaC6_D1)/2
delta_medio_3 = (deltaC3_D2 + deltaC3_D1)/2
delta_medio_1_5 = (deltaC1_5_D2 + deltaC1_5_D1)/2
delta_medio_0_75 = (deltaC0_75_D2 + deltaC0_75_D1)/2
delta_medio_0_375 = (deltaC0_375_D2 + deltaC0_375_D1)/2

sigma_deltaC12 = np.sqrt(1/sum_weights_C12 + 1/(2*sum_weights_B) + 1/(2*sum_weights_R)) * 100
sigma_deltaC6 = np.sqrt(1/sum_weights_C6_D2 + 1/sum_weights_C6_D1 + 1/sum_weights_B + 1/sum_weights_R)/2 * 100
sigma_deltaC3 = np.sqrt(1/sum_weights_C3_D2 + 1/sum_weights_C3_D1 + 1/sum_weights_B + 1/sum_weights_R)/2 * 100
sigma_deltaC1_5 = np.sqrt(1/sum_weights_C1_5_D2 + 1/sum_weights_C1_5_D1 + 1/sum_weights_B + 1/sum_weights_R)/2 * 100
sigma_deltaC0_75 = np.sqrt(1/sum_weights_C0_75_D2 + 1/sum_weights_C0_75_D1 + 1/sum_weights_B + 1/sum_weights_R)/2 * 100
sigma_deltaC0_375 = np.sqrt(1/sum_weights_C0_375_D2 + 1/sum_weights_C0_375_D1 + 1/sum_weights_B + 1/sum_weights_R)/2 * 100

print 'delta(12) =', deltaC12, ' +- ', sigma_deltaC12
print 'delta(6) = ', delta_medio_6, ' +- ', sigma_deltaC6
print 'delta(3) = ', delta_medio_3, ' +- ', sigma_deltaC3
print 'delta(1.5) = ', delta_medio_1_5, ' +- ', sigma_deltaC1_5
print 'delta(0.75) =', delta_medio_0_75, ' +- ', sigma_deltaC0_75
print 'delta(0.375) =', delta_medio_0_375, ' +- ', sigma_deltaC0_375

fig = plt.figure(figsize=(12, 6))
plt.plot(w,t)
plt.axvline(5892.89,c='k')

plt.axvspan(5886.89,5898.89,facecolor='g',alpha=0.3)
plt.axvspan(5874.89,5886.89,facecolor='b',alpha=0.3)
plt.axvspan(5898.89,5910.89,facecolor='r',alpha=0.3)
plt.axvspan(5888.426,5891.386,facecolor='g',alpha=0.4)
plt.axvspan(5894.399,5897.359,facecolor='g',alpha=0.4)
plt.axvspan(5889.156,5890.656,facecolor='g',alpha=0.5)
plt.axvspan(5895.129,5896.629,facecolor='g',alpha=0.5)
plt.axvspan(5889.531,5890.281,facecolor='g',alpha=0.6)
plt.axvspan(5895.504,5896.254,facecolor='g',alpha=0.6)
plt.axvspan(5889.718,5890.094,facecolor='g',alpha=0.7)
plt.axvspan(5895.691,5896.067,facecolor='g',alpha=0.7)
plt.axvspan(5889.812,5890.000,facecolor='g',alpha=0.8)
plt.axvspan(5895.785,5895.973,facecolor='g',alpha=0.8)

'''

plt.axvspan(5874.89,5886.89,facecolor='b',alpha=0.3)
plt.axvspan(5898.89,5910.89,facecolor='r',alpha=0.3)
plt.axvspan(5888.246,5891.246,facecolor='g',alpha=0.4)
plt.axvspan(5894.219,5897.219,facecolor='g',alpha=0.4)
plt.axvspan(5888.996,5890.496,facecolor='g',alpha=0.5)
plt.axvspan(5894.969,5896.469,facecolor='g',alpha=0.5)
plt.axvspan(5889.371,5890.121,facecolor='g',alpha=0.6)
plt.axvspan(5895.344,5896.094,facecolor='g',alpha=0.6)
plt.axvspan(5886.73,5898.73,facecolor='g',alpha=0.3)
plt.axvspan(5889.558,5889.934,facecolor='g',alpha=0.7)
plt.axvspan(5895.531,5895.907,facecolor='g',alpha=0.7)
plt.axvspan(5889.652,5889.840,facecolor='g',alpha=0.8)
plt.axvspan(5895.625,5895.813,facecolor='g',alpha=0.8)
'''
plt.xlabel('$\lambda$ [$\AA$]')
plt.ylabel('R-1')
plt.show()
