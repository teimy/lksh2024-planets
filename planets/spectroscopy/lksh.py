import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [16, 9]
plt.rcParams['font.size'] = 18

def load_data(filename, unpack=True):
    return np.loadtxt(filename, unpack=unpack)

def plot(x, y, xlabel=None, ylabel=None, label=None):
    plt.plot(x, y, label=label)

    if label is not(None):
        plt.legend()

    if xlabel is not(None):
        plt.xlabel(xlabel)
    
    if ylabel is not(None):
        plt.ylabel(ylabel)

def xlim(xmin, xmax):
    plt.xlim([xmin, xmax])

def ylim(ymin, ymax):
    plt.ylim([ymin, ymax])

def gaussian(x, mu, std_dev, amplitude):
    return amplitude*np.exp(-1/2*((x - mu)/std_dev)**2)

def beerlambert(wavelength, k, n, l=100, r=1000):
    l = l # переводим метры в км
    tau = np.exp(-k*n*l)
    if r is None:
        return wavelength, tau
    diff = np.diff(wavelength)
    if np.sum(np.abs((diff - diff[0]))>0.01)>1:
        raise ValueError('Сетка по длине волны должна быть однородная')
    delta_wv = diff[0]
    wv_psf = np.arange(-5000*delta_wv,5001*delta_wv, delta_wv)
    psf = gaussian(wv_psf, mu=0, std_dev=np.mean(wavelength)/r, amplitude=1)
    t = np.convolve(tau, psf, mode='valid')/np.sum(psf)
    wvl_valid = np.convolve(wavelength, psf, mode='valid')/np.sum(psf)
    return wvl_valid, t

def convolve(wavelength, tau, r=1000):
    diff = np.diff(wavelength)
    if np.sum(np.abs((diff - diff[0]))>0.01)>1:
        raise ValueError('Сетка по длине волны должна быть однородная')
    delta_wv = diff[0]
    wv_psf = np.arange(-5000*delta_wv,5001*delta_wv, delta_wv)
    psf = gaussian(wv_psf, mu=0, std_dev=np.mean(wavelength)/r, amplitude=1)
    t = np.convolve(tau, psf, mode='valid')/np.sum(psf)
    wvl_valid = np.convolve(wavelength, psf, mode='valid')/np.sum(psf)
    return wvl_valid, t
def transmission(location, CO2=0.96, H2O=1000*10**(-6), CH4=0.0 , CO=0.0, HDO=0.0, HCl=0.0, O2=0.0, O3=0.0, r=1000):
    if location not in ['Mars2020', 'Zhurong', 'Curiosity']:
        raise ValueError('Проверьте чтовы вы правильно указали локацию. Возможные локации: Mars2020, Zhurong и Curiosity.')
    mixing_ratios = {'CO2': CO2, 'H2O': H2O, 'CH4': CH4, 'CO':CO, 'HDO':HDO, 'HCl':HCl, 'O2':O2, 'O3':O3}
    counter = 0  
    ts = []
    for gas in mixing_ratios.keys():
        if mixing_ratios[gas] != 0.0:
            counter += 1
            wvl, k = np.loadtxt('data/'+location+'/absorption-coefficients/'+gas+'.txt', unpack=True)
            wvl_valid, t_tmp = beerlambert(wvl, k, mixing_ratios[gas], l=100, r=None)
            ts.append(t_tmp)

    t = ts[0]
    for i in range(1,counter):
        t = t*ts[i]

    diff = np.diff(wvl_valid)
    if np.sum(np.abs((diff - diff[0]))>0.01)>1:
        raise ValueError('Сетка по длине волны должна быть однородная')
    delta_wv = diff[0]
    wv_psf = np.arange(-5000*delta_wv,5001*delta_wv, delta_wv)
    psf = gaussian(wv_psf, mu=0, std_dev=np.mean(wvl_valid)/r, amplitude=1)
    t = np.convolve(t, psf, mode='valid')/np.sum(psf)
    wvl_valid = np.convolve(wvl_valid, psf, mode='valid')/np.sum(psf)

    return wvl_valid, t

