# 2023 version
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.modeling.models import BlackBody
from astropy import units as u
from astropy import constants as const

# задаём параметры графиков для функций plot, histogram
plt.rcParams['figure.figsize'] = [14, 9]
plt.rcParams['font.size'] = 18


def load_data(filename, unpack=True):
    if type(filename) != str:
        raise ValueError('Имя файла должно быть указано в формате строки') 
    return np.loadtxt(filename, unpack=unpack)

def plot(x, y, scale='linear', title=None, xlabel=None, ylabel=None, label=None):
    if scale=='linear':
        plt.plot(x, y, label=label) # если нужен обычный график, используем простую функцию plt.plot
    elif scale=='log':
        plt.loglog(x, y, label=label) # если нужен логарифмический масштаб, пользуемся функцией plt.loglog - она логарифмирует обе оси
        plt.xticks([0.5, 1, 2, 5, 10, 20, 50, 100], [0.5, 1, 2, 5, 10, 20, 50, 100]) # передаём числа, которые пишутся по оси X
    else:
        raise ValueError('Неправильное значение "scale". Попробуйте "linear" для обычного графика или "log" для логарифмированных осей.')
    if title is not(None):
        plt.title(title, fontsize=20) # если указано название, выводим его функцией plt.title
    if label is not(None):            # то же самое с легендой и подписями осей
        plt.legend()
    if xlabel is not(None):
        plt.xlabel(xlabel)
    if ylabel is not(None):
        plt.ylabel(ylabel)
        
def histogram(x, bins=10): # здесь подобраны наиболее удобные бины для гистограмм в практике по фотометрии
    plt.hist(x, bins, edgecolor='black') # edgecolor рисует границу у прямоугольников в гистограмме

def spectral_ind(wavelength, flux):
    i1 = np.abs(wavelength-2).argmin()
    i2 = np.abs(wavelength-32).argmin()
    return (np.log10(flux[i2])-np.log10(flux[i1])) / (np.log10(wavelength[i2])-np.log10(wavelength[i1]))

def listdir(path):
    if path[-1]=='/': # если строка path кончается слэшем, к ней нужно просто добавить имена файлов
        return [path+i for i in os.listdir(path)] # возвращаем список, в котором собраны строки с файлами и заданным путём
    else: # если слэша в конце нет, его нужно добавить. Иначе будет что-то типа 'Data/TaurusObj1/txt'
        return [path+'/'+i for i in os.listdir(path)]

def xlim(xmin, xmax):
    plt.xlim([xmin, xmax])

def ylim(ymin, ymax):
    plt.ylim([ymin, ymax])

def acht_star(temp, radius, distance):
    wave = np.arange(0.36, 100, 0.01)*u.um # создаём numpy array (массив) со значениями длин волн, в конце "домножаем" на единице измерения - микрометры
    nu = const.c/wave.to(u.m) # переводим длины волн в частоту (нужно для корректного подсчёта потока)
    bb = BlackBody(temp*u.K) # функция из модуля astropy.modeling.models, считающая спектр чёрного тела (который зависит только от температуры)
    
    rad = radius*const.R_sun # переводим размер объекта в метры
    dist = distance*u.pc # к расстоянию добавляем единицу измерения - парсеки
    coef = np.pi*(rad**2)/((dist.to(u.m))**2) # считаем нормировочный коэффициент - квадрат площади диска объекта на квадрат расстояния до него (расстояние попутно переводим в метры)
    
    flux = bb(wave)*nu*coef # считаем поток и нормируем его на коэффициент
    return wave, flux

def acht_dust(temp, radius, distance): # почти то же самое, как и в прошлой функции
    wave = np.arange(5, 130, 0.01)*u.um
    nu = const.c/wave.to(u.m) 
    bb = BlackBody(temp*u.K)
    
    rad = radius*const.au
    dist = distance*u.pc
    coef = np.pi*(rad**2)/((dist.to(u.m))**2)
    
    flux = bb(wave)*nu*coef
    return wave, flux

def acht(wavel, waver, temp, radius, units, distance=140):
    wave = np.arange(wavel, waver, 0.1)*u.um
    nu = const.c/wave.to(u.m)
    
    bb = BlackBody(temp*u.K)
    dist = distance*u.pc
    
    if units=='rad_sun':
        rad = radius*const.R_sun
        coef = np.pi*(rad**2)/((dist.to(u.m))**2)
    elif units=='au':
        rad = radius*const.au
        coef = np.pi*(rad**2)/((dist.to(u.m))**2)
    
    flux = bb(wave)*nu*coef
    return wave, flux