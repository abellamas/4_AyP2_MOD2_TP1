import pandas as pd
import numpy as np

class Panel():
    def __init__(self, name, density, young_module, loss_factor, poisson_module, lx, ly, thickness):
        self.__name = name
        self.__density = float(density)
        self.__young_module = float(young_module)
        self.__loss_factor = float(loss_factor)
        self.__poisson_module = float(poisson_module)
        self.__lx = float(lx)
        self.__ly = float(lx)
        self.__thickness = float(thickness)
        self.__vel_sound_air = 343
        self.__stiffness = 0
        self.__mass_sup = 0
        self.__freq_critic = 0
        self.__freq_res = 0
    
    @property
    def lx(self):
        return self.__lx
    
    @lx.setter
    def lx(self, lx):
        self.__lx = lx
        
    @property
    def ly(self):
        return self.__ly
    
    @ly.setter
    def ly(self, ly):
        self.__ly = ly
        
    @property
    def thickness(self):
        return self.__thickness
    
    @thickness.setter
    def thickness(self, thickness):
        self.__thickness = thickness
        
    @property
    def vel_sound_air(self):
        return self.__vel_sound_air
    
    @vel_sound_air.setter
    def vel_sound_air(self, velocity):
        self.__vel_sound_air = velocity
    
    @property
    def stiffness(self):
        if self.__thickness != 0:
            self.__stiffness = (self.__young_module/(1-self.__poisson_module**2))*self.__thickness**3/12
            return self.__stiffness
        else:
            raise('Thickness cant be zero')

    @property
    def mass_sup(self):
        if self.thickness != 0:
            self.__mass_sup = self.__density*self.__thickness
            return self.__mass_sup
        else:
            raise('Thickness cant be zero')
    
    @property
    def freq_critic(self):
        if self.__thickness != 0:
            self.__freq_critic = int((self.__vel_sound_air**2/(2*np.pi))*np.sqrt(self.__mass_sup/self.__stiffness))
            return self.__freq_critic
        else:
            raise('Thickness cant be zero')
    
    @property
    def freq_res(self):
        self.__freq_res = (self.__vel_sound_air**2/(4*self.__freq_critic))*(1/self.__lx**2 + 1/self.__ly**2)
        return self.__freq_res
    
    
    def data(self):
        return f'Material: {self.__name}\n Densidad: {self.__density}\n Módulo de Young: {self.__young_module}\n Factor de pérdidas: {self.__loss_factor}\n Módulo de Poisson: {self.__poisson_module}'
    
    def cramer_model(self, frecuencies):
        
        self.mass_sup
        self.stiffness
        self.freq_critic
        self.freq_res
        f_analysis = np.append(frecuencies, self.__freq_critic)
        f_analysis.sort()
        r_vs_freq = []
        
        for f in f_analysis:
            if f < self.__freq_critic:
                r = 20*np.log10(self.__mass_sup*f) - 47 
                r_vs_freq.append(r)
            elif f == self.__freq_critic:
                r = 20*np.log10(self.__mass_sup*f) - 47 - 10*np.log10(np.pi/(4*(self.__loss_factor+self.__mass_sup/(485*np.sqrt(f))))) 
                r_vs_freq.append(r)
            else:
                r = 20*np.log10(self.__mass_sup*f) - 47 - 10*np.log10(np.pi/4*(self.__loss_factor+self.__mass_sup/(485*np.sqrt(f)))) + 10*np.log10(f/self.__freq_critic) + 10*np.log10(1-self.__freq_critic/f) 
                r_vs_freq.append(r)
                
        return f_analysis, r_vs_freq
        