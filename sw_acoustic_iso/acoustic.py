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
        self.__density_air = 1.18
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
    
    def cremer_model(self, frequencies):
        
        self.mass_sup
        self.stiffness
        self.freq_critic
        self.freq_res
        f_analysis = np.append(frequencies, self.__freq_critic)
        f_analysis.sort()
        reduction = []
        
        for f in f_analysis:
            if f < self.__freq_critic:
                r = 20*np.log10(self.__mass_sup*f) - 47 
                reduction.append(r)
            elif f == self.__freq_critic:
                r = 20*np.log10(self.__mass_sup*f) - 47 - 10*np.log10(np.pi/(4*(self.__loss_factor+self.__mass_sup/(485*np.sqrt(f))))) 
                reduction.append(r)
            else:
                r = 20*np.log10(self.__mass_sup*f) - 47 - 10*np.log10(np.pi/4*(self.__loss_factor+self.__mass_sup/(485*np.sqrt(f)))) + 10*np.log10(f/self.__freq_critic) + 10*np.log10(1-self.__freq_critic/f) 
                reduction.append(r)
                
        return f_analysis, reduction
    
    # MODELO DE DAVY
    def __shear(self, f):
        omega = 2 * np.pi * f
        chi = (1 + self.__poisson_module) / (0.87 + 1.12 * self.__poisson_module) 
        chi = chi * chi
        X = self.__thickness**2 / 12 
        QP = self.__young_module / (1 - self.__poisson_module**2) 
        C = -omega**2
        B = C * (1 + 2 * chi / (1 - self.__poisson_module)) * X; 
        A = X * QP / self.__density
        kbcor2 = (-B + np.sqrt(B * B - 4 * A * C)) / (2 * A) 
        kb2 = np.sqrt(-C / A)
        G = self.__young_module / (2 * (1 + self.__poisson_module))
        kT2 = -C * self.__density * chi / G
        kL2 = -C * self.__density / QP
        kS2 = kT2 + kL2
        ASI = 1 + X * (kbcor2 * kT2 / kL2 - kT2)
        ASI = ASI * ASI 
        BSI = 1 - X * kT2 + kbcor2 * kS2 / (kb2 * kb2) 
        CSI = np.sqrt(1 - X * kT2 + kS2 * kS2 / (4 * kb2 * kb2)) 
        
        return ASI / (BSI * CSI)
            
    def __sigma(self, G, freq):
        # Definición de constantes: 
        c0 = self.__vel_sound_air  #Velocidad sonido [m/s] 
        w = 1.3
        beta = 0.234 
        n = 2
        S = self.__lx * self.__ly
        U = 2 * (self.__lx + self.__ly)
        
        twoa = 4 * S / U
        
        k = 2 * np.pi * freq / c0 
        f = w * np.sqrt(np.pi / (k * twoa)) 
        
        if f > 1: 
            f = 1 
        
        h = 1 / (np.sqrt(k * twoa / np.pi) * 2 / 3 - beta)
        q = 2 * np.pi / (k**2 * S)
        qn = q**n
        
        if G < f:
            alpha = h / f - 1 
            xn = (h - alpha * G)**n 
        else:
            xn = G**n 
        
        rad = (xn + qn)**(-1 / n)
        
        return rad
    
    def __single_leaf_davy(self, f):
        
        po = self.__density_air # Densidad del aire [Kg/m3] 
        c0 = self.__vel_sound_air # Velocidad sonido [m/s] 
        cos21Max = 0.9 # Ángulo limite definido en el trabajo de Davy 
        
        critical_frequency = np.sqrt(12 * self.__density * (1 - self.__poisson_module**2) / self.__young_module) * c0**2 / (2 * self.__thickness * np.pi) 
        
        normal = po * c0 / (np.pi * f * self.__mass_sup) 
        
        e = 2 * self.__lx * self.__ly / (self.__lx + self.__ly)
        
        cos2l = c0 / (2 * np.pi * f * e)
        
        if cos2l > cos21Max:
            cos2l = cos21Max 
        
        tau1 = normal**2 * np.log((normal**2 + 1) / (normal**2 + cos2l)) # Con logaritmo en base e (ln)
        ratio = f / critical_frequency
        r = 1 - 1 / ratio
        
        if r < 0: 
            r = 0
        
        G = np.sqrt(r) 

        rad = self.__sigma(G, f)

        netatotal = self.__loss_factor + rad * normal

        z = 2 / netatotal

        y = np.arctan(z) - np.arctan(z * (1 - ratio))
        
        tau2 = normal**2 * rad**2 * y / (netatotal * 2 * ratio)
        tau2 = tau2 * self.__shear(f) 

        if f < critical_frequency: 
            tau = tau1 + tau2
        else: 
            tau = tau2

        single_leaf = -10 * np.log10(tau)

        return single_leaf

    def davy_model(self, filter_oct="third_oct"):
        self.mass_sup
        self.stiffness
        self.freq_critic
        self.freq_res
        
        reduction = []
        averages = 3 # % Promedio definido por Davy
        R = []
            
        if filter_oct == "oct": 
            frequencies = [31.5,63,125,250,500,1000,2000,4000,8000,16000]
            dB = 0.707
            octave = 1 
        
        if filter_oct == "third_oct": 
            frequencies=[20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6000,8000,10000,12500,16000,20000]
            dB = 0.236 
            octave = 3 
        
        f_analysis = np.append(frequencies, self.__freq_critic)
        f_analysis.sort()
        
        for f in f_analysis:  
            n_tot= self.__loss_factor + (self.__mass_sup/(485*np.sqrt(f)))
            ratio = f/self.__freq_critic
            limit = 2**(1/(2*octave))
            
            if (ratio < 1 / limit) or (ratio > limit):
                transmission_lost = self.__single_leaf_davy(f) 
            else:
                av_single_leaf = 0
                for j in range(averages):
                    factor = 2**((2*j-1-averages)/(2*averages*octave))
                    aux=10**(-self.__single_leaf_davy(f*factor)/10)
                    av_single_leaf = av_single_leaf + aux
                
                transmission_lost = -10*np.log10(av_single_leaf/averages)
        
            reduction.append(transmission_lost)
            
        return f_analysis, reduction
        
        # función [single_leaf] = Single_leaf_Davy(frequency, density, Young, Poisson, thickness, lossfactor, length, width) 
        
        
        
        # función [rad] = Sigma(G, freq, width, length) 
         
        
        # función [out] = shear(frequency, density, Young, Poisson, thickness) 
       