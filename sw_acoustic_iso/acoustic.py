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
        self.__freq_density = 0
        self.__freq_res = 0

        #Llamado a calculo de propiedades:
        self.mass_sup
        self.stiffness
        self.freq_critic
        self.freq_res
        self.freq_density
        
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
    def freq_density(self):
        self.__freq_density = (self.__young_module/(2*np.pi*self.__density))*np.sqrt(self.__mass_sup/self.__stiffness)
        return self.__freq_density
    
    @property
    def freq_res(self):
        self.__freq_res = (self.__vel_sound_air**2/(4*self.__freq_critic))*(1/self.__lx**2 + 1/self.__ly**2)
        return self.__freq_res
    
    
    def data(self):
        return f'Material: {self.__name}\n Densidad: {self.__density}\n Módulo de Young: {self.__young_module}\n Factor de pérdidas: {self.__loss_factor}\n Módulo de Poisson: {self.__poisson_module}'
    
    # MODELO DE CREMER
    def cremer_model(self, frequencies):
        
        f_analysis = frequencies
        reduction = []
        
        for f in f_analysis:
            if f < self.__freq_critic or f >= self.__freq_density:
                r = float(round(20*np.log10(self.__mass_sup*f) - 47, 2))
                reduction.append(r)
            elif f == self.__freq_critic:
                n_tot = self.__loss_factor + self.__mass_sup/(485*np.sqrt(f)) 
                r = float(round(20*np.log10(self.__mass_sup*f) - 10*np.log10(np.pi/(4*n_tot)) - 47, 2))
            elif (f >= self.__freq_critic) and (f < self.__freq_density):
                n_tot = self.__loss_factor + self.__mass_sup/(485*np.sqrt(f)) 
                
                r = float(round(20*np.log10(self.__mass_sup*f) - 10*np.log10(np.pi/(4*n_tot)) + 10*np.log10(f/self.__freq_critic) + 10*np.log10(1 - self.__freq_critic/f) - 47,2))
                reduction.append(r)

        return f_analysis, reduction
    
    # MODELO DE SHARP
    
    def sharp_model(self, frequencies):
        
        f_analysis = frequencies
        f_interpolation = [] #guarda frecuencias de interpolación
        reduction = []
        r_1 = [] #guarda los r_1 en f>=fc
        r_2 = [] #guarda los r_2 en f>=fc
        
        for f in f_analysis:
            if f < self.__freq_critic/2:
                r = float(round(10*np.log10(1+((np.pi*self.__mass_sup*f)/(self.__density_air*self.__vel_sound_air))**2) - 5.5,2))
                reduction.append(r)
            elif f>= self.__freq_critic/2 and f<self.__freq_critic:
                f_interpolation.append(f)
            elif f >= self.__freq_critic:
                n_tot = self.__loss_factor + self.__mass_sup/(485*np.sqrt(f))
                
                r_1.append(float(round(10*np.log10(1+((np.pi*self.__mass_sup*f)/(self.__density_air*self.__vel_sound_air))**2) + 10*np.log10((2*f*n_tot)/(np.pi*self.__freq_critic)),2)))
                
                r_2.append(float(round(10*np.log10(1+((np.pi*self.__mass_sup*f)/(self.__density_air*self.__vel_sound_air))**2) - 5.5, 2)))
        
        # print("fc", self.__freq_critic)
        # print(f_analysis[index_stop])
        # print(f_analysis[index_start-1])
        # print("reduction ",reduction)
        if len(reduction) != 0:
            index_start = f_analysis.index(f_interpolation[0]) 
            index_stop = f_analysis.index(f_interpolation[-1]) + 1
            slope = (min(r_1[0], r_2[0]) - reduction[-1])/(f_analysis[index_stop] - f_analysis[index_start-1])
            b = reduction[-1] - slope*f_analysis[index_start-1] 
        
            for f in f_analysis[index_start:index_stop]:
                r = float(round(slope*f + b,2))
                reduction.append(r)

        reduction = reduction + min(r_1, r_2)
        
        return f_analysis, reduction
    
    
    
    def davy_model(self, frequencies):
        ro0 = self.__density_air # Densidad del aire [kg/m^3]
        c0 = self.__vel_sound_air # Velocidad del sonido [m/s]
        espesor = self.__thickness
        alto = self.__lx
        largo = self.__ly
        ro = self.__density # Densidad 
        E = self.__young_module # Módulo de Young
        eta = self.__loss_factor # Factor de pérdidas
        sigma = self.__poisson_module # Módulo de Poisson            
        m = self.__mass_sup # Masa superficial del elemento [kg/m^2]
        B = self.__stiffness  # Rigidez
        fc = self.__freq_critic # Frecuencia crítica [Hz]
        fd = self.__freq_density # Frecuencia de densidad [Hz]
        f11 = self.__freq_res
        R_davy = []
        average = 3 # % promedio definido por Davy
        dB = 0.236
        octava = 3
        
        def Single_leaf_Davy(frecuencia, ro, E, sigma, espesor, eta, alto, largo):
            cos21Max = 0.9 # Ángulo límite definido en el trabajo de Davy
            densidad_superficie = ro*espesor 
            frecuencia_critica = np.sqrt(12*ro*(1-sigma**2)/E)*((c0**2)/(2*np.pi*espesor))
            normal = (ro0*c0)/(np.pi*frecuencia*densidad_superficie)
            normal2 = normal**2
            e = (2*largo*alto)/(largo + alto)
            cos2l = c0/(2*np.pi*frecuencia*e)
            if cos2l > cos21Max:
                cos2l = cos21Max
            tau1 = normal2*np.log((normal2 + 1)/(normal2 + cos2l))
            ratio = frecuencia/frecuencia_critica
            r = 1 - 1/ratio
            if r < 0:
                r = 0
            G = np.sqrt(r)
            rad = Sigma(G, frecuencia, alto, largo)
            rad2= rad**2
            netatotal = eta + rad*normal
            z = 2/netatotal
            y = np.arctan(z) - np.arctan(z*(1-ratio))
            tau2 = (normal2*rad2*y)/(netatotal*2*ratio)
            tau2 = tau2*shear(frecuencia, ro, E, sigma, espesor)
            if frecuencia < frecuencia_critica:
                tau = tau1 + tau2
            else:
                tau = tau2
            single_leaf = -10*np.log10(tau)
            return single_leaf
        
        def Sigma(G, frecuencia, alto, largo):
            w = 1.3
            beta = 0.234
            n = 2
            S = largo*alto
            U = 2*(largo + alto)
            twoa = 4*S/U
            k = (2*np.pi*frecuencia)/c0
            f = w*np.sqrt((np.pi)/(k*twoa))
            if f > 1:
                f = 1
            h = 1/(np.sqrt((k*twoa)/(np.pi))*2/3 - beta)
            q = (2*np.pi)/((k**2)*S)
            qn = q**n
            if G < f:
                alpha = h/f - 1
                xn = (h - alpha*G)**n
            else:
                xn = G**n
            rad = (xn + qn)**(-1/n)
            return rad
        
        def shear(frecuencia, ro, E, sigma, espesor):
            omega = 2*np.pi*frecuencia
            chi = ((1 + sigma)/(0.87 + 1.12*sigma))**2
            X = (espesor**2)/12
            QP = E/(1 - sigma**2)
            C = -(omega)**2
            B = C*(1 + 2*(chi/(1 - sigma)))*X
            A = X*QP/ro
            kbcor2 = (-B + np.sqrt(B**2 - 4*A*C))/(2*A)
            kb2 = np.sqrt(-C/A)
            G = E/(2*(1 + sigma))
            kT2 = -C*ro*chi/G
            kL2 = -C*ro/QP
            kS2 = kT2 + kL2
            ASl = (1 + X*(kbcor2*kT2/kL2 - kT2))**2
            BSl = 1 - X*kT2 + kbcor2*kS2/(kb2**2)  
            CSl = np.sqrt(1 - X*kT2 + (kS2**2)/(4*kb2**2))
            out = ASl/(BSl*CSl)
            return out
        
        for f in frequencies:
            eta_total = eta + m/(485*np.sqrt(f))
            ratio = f/fc
            limit = 2**(1/(2*octava))
            if (ratio < 1/limit) or (ratio > limit):
                TLost = Single_leaf_Davy(f, ro, E, sigma, espesor, eta_total, alto, largo)
                R_davy.append(round(TLost,1))
            else:
                Avsingle_leaf = 0
                for i in range(1, average+1):
                    factor = 2**((2*i - 1 - average)/(2*average*octava))
                    aux = 10**(-(Single_leaf_Davy(f*factor, ro, E, sigma, espesor, eta_total, alto, largo))/10)
                    Avsingle_leaf += aux
                TLost = -10*np.log10(Avsingle_leaf/average)
                R_davy.append(round(TLost,1))
        R_davy = np.array(R_davy)
        
        return frequencies, R_davy
            
    # MODELO DE DAVY
    # def __shear(self, f):
    #     omega = 2 * np.pi * f
    #     chi = (1 + self.__poisson_module) / (0.87 + 1.12 * self.__poisson_module) 
    #     chi = chi * chi
    #     X = self.__thickness**2 / 12 
    #     QP = self.__young_module / (1 - self.__poisson_module**2) 
    #     C = -omega**2
    #     B = C * (1 + 2 *(chi / (1 - self.__poisson_module))) * X; 
    #     A = X * QP / self.__density
    #     kbcor2 = (-B + np.sqrt(B * B - 4 * A * C)) / (2 * A) 
    #     kb2 = np.sqrt(-C / A)
    #     G = self.__young_module / (2 * (1 + self.__poisson_module))
    #     kT2 = -C * self.__density * chi / G
    #     kL2 = -C * self.__density / QP
    #     kS2 = kT2 + kL2
    #     ASI = 1 + X * (kbcor2 * kT2 / kL2 - kT2)
    #     ASI = ASI * ASI 
    #     BSI = 1 - X * kT2 + kbcor2 * kS2 / (kb2 * kb2) 
    #     CSI = np.sqrt(1 - X * kT2 + kS2 * kS2 / (4 * kb2 * kb2)) 
        
    #     return ASI / (BSI * CSI)
            
    # def __sigma(self, G, freq):
    #     # Definición de constantes: 
    #     c_0 = self.__vel_sound_air  #Velocidad sonido [m/s] 
    #     w = 1.3
    #     beta = 0.234 
    #     n = 2
    #     S = self.__lx * self.__ly
    #     U = 2 * (self.__lx + self.__ly)
        
    #     twoa = 4 * S / U
        
    #     k = 2 * np.pi * freq / c_0 
    #     f = w * np.sqrt(np.pi / (k * twoa)) 
        
    #     if f > 1: 
    #         f = 1 
        
    #     h = 1 / (np.sqrt(k * twoa / np.pi) * 2 / 3 - beta)
    #     q = 2 * np.pi / (k**2 * S)
    #     qn = q**n
        
    #     if G < f:
    #         alpha = h / f - 1 
    #         xn = (h - alpha * G)**n 
    #     else:
    #         xn = G**n 
        
    #     rad = (xn + qn)**(-1 / n)
        
    #     return rad
    
    # def __single_leaf_davy(self, f):
        
    #     po = self.__density_air # Densidad del aire [Kg/m3] 
    #     c_0 = self.__vel_sound_air # Velocidad sonido [m/s] 
    #     cos21Max = 0.9 # Ángulo limite definido en el trabajo de Davy 
        
    #     critical_frequency = np.sqrt(12 * self.__density * (1 - self.__poisson_module**2) / self.__young_module) * c_0**2 / (2 * self.__thickness * np.pi) 
        
    #     normal = po * c_0 / (np.pi * f * self.__mass_sup) 
        
    #     e = 2 * self.__lx * self.__ly / (self.__lx + self.__ly)
        
    #     cos2l = c_0 / (2 * np.pi * f * e)
        
    #     if cos2l > cos21Max:
    #         cos2l = cos21Max 
        
    #     tau1 = normal**2 * np.log((normal**2 + 1) / (normal**2 + cos2l)) # Con logaritmo en base e (ln)
    #     ratio = f / critical_frequency
    #     r = 1 - 1 / ratio
        
    #     if r < 0: 
    #         r = 0
        
    #     G = np.sqrt(r) 

    #     rad = self.__sigma(G, f)

    #     netatotal = self.__loss_factor + rad * normal

    #     z = 2 / netatotal

    #     y = np.arctan(z) - np.arctan(z * (1 - ratio))
        
    #     tau2 = normal**2 * rad**2 * y / (netatotal * 2 * ratio)
    #     tau2 = tau2 * self.__shear(f) 

    #     if f < critical_frequency: 
    #         tau = tau1 + tau2
    #     else: 
    #         tau = tau2

    #     single_leaf = -10 * np.log10(tau)

    #     return single_leaf

    # def davy_model(self, filter_oct="third_oct"):
        
    #     reduction = []
    #     averages = 3 # % Promedio definido por Davy
            
    #     if filter_oct == "oct": 
    #         frequencies = [31.5,63,125,250,500,1000,2000,4000,8000,16000]
    #         dB = 0.707
    #         octave = 1 
        
    #     if filter_oct == "third_oct": 
    #         frequencies=[20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6000,8000,10000,12500,16000,20000]
    #         dB = 0.236 
    #         octave = 3 
        
    #     f_analysis = frequencies
        
    #     for f in f_analysis:  
    #         n_tot= self.__loss_factor + (self.__mass_sup/(485*np.sqrt(f)))
    #         ratio = f/self.__freq_critic
    #         limit = 2**(1/(2*octave))
            
    #         if (ratio < 1 / limit) or (ratio > limit):
    #             transmission_lost = self.__single_leaf_davy(f) 
    #             reduction.append(float(round(transmission_lost,2)))
    #         else:
    #             av_single_leaf = 0
    #             for j in range(1, averages+1):
    #                 factor = 2**((2*j-1-averages)/(2*averages*octave))
    #                 aux=10**(-self.__single_leaf_davy(f*factor)/10)
    #                 av_single_leaf = av_single_leaf + aux
                
    #             transmission_lost = -10*np.log10(av_single_leaf/averages)
    #             reduction.append(float(round(transmission_lost,2)))
            
    #     return f_analysis, reduction

    def iso_model(self, frequencies):
        rho_0 = self.__density_air # Densidad del aire [kg/m^3]
        c_0 = self.__vel_sound_air # Velocidad del sonido [m/s]
           
        t = self.__thickness
        height = self.__ly
        width = self.__lx
        f_c = self.__freq_critic
        f_11 = self.__freq_res
        m = self.__mass_sup
        reduction = []

        if self.__lx > self.__ly: # siempre l_1 > l_2 
            l_1 = height
            l_2 = width
        else:
            l_1 = width
            l_2 = height 

        for f in frequencies: 
            def n_total(f):
                n_tot = self.__loss_factor + m/(485*np.sqrt(f))
                return n_tot
            def delta1(f): 
                lamb = np.sqrt(f/f_c)
                delta_1 = (((1 - lamb**2)*np.log((1+lamb)/(1-lamb)) + 2*lamb)/(4*(np.pi**2)*(1-lamb**2)**1.5))
                return delta_1
            
            def delta2(f):
                lamb = np.sqrt(f/f_c)
                delta_2 = (8*(c_0**2)*(1 - 2*lamb**2))/((f_c**2)*(np.pi**4)*l_1*l_2*lamb*np.sqrt(1 - lamb**2))
                return delta_2
            
            def sigma_forced(f): 
                lambda_upper = - 0.964 - (0.5 + l_2/(np.pi*l_1))*np.log(l_2/l_1) + ((5*l_2)/(2*np.pi*l_1)) - (1/(4*np.pi*l_1*l_2*((2*np.pi*f)/c_0)**2))
                sigma = 0.5*(np.log(((2*np.pi*f)/c_0)*np.sqrt(l_1*l_2)) - lambda_upper) # Factor de radiación para transmisión forzadas
                return sigma
            
            def sigma1(f):
                sigma_1 = 1/(np.sqrt(1 - f_c/f)) 
                return sigma_1
            def sigma2(f):
                sigma_2 = 4*l_1*l_2*(f/c_0)**2
                return sigma_2
            def sigma3(f):
                sigma_3 = np.sqrt((2*np.pi*f*(l_1+l_2))/(16*c_0))
                return sigma_3
            
            if f > f_c/2:
                delta_2 = 0
            else:
                delta_2 = delta2(f)
            # Calculamos el factor de radiación por ondas libres  
            if f_11 <= f_c/2: 
                if f >= f_c:
                    sigma = sigma1(f)    
                elif f < f_c:
                    lamb = np.sqrt(f/f_c)
                    delta_1 = (((1 - lamb**2)*np.log((1+lamb)/(1-lamb)) + 2*lamb)/(4*(np.pi**2)*(1-lamb**2)**1.5))
                    delta_2 = (8*(c_0**2)*(1 - 2*lamb**2))/((f_c**2)*(np.pi**4)*l_1*l_2*lamb*np.sqrt(1 - lamb**2))
                    sigma = ((2*(l_1 + l_2)*c_0*delta_1/(l_1*l_2*f_c))) + delta_2
                sigma_2 = sigma2(f)
                if f<f_11 and f<f_c/2 and sigma > sigma_2:
                    sigma = sigma_2
            elif (f_11 > f_c/2):
                sigma_1 = sigma1(f)
                sigma_2 = sigma2(f)
                sigma_3 = sigma3(f)
                if (f < f_c) and (sigma_2 < sigma_3):
                    sigma = sigma_2
                elif (f > f_c) and (sigma_1 < sigma_3):
                    sigma = sigma_1
                else:
                    sigma = sigma_3
            if sigma > 2:
                sigma = 2 
            if f < f_c:
                sigma_f = sigma_forced(f)
                n_tot = n_total(f)
                tau = abs((((2*rho_0*c_0)/(2*np.pi*f*m))**2)*(2*sigma_f + (((l_1 + l_2)**2)/(l_1**2 + l_2**2))*np.sqrt(f_c/f)*(sigma**2)/n_tot))
                reduction.append(float(round(-10*np.log10(tau),2))) 
            elif f == f_c:
                n_tot = n_total(f)
                tau = abs((((2*rho_0*c_0)/(2*np.pi*f*m))**2)*((np.pi*(sigma)**2)/(2*n_tot)))
                reduction.append(float(round(-10*np.log10(tau),2)))
            elif f > f_c:
                n_tot = n_total(f)
                tau = abs((((2*rho_0*c_0)/(2*np.pi*f*m))**2)*((np.pi*f_c*(sigma)**2)/(2*f*n_tot)))
                reduction.append(float(round(-10*np.log10(tau),2)))
        
        return frequencies, reduction
    
