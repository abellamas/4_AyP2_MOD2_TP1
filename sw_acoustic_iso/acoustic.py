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
                r = 20*np.log10(self.__mass_sup*f) - 47 
                reduction.append(r)
            elif f == self.__freq_critic:
                n_tot = self.__loss_factor + self.__mass_sup/(485*np.sqrt(f)) 
                r = 20*np.log10(self.__mass_sup*f) - 10*np.log10(np.pi/(4*n_tot)) - 47
            elif (f >= self.__freq_critic) and (f < self.__freq_density):
                n_tot = self.__loss_factor + self.__mass_sup/(485*np.sqrt(f)) 
                
                r = 20*np.log10(self.__mass_sup*f) - 10*np.log10(np.pi/(4*n_tot)) + 10*np.log10(f/self.__freq_critic) + 10*np.log10(1 - self.__freq_critic/f) - 47
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
                r = 10*np.log10(1+((np.pi*self.__mass_sup*f)/(self.__density_air*self.__vel_sound_air))**2) - 5.5
                # print(f, r)
                reduction.append(r)
            elif f>= self.__freq_critic/2 and f<self.__freq_critic:
                f_interpolation.append(f)
            elif f >= self.__freq_critic:
                n_tot = self.__loss_factor + self.__mass_sup/(485*np.sqrt(f))
                
                r_1.append(10*np.log10(1+((np.pi*self.__mass_sup*f)/(self.__density_air*self.__vel_sound_air))**2) + 10*np.log10((2*f*n_tot)/(np.pi*self.__freq_critic)))
                
                r_2.append(10*np.log10(1+((np.pi*self.__mass_sup*f)/(self.__density_air*self.__vel_sound_air))**2) - 5.5)
        
        print("freqs ",f_interpolation)
        print(reduction)
        index_start = f_analysis.index(f_interpolation[0]) 
        index_stop = f_analysis.index(f_interpolation[-1]) + 1
        print(f_analysis[index_start-1], f_analysis[index_stop])
        # print("start ", index_start)
        slope = (min(r_1[0], r_2[0]) - reduction[-1])/(f_analysis[index_stop] - f_analysis[index_start-1])
        print("slope: ", slope)
        b = reduction[-1] - slope*f_analysis[index_start-1] 
        print("b: ", b)
        for f in f_analysis[index_start:index_stop]:
            print(f)
            r = slope*f + b
            print(f, r)
            reduction.append(r)

        reduction = reduction + min(r_1, r_2)
        
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
        
        f_analysis = frequencies
        
        for f in f_analysis:  
            n_tot= self.__loss_factor + (self.__mass_sup/(485*np.sqrt(f)))
            ratio = f/self.__freq_critic
            limit = 2**(1/(2*octave))
            
            if (ratio <= 1 / limit) or (ratio > limit):
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

    # def iso_12354(self):
    #     ro0 = self.__density_air # Densidad del aire [kg/m^3]
    #     c0 = self.__vel_sound_air # Velocidad del sonido [m/s]
    #     freq_tercio = np.array([20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160,
    #                             200, 250, 315, 400, 500, 630, 800, 1000, 1250,
    #                             1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000,
    #                             10000, 12500, 16000, 20000])    
    #     espesor = float(self.ui.espesor.text())
    #     alto = float(self.ui.alto.text())
    #     largo = float(self.ui.largo.text())
 
    #     R_iso = []
    #     if self.__lx > self.__ly: # siempre l1 > l2 
    #         l1 = largo
    #         l2 = alto
    #     else:
    #         l1 = alto
    #         l2 = largo 
    #     for f in freq_tercio: 
    #         def etatotal(f):
    #             eta_tot = eta + m/(485*np.sqrt(f))
    #             return eta_tot
    #         def delta1(f): 
    #             lamb = np.sqrt(f/fc)
    #             delta_1 = (((1 - lamb**2)*np.log((1+lamb)/(1-lamb)) + 2*lamb)/(4*(np.pi**2)*(1-lamb**2)**1.5))
    #             return delta_1
    #         def delta2(f):
    #             lamb = np.sqrt(f/fc)
    #             delta_2 = (8*(c0**2)*(1 - 2*lamb**2))/((fc**2)*(np.pi**4)*l1*l2*lamb*np.sqrt(1 - lamb**2))
    #             return delta_2
    #         def radforzada(f): 
    #             Lambda = - 0.964 - (0.5 + l2/(np.pi*l1))*np.log(l2/l1) + ((5*l2)/(2*np.pi*l1)) - (1/(4*np.pi*l1*l2*((2*np.pi*f)/c0)**2))
    #             sigma = 0.5*(np.log(((2*np.pi*f)/c0)*np.sqrt(l1*l2)) - Lambda) # Factor de radiación para ondas forzadas
    #             return sigma
    #         def sigma1(f):
    #             sigma_1 = 1/(np.sqrt(1 - fc/f)) 
    #             return sigma_1
    #         def sigma2(f):
    #             sigma_2 = 4*l1*l2*(f/c0)**2
    #             return sigma_2
    #         def sigma3(f):
    #             sigma_3 = np.sqrt((2*np.pi*f*(l1+l2))/(16*c0))
    #             return sigma_3
            
    #         if f > fc/2:
    #             delta_2 = 0
    #         else:
    #             delta_2 = delta2(f)
    #         # Calculamos el factor de radiación por ondas libres  
    #         if f11 <= fc/2: 
    #             if f >= fc:
    #                 rad_libre = sigma1(f)    
    #             elif f < fc:
    #                 lamb = np.sqrt(f/fc)
    #                 delta_1 = (((1 - lamb**2)*np.log((1+lamb)/(1-lamb)) + 2*lamb)/(4*(np.pi**2)*(1-lamb**2)**1.5))
    #                 delta_2 = (8*(c0**2)*(1 - 2*lamb**2))/((fc**2)*(np.pi**4)*l1*l2*lamb*np.sqrt(1 - lamb**2))
    #                 rad_libre = ((2*(l1 + l2)*c0*delta_1/(l1*l2*fc))) + delta_2
    #             sigma_2 = sigma2(f)
    #             if f<f11 and f<fc/2 and rad_libre > sigma_2:
    #                 rad_libre = sigma_2
    #         elif (f11 > fc/2):
    #             sigma_1 = sigma1(f)
    #             sigma_2 = sigma2(f)
    #             sigma_3 = sigma3(f)
    #             if (f < fc) and (sigma_2 < sigma_3):
    #                 rad_libre = sigma_2
    #             elif (f > fc) and (sigma_1 < sigma_3):
    #                 rad_libre = sigma_1
    #             else:
    #                 rad_libre = sigma_3
    #         if rad_libre > 2:
    #             rad_libre = 2 
    #         if f < fc:
    #             rad_forzada = radforzada(f)
    #             eta_total = etatotal(f)
    #             tao = abs((((2*ro0*c0)/(2*np.pi*f*m))**2)*(2*rad_forzada + (((l1 + l2)**2)/(l1**2 + l2**2))*np.sqrt(fc/f)*(rad_libre**2)/eta_total))
    #             R_iso.append(round(float(-10*np.log10(tao)),1)) 
    #         elif f == fc:
    #             eta_total = etatotal(f)
    #             tao = abs((((2*ro0*c0)/(2*np.pi*f*m))**2)*((np.pi*(rad_libre)**2)/(2*eta_total)))
    #             R_iso.append(round(float(-10*np.log10(tao)),1))
    #         elif f > fc:
    #             eta_total = etatotal(f)
    #             tao = abs((((2*ro0*c0)/(2*np.pi*f*m))**2)*((np.pi*fc*(rad_libre)**2)/(2*f*eta_total)))
    #             R_iso.append(round(float(-10*np.log10(tao)),1))
    #     R_iso = np.array(R_iso)
    #     return R_iso
    
    # def sigma_libre():
        
    # def sigma_forzado():
    
    def iso_model(self, frequencies):
        rho_0 = self.__density_air
        c_0 = self.__vel_sound_air
        l_x = self.__lx
        l_y = self.__ly
        t = self.__thickness
        m = self.__mass_sup
        c_l = np.sqrt(self.__young_module/self.__density) #vel long del material m/s
        f_p = c_l/(5.5*self.__thickness)
        f_11 = self.__freq_res
        f_c = c_0/(1.8*c_l*t) 
        f_c_eff = 0
        tau = 0
        reduction = []
        
        
        for f in frequencies:
            n_tot = self.__loss_factor + np.pi/(485*np.sqrt(f))
            k_wave = 2*np.pi*f/c_0
            lambda_f = np.sqrt(f/f_c)
            
            #factor de radiación para transmisión forzada 
            if l_x > l_y:
                lambda_upper = -0.964 - (0.5 + l_y/(np.pi*l_x))*np.log(l_y/l_x)+(5*l_y)/(2*np.pi*l_x)-1/(k_wave**2*4*np.pi*l_x*l_y)
                sigma_f = 0.5*(np.log(k_wave*np.sqrt(l_x*l_y))-lambda_upper)
            elif l_x <= l_y:
                lambda_upper = -0.964 - (0.5 + l_x/(np.pi*l_y))*np.log(l_x/l_y)+(5*l_x)/(2*np.pi*l_y)-1/(k_wave**2*4*np.pi*l_x*l_y)
                sigma_f = 0.5*(np.log(k_wave*np.sqrt(l_x*l_y))-lambda_upper)
                
            #factor de radiación para ondas de flexión libres
            sigma_1 = 1/np.sqrt(1-f_c/f)
            sigma_2 = 4*l_x*l_y*(f/c_0)**2
            sigma_3 = np.sqrt((2*np.pi*f*(l_x+l_y))/(16*c_0))
            
            if f_11<=f_c/2:
                if f>=f_c:
                    if f<f_p:
                        f_c_eff = f_c*(4.05*(t*f/c_l)+np.sqrt(1+(4.05*t*f/c_l)))
                    else:
                        f_c_eff = 2*f_c*(f/f_c)**3
                        
                    sigma_1 = 1/np.sqrt(1-f_c/f)
                    sigma = sigma_1
                elif f<f_c:
                    delta_1 = ((1-lambda_f**2)*np.log((1+lambda_f)/(1-lambda_f))+2*lambda_f)/(4*np.pi**2*(1-lambda_f**2)**1.5)
                    
                    if f>f_c/2:
                        delta_2 = 0
                    else:
                        delta_2 = (8*c_0**2*(1-2*lambda_f**2))/(f_c**2*np.pi**4*l_x*l_y*lambda_f*np.sqrt(1-lambda_f**2))
                    
                    sigma = ((2*(l_x+l_y)*c_0)/(l_x*l_y*f_c))*delta_1 + delta_2
                    
                    if f_11<f_c and f<f_11 and sigma > sigma_2:
                        sigma = sigma_2
            else:
                if f<f_c and sigma_2<sigma3:
                    sigma = sigma_2
                elif f>f_c and sigma_1<sigma_3:
                    sigma = sigma_1
                else:
                    sigma = sigma_3
                
                    
            
            if f>f_c:
                tau = ((2*rho_0*c_0)/(2*np.pi*f*m))**2 * (np.pi*f_c*sigma**2)/(2*f*n_tot)
            elif f==round(f_c,0):
                tau = ((2*rho_0*c_0)/(2*np.pi*f*m))**2 * (np.pi*sigma**2)/(2*n_tot)
            else:
                tau = ((2*rho_0*c_0)/(2*np.pi*f*m))**2 * (2*sigma_f + ((l_x+l_y)**2/(l_x**2+l_y**2))*sqrt(f_c/f)*sigma**2/n_tot)
            
            r = -10*np.log10(tau)
            reduction.append(r)
        
        return frequencies, reduction