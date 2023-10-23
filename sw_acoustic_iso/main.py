import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import openpyxl
import excel
import acustica

def main():
    df_material = excel.ExcelToDataFrame('materiales_db.xlsx', 'Hoja1').get_dataframe()
    print(df_material)
    id_material = int(input('Ingrese el material a utilizar: '))-1
    # print(df_material.iloc[id_material])
    material = acustica.Material(df_material.iloc[id_material])
    material.lx = 4
    material.ly = 3
    material.thickness = 0.1
    print(material.data())
    print(f'Rigidez: {material.stiffness}')
    print(f'Masa superficial: {material.mass_sup}')
    print(f'Frecuencia crítica: {material.freq_critic}')
    print(f'Frecuencia de resonancia: {material.freq_res}')
    f_per_thirds = [20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000]
    f_cramer, r_cramer = material.cramer_model(f_per_thirds)
    print(r_cramer)
    
    
    plt.figure(figsize=(12,5))
    plt.semilogx(f_cramer, r_cramer, label="Modelo de Cramer")
    plt.xlim(f_cramer[0], f_cramer[-1])
    plt.ylim(0,140)
    plt.xlabel('Frecuencia [Hz]')
    plt.ylabel('R [dB]')
    plt.xticks(f_per_thirds, f_per_thirds, rotation=45)
    plt.legend(loc="lower right")
    # plt.ylim(-10,2)
    plt.grid()
    plt.show()
    
    
if __name__ == '__main__':
  main()
