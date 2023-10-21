import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import openpyxl
import excel

def main():
    excel_file =  excel.ExcelDataframe('materiales_db.xlsx', 'Hoja1')
    print(excel_file.get_dataframe())

if __name__ == '__main__':
  main()
