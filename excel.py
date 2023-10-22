import pandas as pd
import openpyxl

class ExcelToDataFrame:
    def __init__(self, excel, sheet):
        self.__excel = excel  # archivo
        self.__sheet = sheet  # hoja
        self.__dataframe = pd.read_excel(self.__excel, self.__sheet, header=0).set_index('Id.')

    def get_dataframe(self):
        return self.__dataframe