import pandas as pd
import openpyxl

class ExcelDataframe:
    def __init__(self, excel, sheet):
        self.excel = excel  # archivo
        self.sheet = sheet  # hoja
        self.dataframe = pd.read_excel(self.excel, self.sheet, header=1)

    def get_dataframe(self):
        
        return self.dataframe