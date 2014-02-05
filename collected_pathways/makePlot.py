import os
import csv


def openCSV():
    with open("knock_TSS_change_edit_for_latex.csv","rb") as csvfile:
        sr = csv.reader(csvfile,delimiter = ',')
        for row in sr:
            print row
            break
    

if __name__ == '__main__':
    openCSV()
