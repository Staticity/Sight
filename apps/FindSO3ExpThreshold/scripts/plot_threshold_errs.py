import os
import sys
import pandas as pd
import matplotlib.pyplot as plt

def RootDir():
    return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

def GetCSVPath(precision, filename):
    return os.path.join(RootDir(), "outputs", precision, filename)

def LoadCSV(precision, filename):
    filepath = GetCSVPath(precision, filename)
    if not os.path.exists(filepath):
        return pd.DataFrame()
    return pd.read_csv(filepath)

if __name__ == "__main__":
    float_results = LoadCSV('float', 'result.csv')
    double_results = LoadCSV('double', 'result.csv')
    
    float_results.plot(x='Threshold', y='TotalErr', title='float summary')
    double_results.plot(x='Threshold', y='TotalErr', title='double summary')

    verbose_float_results = LoadCSV('float', 'verbose.csv')
    verbose_double_results = LoadCSV('double', 'verbose.csv')

    thresholds = list(double_results['Threshold'])

    if not verbose_float_results.empty or not verbose_double_results.empty:
        # Render the existing plots
        plt.draw()
        plt.pause(.001)

        while True:
            idxStr = input("Please enter your desired threshold index, or enter to quit:")
            if idxStr == "":
                break
                
            threshIdx = int(idxStr)
            if threshIdx < 0 or threshIdx >= len(thresholds):
                print(f"{threshIdx} is not in the range [0, {len(thresholds) - 1}]")
                continue
            threshold = thresholds[threshIdx]

            if not verbose_float_results.empty:
                fdata = verbose_float_results[verbose_float_results['threshIdx'] == threshIdx]
                ax = fdata.plot(x='theta', y='err', title=f'float [threshold = {threshold}]')
                ax.set_ylim([0, .001])

            if not verbose_double_results.empty:
                ddata = verbose_double_results[verbose_double_results['threshIdx'] == threshIdx]
                ax = ddata.plot(x='theta', y='err', title=f'double [threshold = {threshold}]')
                ax.set_ylim([0, .001])
            
            plt.draw()
            plt.pause(.001)
            

    else:
        plt.show()
