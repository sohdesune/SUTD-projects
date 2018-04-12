# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 00:58:08 2018

@author: sohdesune
"""
'''
ln|T_w - T| - ln|T_w - T_amb| = -(1/tau) * t

1. Extract raw data for first 20s from csv
2. Plot complicated ln function vs t
3. Compute tau
4. [Remove outliers]
5. Write tau vs T_w to txt                                                  '''

from math import log as ln
from numpy import polyfit

'''============================================== Reading raw data from csv '''
T_amb = [26.375, 26.687, 28.312, 25.562, 28.312, 27.062, 31.125, 31.750, 
         28.625, 29.687, 26.375, 28.062, 30.125, 25.625, 27.250, 29.687, 
         31.125, 26.125, 33.437, 30.000, 27.000, 24.687, 31.000, 33.500, 
         33.187, 32.937, 29.500, 29.062, 28.062, 30.375, 30.437, 26.687, 
         32.312, 30.937, 23.937, 27.500, 32.125, 31.125, 32.250, 31.875, 
         25.250, 29.375, 34.312, 24.250, 31.750, 30.875, 29.687, 31.812, 
         30.875, 32.562, 30.812, 26.875, 33.187, 31.062, 25.062, 31.312]
T_w = [11.8, 11.8, 12.6, 12.6, 12.6, 13.2, 13.3, 13.3, 14.3, 14.3, 14.3, 16.0, 
       16.1, 17.4, 17.4, 18.8, 18.9, 20.0, 21.3, 21.3, 21.3, 22.5, 22.9, 29.5, 
       29.6, 29.8, 35.0, 35.3, 35.7, 37.4, 37.9, 38.5, 40.8, 41.3, 41.9, 43.6,
       43.9, 44.3, 46.3, 46.6, 47.0, 48.5, 48.8, 49.1, 50.1, 50.4, 50.9, 51.1,
       51.4, 51.7, 51.9, 52.3, 56.2, 56.7, 56.9, 57.4]

data = 'directory to csv file with temp vs time data'
f = open(data, 'r')
print('\nReading data from csv file.')
print('Directory:\n{}\n'.format(data))

line = f.readline()
i = 0
all_results = []
#each entry: [T_w, T_amb, list of x values, list of y values]

#each while loop reads the time and temp data for one T_w set
while line != '':
    x_val = []
    y_val = []
    
    time = line.strip().split(';')
    
    for elem in time:
        x_val.append(float(elem))
    
    line = f.readline()
    temp = line.strip().split(';')
    
    #compute ln values
    for elem in temp:
        try:
            value = ln(abs(T_w[i] - float(elem))) - ln(abs(T_w[i] - T_amb[i]))
        except ValueError:
            print('ValueError at T = {}'.format(elem))
            print('Occurred for T_w = {}, T_amb = {}'.format(T_w[i], T_amb[i]))
            value = ln(0.001)
        y_val.append(value)
    
    dataset = [x_val, y_val]
    all_results.append([T_w[i], T_amb[i], dataset])
    
    line = f.readline() #skip blank row
    line = f.readline()
    i += 1

f.close()
print('\nData compiled and modified into the complicated logarithm.')

'''====================================================== Performing linreg '''
linreg_results = []
#each entry: [T_w, gradient, y-intercept]

for result in all_results:
    grad, y_int = polyfit(result[2][0], result[2][1], 1)
    #print('T_w = {}: gradient {:+.3f}, y-intercept {:+.3f}'.format(result[0], grad, y_int))
    linreg_results.append([result[0], grad, y_int])

print('\nLinear regression performed for abovementioned logarithm vs time.')

'''========================================================= Determining tau'''
tau = [(-1/item[1]) for item in linreg_results]
twater = [item[0] for item in linreg_results]

print('\nTau values computed.')

'''=================================================== Plot regression line '''
grad, y_int = polyfit(twater, tau, 1)
print('\nRegression line calculated for full data set of tau against T_water.')
print('Gradient: {:.3f}   y-intercept: {:.3f}'.format(grad, y_int))

'''=========================== Remove anomalies and re-plot regression line '''
'''     dist from regr line = sqrt( vector^2 - projection^2 )

        projection, p = proj matrix, P * vector, b
        P = [1 grad]^T * [1 grad] / [1 grad] * [1 grad]^T
          = [ 1  g  ]
            [ g g^2 ] / (g^2 + 1)
        b = [x y+c]^T
        Pb = [x+gy+gc   g(x+gy+gc)]^T   /   (g^2+1)              '''

def dist_from_regr(g, c, x, y):
    x_proj = (x + g*y) / (1 + g**2)
    y_proj = g * x_proj + c
    distance = ((x-x_proj)**2 + (y-y_proj)**2)**0.5
    return distance

num_outliers = 0      #number of outliers you wish to remove

removed = 0
while removed < num_outliers:
    dist_list = []
    for i in range(len(twater)):
        dist_list.append(dist_from_regr(grad, y_int, twater[i], tau[i]))
    m = dist_list.index(max(dist_list))
    m_twater = twater.pop(m)
    m_tau = tau.pop(m)
    print('\n{:.1f},{:.1f} removed for being {:.1f} away from regression line.'.format(
            m_twater, m_tau, dist_list[m]))
    grad, y_int = polyfit(twater, tau, 1)
    print('New regression line plotted after removing outlier.')
    removed += 1

print('\n========================================================\n\nRESULT\n')
print('{} outliers removed from original data.'.format(num_outliers))
grad, y_int = polyfit(twater, tau, 1)
print('Final regression line plotted from {} pairs of values.'.format(len(twater)))
print('Gradient: {:.3f}   y-intercept: {:.3f}'.format(grad, y_int))

'''============================================= Write cleaned data to file '''
#send data to txt file to settle the remaining manipulations in Excel
def send_data():
    sendto = 'txt for writing to'
    f2 = open(sendto, 'a')
    for i in range(len(twater)):
        f2.write('{},{}\n'.format(twater[i], tau[i]))
    
    f2.close()
    print('\nCleaned data set written to text file for further processing.')
    print('Destination:\n{}'.format(sendto))

#checkpoint to ensure intentional writing
answer = input('Are you sure you want to write the results to txt? Y/N: ')
if answer == 'Y' or answer == 'y':
    send_data()
else:
    print('Data not written.')