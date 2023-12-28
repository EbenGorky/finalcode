from tkinter import *
from matplotlib import pyplot as plt
from math import *
import numpy

def damped_oscillations(m,b,k):
    try:
        x = numpy.linspace(0,100,1000)
        y = []
        for i in x:
            y.append(e**(-b*i/2*m)*cos(((k/m - (b**2/4*m**2))**(1/2))*i))
        plt.plot(x,y)
        plt.show()
    except:
        print('hello world!')
        x = numpy.linspace(0,100,1000)
        y = []
        for i in x:
            y.append(e**(-b*i/2*m)*cos(((k/m + (b**2/4*m**2))**0.5)*i))
        plt.plot(x,y)
        plt.show()

def sum_of_sin_wav(Numner_of_sin_wave):
    def abs():
        x = numpy.linspace(0,20,1000)
        y = list()
        y_ = 0
        for i in x:
            #y.append(sin(i))
            y_ = 0
            for j in a:
                y_ += sin(float(j.get())*i)
            y.append(y_)
        plt.plot(x,y)
        plt.show()

    a = list()
    w = Toplevel()
    for i in range(0,2*Numner_of_sin_wave,2):
        Label(w,text = f'enter the function of {int(i/2+1)}st function').grid(row = i,column = 0)
        e = Entry(w)
        e.grid(row = i+1,column = 0)
        a.append(e)
    b = Button(w,text='draw graph',command = abs).grid(row = len(a)*2,column = 0)
    print(a)
