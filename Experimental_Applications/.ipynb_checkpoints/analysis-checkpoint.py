import numpy as np
import time
from plotly import graph_objs as go
import scipy
from scipy.optimize import curve_fit
import os
import scipy.stats as stats
import tkinter as tk
from tkinter import *
from tkinter.messagebox import askyesno

# # Tkinter Messagebox
# def tkinter_permission(file_type):
#     root = tk.Tk()
#     root.title('Permission accesss')
#     root.geometry('300x150')
#     root.eval('tk::PlaceWindow . center')
#     def confirmation():
#         globals()['permission'] = askyesno(title = None,message='Confirmation of Overwriting!')
#         if globals()['permission']:
#             top = Toplevel(root)     # using Toplevel
#             top.update_idletasks()
#             screen_width = top.winfo_screenwidth()
#             screen_height = top.winfo_screenheight()
#             size = tuple(int(_) for _ in top.geometry().split('+')[0].split('x'))
#             x = screen_width/2 - size[0]/2
#             y = screen_height/2 - size[1]/2
#             top.geometry("+%d+%d" % (x, y))

#             top.title(None)
#             Message(top, text=f'{file_type} is going to be overwritten!', padx=100, pady=100).pack()
#             top.after(1500, top.destroy)       
#             root.after(1500, root.destroy)
#         else:
#             root.destroy()
#     root_button = Button(root,text=f'Do you want to overwrite\n the exiting {file_type}?',command=confirmation)
#     root_button.pack(side = TOP, expand=0.5)
#     root.mainloop()


# Replacing mistakes in file naming
def replace_space(name):
    name = name.replace(' ', '_').replace('.','_').replace('__','_').replace('___','_').replace(':','_')
    if name[-1]=='_':  
        name=name[:-1]
    return name

# Saving file in given directory
def npz_save(folder_path,file_name,**dict_args):

    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    total_path = os.path.join(folder_path, file_name)
    np.savez(total_path,**dict_args)

    if os.path.exists(total_path)==False:
        raise Exception('Saved file does not exist!\n')
    elif os.stat(total_path).st_size == False:
        raise Exception('Saved file is empty!\n')
    else:
        print(f"saving data_file '{file_name}' is successful!\n")

    return total_path

# Curve_fitting Function
def exponential(x,y0,y_max,tau):
    return y0+y_max*np.exp(-x/tau)
def sigmoid(x,x0):
    return 1/(1+np.exp(-(x-x0)))
def inverse_sigmoid(x,x0):
    return 1/(1+np.exp(+(x-x0)))

# Function to calculate the Lifetime
def fit_func(x_old,y_old,fit_range=False,func='exp',guess_params=np.array([0.7,0.9,0.5e6])):
    indices = np.where(y_old!=0)
    yOld = y_old[indices]; xOld = x_old[indices]
    
    if type(fit_range)==np.ndarray:
        range_indicies = np.where(np.logical_and(xOld>=fit_range[0],xOld<=fit_range[1]))
        xOld = xOld[range_indicies]; yOld = yOld[range_indicies]
        
    if func.lower()=='exp':
        coefficient, covariance_matrix = curve_fit(exponential,xOld,yOld,p0=guess_params,absolute_sigma=False)
        x_new = xOld; y_new = exponential(x_new,*coefficient)
    if func.lower()=='sigmoid':
        coefficient, covariance_matrix = curve_fit(sigmoid,xOld,yOld,p0=guess_params,absolute_sigma=False)
        x_new = xOld; y_new = exponential(x_new,*coefficient)
    if func.lower()=='inverse_sigmoid':
        coefficient, covariance_matrix = curve_fit(inverse_sigmoid,xOld,yOld,p0=guess_params,absolute_sigma=False)
        x_new = xOld; y_new = exponential(x_new,*coefficient)

    error_bars = np.sqrt(np.diag(covariance_matrix))
    condition_number =  np.format_float_scientific(np.linalg.cond(covariance_matrix),precision=2)

#     # Different ways of 'Goodness of Fit' Test
#     chi_square_test, p_value = stats.chisquare(yOld, y_new)
#     ss_res = np.sum(np.square(yOld-y_new )); ss_total = np.sum(np.square(yOld-np.mean(yOld)))
#     r_squared = 1-(ss_res/ss_total)
#     mean_squared_error = np.square(np.subtract(y_new,yOld)).mean()

    print(f'Lifetime in nano_second is : {np.round(coefficient[2],1)} \u00B1 {np.round(error_bars[2],1)}\n')
#     print(f'Chi_square, p-value, R_squared,MeanSquaredError and Condition Number are : {np.round(chi_square_test,3)}\t{np.round(p_value,3)}\
#     \t{np.round(r_squared,3)}\t{np.round(mean_squared_error,5)}\t{condition_number}.\n')
#     if p_value<=0.05:
#         print('The p_value of fitting is low. Please check the fitting!')
    return x_new,y_new,coefficient,error_bars





#fig template 
fig_template = go.layout.Template()
fig_template.layout = {
    'template': 'simple_white+presentation',
    'autosize': False,
    'width': 800,
    'height': 600,
    # 'opacity': 0.2,
    'xaxis': {
        'ticks': 'inside',
        'mirror': 'ticks',
        'linewidth': 1.5+0.5,
        'tickwidth': 1.5+0.5,
        'ticklen': 6,
        'showline': True,
        'showgrid': False,
        'zerolinecolor': 'white',
        },
    'yaxis': {
        'ticks': 'inside',
        'mirror': 'ticks',
        'linewidth': 1.5+0.5,
        'tickwidth': 1.5+0.5,
        'ticklen': 6,
        'showline': True,
        'showgrid': False,
        'zerolinecolor': 'white'
        },
    'font':{'family':'mathjax',
            'size': 22,
            }
}

date = time.ctime()[4:10].replace(' ','_')
def simple_plot(x,y,show,x_name = "Time (&mu;s)",y_name = "Counts (T<sub>1</sub>)",mode='markers',title=date):
    
    fig = go.Figure()
    fig.add_scatter(x=x,y=y,mode=mode)
    
    if show==True:
        fig.update_layout(template = fig_template,width=800,height=600,
                         title=title,
                          xaxis_title=x_name,yaxis_title=y_name,
                         ) 
        fig.show()
    return fig.data[0]

def add_figures(figs,show,x_name = "Time (&mu;s)",y_name = "Counts (T<sub>1</sub>)",mode='markers',title=date):
    fig=go.Figure()
    for i in figs:
        fig.add_traces(i)
    if show==True:
        fig.update_layout(template = fig_template,width=800,height=600,
                         title=title,
                          xaxis_title=x_name,yaxis_title=y_name,
                         ) 
        fig.show()