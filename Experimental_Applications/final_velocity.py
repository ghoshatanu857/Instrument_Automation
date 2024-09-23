import math
from PIL import Image
from PIL import ImageFilter
import numpy as np
import matplotlib.pyplot as plt

import tkinter as tk
import tkinter.simpledialog as inp
import tkinter.messagebox as out
from tkinter import ttk


def get_peak(fname, area):
    img = Image.open(fname)
    crop_img = img.crop(area)
    detail_img = crop_img.filter(ImageFilter.DETAIL)
    bw = detail_img.convert('L')  # convert image to monochrome
    pix = bw.load()
    size = bw.size
    xmax = int(size[0])
    ymax = int(size[1])
    xx = np.zeros(xmax)
    for i in range(0, xmax):
        intensity = []
        for j in range(0, ymax):
            intensity.append(int(pix[i, j]))

        imax = max(intensity)
        imin = min(intensity)
        tol = (imin + imax) / 2.0
        for j in range(ymax-1, -1, -1):
            piv = int(pix[i, j])
            if piv <= tol:
                xx[i] = j
                break

    ipeak = min(xx)
    return ipeak


def get_coord(fname, area, vpt):
    img = Image.open(fname)
    crop_img = img.crop(area)
    detail_img = crop_img.filter(ImageFilter.DETAIL)
#    edge_img = crop_img.filter(ImageFilter.EDGE_ENHANCE_MORE)
#    edgex_img = crop_img.filter(ImageFilter.FIND_EDGES)
#    crop_img.show()
#    detail_img.show()
#    edge_img.show()
#    edgex_img.show()
#    crop_img.save("crop_img.jpg")
#    detail_img.save("detail_img.jpg")
#    edge_img.save("edge_img.jpg")
#    edgex_img.save("edgex_img.jpg")

    bw = detail_img.convert('L')  # convert image to monochrome
#    bw = crop_img.convert('L')  # convert image to monochrome
#    bw.show()
#    bw.save("out.jpg")
#    r, g, b = crop_img.split()
#    r.show()
    pix = bw.load()
    size = bw.size
#    xmax = int(size[0])
    ymax = int(size[1])
#    print(xmax, ymax)
#    print(pix[0, 0])
#    print(pix[25, 25])
#    print(pix[131,86])
#    print(pix[251,86])
#    print(pix[257, 165])
#    print(pix[16, 165])
    edge_pt = 0
    intensity = []
    for i in range(0, ymax):
        intensity.append(int(pix[vpt, i]))

    imax = max(intensity)
    imin = min(intensity)
    tol = (imin + imax) / 2.0
    for i in range(ymax-1, -1, -1):
        piv = int(pix[vpt, i])
        if piv >= tol:
            edge_pt = i
            break

    return edge_pt
#    print(edge_pt)


def get_minima(imax, amax, fy):
    amin = fy[0]
    aval = 0
    for i in range(1, imax):
        if fy[i] < amin:
            amin = fy[i]
            aval = i
#    print(aval, amin)
    fxy = np.zeros(imax)
    for i in range(0, imax):
        if fy[i] > amax or fy[i] < amin:
            fxy[i] = 0
        else:
            fxy[i] = fy[i]

    return aval, amin, fxy


def rescale(imax, aval, ndim, fy):
    fxy = np.zeros(ndim)
    for i in range(aval, imax):
        fxy[i-aval] = fy[i]

    return fxy


def interpolation_linear(xm, fx):
    det = 0
    gx = np.zeros(xm)
    gx[0] = fx[0]
    x0 = 0
    y0 = fx[0]
    for i in range(1, xm):
        if fx[i] == 0:
            det = 1
        else:
            gx[i] = fx[i]
            x1 = i
            y1 = fx[i]
            if det == 1:
                dyx = (y1 - y0)/(x1 - x0)
                for j in range(x0+1, x1):
                    gx[j] = y0 + dyx*(j-x0)
                det = 0
            x0 = x1
            y0 = y1

    return gx


def interpolation(xm, fx):
    gx = fx
    x0 = 0
    y0 = fx[0]
    npt = 1
    det = 0

    x1 = 0      # just for initiation purpose
    y1 = 0      # just for initiation purpose

    for i in range(1, xm):
        if fx[i] != 0:
            npt += 1
            if npt == 2:
                x1 = i
                y1 = fx[i]
            elif npt >= 3:
                x2 = i
                y2 = fx[i]
                if det == 1:
                    f01 = (y1 - y0)/(x1 - x0)
                    f12 = (y2 - y1)/(x2 - x1)
                    f012 = (f12 - f01)/(x2 - x0)
                    for j in range(x0+1, x2):
                        gx[j] = y0 + f01*(j-x0) + f012*(j-x1)*(j-x0)

                    det = 0
                    x0 = x1
                    y0 = y1
                    x1 = x2
                    y1 = y2
                else:
                    x0 = x1
                    y0 = y1
                    x1 = i
                    y1 = fx[i]

            else:
                x1 = i
                y1 = fx[i]
                gx[i] = 0
        else:
            det = 1

    return gx


def get_velocity(xn, fx):
    dfx = np.zeros(xn)
    dfx[0] = 0
    for i in range(1, xn):
        dfx[i] = (fx[i+1] - fx[i-1])/2

    return dfx

# #####################################################################
# #####################################################################
# root = tk.Tk()

# out.showinfo("Welcome", "Hello, Sudip.")

# var = 1

# while var == 1:
#     password = inp.askstring("Password", "Give your password", show="*")
# #    print(password)
#     if password == "conchita":
#         break

# topFrame = tk.Frame(root, height=150, width=150, relief=tk.SUNKEN)
# # topFrame.pack()
# topFrame.grid()

# middleFrame = tk.Frame(root, height=150, width=150, relief=tk.SUNKEN)
# # topFrame.pack()
# middleFrame.grid()

# middleFrame1 = tk.Frame(root, height=150, width=150, relief=tk.SUNKEN)
# # topFrame.pack()
# middleFrame1.grid()

# middleFrame2 = tk.Frame(root, height=150, width=150, relief=tk.SUNKEN)
# # topFrame.pack()
# middleFrame2.grid()

# bottomFrame = tk.Frame(root, height=150, width=150, relief=tk.SUNKEN)
# # topFrame.pack()
# bottomFrame.grid()

# xfile = tk.Label(topFrame, text="File Information", bg="blue", fg="white", font="-weight bold")
# xfile.grid(row=0, columnspan=5, sticky=tk.N)

# label_00 = tk.Label(topFrame, text="File Name", bg="yellow green", fg="black")
# label_01 = tk.Label(topFrame, text="File Extension", bg="yellow green", fg="black")
# label_1 = tk.Label(topFrame, text="File No", bg="yellow green", fg="black")
# label_2 = tk.Label(topFrame, text="From", bg="yellow green", fg="black")
# label_3 = tk.Label(topFrame, text="To", bg="yellow green", fg="black")

# fn = tk.StringVar()
# fe = tk.StringVar()
# f1 = tk.StringVar()
# f2 = tk.StringVar()

# entry_00 = tk.Entry(topFrame, textvariable=fn, width=10)
# entry_00.insert(tk.END, 'capture-')
# entry_01 = tk.Entry(topFrame, textvariable=fe, width=5)
# entry_01.insert(tk.END, 'jpg')
# entry_1 = tk.Entry(topFrame, textvariable=f1, width=4)
# entry_2 = tk.Entry(topFrame, textvariable=f2, width=4)

# label_00.grid(row=1, sticky=tk.E, pady=10)
# label_01.grid(row=1, column=3, sticky=tk.E, pady=10)
# label_1.grid(row=2, sticky=tk.W, pady=5)
# label_2.grid(row=2, column=1, sticky=tk.E, pady=5)
# label_3.grid(row=2, column=3, sticky=tk.E, pady=5)

# entry_00.grid(row=1, column=1, sticky=tk.W, pady=5)
# entry_01.grid(row=1, column=4, sticky=tk.W, pady=5)

# entry_1.grid(row=2, column=2, sticky=tk.W, pady=10)
# entry_2.grid(row=2, column=4, sticky=tk.W, pady=10)

# label_4 = tk.Label(middleFrame, text="Crop Area Information", bg="blue", fg="white", font="-weight bold")
# label_4.grid(row=0, columnspan=5, sticky=tk.N, pady=10)

# ca1 = tk.StringVar()
# ca2 = tk.StringVar()
# ca3 = tk.StringVar()
# ca4 = tk.StringVar()

# entry_3 = tk.Entry(middleFrame, textvariable=ca1, width=4)
# entry_4 = tk.Entry(middleFrame, textvariable=ca2, width=4)
# entry_5 = tk.Entry(middleFrame, textvariable=ca3, width=4)
# entry_6 = tk.Entry(middleFrame, textvariable=ca4, width=4)

# entry_3.grid(row=1, column=0, pady=5)
# entry_4.grid(row=1, column=1, pady=5)
# entry_5.grid(row=2, column=3, pady=5)
# entry_6.grid(row=2, column=4, pady=5)

# label_05 = tk.Label(middleFrame1, text="Tip Value of the Scaffold: ", fg="black", font="-weight bold")
# label_05.grid(row=0, column=0, sticky=tk.N, pady=10)

# tip = tk.StringVar()
# entry_05 = tk.Entry(middleFrame1, textvariable=tip, width=4)
# entry_05.grid(row=0, column=1, pady=5)

# xfile = tk.Label(middleFrame2, text="Coordinate Scaling Information", bg="blue", fg="white", font="-weight bold")
# xfile.grid(row=0, columnspan=5, sticky=tk.N, pady=10)

# label_5 = tk.Label(middleFrame2, text="X axis scaling", bg="yellow green", fg="black")
# label_6 = tk.Label(middleFrame2, text="Y axis scaling", bg="yellow green", fg="black")
# label_5.grid(row=1, sticky=tk.E, pady=10)
# label_6.grid(row=1, column=3, sticky=tk.E, pady=10)

# label_5 = tk.Label(middleFrame2, text=":")
# label_5.grid(row=1, column=2, sticky=tk.E, padx=10, pady=10)

# fax = tk.StringVar()
# fay = tk.StringVar()

# entry_7 = tk.Entry(middleFrame2, textvariable=fax, width=5)
# entry_7.insert(tk.END, '1')
# entry_8 = tk.Entry(middleFrame2, textvariable=fay, width=5)
# entry_8.insert(tk.END, '1')

# entry_7.grid(row=1, column=1, sticky=tk.W, pady=10)
# entry_8.grid(row=1, column=4, sticky=tk.W, pady=10)

# printButton = tk.Button(bottomFrame, text="Next", command=root.quit, bg="green2", fg="blue4", font="-weight bold")
# printButton.grid(row=0, columnspan=5, sticky=tk.S)

# root.mainloop()

# label_06 = tk.Label(bottomFrame, text="Work in Progress...", fg="blue2")
# label_06.grid(row=1, columnspan=5, sticky=tk.N)


# strVar = tk.StringVar()
# strVar2 = tk.StringVar()

# label_07 = tk.Label(bottomFrame, textvariable=strVar, fg="blue2")
# label_07.grid(row=2, columnspan=5, sticky=tk.S)

# label_08 = tk.Label(bottomFrame, textvariable=strVar2, fg="blue2")
# label_08.grid(row=3, columnspan=5, sticky=tk.S)

# # s = ttk.Style()
# # s.theme_use('alt')

# progress = ttk.Progressbar(bottomFrame, length=350)
# progress.grid(row=4)

# stepmax = 100

# tini = int(f1.get())
# tfin = int(f2.get())
# crop_area = (int(ca1.get()), int(ca2.get()), int(ca3.get()), int(ca4.get()))

# fnam = str(fn.get())
# fext = '.' + str(fe.get())

# xscale = float(fax.get())
# yscale = float(fay.get())

# #####################################################################
# #####################################################################
# # ouf = open('velocity_vs_time.txt', 'w')
# ouf2 = open('position_vs_time.txt', 'w')

# # fnam = "capture-"
# # fext = ".jpg"
# # midd = str('{0:04}'.format(895))
# # name = fnam + midd + fext
# # name = "capture-0895.jpg"
# # crop_area = (618, 484, 932, 762)
# # print(crop_area[0])
# # get_coord(name, crop_area, 130)

# # peak = 768
# xpt = int(tip.get()) - crop_area[0]
# # print(xpt)

# # tini = 300
# # tfin = 555
# #####################################################################
# tnum = tfin-tini+1
# rmax = tnum+1
# step = (stepmax-0.1)/(rmax-1)

# midd = str('{0:04}'.format(int(tini)))
# name = fnam + midd + fext

# # xpt = get_peak(name, crop_area)
# # print(xpt)

# coord = np.zeros(tnum)
# for index in range(tini, tfin+1):
#     midd = str('{0:04}'.format(int(index)))
#     name = fnam + midd + fext
#     ed_pt = get_coord(name, crop_area, xpt)
#     coord[index-tini] = ed_pt
#     strVar.set("Image named " + name + " is processed successfully")
#     strVar2.set("Process Completed " + str(index-tini+1) + " out of " + str(tnum))
# #    strVar.set("New Text!" + str(index))
#     root.update_idletasks()

#     progress.step(step)
#     progress.update_idletasks()

# # tax = np.arange(tnum)
# # plt.plot(tax, coord, 'ro')
# # plt.show()

# up_l = (crop_area[3]-crop_area[1])*2.0/3.0

# lv, lo_l, coord2 = get_minima(tnum, up_l, coord)

# # print(up_l)
# # xcor = interpolation(lv, lo_l, up_l, tnum, coord)
# # plt.plot(tax, coord2, 'ro')# plt.show()

# xdim = tnum - lv

# coord3 = rescale(tnum, lv, xdim, coord2)

# # tax = np.arange(tnum)
# # plt.plot(tax, coord2, 'ro')
# # plt.show()

# # xcor = interpolation_linear(xdim, coord3)
# xcor = interpolation(xdim, coord3)

# # xdim = tnum - lv
# xax = np.zeros(xdim)
# for index in range(0, xdim):
#     xax[index] = float(index)

# # print(xscale, yscale)

# xval0 = xcor[0]
# for index in range(0, xdim):
#     xax[index] *= xscale
#     xcor[index] = (xcor[index]-xval0)*yscale
# #    print(xax[index], xcor[index])

# # xax = np.arange(xdim)
# ##############################################################################
# #     writing in the file  ==> position vs time plot
# ##############################################################################
# for ij in range(0, xdim):
#     ouf2.write('{:f} {:f}\n'.format(xax[ij], xcor[ij]))
# ##############################################################################
# #       showing the position vs time plot
# ##############################################################################
# plt.plot(xax, xcor, marker='o', linestyle='--', color='r')
# plt.xlabel('Time')
# plt.ylabel('Position')
# plt.title('Position vs. Time Plot')
# plt.show()
# ##############################################################################
# #                     log log plot
# ##############################################################################
# '''
# zax = np.zeros(xdim-1)
# zval = np.zeros(xdim-1)

# for index in range(1, xdim):
#     zax[index-1] = math.log10(xax[index])
#     zval[index-1] = math.log10(xcor[index])

# plt.plot(xax, xcor, marker='o', linestyle='--', color='r')
# plt.xlabel('Time')
# plt.ylabel('Position')
# plt.title('Position vs. Time Plot')
# plt.show()

# plt.plot(zax, zval, marker='o', linestyle='--', color='b')
# plt.xlabel('Log(Time)')
# plt.ylabel('Log(Position)')
# plt.title('Log(Position) vs. Log(Time) Plot')
# plt.show()

# for ij in range(0, xdim):
#     ouf2.write('{:f} {:f}\n'.format(xax[ij], xcor[ij]))
# '''
# ##############################################################################
# #                     Velocity Calculation
# ##############################################################################
# '''
# vdim = xdim - 1
# vel = get_velocity(vdim, xcor)

# vax = np.arange(vdim)
# plt.plot(vax, vel, 'ro')
# plt.show()

# vax = np.arange(vdim)
# for index in range(0, vdim):
#     if vel[index] < 0:
#         vel[index] = 0

# for ij in range(0, vdim):
#     ouf.write('{:f} {:f}\n'.format(vax[ij], vel[ij]))

# plt.plot(vax, vel, marker='o', linestyle='--')
# plt.xlabel('Time')
# plt.ylabel('Velocity')
# plt.title('Velocity vs. Time Plot')
# plt.show()
# '''
# ##############################################################################

# label_06 = tk.Label(bottomFrame, text="The Task is Finished", fg="blue4")
# label_06.grid(row=5, columnspan=5, sticky=tk.N)

# printButton = tk.Button(bottomFrame, text="Finish", command=root.quit, bg="green2", fg="blue4", font="-weight bold")
# printButton.grid(row=6, columnspan=5, sticky=tk.S)

# root.mainloop()
