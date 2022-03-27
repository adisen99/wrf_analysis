# Python  modulde containing functions for wrf analysis in python
# for the analysis of Uttarakhand 2013 extreme rainfall event

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import matplotlib as mpl
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from wrf import (to_np, interplevel, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, ALL_TIMES, CoordPair, vertcross)

# utility functions

dom1 = Dataset("./data/wrfout_d01.nc")
dom2 = Dataset("./data/wrfout_d02.nc")


test1 = getvar(dom1, "RAINC") # important keep it
test2 = getvar(dom2, "RAINC") # important keep it

# Get the cartopy mapping object
# Get the latitude and longitude points
lats1, lons1 = latlon_coords(test1)
lats2, lons2 = latlon_coords(test2)
cart_proj1 = get_cartopy(test2)
cart_proj2 = get_cartopy(test2)

def plot_background1(ax):
    ax.add_feature(cfeature.COASTLINE)
    # ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    gridliner = ax.gridlines(draw_labels=True, linewidth=1.0, color="gray", linestyle="--", alpha=0.5)
    gridliner.top_labels=False
    gridliner.right_labels=False
    gridliner.xlines=False
    gridliner.ylines=False
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(test1))
    ax.set_ylim(cartopy_ylim(test1))
    return ax

def plot_background2(ax):
    ax.add_feature(cfeature.COASTLINE)
    # ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    gridliner = ax.gridlines(draw_labels=True, linewidth=1.0, color="gray", linestyle="--", alpha=0.5)
    gridliner.top_labels=False
    gridliner.right_labels=False
    gridliner.xlines=False
    gridliner.ylines=False
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(test2))
    ax.set_ylim(cartopy_ylim(test2))
    return ax

# Using function

def plot_single(var, day, cmap, levels, cbar_label, domain=1):
    if domain == 1:
        ax = plt.axes(projection = cart_proj1)
        plot_background1(ax)
        ax.plot([77,82], [28, 28], color='k', transform=ccrs.PlateCarree())
        ax.plot([77,82], [32, 32], color='k', transform=ccrs.PlateCarree())
        ax.plot([77,77], [28, 32], color='k', transform=ccrs.PlateCarree())
        ax.plot([82,82], [28, 32], color='k', transform=ccrs.PlateCarree())
        ax.scatter(78.1, 30.1, s=20.0, color='k', marker='o', transform=ccrs.PlateCarree())
        cf = ax.contourf(to_np(lons1), to_np(lats1), to_np(var[day]), levels,
                         transform=ccrs.PlateCarree(),
                         cmap=get_cmap(cmap), extend='both')

        cb = plt.colorbar(cf, ax=ax, shrink=.90)
        cb.set_label(cbar_label)
    if domain == 2:
        ax = plt.axes(projection = cart_proj2)
        plot_background1(ax)
        ax.plot([77,82], [28, 28], color='k', transform=ccrs.PlateCarree())
        ax.plot([77,82], [32, 32], color='k', transform=ccrs.PlateCarree())
        ax.plot([77,77], [28, 32], color='k', transform=ccrs.PlateCarree())
        ax.plot([82,82], [28, 32], color='k', transform=ccrs.PlateCarree())
        ax.scatter(78.1, 30.1, s=40.0, color='k', marker='o', transform=ccrs.PlateCarree())
        cf = ax.contourf(to_np(lons2), to_np(lats2), to_np(var[day]), levels,
                         transform=ccrs.PlateCarree(),
                         cmap=get_cmap(cmap), extend='both')

        cb = plt.colorbar(cf, ax=ax, shrink=.90)
        cb.set_label(cbar_label)
    else:
        print("wrong value of domain input")

def plot_figure(var, cmap, levels, cbar_label, domain=1, vmin=None, vmax=None):
    if domain == 1:
        fig, axarr = plt.subplots(nrows=3, ncols=4, figsize=(18, 15), constrained_layout=True, sharex=True, sharey=True, subplot_kw={'projection': cart_proj1})
        axlist = axarr.flatten()
        for ax in axlist[0:-1]:
            plot_background1(ax)
        for i in range(0, 11):
            # c = axlist[i].contour(to_np(lons), to_np(lats), to_np(var[i]), colors="black", linewidths=0.005,
            #              transform=ccrs.PlateCarree())
            axlist[i].plot([77,82], [28, 28], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([77,82], [32, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([77,77], [28, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([82,82], [28, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].scatter(78.1, 30.1, s=20.0, color='k', marker='o', transform=ccrs.PlateCarree())
            if vmin == None and vmax == None:
                cf = axlist[i].contourf(to_np(lons1), to_np(lats1), to_np(var[i]), levels,
                         transform=ccrs.PlateCarree(),
                         cmap=get_cmap(cmap), vmin=np.round(var.min(), -1), vmax=np.round(var.max(), -1), extend='both')
            else:
                cf = axlist[i].contourf(to_np(lons1), to_np(lats1), to_np(var[i]), levels,
                         transform=ccrs.PlateCarree(),
                         cmap=get_cmap(cmap), vmin=vmin, vmax=vmax, extend='both')

            # Add a color bar
            # cb = plt.colorbar(cf, ax=axlist[i], shrink=.90)
            # cb.set_label('SAT (K)', size='x-large')
            axlist[i].set_title(str(10+i) + "th June 2013")

        # fig.colorbar(cf, ax = axlist, orientation='horizontal', shrink=0.80, aspect=30, pad=0.03, label=cbar_label)
        if vmin == None and vmax == None:
            bounds = np.linspace(np.round(var.min(), -1), np.round(var.max(), -1), levels)
        else:
            bounds = np.linspace(vmin, vmax, levels)
        norm = mpl.colors.BoundaryNorm(bounds, get_cmap(cmap).N, extend='both')
        # norm = mpl.colors.Normalize(vmin=var.min(), vmax=var.max())
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=get_cmap(cmap)), ax=axlist, orientation='horizontal', shrink=0.80, aspect=30, pad=0.03, label=cbar_label)
        fig.delaxes(axlist[-1])
        fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.08, hspace=0., wspace=0.)
    elif domain == 2:
        fig, axarr = plt.subplots(nrows=3, ncols=4, figsize=(18, 15), constrained_layout=True, sharex=True, sharey=True, subplot_kw={'projection': cart_proj2})
        axlist = axarr.flatten()
        for ax in axlist[0:-1]:
            plot_background2(ax)
        for i in range(0, 11):
            # c = axlist[i].contour(to_np(lons), to_np(lats), to_np(var[i]), colors="black", linewidths=0.005,
            #              transform=ccrs.PlateCarree())
            axlist[i].plot([77,82], [28, 28], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([77,82], [32, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([77,77], [28, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([82,82], [28, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].scatter(78.1, 30.1, s=40.0, color='k', marker='o', transform=ccrs.PlateCarree())
            if vmin == None and vmax == None:
                cf = axlist[i].contourf(to_np(lons2), to_np(lats2), to_np(var[i]), levels,
                         transform=ccrs.PlateCarree(),
                         cmap=get_cmap(cmap), vmin=np.round(var.min(), -1), vmax=np.round(var.max(), -1), extend='both')
            else:
                cf = axlist[i].contourf(to_np(lons2), to_np(lats2), to_np(var[i]), levels,
                         transform=ccrs.PlateCarree(),
                         cmap=get_cmap(cmap), vmin=vmin, vmax=vmax, extend='both')

            # Add a color bar
            # cb = plt.colorbar(cf, ax=axlist[i], shrink=.90)
            # cb.set_label('SAT (K)', size='x-large')
            axlist[i].set_title(str(10+i) + "th June 2013")

        # fig.colorbar(cf, ax = axlist, orientation='horizontal', shrink=0.80, aspect=30, pad=0.03, label=cbar_label)
        if vmin == None and vmax == None:
            bounds = np.linspace(np.round(var.min(), -1), np.round(var.max(), -1), levels)
        else:
            bounds = np.linspace(vmin, vmax, levels)
        norm = mpl.colors.BoundaryNorm(bounds, get_cmap(cmap).N, extend='both')
        # norm = mpl.colors.Normalize(vmin=var.min(), vmax=var.max())
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=get_cmap(cmap)), ax=axlist, orientation='horizontal', shrink=0.80, aspect=30, pad=0.03, label=cbar_label)
        fig.delaxes(axlist[-1])
        fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.08, hspace=0., wspace=0.)
    else:
        print("wrong value of domain input")

def plot_wind(var, uwnd, vwnd, cmap, levels, cbar_label, method=None, domain=1, vmin=None, vmax=None):
    if domain == 1:
        fig, axarr = plt.subplots(nrows=3, ncols=4, figsize=(18, 15), constrained_layout=True, sharex=True, sharey=True, subplot_kw={'projection': cart_proj1})
        axlist = axarr.flatten()
        for ax in axlist[0:-1]:
            plot_background1(ax)
        for i in range(0, 11):
#             c = axlist[i].contour(to_np(lons), to_np(lats), to_np(var[i]), \
#                                   colors="black", linewidths=0.005, \
#                                   transform=ccrs.PlateCarree())
            axlist[i].plot([77,82], [28, 28], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([77,82], [32, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([77,77], [28, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([82,82], [28, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].scatter(78.1, 30.1, s=20.0, color='k', marker='o', transform=ccrs.PlateCarree())
            if vmin == None and vmax == None:
                cf = axlist[i].contourf(to_np(lons1), to_np(lats1), to_np(var[i]), levels, \
                                    transform=ccrs.PlateCarree(),cmap=get_cmap(cmap), \
                                    vmin=np.round(var.min(),-1), vmax=np.round(var.max(),-1), \
                                    extend='both')
            else:
                cf = axlist[i].contourf(to_np(lons1), to_np(lats1), to_np(var[i]), levels, \
                                    transform=ccrs.PlateCarree(),cmap=get_cmap(cmap), \
                                    vmin=vmin, vmax=vmax, \
                                    extend='both')

            if method == 'stream':
                axlist[i].streamplot(to_np(lons1), to_np(lats1), to_np(uwnd[i]), \
                                     to_np(vwnd[i]), transform=ccrs.PlateCarree(), \
                                     linewidth=0.5, color='k')
            else:
                q = axlist[i].quiver(to_np(lons1)[::7, ::7], to_np(lats1)[::7, ::7], \
                                 to_np(uwnd[i])[::7, ::7], to_np(vwnd[i])[::7, ::7], \
                                 pivot='middle', transform=ccrs.PlateCarree(), scale=200, \
                                    units='width')
                qk = axlist[i].quiverkey(q, 0.05, 0.05, 20, '20 m/s', labelpos='W', \
                                    coordinates='figure', fontproperties={'family':'serif', 'size':13})
            # Add a color bar
            # cb = plt.colorbar(cf, ax=axlist[i], shrink=.90)
            # cb.set_label('SAT (K)', size='x-large')
            axlist[i].set_title(str(10+i) + "th June 2013")

        if vmin == None and vmax == None:
            bounds = np.linspace(np.round(var.min(), -1), np.round(var.max(), -1), levels)
        else:
            bounds = np.linspace(vmin, vmax, levels)
        norm = mpl.colors.BoundaryNorm(bounds, get_cmap(cmap).N, extend='both')
        # norm = mpl.colors.Normalize(vmin=var.min(), vmax=var.max())
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=get_cmap(cmap)), ax=axlist, \
                     orientation='horizontal', shrink=0.80, aspect=30, pad=0.03, label=cbar_label)
        fig.delaxes(axlist[-1])
        fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.08, hspace=0., wspace=0.)
    elif domain == 2:
        fig, axarr = plt.subplots(nrows=3, ncols=4, figsize=(18, 15), constrained_layout=True, sharex=True, sharey=True, subplot_kw={'projection': cart_proj2})
        axlist = axarr.flatten()
        for ax in axlist[0:-1]:
            plot_background2(ax)
        for i in range(0, 11):
#             c = axlist[i].contour(to_np(lons), to_np(lats), to_np(var[i]), \
#                                   colors="black", linewidths=0.005, \
#                                   transform=ccrs.PlateCarree())
            axlist[i].plot([77,82], [28, 28], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([77,82], [32, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([77,77], [28, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].plot([82,82], [28, 32], color='k', transform=ccrs.PlateCarree())
            axlist[i].scatter(78.1, 30.1, s=40.0, color='k', marker='o', transform=ccrs.PlateCarree())
            if vmin == None and vmax == None:
                cf = axlist[i].contourf(to_np(lons2), to_np(lats2), to_np(var[i]), levels, \
                                    transform=ccrs.PlateCarree(),cmap=get_cmap(cmap), \
                                    vmin=np.round(var.min(),-1), vmax=np.round(var.max(),-1), \
                                    extend='both')
            else:
                cf = axlist[i].contourf(to_np(lons2), to_np(lats2), to_np(var[i]), levels, \
                                    transform=ccrs.PlateCarree(),cmap=get_cmap(cmap), \
                                    vmin=vmin, vmax=vmax, \
                                    extend='both')

            if method == 'stream':
                axlist[i].streamplot(to_np(lons2), to_np(lats2), to_np(uwnd[i]), \
                                     to_np(vwnd[i]), transform=ccrs.PlateCarree(), \
                                     linewidth=0.5, color='k')
            else:
                q = axlist[i].quiver(to_np(lons2)[::7, ::7], to_np(lats2)[::7, ::7], \
                                 to_np(uwnd[i])[::7, ::7], to_np(vwnd[i])[::7, ::7], \
                                 pivot='middle', transform=ccrs.PlateCarree(), scale=200, \
                                    units='width')
                qk = axlist[i].quiverkey(q, 0.05, 0.05, 20, '20 m/s', labelpos='W', \
                                    coordinates='figure', fontproperties={'family':'serif', 'size':13})
            # Add a color bar
            # cb = plt.colorbar(cf, ax=axlist[i], shrink=.90)
            # cb.set_label('SAT (K)', size='x-large')
            axlist[i].set_title(str(10+i) + "th June 2013")

        if vmin == None and vmax == None:
            bounds = np.linspace(np.round(var.min(), -1), np.round(var.max(), -1), levels)
        else:
            bounds = np.linspace(vmin, vmax, levels)
        norm = mpl.colors.BoundaryNorm(bounds, get_cmap(cmap).N, extend='both')
        # norm = mpl.colors.Normalize(vmin=var.min(), vmax=var.max())
        fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=get_cmap(cmap)), ax=axlist, \
                     orientation='horizontal', shrink=0.80, aspect=30, pad=0.03, label=cbar_label)
        fig.delaxes(axlist[-1])
        fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.08, hspace=0., wspace=0.)
    else:
        print("wrong value of domain input")
