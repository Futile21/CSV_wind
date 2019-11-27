import os
import xarray as xr
from netCDF4 import Dataset  # dataset allows us to read the netcdf4 files
import numpy as np
import datetime
import pandas as pd

import bisect

import os
from pathlib import Path
import dash
import plotly.graph_objs as go
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output,State
import pandas as pd
import base64
import pandas as pd
# import lorem
import numpy as np
from scipy.stats import norm
import json
from scipy.stats import norm
import copy
import xarray as xr



lat =[]
lon=[]
data=[]
time=[]

HGT=[50,80,100,150]
speed=['A<','B<','C<','D<','E<']

DATA = []


DATA_DIR = r"C:\Users\futil\OneDrive\GIZ\Internship\seb_test_ring\CVSmaker\smallData"
filelist = os.listdir(DATA_DIR)
files_nc=[DATA_DIR+"\\"+i for i in filelist if i.endswith('.nc')]

print (len(files_nc))


for Net in files_nc:
    root = Dataset(Net)
    lat.append(float(root.variables["XLAT"][0]))
    lon.append(float(root.variables["XLONG"][0]))

    SPDarray = []

    for i in np.arange(0, 4):
        WindDri = np.array(root["WDRO"][:, i:i + 1]).flatten()
        WindSpeed = np.array(root["WSPD"][:, i:i + 1]).flatten()
        windLen = len(WindDri)
#         print (len(WindDri))
#         print (len(WindSpeed))
#         print(f'the cont g is {contg}')
#         print(f'the cont l is {i+1}')


        N = []; NNW = []; NW = []; WNW = []
        W = []; WSW = []; SW = []; SSW = []
        S = []; SSE = []; SE = []; ESE = []
        E = []; ENE = []; NE = []; NNE = []

        cont = 0

        for i in WindDri:
            if 360 - 11.25 <= i or i < 11.25 + 22.5 * 0:
                N.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 1:
                NNW.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 2:
                NW.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 3:
                WNW.append(WindSpeed[cont])

            elif i < 11.25 + 22.5 * 4:
                W.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 5:
                WSW.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 6:
                SW.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 7:
                SSW.append(WindSpeed[cont])


            elif i < 11.25 + 22.5 * 8:
                S.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 9:
                SSE.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 10:
                SE.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 11:
                ESE.append(WindSpeed[cont])


            elif i < 11.25 + 22.5 * 12:
                E.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 13:
                ENE.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 14:
                NE.append(WindSpeed[cont])
            elif i < 11.25 + 22.5 * 15:
                NNE.append(WindSpeed[cont])

            cont += 1
        N = sorted(N)
        NNW = sorted(NNW)
        NW = sorted(NW)
        WNW = sorted(WNW)

        W = sorted(W)
        WSW = sorted(WSW)
        SW = sorted(SW)
        SSW = sorted(SSW)

        S = sorted(S)
        SSE = sorted(SSE)
        SE = sorted(SE)
        ESE = sorted(ESE)

        E = sorted(E)
        ENE = sorted(ENE)
        NE = sorted(NE)
        NNE = sorted(NNE)

        a = 2.5    ###### Change
        b = 5      ###### Change
        c = 7.5    ###### Change
        d = 10     ###### Change
        e = 12.5   ###### Change

        #print(f"The len of N is {len(N)}")
        NspdA = len(N) - 1 - bisect.bisect_left(N, a)  # length start at 1 but bisect returns index
        NspdB = len(N) - 1 - bisect.bisect_left(N, b)
        NspdC = len(N) - 1 - bisect.bisect_left(N, c)
        NspdD = len(N) - 1 - bisect.bisect_left(N, d)
        NspdE = len(N) - 1 - bisect.bisect_left(N, e)

        #print(f"The len of N is {len(NNW)}")
        NNWspdA = len(NNW) - 1 - bisect.bisect_left(NNW, a)  # length start at 1 but bisect returns index
        NNWspdB = len(NNW) - 1 - bisect.bisect_left(NNW, b)
        NNWspdC = len(NNW) - 1 - bisect.bisect_left(NNW, c)
        NNWspdD = len(NNW) - 1 - bisect.bisect_left(NNW, d)
        NNWspdE = len(NNW) - 1 - bisect.bisect_left(NNW, e)

        #print(f"The len of N is {len(NW)}")
        NWspdA = len(NW) - 1 - bisect.bisect_left(NW, a)  # length start at 1 but bisect returns index
        NWspdB = len(NW) - 1 - bisect.bisect_left(NW, b)
        NWspdC = len(NW) - 1 - bisect.bisect_left(NW, c)
        NWspdD = len(NW) - 1 - bisect.bisect_left(NW, d)
        NWspdE = len(NW) - 1 - bisect.bisect_left(NW, e)

        #print(f"The len of N is {len(WNW)}")
        WNWspdA = len(WNW) - 1 - bisect.bisect_left(WNW, a)  # length start at 1 but bisect returns index
        WNWspdB = len(WNW) - 1 - bisect.bisect_left(WNW, b)
        WNWspdC = len(WNW) - 1 - bisect.bisect_left(WNW, c)
        WNWspdD = len(WNW) - 1 - bisect.bisect_left(WNW, d)
        WNWspdE = len(WNW) - 1 - bisect.bisect_left(WNW, e)

        #print(f"The len of W is {len(W)}")
        WspdA = len(W) - 1 - bisect.bisect_left(W, a)  # length start at 1 but bisect returns index
        WspdB = len(W) - 1 - bisect.bisect_left(W, b)
        WspdC = len(W) - 1 - bisect.bisect_left(W, c)
        WspdD = len(W) - 1 - bisect.bisect_left(W, d)
        WspdE = len(W) - 1 - bisect.bisect_left(W, e)

        #print(f"The len of W is {len(WSW)}")
        WSWspdA = len(WSW) - 1 - bisect.bisect_left(WSW, a)  # length start at 1 but bisect returns index
        WSWspdB = len(WSW) - 1 - bisect.bisect_left(WSW, b)
        WSWspdC = len(WSW) - 1 - bisect.bisect_left(WSW, c)
        WSWspdD = len(WSW) - 1 - bisect.bisect_left(WSW, d)
        WSWspdE = len(WSW) - 1 - bisect.bisect_left(WSW, e)

        #print(f"The len of W is {len(SW)}")
        SWspdA = len(SW) - 1 - bisect.bisect_left(SW, a)  # length start at 1 but bisect returns index
        SWspdB = len(SW) - 1 - bisect.bisect_left(SW, b)
        SWspdC = len(SW) - 1 - bisect.bisect_left(SW, c)
        SWspdD = len(SW) - 1 - bisect.bisect_left(SW, d)
        SWspdE = len(SW) - 1 - bisect.bisect_left(SW, e)

        #print(f"The len of W is {len(SSW)}")
        SSWspdA = len(SSW) - 1 - bisect.bisect_left(SSW, a)  # length start at 1 but bisect returns index
        SSWspdB = len(SSW) - 1 - bisect.bisect_left(SSW, b)
        SSWspdC = len(SSW) - 1 - bisect.bisect_left(SSW, c)
        SSWspdD = len(SSW) - 1 - bisect.bisect_left(SSW, d)
        SSWspdE = len(SSW) - 1 - bisect.bisect_left(SSW, e)

        #print(f"The len of S is {len(S)}")
        SspdA = len(S) - 1 - bisect.bisect_left(S, a)  # length start at 1 but bisect returns index
        SspdB = len(S) - 1 - bisect.bisect_left(S, b)
        SspdC = len(S) - 1 - bisect.bisect_left(S, c)
        SspdD = len(S) - 1 - bisect.bisect_left(S, d)
        SspdE = len(S) - 1 - bisect.bisect_left(S, e)

        #print(f"The len of S is {len(SSE)}")
        SSEspdA = len(SSE) - 1 - bisect.bisect_left(SSE, a)  # length start at 1 but bisect returns index
        SSEspdB = len(SSE) - 1 - bisect.bisect_left(SSE, b)
        SSEspdC = len(SSE) - 1 - bisect.bisect_left(SSE, c)
        SSEspdD = len(SSE) - 1 - bisect.bisect_left(SSE, d)
        SSEspdE = len(SSE) - 1 - bisect.bisect_left(SSE, e)

        #print(f"The len of S is {len(SE)}")
        SEspdA = len(SE) - 1 - bisect.bisect_left(SE, a)  # length start at 1 but bisect returns index
        SEspdB = len(SE) - 1 - bisect.bisect_left(SE, b)
        SEspdC = len(SE) - 1 - bisect.bisect_left(SE, c)
        SEspdD = len(SE) - 1 - bisect.bisect_left(SE, d)
        SEspdE = len(SE) - 1 - bisect.bisect_left(SE, e)

        #print(f"The len of S is {len(ESE)}")
        ESEspdA = len(ESE) - 1 - bisect.bisect_left(ESE, a)  # length start at 1 but bisect returns index
        ESEspdB = len(ESE) - 1 - bisect.bisect_left(ESE, b)
        ESEspdC = len(ESE) - 1 - bisect.bisect_left(ESE, c)
        ESEspdD = len(ESE) - 1 - bisect.bisect_left(ESE, d)
        ESEspdE = len(ESE) - 1 - bisect.bisect_left(ESE, e)

        #print(f"The len of E is {len(E)}")
        EspdA = len(E) - 1 - bisect.bisect_left(E, a)  # length start at 1 but bisect returns index
        EspdB = len(E) - 1 - bisect.bisect_left(E, b)
        EspdC = len(E) - 1 - bisect.bisect_left(E, c)
        EspdD = len(E) - 1 - bisect.bisect_left(E, d)
        EspdE = len(E) - 1 - bisect.bisect_left(E, e)

        #print(f"The len of E is {len(ENE)}")
        ENEspdA = len(ENE) - 1 - bisect.bisect_left(ENE, a)  # length start at 1 but bisect returns index
        ENEspdB = len(ENE) - 1 - bisect.bisect_left(ENE, b)
        ENEspdC = len(ENE) - 1 - bisect.bisect_left(ENE, c)
        ENEspdD = len(ENE) - 1 - bisect.bisect_left(ENE, d)
        ENEspdE = len(ENE) - 1 - bisect.bisect_left(ENE, e)

        #print(f"The len of E is {len(NE)}")
        NEspdA = len(NE) - 1 - bisect.bisect_left(NE, a)  # length start at 1 but bisect returns index
        NEspdB = len(NE) - 1 - bisect.bisect_left(NE, b)
        NEspdC = len(NE) - 1 - bisect.bisect_left(NE, c)
        NEspdD = len(NE) - 1 - bisect.bisect_left(NE, d)
        NEspdE = len(NE) - 1 - bisect.bisect_left(NE, e)

        #print(f"The len of E is {len(NNE)}")
        NNEspdA = len(NNE) - 1 - bisect.bisect_left(NNE, a)  # length start at 1 but bisect returns index
        NNEspdB = len(NNE) - 1 - bisect.bisect_left(NNE, b)
        NNEspdC = len(NNE) - 1 - bisect.bisect_left(NNE, c)
        NNEspdD = len(NNE) - 1 - bisect.bisect_left(NNE, d)
        NNEspdE = len(NNE) - 1 - bisect.bisect_left(NNE, e)

        SPDA = np.array(
            [NspdA, NNWspdA, NWspdA, WNWspdA, WspdA, WSWspdA, SWspdA, SSWspdA, SspdA, SSEspdA, SEspdA, ESEspdA, EspdA,
             ENEspdA, NEspdA, NNEspdA, ]) / windLen * 100
        SPDB = np.array(
            [NspdB, NNWspdB, NWspdB, WNWspdB, WspdB, WSWspdB, SWspdB, SSWspdB, SspdB, SSEspdB, SEspdB, ESEspdB, EspdB,
             ENEspdB, NEspdB, NNEspdB, ]) / windLen * 100
        SPDC = np.array(
            [NspdC, NNWspdC, NWspdC, WNWspdC, WspdC, WSWspdC, SWspdC, SSWspdC, SspdC, SSEspdC, SEspdC, ESEspdC, EspdC,
             ENEspdC, NEspdC, NNEspdC, ]) / windLen * 100
        SPDD = np.array(
            [NspdD, NNWspdD, NWspdD, WNWspdD, WspdD, WSWspdD, SWspdD, SSWspdD, SspdD, SSEspdD, SEspdD, ESEspdD, EspdD,
             ENEspdD, NEspdD, NNEspdD, ]) / windLen * 100
        SPDE = np.array(
            [NspdE, NNWspdE, NWspdE, WNWspdE, WspdE, WSWspdE, SWspdE, SSWspdE, SspdE, SSEspdE, SEspdE, ESEspdE, EspdE,
             ENEspdE, NEspdE, NNEspdE, ]) / windLen * 100

        SPD=np.array([SPDA,SPDB,SPDC,SPDD,SPDE])
        SPDarray.append(SPD)


    DATA.append(SPDarray)





print(np.shape(np.array(DATA)))












