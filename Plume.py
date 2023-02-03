def xy_latlon(img_da, area, lat, lon): 
    """Returns the image pixel coordinates (x,y) from (lat,lon) coordinates. 
    
    Parameters
    ----------
    img_da: an image band as an xarray dataarray from the satpy
            Scene.to_xarray_dataset method.
    area: stapy area extent definition
    lat, lon: coordinates of the point of which data coordinates x,y must 
    be calculated
    
    Returns
    -------
    x,y: coordinates of the point in data coordinates (pixel_x, pixel_y)
    Usage
    -----
    plt_xy_latlon(img.IR_108, 'area', lat, lon)  or
    plt_xy_latlon(img['IR_108'], 'area', lat, lon)   
    Author
    ------
    LM 20210317"""
    
    x,y = img_da.attrs[area].get_xy_from_lonlat(lon, lat)
    return x,y

def get_AOI_from_xy(img_ds, x, y, radius=1): 
    """Returns the Area Of Interest xarray dataset as a square image 
    centered on image pixel of coordinates (x,y) and side=(2*radius+1). 
    
    AOI point (x,y) can be calculated from (lat,lon) coordinates with 
    xy_latlon method.
    
    Parameters
    ----------
    img_ds: an image band as an xarray dataarray from the satpy
            Scene.to_xarray_dataset method.
    lat, lon: coordinates of the point of which data coordinates x,y 
    must be calculated
    
    Returns
    -------
    AOI: xarray dataset subset of img_ds
    
    Usage
    -----
    get_AOI_from_xy(img, x, y [,radius=radius]) 
    
    Author
    ------
    LM 20210317"""

    import numpy as np

    x_range = np.arange(x-radius, x+radius+1)
    y_range = np.arange(y-radius, y+radius+1)
    #print('x_range =', x_range,'y_range =', y_range)
    AOI_ds = img_ds.isel(x=x_range, y=y_range)
    return AOI_ds

def   get_mintemp_from  (img,chn) : 
    """Returns the minimum temperature in the Area Of Interest xarray dataset centered on image img in channel chn.
    
    Parameters
    ----------    
    img: xarray dataset from the selected subset of satpy 
        Scene.to_xarray_dataset method or any satpy xarray dataarray.
    chn: str, channel on which is computed the minimum temperature
    
    Returns
    -------
    AOI: xarray dataset subset of img_ds
    
    Usage
    -----
    Output: minvalue (TODO: and index, if needed)
    
    Usage
    -----
    get_mintemp_from(img, chn) 
    
    Author
    ------
    LM 20210324"""
    
    import numpy as np
    
    return np.nanmin(img[chn].values)
    
class DarkestPixel:
    
    import numpy as np        
    from scipy import interpolate

    
    def __init__(self, BT, Z, T, hmax, dT=0.0, int_flag='first', 
                 ext_flag='near'):
        self.BT = BT
        self.Z = Z
        self.T = T
        self.hmax = hmax
        self.dT = dT
        self.int_flag = int_flag
        self.ext_flag = ext_flag
            
    def darkest_pixel(self):

        """
        INPUT:
            L = Radiance (scalare o anche vettore)
            WL = Wavelentgh
            Zpro = Profile altitude (m)
            Tpro = Profile Temperature (C)
            int_flag = 'first' or 'last' (quale quota prendo se c'e' più di un
            punto di intersezione)
            ext_flag = flag per punti dove non interseca
            if = 'extrap' estrapola linearmente
            if = 'near' prende il più vicino
            altrimenti mette NaN¶
            dT = delta di temperatura da aggiungere
            hmax = quota massima possibile in km (tropopausa); 
            se H>hmax --> H=hmax;
            
        OUTPUT:
            H = Quota TOP in km (scalare o vettore)
            BT = Brightness Temperature (C) + dT
            
        Inizializzo
            H=L*0; [Tprou,ia,~] = unique(Tpro,'stable')
            Planck
            BT=planck(WL,L,1)-273.15 + dT
            
        """
    
        
        """Inizializzo questo vettore Tprou (Tprofilounico) di valori unici di
        T che Lorenzo usa per l'estrapolazione lineare (occhio che uniq in
        matlab ha la keyword 'stable' che significa che i valori unici non sono
        sorted come invece è per np.unique)"""
        
        import numpy as np
        from scipy import interpolate
        Tprou, ia, ic = np.unique(self.T,return_index=True,return_inverse=True) 
                                  #, axis=0)
        # matrice BT di 2dim nxn trasformata in un vettore 1dim di nxn elementi
        BT_flat = self.BT.ravel()    
        # inizializzo con zeri un vettore di altezze H di nxn elementi
        H = np.zeros(len(BT_flat), dtype=np.float64)  
        
    #    print('T ==', T)
    #    print('Z ==', Z)
        #print('H ==', H)
        #print('int_flag, ext_flag, dT, len(BT_flat) = ', self.int_flag,  
              #self.ext_flag, self.dT, len(BT_flat))
    
        for i in range(len(BT_flat)):
            # intersezione di BT_flat e T, con indici di elementi uguali
            dp_eq, i_BT_flat, i_T = np.intersect1d(BT_flat[i], self.T, 
                                                   return_indices=True)      
            #print('dp_eq ==',dp_eq)
            # se BT_flat e T non hanno valori uguali --> bisogna interpolare
            if dp_eq.size == 0:
                # init lista z dark pixel che per ogni elemento BT_flat 
                #conterrà i valori interpolati su T
                zdp = []                              
                for p in range(self.Z.size-1):
                    t = self.T[p:p+2]
                    z = self.Z[p:p+2]
                     # np.interp vuole x crescenti
                    if self.T[p] > self.T[p+1]:                 
                        t = np.flip(self.T[p:p+2])
                        z = np.flip(self.Z[p:p+2])
                    if BT_flat[i] > t[0] and BT_flat[i] < t[1]:
                        # qui li prendo tutti TOP
                        zdp.append(np.interp(BT_flat[i], t, z)/1000)  
                        #print('i, p, zdp, BT[i], t, z ==', i, p, zdp[-1], 
                             # BT_flat[i], t, z)               
                if len(zdp)>0:
                    #print('      i, p, zdp, len(zdp)  ==', i, p, zdp, len(zdp))
                    if self.int_flag == 'first': 
                        H[i] = zdp[0]    
                    else: 
                        H[i] = zdp[-1]
                # 'bisogna estrapolare'
                else:                                    
                    #print('  --> bisogna estrapolare')
                    #print('i, p ==', i, p)
                    if self.ext_flag == 'extrap':
                        f = interpolate.interp1d(Tprou, self.Z[ia], 
                                                 kind='linear', fill_value='extrapolate')
                        #print('Tprou, Z[ia], BT[i], f(BT[i]) ==', Tprou, self.Z[ia], BT_flat[i], f(BT_flat[i]))
                        H[i] = f(BT_flat[i])/1000          # km TOP
                    elif self.ext_flag == 'near':
                        f = interpolate.interp1d(Tprou, self.Z[ia], kind='nearest', fill_value='extrapolate')    #km TOP
                        H[i] = f(BT_flat[i])/1000          # km TOP
                        #print('Tprou, Z[ia], BT[i], f(BT[i]) ==', Tprou, self.Z[ia], BT_flat[i], f(BT_flat[i]))
                    else:
                        H[i] = NaN
            else:                # BT_flat e T hanno valori uguali !!!
                if self.int_flag == 'first':
                    #print('i, H[i], hmax ==', i, H[i], self.hmax)
                    #print('i_BT_flat ==', i_BT_flat)
                    #print('len(i_T) ==', len(i_T))
                    #print('Z[i_T[0]] ==', self.Z[i_T[0]])
                    #print('T[80] ==', self.T[80])
                    H[i] = self.Z[i_T[0]]/1000
                else:
                    H[i] = self.Z[i_T[-1]]/1000
            if H[i]>self.hmax: H[i] = self.hmax         
            
            #print('      H[i] ==', H[i])
        return H