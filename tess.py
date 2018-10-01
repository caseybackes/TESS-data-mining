import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits 
from requests import get
from requests.exceptions import RequestException
from contextlib  import closing
from bs4 import BeautifulSoup
import argparse
import requests


## Web scrapping Routines to download data from URLs
def simple_get(url):
    """
    Attempts to get the content at `url` by making an HTTP GET request.
    If the content-type of response is some kind of HTML/XML, return the
    text content, otherwise return None
    """
    try:
        with closing(get(url, stream=True)) as resp:
            if is_good_response(resp):
                return resp.content
            else:
                return None
    except RequestException as e:
        log_error('Error during requests to {0} : {1}'.format(url, str(e)))
        return None


def is_good_response(resp):
    """
    Returns true if the response seems to be HTML, false otherwise
    """
    content_type = resp.headers['Content-Type'].lower()
    return (resp.status_code == 200 
            and content_type is not None 
            and content_type.find('html') > -1)


def log_error(e):
    """
    It is always a good idea to log errors. 
    This function just prints them, but you can
    make it do anything.
    """
    print(e)

##  Construct a URL that exists on MAST database and use BeautifulSoup4 to download the HTML
##  content at that URL and restructure it into array objects.  
def fetch_lc_list(tid):
    """
    args = ()
    MAST Directory Structure:
        <TID[0:1]>/<TID[2:4]>/<TID[5:7]>/<TID[8:10]>/
    Product:
        Light Curve file
    File Name Convention:
        tess{date-time}-{tid}-{scid}-{cr}_lc.fits	
        {date-time} = The timestamp associated with this file, in the yyyydddhhmmss format.
        {tid} = A zero-padded, 16-digit target identifier that refers to an object in the TESS Input Catalog.
        {scid} = A zero-padded, four-digit identifier of the spacecraft configuration map used to process this data.
        {cr} = A string character that denotes the cosmic ray mitigation procedure. Possible values are:
            'x': No mitigation performed at the SPOC.
            's': Mitigation performed on the spacecraft.
            'a': A SPOC mitigation algorithm was used.
            'b': Both a SPOC and onboard spececraft algorithm was used.
    Description:
        Calibrated, extracted light curves.
    Example File URL:
        https://archive.stsci.edu/missions/tess/ete-6/tid/00/000/003/105/tess2019128220341-0000000310506508-0016-s_lc.fits

    This has an example TID of 3105 ->                                      ^^^^^^
    """
    tid = '0'*(11-len(str(tid)))+str(tid)
    tidargs =  (tid[0:2], tid[2:5] , tid[5:8] ,tid[8:12])
    tid_path = 'https://archive.stsci.edu/missions/tess/ete-6/tid/%s/%s/%s/%s/' % (tidargs[0],tidargs[1],tidargs[2],tidargs[3])
    lc_array=list()
    request = requests.get(tid_path)

    if request.status_code == 200:
        html = BeautifulSoup(simple_get(tid_path), 'html.parser')
        print "Defining unique URL paths for TID %s..." %(tid)
        for a in html.select('a'):
            if "_lc.fits" in str(a):
                lc_array.append(tid_path +str(a).split('"')[1])
        print "Successfully defined %s URL paths."% len(lc_array)
    if request.status_code == 404: 
        print "URL does not exist for TID %s" % (tid)
        exit()
    return lc_array


## Define a class of LightCurve object wth observation data and metadata 
class LightCurve:
    def __init__(self, urlpath):
        with fits.open(urlpath) as f:
            self.hdu_data = f[1].data 
            self.hdu_header = f[1].header
            self.ticid =self.hdu_header['OBJECT']
            self.ra = self.hdu_header['RA_OBJ']
            self.dec = self.hdu_header['DEC_OBJ']
            self.date_obs = self.hdu_header['DATE-OBS']
            self.time=np.array(self.hdu_data["TIME"])
            self.flux=np.array(self.hdu_data["PDCSAP_FLUX"])
            self.fluxmean  =np.array( self.hdu_data['PDCSAP_FLUX'])/(np.nanmean(self.hdu_data["PDCSAP_FLUX"]))
            self.fluxmean2=np.array(self.hdu_data['PDCSAP_FLUX'])/(np.nanmean(self.hdu_data["PDCSAP_FLUX"]))
            self.data    = np.array( self.hdu_data['PDCSAP_FLUX'])/(np.nanmean(self.hdu_data["PDCSAP_FLUX"]))
            #self.tid =lc.split('-')[2])


# def linear_model(a_p,R_planet,R_star,dz):
def linear_model():
    '''
    NON LIMB DARKENING MODEL
    a_p is the semimajor axis at periastron in km
    R_star is the stellar radius in km
    R_planet is the planet's radius in km
    dz is the step size in impact parameter (km)
    '''
    a_p = 1.5e8 # 1 AU periastron orbital distance
    dz = 100

    R_star = 695508.0 #R_sun
    R_planet = 69911.0 # R_jupiter
    d_sep = abs(np.arange(-a_p/10, a_p/10, dz))#np.arange(-a_p/10.0,a_p/10.0,dz+1.0) # the dz+1 avoids a zero in this array
    impact_param = d_sep / R_star
    flux = []
    p =R_planet/R_star 
    # for earth-sun system, this is about 0.00935
    for z in impact_param:

       #OUT OF TRANSIT (NOMINAL) FLUX OF THE SYSTEM (planet is not transiting)
        if (1.0+p) < z:  # 1 + Rp/R* less than the distance from 
            #print "(1 + p) < ",z
            obscured_flux = 0.0

        # INGRESS/EGRESS OF PLANET ENTERING/EXITING THE TRANSIT
        if abs(1.0-p) < z and z <= (1.0+p):
            k_1 = np.arccos((1.0-p**2.0 + z**2.0)/(2.0*z))
            k_0 = np.arccos((p**2.0 + z**2.0 -1.0)/(2*p*z))
            obscured_flux = 1.0/np.pi * (p**2.0*k_0 + k_1 - np.sqrt(   (4.0*z**2.0-(1.0+z**2.0-p**2.0)**2.0) /4.00 ))
            #print "obscured flux is %s" % obscured_flux

        #  FLUX DURING FULL ECLIPSE
        if z <= (1.0-p):
            obscured_flux = p**2.0

        # ADD THIS MODEL POINT DURING DZ TO THE FLUX ARRAY
        flux.append(1.0- obscured_flux)

    plt.plot( flux)
    plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-tid', type= int, help="TID of TESS object")
    parser.add_argument('-a_p', type = int, help="Exoplanet periastron dist in km")
    parser.add_argument('-fl', type = bool, help = 'fit linear model to TID' )
    args = parser.parse_args()
    
    tid =args.tid
    
      
    # pull a TID data set and measure the fit of a linear model to the simulated data
    TID_urllist = fetch_lc_list(tid)
    
    lc1 = LightCurve(TID_urllist[0])
    plt.plot(lc1.data)
    plt.title(str('TID '  +str(tid) ) )
    plt.show()


if __name__ == '__main__':
    main()