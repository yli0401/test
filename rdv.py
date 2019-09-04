# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.
"""

import time
from datetime import datetime
from threading import Timer 
from tkinter import *
import winsound

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.common.exceptions import NoSuchElementException
from selenium.common.exceptions import TimeoutException
 
chromeOptions = Options()
chromeOptions.add_experimental_option('useAutomationExtension', False)
# driver = webdriver.Chrome(executable_path='C:/Users/yli/anaconda3/Scripts/chromedriver.exe')
# driver = webdriver.Chrome('http://www.essonne.gouv.fr/booking/create/14056/0',chrome_options=chrome_options) 

def play_alert():
    count = 0
    while count<20:
        winsound.PlaySound('alert', winsound.SND_ASYNC)
        count += 1
        time.sleep(3)

def google_test():
    driver = webdriver.Chrome(chrome_options=chromeOptions, desired_capabilities=chromeOptions.to_capabilities())
    driver.get('http://www.google.fr')
    driver.quit()
    

def possible():
    print(datetime.now())    
    driver = webdriver.Chrome(chrome_options=chromeOptions, desired_capabilities=chromeOptions.to_capabilities())
    try:
        driver.get('http://www.essonne.gouv.fr/booking/create/14056')
        # driver.get('http://www.google.fr')
        time.sleep(1)
        last_height = driver.execute_script("return document.body.scrollHeight")
        driver.execute_script("window.scrollTo(0, document.body.scrollHeight);")
        time.sleep(1)
        driver.find_element_by_name("condition").click()
        time.sleep(1)
        driver.find_element_by_name("nextButton").click()
        time.sleep(2)
        driver.find_element_by_id("planning21029").click()
        time.sleep(2)
        driver.find_element_by_name("nextButton").click()
        driver.find_element_by_id("FormBookingCreate")
        if "Il n'existe plus de plage horaire libre pour votre demande de rendez-vous." in driver.page_source:
            print("No available place")
            time.sleep(2)
            driver.quit()
        else:
            play_alert()
            root = Tk()
            root.geometry("800x100+500+500")
            label = Label(root,
                         text = "Yes, we get the RDV!",
                         fg = "black",
                         font = ("times",20),
                         width = 800,
                         height = 100,
                         anchor = "c",
                         justify = "right"
                         )
            label.pack()
            root.attributes("-topmost", True)
            root.mainloop() 
    except (NoSuchElementException, TimeoutException) as er:
        print("Find the error.", er)
        driver.quit()



class MyTimer( object ):
    def __init__( self, start_time, interval, callback_proc, args=None, kwargs=None ):
        self.__timer = None
        self.__start_time = start_time
        self.__interval = interval
        self.__callback_pro = callback_proc
        self.__args = args if args is not None else []
        self.__kwargs = kwargs if kwargs is not None else {}
        
    def exec_callback( self, args=None, kwargs=None ):
        self.__callback_pro( *self.__args, **self.__kwargs )
        self.__timer = Timer( self.__interval, self.exec_callback )
        self.__timer.start()
        
    def start( self ):
        interval = self.__interval - ( datetime.now().timestamp() - self.__start_time.timestamp() )
        self.__timer = Timer( interval, self.exec_callback )
        self.__timer.start()
        
    def cancel( self ):
        self.__timer.cancel() 
        self.__timer = None



if __name__=="__main__":
    start = datetime.now()
    possible()
    tmr = MyTimer( start, 60*5, possible)
    tmr.start()


