# -*- coding: utf-8 -*-
from __future__ import unicode_literals

# Import PyQt5 classes
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWebKit import *
from PyQt5.QtNetwork import *
from PyQt5.QtWidgets import *
from PyQt5.QtWebKitWidgets import *
from PyQt5.QtPrintSupport import *

from plugins import VisualisationPlugin

from collections import defaultdict, OrderedDict

import os, re, time
from copy import copy

import numpy as np
from sklearn.decomposition import PCA

import ui, db, utils, threads
from data import DataSet, DataDefinition
from views import D3PrerenderedView

from poster.encode import multipart_encode

import requests


class ZeitgeistApp( ui.AnalysisApp ):
    def __init__(self, **kwargs):
        super(ZeitgeistApp, self).__init__(**kwargs)

        # Define automatic mapping (settings will determine the route; allow manual tweaks later)
        
        self.addDataToolBar()
        self.addFigureToolBar()
        
        #self.data.add_input('input') # Add input slot
        self.views.addView( D3PrerenderedView(self), 'View')
        
        #self.table = QTableApp()
        #self.table.setModel(self.data.o['output'].as_table)
        #self.tabs.addTab(self.table,'Table')
        
        self.predefined_geist = {
            'Omics':['metabolomics','genomics','transcriptomics','proteomics','metabonomics'],
            'Arthritis':['RA','SLE','OA','JA','REA'],
            'Cytokines':['IL1','IL1alpha','IL1beta','IL2','IL3','IL4','IL5','IL6','IL10','IL10R','IL17','IL30','IFNgamma','TNFalpha'],
            'Viruses':['H1N1','H5N1', 'EBV', 'HIV','HPV','HSV','HPV5','HPV8']
        }

        self.config.set_defaults({
            'predefined_geist':'Omics',
        })
        
        
        t = self.addToolBar('Zeitgeist')
        t.zg_control = QComboBox()
        t.zg_control.addItems( [h for h in self.predefined_geist.keys()] )
        self.config.add_handler('predefined_geist', t.zg_control)
        t.addWidget(t.zg_control)
        
        self.toolbars['zeitgeist'] = t

        
        # Setup data consumer options
        self.data.consumer_defs.append( 
            DataDefinition('input', {
                'entities_t':(None, ['Metabolite','Pathway','Gene','Protein']),
            })
        )

        self.finalise()


    def batches(self, data, batch_size):
        for i in range(0, len(data), batch_size):
                yield data[i:i+batch_size]    
     
    # Generate Zeitgeist
    def generate(self, input=None):    
        dso = input
        # Build a list object of class, x, y
        #terms = ['IL1','IL-6','TNF-a','TNF-alpha','IL-10','IL-21','IL-12']
        #terms = ['glucose','lactate','pyruvate','succinate','maleate']
        #terms = ['metabolomics','genomics','transcriptomics','proteomics','metabonomics']
        terms = self.predefined_geist[ self.config.get('predefined_geist') ]
                
        figure_data = OrderedDict()
        
        # Pubmed alternative (one at a time; multi-request for each year?)
        # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&rettype=count&term=maltose%20[tiab]%20and%202011:2013%20[dp]&retmode=xml
        values = {
            # Raw data download doesn't pay attention to scaling, etc. boo!
            'db':'pubmed', # and abstracts',
            'retmode':'xml',
            'rettype':'count',
            }
        
        year_range = 5
        years = range(1960,2012,year_range)

        for year in years:
            figure_data[year] = []
        
        terms_n = float(len(terms))
        
        for n, term in enumerate(terms): # Iterate in 5s using a generator (query limit)
            previous_n = 0
            self.updateProgress( (n/terms_n)*100 )
            for year in years:
                values['term'] = "%s [tiab] %d:%d [ppdat]" % (term, year, year+year_range-1)
                
                n = self.plugin.get_cache_item(values)
                
                if not n:
                    time.sleep(0.5)
                    #?db=pubmed&retmode=xml&rettype=count&term=maltose%20[tiab]%20and%20glucose%20[tiab]%20and%202011:2013%20[dp]
                    r = requests.get('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi', params=values )
                    match = re.search('<Count>(\d+)</Count>', r.text, re.MULTILINE)
                    if match:
                        n = match.group(1)
                        self.plugin.put_cache_item(values, n)
                    else:
                        continue # next
                
                # We have an n!
                n = int(n)        
                n_delta_rel = (n-previous_n)/float(n) if n >0 else 0
                figure_data[year].append( (term, n) )  
                previous_n = n
                    
        metadata = {
            'figure':{
                'data':figure_data,
                },
        }

        return {'metadata': metadata }       

    def prerender(self, metadata):
        return {'View':{'metadata':metadata, 'template':'d3/stacked-area.svg'}}                
     
        

class Zeitgeist(VisualisationPlugin):

    def __init__(self, **kwargs):
        super(Zeitgeist, self).__init__(**kwargs)
        self.register_app_launcher( ZeitgeistApp )
        
        
