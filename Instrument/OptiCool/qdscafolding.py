#!/usr/bin/env python3
#
# qdsimulate.py
#

from enum import Enum, auto
from parse_inputs import instrumentType
    
# simple simulator, with no time dependence or dynamics
class QDInstrumentSim:
    def __init__(self, instrument_type):
        if instrument_type == instrumentType.DYNACOOL.name:
            self._class_id = 'QD.MULTIVU.DYNACOOL.1'
        elif instrument_type == instrumentType.PPMS.name:
            self._class_id = 'QD.MULTIVU.PPMS.1'
        elif instrument_type == instrumentType.VERSALAB.name:
            self._class_id = 'QD.MULTIVU.VERSALAB.1'
        elif instrument_type == instrumentType.MPMS3.name:
            self._class_id = 'QD.MULTIVU.MPMS3.1'
        elif instrument_type == instrumentType.OPTICOOL.name:
            self._class_id = 'QD.MULTIVU.OPTICOOL.1'
        else:
            raise Exception('Unrecognized instrument type: {0}.'.format(instrument_type))
        self.simCode = 1
        self.simMode = 1
        self.simTemp = 1
        self.simRate = 1
        self.simMode = 1
        self.simField = 1
        self.simRate = 1
        self.simApproach = 1
        self.simMode = 1  
        self.simCode = 1  
        
    def set_temperature(self, temperature, rate, mode):
        """Sets temperature and returns MultiVu error code"""
        self.simTemp = temperature
        self.simRate = rate
        self.simMode = mode
        return 0

    def get_temperature(self):
        """Gets and returns temperature info as (MultiVu error, temperature, status)"""
        return 0, self.simTemp, self.simMode

    def set_field(self, field, rate, approach, mode):
        """Sets field and returns MultiVu error code"""
        self.simField = field
        self.simRate = rate
        self.simApproach = approach
        self.simMode = mode
        return 0

    def get_field(self):
        """Gets and returns field info as (MultiVu error, field, status)"""
        return 0, self.simField, self.simMode

    def set_chamber(self, code):
        """Sets chamber and returns MultiVu error code"""
        self.simCode = code
        return 1

    def get_chamber(self):
        """Gets chamber status and returns (MultiVu error, status)"""
        return 0, self.simCode
