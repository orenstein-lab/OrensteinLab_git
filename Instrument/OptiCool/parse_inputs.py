# -*- coding: utf-8 -*-
"""
This class is a tool to parse the command-line inputs

"""

import sys
import re
from os import path
import ntpath
from enum import Enum, auto


class instrumentType(Enum):
    DYNACOOL = auto()
    PPMS = auto()
    VERSALAB = auto()
    MPMS3 = auto()
    OPTICOOL = auto()
    na = auto()

separator=" "

class inputs():
    def __init__(self, instrumentRequired=True, sep=separator):
        self.instrumentRequired = instrumentRequired
        self.separator = sep
            
    def path_leaf(self, path):
        head, tail = ntpath.split(path)
        return tail or ntpath.basename(head)
    
    def parseInput(self, inputArgs, sep=separator):
        """Arguments flags are:
            --help to display the help text
            --s to simulate data
            --ip=<host address> to specify the host IP address (default = 'localhost')
            
            An argument without a flag is the instrument
            
            This function returns a list of [instrument, simulate, host]"""
            
        
        # convert the inputArgs into a list separated by 'sep' and remove
        # leading or trailing white space.
        argLine = " ".join(inputArgs)
        args = argLine.split(sep)
        inputList = [reg.lstrip() for reg in args]
        # remove any empty items
        inputList = list(filter(None, inputList))
        
        # reg-ex for finding flags in the input
        helpArgs = re.compile("-[-]?(?i)(h)")
        simArgs = re.compile("-[-]?(?i)(s)")
        ipArgs = re.compile("-[-]?(?i)(ip)[ =]?")
        
        instrument = ''
        simulateMode = False
        host = "localhost"
        showHelp = False
        for item in inputList:
            # check for help string
            if helpArgs.search(item):
                showHelp = True
            # check if the filename input flag is used, and get the filename
            elif simArgs.search(item):
                simulateMode=True
            elif ipArgs.search(item):
                # remove the flag
                host = ipArgs.sub("",item)
                # remove any white space
                host = host.lstrip()
                
            else:
                # Add the input to the instrument name.  If user enters
                # more than one instrument, then this will cause an error
                # once outside of the for-loop
                instrument += item.upper()
        
        if self.instrumentRequired:
            # Check to see if the input is a valid instrument
            if instrument not in instrumentType._member_names_ \
                or instrument == instrumentType.na: 
                showHelp = True
        else:
            # Check to see if any input without a flag was entered
            if instrument != '':
                showHelp = True
        
        if showHelp:
                print(self.helpText())
                sys.exit()
                
        return [instrument, simulateMode, host]
    
    def helpText(self) -> str:
        programName = self.path_leaf(sys.argv[0])
        helpText = ("Input Options: \n")
        helpText += ("To display this help text:\n")
        helpText += (f"\tpython {programName} -h\n")
        helpText += ("To run the scaffolding and test the server:\n")
        helpText += (f"\tpython {programName} -s\n")
        helpText += ("To specify the host IP address:\n")
        helpText += (f"\tpython {programName} -ip=<host IP address> (default = 'localhost')\n")
        if self.instrumentRequired:
            helpText += (f"\nArguments without a flag are the instrument, and only one instrument is allowed.\n")
            helpText += (f"Recognized instruments are:\n")
            for i in instrumentType.__members__.values():
                if i.name != instrumentType.na.name:
                    helpText += (f"\t{i.name}\n")
            helpText += (f"For example, to start the server for the PPMS flavor of MultiVu:")
            helpText += (f"\n\tpython {programName} PPMS\n")
            helpText += ("Note that MultiVu must be running before starting the server.")
        helpText += ("\n\nCOMMAND OPTIONS:\n")
        helpText += ("\tTEMP? - returns the temperature, unit of measurement, and status.\n")
        helpText += ("\tTEMP target,rate,mode - sets the target temperature (K), rate (K/min), "
                     "and mode (0: Fast Settle, 1: No Overshoot)\n")
        helpText += ("\n\tFIELD? - returns the field (oe) and the state.\n")
        helpText += ("\tFIELD field, rate, approach, mode - sets the field set point (oe), the rate to reach field (oe/s), the Approach mode, and Field mode\n")
        helpText += ("\tAPPROACH MODE:\n\t0: Linear\n\t1: No Overshoot\n\t2: Ocsillate\n")
        helpText += ("\tFIELD MODE:\n\t0: Persistent (PPMS and MPMS3 only)\n\t1: Driven\n")
        helpText += ("\n\tCHAMBER? - returns the chamber state\n")
        helpText += ("\tCHAMBER code - sets the chamber mode according to the Code:\n")
        helpText += ("\t0: Seal\n\t1: Purge/Seal\n\t2: Vent/Seal\n\t3: Pump continuous\n\t4: Vent continuous\n\t5: High vacuum\n")
        return helpText
    


# if __name__ == '__main__':
#     main()
