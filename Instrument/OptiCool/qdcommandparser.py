#!/usr/bin/env python3
#
# qdcommandparser.py
#

import sys
from parse_inputs import inputs
from qdscafolding import QDInstrumentSim
from qdinstrument import QDInstrument

class QdCommandParser:

    def __init__(self, instrument_type, line_term = '\r\n', simulateMode=False):
        self.cmd_list = {'TEMP': (self.set_temperature, self.get_temperature),
                         'FIELD': (self.set_field, self.get_field),
                         'CHAMBER': (self.set_chamber, self.get_chamber),
                         'CHAMB': (self.set_chamber, self.get_chamber)}
        if simulateMode:
            self._instrument = QDInstrumentSim(instrument_type)
        else:
            if sys.platform == 'win32':
                try:
                    import win32com.client
                    import pythoncom
                except ImportError:
                    print("Must import the pywin32 module.  Use:  ")
                    print(f"\tconda install -c anaconda pywin32")
                    print("   or")
                    print(f"\tpip install pywin32")
                    exit()

                self._instrument = QDInstrument(instrument_type)
            else:
                print('The server only works on a Windows machine. However, the server can be tested using the -s flag.\n')
                instrumentInfo = inputs()
                instrumentInfo.parseInput(['-h'])

        self._line_term = line_term

    def parse_cmd(self, data):
        cmd = data.split(' ')[0]
        for test_cmd in self.cmd_list:
            if cmd.find(test_cmd) == 0:
                if cmd.find(test_cmd + '?') == 0:
                    return str(self.cmd_list[test_cmd][1]()) + self._line_term
                else:
                    try:
                        cmd, arg_string = data.split(' ', 1)
                    except:
                        return 'No argument(s) given for command {0}.'.format(test_cmd) + self._line_term
                    return str(self.cmd_list[test_cmd][0](arg_string)) + self._line_term
        return 'Unknown command: {0}.'.format(data) + self._line_term

    def get_temperature(self):
        ret = self._instrument.get_temperature()
        # print(ret)

        # State code dictionary
        TempStates = {
            "1": "Stable",
            "2": "Tracking",
            "5": "Near",
            "6": "Chasing",
            "7": "Pot Operation",
            "10": "Standby",
            "13": "Diagnostic",
            "14": "Impedance Control Error",
            "15": "General Failure",
        }
        # Rebuild the return tuple with status code translated
        ret = (ret[2], ret[0], TempStates[str(ret[1])])

        # Suppressing the error state; not sure what users do with that info
        return '"TEMP?",{0:7.3f},"K","{2}"'.format(*ret)

    def set_temperature(self, arg_string):
        try:
            temperature, rate, mode = arg_string.split(',')
            temperature = float(temperature)
            rate = float(rate)
            mode = int(mode)
            err = self._instrument.set_temperature(temperature, rate, mode)
            if err == 0:
                return '"TEMP Command Received"'
            else:
                return err
        except:
            return '"Argument error in TEMP command"'

    def get_field(self):
        ret = self._instrument.get_field()

        # State code dictionary
        MagStates = {
            "1": "Stable",
            "2": "Switch Warming",
            "3": "Switch Cooling",
            "4": "Holding (Driven)",
            "5": "Iterate",
            "6": "Ramping",
            "7": "Ramping",
            "8": "Resetting",
            "9": "Current Error",
            "10": "Switch Error",
            "11": "Quenching",
            "12": "Charging Error",
            "14": "PSU Error",
            "15": "General Failure",
        }
        # Rebuild the return tuple with status code translated
        ret = (ret[2], ret[0], MagStates[str(ret[1])])

        # Suppressing the error state; not sure what users do with that info
        return '"FIELD?",{0:8.2f},"Oe","{2}"'.format(*ret)


    def set_field(self, arg_string):
        try:
            field, rate, approach, mode = arg_string.split(',')
            field = float(field)
            rate = float(rate)
            approach = int(approach)
            mode = int(mode)
            err = self._instrument.set_field(field, rate, approach, mode)
            if err == 0:
                return '"FIELD Command Received"'
            else:
                return err
        except:
            return '"Argument error in FIELD command"'

    def get_chamber(self):
        ret = self._instrument.get_chamber()
        # print(ret)

        ChamberStates = {
            "0": "Sealed",
            "1": "Purged and Sealed",
            "2": "Vented and Sealed",
            "3": "Sealed",
            "4": "Performing Purge/Seal",
            "5": "Performing Vent/Seal",
            "6": "Pre-HiVac",
            "7": "HiVac",
            "8": "Pumping Coninuously",
            "9": "Flooding Continuously",
            "14": "HiVac Error",
            "15": "General Failure",
        }

        # Rebuild the return tuple with status code translated
        ret = (ret[0], ChamberStates[str(ret[1])])

        # Suppressing the error state; not sure what users do with that info
        return '"CHAMBER?",,,"{1}"'.format(*ret)

    def set_chamber(self, arg_string):
        try:
            code = arg_string
            code = int(code)
            err = self._instrument.set_chamber(code)
            # returned err values seem different for the chamber compared with temp/field
            # 1 indicates success here, whereas it's zero for the others?
            if err == 1:
                return '"CHAMBER Command Sent"'
            else:
                return err
        except:
            return '"Argument error in CHAMBER command"'
