#!/usr/bin/env python3
"""

This is a custom task for finding a patient's P/F ratio, which is relevant
to assessing the amount of oxygen in the blood. The P/F ratio is defined as:

    (partial pressure of O2 in arterial blood) / (fraction of inspired oxygen)

Using common symbols, this can be written as:

    P/F = PaO2 / FiO2

Generally a P/F ratio less than 300 means the patient has hypoxemia.

The P/F ratio can be estimated from pulse ox measurements and O2 flow rates.

Conversion of pulse ox value (SpO2 %) to O2 partial pressure (mmHg):

    SpO2 (%)    PaO2 (mmHg)
    --------    -----------
      80            44
      81            45
      82            46
      83            47
      84            49
      85            50
      86            52
      87            53
      88            55
      89            57
      90            60
      91            62
      92            65
      93            69
      94            73
      95            79
      96            86
      97            96
      98           112
      99           145

The FiO2 percentage can be estimated from O2 flow rates of different breathing
devices. Normal breathing in room air produces an FiO2 value of 20%.

For a nasal cannula, each L/min of O2 adds 4% to the FiO2 value:
FiO2 = 0.2 + (0.04) * flow_rate_l_per_min, 1 <= flow_rate <= 10


       flow rate (L/min)    FiO2 (%)
       ----------------     --------
           1                  24
           2                  28
           3                  32
           4                  36
           5                  40
           6                  44
           7                  48
           8                  52
           9                  56
          10                  60

For a nasopharyngeal catheter:
FiO2 = (0.1) * flow_rate_l_per_min, 4 <= flow_rate <= 6

       flow rate (L/min)    FiO2 (%)
       ----------------     --------
           4                  40
           5                  50
           6                  60

For a face mask with no reservoir (rebreather mask):
FiO2 = 0.35 + (0.04) * (flow_rate_l_per_min - 5), 5 <= flow_rate <= 10

       flow rate (L/min)    FiO2 (%)
       ----------------     --------
           5                  35
           6                  39
           7                  43
           8                  47
           9                  51
          10                  55

For a face mask with reservoir (non-rebreather):
FiO2 = 0.6 + (0.10) * (flow_rate_l_per_min - 6), 6 <= flow_rate <= 9
FiO2 = 0.9 + (0.05) * (flow_rate_l_per_min - 9), 9 <= flow_rate <= 10

       flow rate (L/min)    FiO2 (%)
       ----------------     --------
           6                  60
           7                  70
           8                  80
           9                  90
          10                  95

Venturi Mask:

FiO2 = 0.24 + (0.020)(flow_rate_l_per_m -  2),  2 <= flow_rate <= 4
       0.28 + (0.015)(flow_rate_l_per_m -  4),  4 <= flow_rate <= 6
       0.31 + (0.020)(flow_rate_l_per_m -  6),  6 <= flow_rate <= 8
       0.35 + (0.025)(flow_rate_l_per_m -  8),  8 <= flow_rate <= 10
       0.40 + (0.040)(flow_rate_l_per_m - 10), 10 <= flow_rate <= 15

       flow rate (L/min)    FiO2 (%)
       ----------------     --------
           2                  24
           4                  28
           6                  31
           8                  35
          10                  40
          15                  60

"""

import os
import re
import sys
import json
import argparse
from pymongo import MongoClient
from collections import namedtuple

# try:
#     #from . import finder_overlap as overlap
#     import finder_overlap as overlap
# except:
#     from algorithms.finder import finder_overlap as overlap

# default value for all fields
EMPTY_FIELD = None

O2_VALUE_FIELDS = [
    'text',
    'start',
    'end',
    'spo2',
    'pao2',
    'pf_ratio'
    'flow_rate',
    'device'
]
O2Value = namedtuple('O2Value', O2_VALUE_FIELDS)


###############################################################################

_VERSION_MAJOR = 0
_VERSION_MINOR = 1

# set to True to enable debug output
_TRACE = False


#_str_connector = r'\s?([-/:=\s]|of|on)\s?'
_str_connector = r'([-/:=\s]|of|on a|on|to|with)+'

# O2 saturation header
_str_o2_sat_hdr = r'\b(spo2|sao2|pox|so2|(o2|oxygen)[-\s]?saturation|'       +\
    r'o2[-\s]sat\.?s?|satting|o2sats?|sat\.?s?|pulse ox|o2|'                 +\
    r'desatt?ing|desat)%?'

_str_units = r'\(?(percent|pct\.?|%|cmH2O|mmHg)\)?'

# o2 saturation value
_str_o2_val_range = r'\b(was|from)?\s?(?P<val1>\d+)(\s?' + _str_units + r')?'+\
    r'(\s?(\-|to)\s?)(?P<val2>\d+)(\s?' + _str_units + r')?'
_str_o2_value = r'(?P<value>\d+)(\s?' + _str_units + r')?'
_str_o2_val = r'(' + _str_o2_val_range + r'|' + _str_o2_value + r')'

# O2 flow rate in L/min
_str_flow_rate = r'\d+\s?(Liters|L)(/min\.?)?'

# devices and related acronyms
#         RA : radial artery, right atrial
#         FM : face mask
#        BVM : bag valve mask
#  venti mask: Venturi mask
#       bipap: bilevel positive airway pressure
# _str_device = r'\(?(bipap|non[-\s]?rebreather(\smask)?|nasal\s?cannula|'     +\
#     r'room air|([a-z]+)?(-|\s)?mask|cannula|psv|nrb|ra|fm|nc|air|bvm|'       +\
#     r'rebreather(\smask)?|on\svent(ilator)?)\)?'

_str_device = r'\(?(nc|nrb|bvm|ra|fm|room air|air|nasal\s?cannula|'          +\
    r'cannula|non[-\s]?rebreather(\smask)?|[a-z\s]+[-\s]?mask|mask|'         +\
    r'vent(ilator)?|\d+%\s?[a-z]+[-\s]?mask|\d+%\s[a-z]+[-\s]?tent|'         +\
    r'bipap\s\d+/\d+\s?(with\s\d+L|\d+%))\)?'

# finds 98% RA" and similar (value to the left)
_str_o2_0 = r'(?<![-:=/])(?<=\s)' + _str_o2_val + r'(' + _str_connector + r')?' +\
    r'(?P<device>' + _str_device +r')'
_regex_o2_0 = re.compile(_str_o2_0, re.IGNORECASE)

# finds "spo2: 98% on 2L NC" and similar
_str_o2_1 = _str_o2_sat_hdr + r'(' + _str_connector + r')?' + _str_o2_val    +\
    r'(' + _str_connector + r')?' + _str_flow_rate                           +\
    r'(' + _str_connector + r')?' + r'(?P<device>' + _str_device + r')'
_regex_o2_1 = re.compile(_str_o2_1, re.IGNORECASE)

# finds "spo2: 98%/NC" and similar
_str_o2_2 = _str_o2_sat_hdr + r'(' + _str_connector + r')?' + _str_o2_val    +\
    r'(' + _str_connector + r')?' + r'(?P<device>' + _str_device + r')'
_regex_o2_2 = re.compile(_str_o2_2, re.IGNORECASE)

# finds "spo2: 98%/3L" and similar
_str_o2_3 = _str_o2_sat_hdr + r'(' + _str_connector + r')?' + _str_o2_val    +\
    r'(' + _str_connector + r')?' + _str_flow_rate
_regex_o2_3 = re.compile(_str_o2_3, re.IGNORECASE)

# finds "spo2=98%" and similar
_str_o2_4 = _str_o2_sat_hdr + r'(' + _str_connector + r')?' + _str_o2_val
_regex_o2_4 = re.compile(_str_o2_4, re.IGNORECASE)

_REGEXES = [
    _regex_o2_0,
    _regex_o2_1,
    _regex_o2_2,
    _regex_o2_3,
    _regex_o2_4
]

# # o2 partial pressure
# _str_pao2 = r'\b(pao2|partial pressure of (oxygen|o2))'

# # fraction of inspired oxygen
# _str_fio2 = r'\bfio2'



# connector of O2 sat to device:
#     on or symbol
# flow rate:
#     \dL or \dL/min
# device:
#     \d\d%\s[a-z]+(-\s)?mask

# on \dL(/min)? device
# \dLdevice

# on device
# on 55% [a-z]+(-\s)?mask
# /device
# device
# \dL and no device


# # O2 flow rates and percentages for devices
# _str_flow_rate = r'(' + _str_connector + r')' + r'\d\s?(l|L)' +\
#     _str_connector + r'(' + _str_device + r')?'
# _regex_flow = re.compile(_str_flow_rate, re.IGNORECASE)



# _str_pulse_ox_1 = _str_spo2 + r'(' + _str_connector + r')?' + r'\d+' +\
#     r'(' + _str_units + r'\s?)?' + r'(on\s?|[-/:=\s])?'              +\
#     r'(\d+\s?(liters|l)\s?)?' + r'(' + _str_connector + r')?'        +\
#     r'(' + _str_device + r')?'

# _str_pulse_ox_2 = _str_spo2 + r'%(' + _str_connector + r')' +\
#     r'\d+(' + _str_connector + r')' + _str_device

# _str_pulse_ox = r'(' + _str_pulse_ox_1 + r'|' + _str_pulse_ox_2 + r')'
# _regex_pulse_ox = re.compile(_str_pulse_ox, re.IGNORECASE)



###############################################################################
def _enable_debug():

    global _TRACE
    _TRACE = True
    

###############################################################################
def _cleanup(sentence):
    """
    Apply some cleanup operations to the sentence and return the
    cleaned sentence.
    """

    # collapse repeated whitespace
    sentence = re.sub(r'\s+', r' ', sentence)

    return sentence

    
###############################################################################
def run(sentence):
    """
    Find values related to oxygen saturation, flow rates, etc. Compute values
    such as P/F ratio when possible. Returns a JSON array containing info
    on all values extracted or computed.
    """

    cleaned_sentence = _cleanup(sentence)

    for i, regex in enumerate(_REGEXES):
        match = regex.search(cleaned_sentence)
        if match:
            match_text = match.group()
            print('\t[{0}]: {1}'.format(i, match_text))
    
    # match = _regex_pulse_ox.search(cleaned_sentence)
    # if match:
    #     print('\t{0}'.format(match.group()))

    # match = _regex_flow.search(cleaned_sentence)
    # if match:
    #     print('\t{0}'.format(match.group()))
    

###############################################################################
def get_version():
    return '{0} {1}.{2}'.format(_MODULE_NAME, _VERSION_MAJOR, _VERSION_MINOR)

                        
###############################################################################
if __name__ == '__main__':

    # for command-line testing only

    parser = argparse.ArgumentParser(
        description='test P/F ratio task locally')

    parser.add_argument('--debug',
                        action='store_true',
                        help='print debugging information')

    args = parser.parse_args()

    if 'debug' in args and args.debug:
        enable_debug()

    SENTENCES = [
        'Vitals were HR=120, BP=109/44, RR=29, POx=93% on 8L FM',
        'Vitals: T: 96.0  BP: 90/54 P: 88 R: 16 18 O2:88/NRB',
        'Vitals: T 98.9 F BP 138/56 P 89 RR 28 SaO2 100% on NRB',
        'Vitals were T 98 BP 163/64 HR 73 O2 95% on 55% venti mask',
        'VS: T 95.6 HR 45 BP 75/30 RR 17 98% RA.',
        'VS T97.3 P84 BP120/56 RR16 O2Sat98 2LNC',
        'Vitals: T: 99 BP: 115/68 P: 79 R:21 O2: 97',
        'Vitals - T 95.5 BP 132/65 HR 78 RR 20 SpO2 98%/3L',
        'VS: T=98 BP= 122/58  HR= 7 RR= 20  O2 sat= 100% 2L NC',
        'Vitals: T: 97.7 P:100 R:16 BP:126/95 SaO2:100 Ra',
        'VS:  T-100.6, HR-105, BP-93/46, RR-16, Sats-98% 3L/NC',
        'VS - Temp. 98.5F, BP115/65 , HR103 , R16 , 96O2-sat % RA',
        'Vitals: Temp 100.2 HR 72 BP 184/56 RR 16 sats 96% on RA',
        'PHYSICAL EXAM: O: T: 98.8 BP: 123/60   HR:97    R 16  O2Sats100%',
        'VS before transfer were 85 BP 99/34 RR 20 SpO2% 99/bipap 10/5 50%.',
        'Initial vs were: T 98 P 91 BP 122/63 R 20 O2 sat 95%RA.',
        'Initial vitals were HR 106 BP 88/56 RR 20 O2 Sat 85% 3L.',
        'Initial vs were: T=99.3 P=120 BP=111/57 RR=24 POx=100%.',        
        'At transfer vitals were HR=120 BP=109/44 RR=29 POx=93% on 8L FM.',
        "Vitals as follows: BP 120/80 HR 60-80's RR  SaO2 96% 6L NC.",
        'Vital signs were T 97.5 HR 62 BP 168/60 RR 18 95% RA.',
        'T 99.4 P 160 R 56 BP 60/36 mean 44 O2 sat 97% Wt 3025 grams ',
        'Vital signs were T 97.0 BP 85/44 HR 107 RR 28 and SpO2 91% on NRB.',
        'Vitals were T 95.6 HR 67 BP 143/79 RR 16 and O2 sat 92% on room air and 100% on 3 L/min nc',
        'Vitals were Temp. 100.8 Pulse: 103 RR: 28 BP: 84/43 O2Sat: 88 O2 Flow: 100 (Non-Rebreather).',
        'Vitals were T 97.1 HR 76 BP 148/80 RR 25 SpO2 92%/RA.',
        'Tm 96.4, BP= 90-109/49-82, HR= paced at 70, RR= 24, O2 sat= 96% on 4L',
        'Vitals were T 97.1 BP 80/70 AR 80 RR 24 O2 sat 70% on 50% flowmask',
        'Vitals: T: 98.9 degrees Farenheit BP: 120/49 mmHg supine ' +\
        'HR 84 bpm RR 13 bpm O2: 100% PS 18/10 FiO2 40%',
        'Vitals T 99.2 (baseline ~96-97), BP 91/50, HR 63, RR 12, ' +\
        'satting 95% on trach mask',
        'O2 sats 98-100%',
        'Pt. desating to 88%',
        "O2 sats increasing back up to low to mid 90's on Bipap",
        'spo2 difficult to monitor but appeared to remain ~ 96-100% on bipap 8/5',
        'using BVM w/ o2 sats 74%',
        'Placed on BiPAP 10/5 with 100% FiO2',
        'Continues to need bipap 12/10 with 10L bleed in',
        'desat to 83 with 100% face tent and 4 l n.c.',

        'Ventilator mode: CMV/ASSIST/AutoFlow   Vt (Set): 550 (550 - 550) mL ' +\
        'Vt (Spontaneous): 234 (234 - 234) mL   RR (Set): 16 ' +\
        'RR (Spontaneous): 0   PEEP: 5 cmH2O   FiO2: 70%   RSBI: 140 ' +\
        'PIP: 25 cmH2O   SpO2: 98%   Ve: 14.6 L/min',

        'Vt (Spontaneous): 608 (565 - 793) mL   PS : 15 cmH2O   ' +\
        'RR (Spontaneous): 27   PEEP: 10 cmH2O   FiO2: 50%   '    +\
        'RSBI Deferred: PEEP > 10   PIP: 26 cmH2O   SpO2: 99%   ' +\
        'ABG: 7.41/39/81/21/0   Ve: 17.4 L/min   PaO2 / FiO2: 164',
        
        'Respiratory: Vt (Set): 600 (600 - 600) mL   Vt (Spontaneous): 743 ' +\
        '(464 - 816) mL  PS : 5 cmH2O   RR (Set): 14   RR (Spontaneous): 19' +\
        ' PEEP: 5 cmH2O   FiO2: 50%   RSBI: 49   PIP: 11 cmH2O   '           +\
        'Plateau: 20 cmH2O   SPO2: 99%   ABG: 7.34/51/109/25/0   '           +\
        'Ve: 10.3 L/min   PaO2 / FiO2: 218',
        
        'an oxygen saturation of 96% on 2 liters',

        'the respiratory rate was 21,\nand the oxygen saturation was 80% ' +\
        'to 92% on a 100% nonrebreather mask',

        'temperature 100 F., orally.  O2 saturation 98% on room air'
    ]

    for i, sentence in enumerate(SENTENCES):
        print('[{0:2d}]: {1}'.format(i, sentence))
        result = run(sentence)

    
