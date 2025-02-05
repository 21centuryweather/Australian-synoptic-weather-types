import pickle

colors=[(134/255,0/255,34/255),(241/255,0/255,241/255),(255/255,134/255,255/255),(255/255,241/255,255/255), #Pink WH
        (255/255,204/255,51/255),(255/255,245/255,204/255), #Yellow CH
        (153/255,15/255,15/255),(178/255,44/255,44/255),(204/255,81/255,81/255),(229/255,126/255,126/255),(255/255,178/255,178/255), #Red EH
        (153/255,84/255,15/255),(204/255,142/255,81/255),(255/255,216/255, 178/255), #Brown TH
        (107/255,153/255,15/255),(163/255,204/255,81/255),(195/255,229/255,126/255), #Green FH (133/255,178/255,44/255),
        (66/255,44/255,178/255),(143/255,126/255,229/255), #Purple WCT
        (5/255,67/255,113/255),(15/255,107/255,153/255),(44/255,133/255,178/255),(81/255,163/255,204/255),(126/255,195/255,229/255),(178/255,229/255,255/255), #Blue COL
        (0/255,60/255,48/255),(1/255,102/255,95/255),(53/255,151/255,143/255),(128/255,205/255,193/255),(199/255,234/255,229/255)
       ]

SWTs = {
    23: {"WR": "WH"   ,"SWT": "A" , "Description": "", "color": colors[0]},
    5:  {"WR": "WH"   ,"SWT": "B" , "Description": "", "color": colors[1]},
    19: {"WR": "WH"   ,"SWT": "C" , "Description": "", "color": colors[2]},
    1:  {"WR": "WH"   ,"SWT": "D" , "Description": "", "color": colors[3]},
    17: {"WR": "CH"   ,"SWT": "A" , "Description": "", "color": colors[4]},
    6:  {"WR": "CH"   ,"SWT": "B" , "Description": "", "color": colors[5]},
    15: {"WR": "EH"   ,"SWT": "A" , "Description": "", "color": colors[6]},
    10: {"WR": "EH"   ,"SWT": "B" , "Description": "", "color": colors[7]},
    29: {"WR": "EH"   ,"SWT": "C" , "Description": "", "color": colors[8]},
    11: {"WR": "EH"   ,"SWT": "D" , "Description": "", "color": colors[9]},
    20: {"WR": "EH"   ,"SWT": "E" , "Description": "", "color": colors[10]},
    16: {"WR": "TH"   ,"SWT": "A" , "Description": "", "color": colors[11]},
    28: {"WR": "TH"   ,"SWT": "B" , "Description": "", "color": colors[12]},
    4:  {"WR": "TH"   ,"SWT": "C" , "Description": "", "color": colors[13]},
    18: {"WR": "FH"   ,"SWT": "A" , "Description": "", "color": colors[14]},
    22: {"WR": "FH"   ,"SWT": "B" , "Description": "", "color": colors[15]},
    12: {"WR": "FH"   ,"SWT": "C" , "Description": "", "color": colors[16]},
    21: {"WR": "WCT"  ,"SWT": "A" , "Description": "", "color": colors[17]},
    25: {"WR": "WCT"  ,"SWT": "B" , "Description": "", "color": colors[18]},
    14: {"WR": "COL"  ,"SWT": "A" , "Description": "", "color": colors[19]},
    7:  {"WR": "COL"  ,"SWT": "B" , "Description": "", "color": colors[20]},
    9:  {"WR": "COL"  ,"SWT": "C" , "Description": "", "color": colors[21]},
    30: {"WR": "COL"  ,"SWT": "D" , "Description": "", "color": colors[22]},
    27: {"WR": "COL"  ,"SWT": "E" , "Description": "", "color": colors[23]},
    13: {"WR": "COL"  ,"SWT": "F" , "Description": "", "color": colors[24]},
    24: {"WR": "AM"   ,"SWT": "A" , "Description": "", "color": colors[25]},#colors[0], "x": True}, 
    26: {"WR": "AM"   ,"SWT": "B" , "Description": "", "color": colors[26]},#colors[4], "x": True}, 
    2:  {"WR": "AM"   ,"SWT": "C" , "Description": "", "color": colors[27]},#colors[4], "x": True}, 
    8:  {"WR": "AM"   ,"SWT": "D" , "Description": "", "color": colors[28]},#colors[6], "x": True}, 
    3:  {"WR": "AM"   ,"SWT": "E" , "Description": "", "color": colors[29]},#colors[14], "x": True}
}

# save dictionary to person_data.pkl file
with open('SWT_WR_definition.pkl', 'wb') as fp:
    pickle.dump(SWTs, fp)
    print('dictionary saved successfully to file')

import pickle
with open(outpath+'SWT_WR_definition_v2.pkl', 'rb') as fp:
    defs = pickle.load(fp)
    #print('Person dictionary')
    #print(person)
