import sys
import argparse
from numpy import round
from Utils.tools import *
str_reg = '# Region file format: DS9 version 4.1\n\
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n\
image\n'

if __name__=="__main__":

    parser = argparse.ArgumentParser(description="Produce the mask.reg to import into ds9 images")
    parser.add_argument("-i","--invert_mask",dest="invert",action="store_true",default=False,help="if true, find mask.reg and print the corresponding masks to input in the setting file")
    parser.add_argument('SETTING_FILE',nargs="+", help="setting file to consider")
    present_program(sys.argv[0])

    args      = parser.parse_args()
    invert    = args.invert
    settings  = get_setting_module(args.SETTING_FILE,1)
    for setting in settings:
        print_setting(setting)
        reg_file = get_savefigpath(setting=setting)+"/mask.reg"
        if invert:
            with open(reg_file, 'r') as f:    
                reg_lines = f.readlines()
            xs,ys,rs  = [],[],[]
            for rl  in  reg_lines:
                if rl[:7]== "circle(":
                    xyr = rl.split("(")[1].split(")")[0]
                    x,y,r = xyr.split(",")
                    xs.append(float(x))
                    ys.append(float(y))
                    rs.append(float(r))
            #must not be corrected by -1 bc it's done in the next step in the setting file 
            print(f"x_mask = {[ round(x,6) for x in xs] }" )
            print(f"y_mask = {[ round(y,6) for y in ys] }" )
            print(f"r_mask = {[ round(r,6) for r in rs] }" )
        else:
            x_mask = setting.x_mask+1
            y_mask = setting.y_mask+1
            r_mask = setting.r_mask

            for x,y,r in zip(x_mask,y_mask,r_mask):
                str_reg+=f"circle({x},{y},{r}) # color=chocolate dash=1\n"

            with open(reg_file, "w") as text_file:
                    print(str_reg, file=text_file)
            print(f"Saved file {reg_file}")
    success(sys.argv[0])

