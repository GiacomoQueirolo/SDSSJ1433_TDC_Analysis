from H0.H0_Data import *
import argparse as ap

from Plots.plotting_tools import get_median_with_error
from Utils.statistical_tools import marginalise_prob,simp,simp_wo_sig


def print_table_Dphi(Df_pdf,Df_bins,labels):
    bin_densities = marginalise_prob(Df_pdf,Df_bins) 
    Dphi_res = []
    for i in range(3):
        for j in range(3):
            if [i,j]!=[0,1] and [i,j]!=[0,2] and [i,j]!=[1,2]:
                if i==j:
                    Dphi_res.append(get_median_with_error(bin_densities[i],Df_bins[i],ret_str=True))

    newline = "\n"
    str_tbl  = r"\begin{table}\centering"
    str_tbl += newline
    str_tbl += r"\resizebox{.9\columnwidth}{!}{"
    str_tbl += newline

    str_tbl +=r"\begin{tabular}{ccc}"
    str_tbl += newline
    str_tbl += r"\hline"
    str_tbl += newline
    str_tbl +="  "
    for i,l in enumerate(labels):
        if i!=len(labels)-1:
            fin_and=" & "
        else:
            fin_and=""
        str_tbl += f"$\Delta \phi$ {l} [arcsec$^2$]"+fin_and
        
    str_tbl += r"\\"
    str_tbl += newline  
    str_tbl += r"\hline"
    str_tbl += newline  
    str_tbl +="  "
    for i,r in enumerate(Dphi_res):
        if i!=len(Dphi_res)-1:
            fin_and=" & "
        else:
            fin_and=""        
        str_tbl+= r+fin_and

    str_tbl +=r"\\"
    str_tbl += newline  

    str_tbl += r"\hline"
    str_tbl +=newline
    str_tbl +=r"\end{tabular} }"
    str_tbl +=newline
    str_tbl +=r"\caption{Combined posterior for the difference of Fermat potential at the image positions.}"
    str_tbl +=newline
    str_tbl +=r"\label{tab:comb_df_res}"
    str_tbl +=newline
    str_tbl +=r"\end{table}"
    print(str_tbl)



if __name__ == '__main__':

    # Adapted copy from print_Df_res.p
    help_fermatpot = "Name of the posterior of the difference of Fermat Potential (name of combined setting file)"

    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Print latex table for Dphi results from fermat potential")
    parser.add_argument("-dfn","--Df_name",dest='df_name', type=str,default=".",
                        metavar='df_name', action='store',
                        help=help_fermatpot)
     


    args      = parser.parse_args()
    dphi_name = args.df_name
    cmb,pdf,bins = get_Df_post(dphi_name,link=False) 
    if cmb.BC:
        labels = ["AB","AC","BC"]
    else:
        labels = ["AB","AC","AD"]
    print_table_Dphi(Df_pdf=pdf,Df_bins=bins,labels=labels)
