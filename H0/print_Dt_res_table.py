from H0.H0_Data import *
import argparse as ap


def print_table_Dt(Dt_res_cl):
   
    labels = Dt_res_cl.labels
    Dt_res = Dt_res_cl.results
    Dt_err = Dt_res_cl.error.tot
    

    newline = "\n"
    str_tbl  = r"\begin{table}\centering"
    str_tbl += newline
    str_tbl += r"\resizebox{.8\columnwidth}{!}{"
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
        str_tbl += f"$\Delta t$ {l} [d]"+fin_and
        
    str_tbl += r"\\"
    str_tbl += newline  
    str_tbl += r"\hline"
    str_tbl += newline  
    str_tbl +="  "
    for i,(r,er) in enumerate(zip(Dt_res,Dt_err)):
        if i!=len(Dt_res)-1:
            fin_and=" & "
        else:
            fin_and=""
        str_tbl+= str(np.round(r,1))+r" \pm "+str(np.round(er,1))+fin_and
    str_tbl +=r"\\"
    str_tbl += newline  

    str_tbl += r"\hline"
    str_tbl +=newline
    str_tbl +=r"\end{tabular} }"
    str_tbl +=newline
    str_tbl +=r"\caption{Combined results for the time delay between the various couple of images.}"
    str_tbl +=newline
    str_tbl +=r"\label{tab:timedelay_results}"
    str_tbl +=newline
    str_tbl +=r"\end{table}"
    print(str_tbl)



if __name__ == '__main__':

    # Adapted copy from print_Df_res.p
    help_timedelay = "Name of the posterior of the difference of Time Delay (name of config file)"

    parser = ap.ArgumentParser(prog="python {}".format(os.path.basename(__file__)),
                               description="Print latex table for Dt results from time delay")
    parser.add_argument("-dtn","--Dt_name",dest='dt_name', type=str,default=".",
                        metavar='dt_name', action='store',
                        help=help_timedelay)  

    args      = parser.parse_args()
    dt_name   = args.dt_name
    cnf,Dt_res = get_Dt_post(dt_name,link=False)
    print_table_Dt(Dt_res_cl=Dt_res)
