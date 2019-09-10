# # Remove RFI
python clean.py files
# (1)adding "-plot" can plot the dynamic have been remove RFI. eg:
python clean.py -plot files 

# #scintillation
python scin.py files
# (1) -h, --help           show this help message and exit
# (2) -fill, --fillempty   fill the empty place where the RFI had been removed
# (3) -end , --endmethod   There have two methods to find the best end, one of the them is from lars(0), the other is useing the slope rate of ACF(1)
# (4) -plot, --plot        plot the dynamic spectrum, ACF and secendary. Options: all, dysp, acf, sec
# (5) -s, --savefigure     save the dysp
# (6) -q, --quiet          Do not print text information.


# # combine files into one dynamic spectrum. there have to use txt files or exsist txt files
python combine.py files
