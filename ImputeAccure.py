# -*- coding: utf-8 -*-
# Created on Tue May 18 15:17:37 2021
# @author: arosenb

 
import os, sys, getopt, platform
from pathlib import Path

import gzip

import numpy as np
from colorama import Fore, Style
import datetime
from datetime import time
import time


########################## FUNKTIONEN def                   #######################################

# Rosenberger: Bestimmen des best_guess Genotypen, auch bei gleichgrßen Dosages
def bestguess(list):
    # find number of best-guess genotypes in dosages i,i
    winner=np.flatnonzero(list == np.max(list))
    zmean = np.mean(winner)   
    return zmean

def bestguess_array(list):
    # find number of best-guess genotypes in dosages i,i
    winner = np.floor(list/np.max(list)) 
    z_array = winner / sum(winner)
    return z_array


####### exponentially weighted moving average ###############
import numpy as np

def ewma(x, alpha):
    '''
    Returns the exponentially weighted moving average of x.

    Parameters:
    -----------
    x : array-like
    alpha : float {0 <= alpha <= 1}

    Returns:
    --------
    ewma: numpy array
          the exponentially weighted moving average
    '''
    # Coerce x to an array
    x = np.array(x)
    n = x.size

    # Create an initial weight matrix of (1-alpha), and a matrix of powers
    # to raise the weights by
    w0 = np.ones(shape=(n,n)) * (1-alpha)
    p = np.vstack([np.arange(i,i-n,-1) for i in range(n)])

    # Create the weight matrix
    w = np.tril(w0**p,0)

    # Calculate the ewma
    return np.dot(w, x[::np.newaxis]) / w.sum(axis=1)


# Thormann: Basisgerüst, Bestimmung Iam hiQ 
# Rosenberger: Erweiterung durch info, r2_BEAGLE, r2_MACH
# Rosenberger: Anpassen der des best_guess Genotypen. siehe def oben 
def process_row_data(dist_data):
#!! Rosenberger: f_A zur Ausgaben hinzugefügt
   
    n_indiv = int(0)
    f_A = 0.0
    q_mean = 0.0
    dosage_mean = np.zeros(3)
    best_guess_mean = np.zeros(3)
    e = 0.0
    e2 = 0.0
    e2_ = 0.0
    f = 0.0
    fe2 = 0.0
    z = 0.0
    z_ = 0.0
    z2 = 0.0
    ze = 0.0
#!! Rosenberger: Rechengrößen e, f und e² gemäß Marchini 2010        
        
    for j in range(len(dist_data)):  # samples
           
        if sum(dist_data[j,:]) != 1 and sum(dist_data[j,:]) > 0:
            dist_data[j,:]=dist_data[j,:]/sum(dist_data[j,:])
            dist_data[j,:]=dist_data[j,:]/sum(dist_data[j,:])
            dist_data[j,:]=dist_data[j,:]/sum(dist_data[j,:])
#!! Rosenberger: Sicheren, dass Dosage-Summe =1   
               
        for k in range(3):
            q_mean += dist_data[j, k] * (1.0 - dist_data[j, k])
            f_A += k * dist_data[j, k] 
#!! Rosenberger: Hier f_A=Allel-count               
        dosage_mean += dist_data[j]
            
        if sum(dist_data[j]) > 0: best_guess_mean = best_guess_mean + bestguess_array(dist_data[j])
#!! Rosenberger:  best_guess_maen durch Sum best_guess_array ersetzt. multile best-guess genotypen adeq. berücksichtigt
        z_ = bestguess(dist_data[j]) 
        z += z_
        z2 += z_**2
        e_ = dist_data[j,1] + 2 * dist_data[j,2]
        e += e_
        e2_ = (dist_data[j,1] + 2 * dist_data[j,2]) ** 2
        e2 += e2_
        f_ = dist_data[j,1] + 4 * dist_data[j,2]
        f += f_
        fe2 += (f_ - e2_)
        ze += z_*e_
#!! Rosenberger: Rechengrößen e, e², z, z² und f gemäß Marchini 2010. Nature Reviews Genetics. Vol 11, pages 499–511  
#!! Rosenebrger: x_ nicht auchsummiert, x aufsummiert über j=Personen/Proben 

    n_indiv = int(round(sum(sum(dist_data)),0))
#!! Rosenberger: N Personen entspricht Summe aller Dosages, ausgeschlossen Personen/Proben Dosages=Tripel-Zero     
    if n_indiv>0: f_A /= (2*n_indiv)
    else: f_A=0
#!! Rosenberger: Hier f_A=Allel-count / Anz. Haplotypen (2xPersonen)
    if f_A == 0:
        # f_A = (2*max(1000,n_indiv)+1)/(2*max(1000,n_indiv)+2)
        f_A = np.finfo(np.float64).eps
    if f_A == 1:
        # f_A = 1-(2*max(1000,n_indiv)+1)/(2*max(1000,n_indiv)+2)   
        f_A = 1 - np.finfo(np.float64).eps
#!! Rosenberger: Berechungsprobelem bei f_A 0 oder 1 vermeiden (Bayes-Schätzer x+1/x+2)           
     
    q_hwe = -2.0 * f_A * (f_A - 1.0) * (3.0 * f_A ** 2 - 3.0 * f_A + 2.0)
    if n_indiv>0: 
        q_mean /= n_indiv
        dosage_mean /= n_indiv
        best_guess_mean /= n_indiv   
#!! Rosenberger: Nenner auf Anzahl valide Personen geändert, ohne dossage [0,0,0] 
  
    # calculate metrics
    if n_indiv > 0: 
        iam_chance = 1.0 - q_mean * 3.0 / 2.0
        iam_hwe = 1.0 - q_mean / q_hwe
        hiq = 1.0 - np.sqrt(1.0 - min(1,np.sum(np.sqrt(best_guess_mean)*np.sqrt(dosage_mean))))
        MAF = f_A * 100
        if ( 2*f_A*(1-f_A) ) >0:
            r2_MACH = ( (e2/n_indiv) - (e/n_indiv)**2 ) / ( 2*f_A*(1-f_A) ) 
        else: r2_MACH = -9
        if (f - (e**2/n_indiv)) * (z2 - z**2/n_indiv) >0:
            r2_BEAGLE = (ze - z*e/n_indiv)**2 / (f - (e**2/n_indiv)) / (z2 - z**2/n_indiv)
        else: r2_BEAGLE = -9
        if ( 2*n_indiv*f_A*(1-f_A) ) >0:
            info_IMPUTE = 1 - fe2 / ( 2*n_indiv*f_A*(1-f_A) )
        else: info_IMPUTE = -9                    
    else:   
        iam_chance = -9 
        iam_hwe = -9
        hiq = -9
        MAF = -9
        r2_MACH = -9
        r2_BEAGLE = -9
        info_IMPUTE = -9
                     
#!! Rosenberger: Formel für hiq geändert   
#!! Rosenebrger: iam_hiq zu hiq geändert (im ganzen Programm)  
#    print ("ICH BIN HIER")
#    print ("dist_data[0,1]: ", dist_data[0,1])
#    print (dist_data[0,1],dist_data[0,2],dist_data[0,3],"...")
#    print ("MAF (f_A)", f_A )
#    print ("N  wird summiert auf ", n_indiv )
#        print ("e  wird summiert auf ", e )
#        print ("e² wird summiert auf ",e2)   
#        print ("f wird summiert auf ",f)  
#        print ("f-e² wird summiert auf ",fe2)
#        print ("z wird summiert auf ",z)        
#        print ("z2 wird summiert auf ",z2)         
#        print ("ze wird summiert auf ",ze)         
         
#!! Rosenberger: Berechen von info und r² gemäß Marchini 2010. Marchini 2010. Nature Reviews Genetics. Vol 11, pages 499–511
#!! Rosenberger: r2_BEAGLE enstpriccht Brownings Allelic-R² (Browning 2009.The American Journal of Human Genetics 84, 210–223

    # Ausgabewerte runden
    iam_chance=round(iam_chance,3)
    iam_hwe=round(iam_hwe,3)
    hiq=round(hiq,3)
    N=round(n_indiv,0)
    MAF=round(MAF,1)
    info_IMPUTE=round(info_IMPUTE,3)
    r2_MACH=round(r2_MACH,3)      
    r2_BEAGLE=round(r2_BEAGLE,3)  
              
#    print(Style.RESET_ALL)        
    return iam_chance, iam_hwe, N, MAF, hiq, info_IMPUTE, r2_MACH, r2_BEAGLE
#!! Rosenberger: N, MAF, r2_MACH undr2_BEAGLE  zur Ausgaben hinzugefügt

########################## BASIS-SPEZIFIKATIONEN    #######################################
########################## Folder File              #######################################
if __name__ == "__main__":
#    print("---- start of main program ----",__name__)
    code_name = sys.argv[0]
    argv = sys.argv[1:]

    if getattr(sys, 'frozen', False):
        file_dir = os.path.dirname(sys.executable)
    #elif __file__:
    #    file_dir = os.path.dirname(__file__)
    #    print(os.path.dirname(__file__))
    else:
        file_dir = '.'

    # default values
    leading_columns = 5  # number of leading columns
    calculate_third_column = False  # if only 2 columns per person, the third has to be calculated
    exclude_ppl = np.zeros(0)  # ids of samples to exclude, order by their occurence in the input file
    exclude_markers = np.zeros(0)  # id of markers to exclude, refers to column marked by "marker_pos"
    marker_pos = 1  # column with marker id
    unpack_file_form = ''  # either '' or .gz
    columns = np.zeros(0, dtype='O')

    data_path = ''
    data_name = ''
    file_name = ''

    try:
        opts, args = getopt.getopt(argv, 'hf:i:l:c:p:m:n:', ['specificationFile=', 'input=', 'leadingColumns=', 'calcThirdColumn=', 'pplToExclude=', 'markersToExclude=', 'columnNames='])
    except getopt.GetoptError:
        print(code_name + ' -f <Specificationfile>')
        print('or')
        print(code_name + ' -i <Inputfile> -l <leading columns>=%i -c <calculate third column>=%r -p <list with samples to exclude>=[] -m <list with markers to exclude>=[] -n <column names separated by comma>' % (leading_columns, calculate_third_column))
        input('Press enter to exit')
        sys.exit(1)
    
    if len(opts) == 0:
        # ask for arguments via command line
        opts = []
        arg = input('Please specify relative path to input file:')
        opts.append(('-i', arg))
        arg = input('Please specify leading columns (default 5):')
        if arg == '':
            arg = '5'
        opts.append(('-l', arg))
        arg = input('Please specify whether third column needs to be calculated (default false):')
        if arg == '':
            arg = 'false'
        opts.append(('-c', arg))
        arg = input('Please specify relative path to file with samples to exclude (default empty):')
        if arg != '':
            opts.append(('-p', arg))
        arg = input('Please specify relative path to file with markers to exclude (default empty):')
        if arg != '':
            opts.append(('-m', arg))
        arg = input('Please specify names of columns separated by comma (default col[column number]):')
        if arg != '':
            opts.append(('-n', arg))
    else:
        # check for specification file argument
        spec_file_path = ''
        for opt, arg in opts:
            if opt in ('-f', '--specificationFile'):
                spec_file_path = arg
        
        if spec_file_path == '':
            spec_file = None
        else:
            spec_file = open(spec_file_path)
            opts = []  # create new array based on input file
            for row in spec_file:
                new_opt, new_arg = row.split(' ')
                new_arg = new_arg[:-1]  # get rid of line break
                if len(new_arg) > 0:
                    opts.append((new_opt, new_arg))

    # handle input arguments
    for opt, arg in opts:
        if opt == '-h':
            print(code_name + ' -f <Specificationfile>')
            print('or')
            print(code_name + ' -i <Inputfile> -l <leading columns>=%i -c <calculate third column>=%r -p <list with samples to exclude>=[] -m <list with markers to exclude>=[] -n <column names separated by comma>' % (leading_columns, calculate_third_column))
            input('Press enter to exit')
            sys.exit()
        elif opt in ('-i', '--input'):
            path_parts = Path(arg)
            data_path = path_parts.parent
            data_name = path_parts.name
            file_name_parts = path_parts.suffix
            if file_name_parts[-1] == '.gz':
                file_name = ''.join(path_parts.name.split('.')[:-2])
                unpack_file_form = '.gz'
            else:
                file_name = ''.join(path_parts.name.split('.')[:-1])
        elif opt in ('-l', '--leadingColumns'):
            leading_columns = int(arg)
        elif opt in ('-c', '--calcThirdColumn'):
            if (arg == '0') | (arg == 'False') | (arg == 'false'):
                calculate_third_column = False
            elif (arg == '1') | (arg == 'True') | (arg == 'true'):
                calculate_third_column = True
            else:
                print('Invalid input for <calculate third column>, please use either 0, False, or false or 1, True, or true')
                input('Press enter to exit')
                sys.exit(1)
        elif opt in ('-p', '--ppltToExclude'):
            if not os.path.isfile(arg):
                print('File with samples to exclude DOES NOT exists!')
                input('Press enter to exit')
                sys.exit()
            ppl_file = open(arg)
            exclude_ppl = np.array(ppl_file.read().split(','), int) - 1  # first ID is 1
            ppl_file.close()
        elif opt in ('-m', '--markersToExclude'):
            if not os.path.isfile(arg):
                print('File with markers to exclude DOES NOT exists!')
                input('Press enter to exit')
                sys.exit()
            marker_file = open(arg)
            exclude_markers = np.array(marker_file.read().split(','))
            if '\n' in exclude_markers[-1]:
                exclude_markers[-1] = exclude_markers[-1][:-1]
            marker_file.close()
        elif opt in ('-n', '--columnNames'):
            columns = arg.split(',')

    if leading_columns < 3:
        print('There needs to be at least 3 leading columns!')
        input('Press enter to exit')
        sys.exit()

    # prepare column names
    if len(columns) != leading_columns:
        if len(columns) != 0:  # not default, wrong input has been provided, write error message
            print('Number of column headers does not correspond to specified number of leading columns! Switching to default names!')
        columns = np.array(['col%i' % i for i in range(leading_columns)], 'O')
    if columns[1] != 'SNP':
        columns[1] = 'SNP'
        print('Renamed second column header to SNP!')
    if columns[2] != 'position':
        columns[2] = 'position'
        print('Renamed third column header to position!')
    header = ''
    for i in range(leading_columns):
        header += columns[i] + ' '
    header += 'N MAF Iam_chance Iam_hwe hiQ info_IMPUTE r2_MACH r2_BEAGLE\n'
            
    print("--- Specifications        ---------------------------")
    print('  leading columns:                  ', leading_columns)
    print('  calculate third column/probability:',calculate_third_column)
    print('  exclude marker:                   ',len(exclude_markers), "entries:",exclude_markers)
    print('  exclude probes:                   ',len(exclude_ppl), "entries:",exclude_ppl)
    print('  column with marker id:            ',marker_pos) 
    print('---- Data and Result file: -------------------------')
    print('  data_path:     ',data_path)

    if data_path == '':
        print('Please specify an input file!')
        input('Press enter to exit')
        sys.exit()

    data_file = Path(file_dir, data_path, data_name + unpack_file_form)
    result_file = Path(file_dir, data_path, file_name + '.accuracy')
    print('  file_name:     ',file_name)
    print('  input-file:    ',data_file)
    if (not unpack_file_form):
        print('  packed/unpacked: unpacked',)
    else:
        print('  packed/unpacked',unpack_file_form)
    print('  output-file:   ',result_file)
    print('--- ---------------------- -------------------------')

    if os.path.isfile(data_file):
        print('Input file exists!')
    else:
        print('Input file DOES NOT exists!')
        input('Press enter to exit')
        sys.exit()

    if os.path.isfile(result_file):
        print('Output file already exists!')
    else:
        print('Output file DOES NOT exist yet!')



    ########################## 1.+2. Zeile in Ausgabefile erzeugen ################################
    #!! Rosenberger: 1.Zeile Datum und anz-Zeile przesiert in parallel
    #!! Rosenberger: 2.Zeile Variablennamen in Ausgabefile schreiben  
    print(result_file)
    results = open(str(result_file), 'a')       
    with data_file.open() as file:
        results.truncate(0)
     #!! Rosenberger: Alle Zeile in Ausgebedatei löschen   
        results.write("created at ")  
        results.write(datetime.datetime.now().strftime('%d.%m.%Y %H:%M:%S'))
        results.write('\n')
        results.write(header)
    results.close() 




    ########################## ABLAUF                   #######################################
    num_samples = None
            
    #results = result_file.open('a')      

    row_counter = 0
    prep_time = 0.0
    calc_time = 0.0
    write_time = 0.0

    open_func = open if unpack_file_form == '' else gzip.open

    with open(result_file, 'a') as results:    
        with open_func(data_file, 'rt') as file:
            tic = time.perf_counter()

            for row in file:
                tic1 = time.perf_counter()
                # read data
                row_data = row.split(' ')
                if num_samples is None:
                    num_samples = (len(row_data) - leading_columns) / 2.0 if calculate_third_column else (len(row_data) - leading_columns) / 3.0
                    if num_samples % 1.0 != 0:
                        print('Invalid column specifications!')
                        input('Press enter to exit')
                        sys.exit(2)
                    num_samples = int(num_samples)
                    print('Number of samples: %i' % num_samples)
                # skip marker
                if len(exclude_markers) > 0:
                    if row_data[marker_pos] in exclude_markers:
                        continue

                # collect rows for calculation as batch
                names = row_data[:leading_columns]
                dist_data = row_data[leading_columns:]

                # preprocess data
                if calculate_third_column:
                    input_data = np.zeros((num_samples, 3))
                    input_data[:, :2] = np.vstack(np.asarray(dist_data, float)).reshape((num_samples, 2))
                    input_data[:, 2] = 1.0 - (input_data[:, 0] + input_data[:, 1])
                else:
                    input_data = np.vstack(np.asarray(dist_data, float)).reshape((num_samples, 3))
                if len(exclude_ppl) > 0:
                    if any(exclude_ppl >= len(input_data)) | any(exclude_ppl < 0):
                        print('Invalid ID(s) among excluded samples detected! They will be ignored!')
                        exclude_ppl = exclude_ppl[(exclude_ppl >= 0) & (exclude_ppl < len(input_data))]
                    input_data[exclude_ppl] *= 0.0
                # exclude samples with negative values
                input_data[input_data[:, 0] < 0.0] *= 0.0
                prep_time += time.perf_counter() - tic1

                # calculate metrics
                tic1 = time.perf_counter()
                iam_chance, iam_hwe, N, MAF, hiq, info_IMPUTE, r2_MACH, r2_BEAGLE = process_row_data(input_data)
                calc_time += time.perf_counter() - tic1

                # write results
                tic1 = time.perf_counter()
                row_result = names[0] + ' '
                for i in range(1, leading_columns):
                    row_result += names[i] + ' '
                row_result += str(int(N)) + ' '                
                row_result += str(MAF) + '% '  
                row_result += str(iam_chance) + ' '
                row_result += str(iam_hwe) + ' '
                row_result += str(hiq) + ' '
                row_result += str(info_IMPUTE) + ' '
                row_result += str(r2_MACH) + ' '
                row_result += str(r2_BEAGLE) + '\n'
               
                #!! Rosenberger: N, MAF, info_IMPUTE, r2_MACH, r2_BEAGLE  zur Ausgaben hinzugefügt

                row_counter += 1
                results.write(row_result)
                write_time += time.perf_counter() - tic1

            toc = time.perf_counter()

    # print runtime
    rt = toc - tic
    print('\nRuntime: %f' % rt)
    print('Preparation time: %f' % prep_time)
    #   print(Fore.RED + Style.BRIGHT + '\nRuntime: %f' % rt)
    #   print(Fore.WHITE + Style.BRIGHT + 'Preparation time: %f' % prep_time)
    print('Calculation time: %f' % calc_time)
    print('Writing time: %f' % write_time)
    print('Average runtime per row: %f' % (rt / row_counter))
 #   print(Style.RESET_ALL)  
    results.close()

######   EWMA berechnen         #####################################################

    #!! Starke: Ergebnisfile einlesen
    lines = open(result_file, "r")
    lines = lines.readlines()
    #!! Starke: str zu Liste von Listen aufsplitten
    for x in range(1, len(lines)):
        lines[x] = lines[x].split(' ')
   
    #!! Starke: Marker nach Position sortieren
    # 'position'-Spalte in Integer umformen:
    for i in range(2, len(lines)):
        lines[i][2] = int(lines[i][2])
    # nach 'position' sortieren:
    from operator import itemgetter
    lines[2:] = sorted(lines[2:], key=itemgetter(2))
   
    # 'position' wieder zu str umformen
    for i in range(2, len(lines)):
        lines[i][2] = str(lines[i][2])
   
     #!! Starke: Spalte mit Iam_hwe bzw hiq-Werten extrahieren
    hwe = []
    hiq = []
    for x in lines[2:len(lines)]:
        hwe.append(x[3 + leading_columns])
        hiq.append(x[4 + leading_columns])
    hwe = [float(x) for x in hwe]
    hiq = [float(x) for x in hiq]
   
    #!! Starke: temporäre Listen für hiq und hwe erzeugen 
    #in denen für die Berechnung des ewma die fehlenden Werte mit den Schwellenwerten ersetzt werden
    hwe_temp=hwe
    hiq_temp=hiq
    
    for i in range(len(hwe_temp)):
        if hwe_temp[i] == -9.0:
            hwe_temp[i] = 0.47
        if hiq_temp[i]== -0.9:
            hiq_temp[i] = 0.97
    
     #!! Starke: EWMA von Iam_hwe bzw. hiq berechnen (fehlende Werte sind als Schwellenwerte einbezogen)
    alpha = 0.1
    hwe_ewma=[round(x,3) for x in ewma(hwe_temp,alpha)]
    hiq_ewma=[round(x,3) for x in ewma(hiq_temp,alpha)]
   
   
    #!! Starke: ewmas klassifizieren und accuracy scale in neue Spalte (Liste) schreiben
    accuracy = ['accuracy']
    for i in range(len(hwe_ewma)):
        if hwe_ewma[i] > 0.47 and hiq_ewma[i] > 0.97:
            accuracy.append("cold")
        elif hwe_ewma[i] < 0.47/2 and hiq_ewma[i] < 0.97/2:
            accuracy.append("very hot")
        elif hwe_ewma[i] < 0.47 and hiq_ewma[i] < 0.97:
            accuracy.append("hot")
        else:  # hwe_ewma[i]<0.47 or hiq_ewma[i]<0.97:
            accuracy.append("tepid")
   
    #!! Starke: accuracy scale als neue Spalte in Ergebnisfile einfügen
    for x in range(1, len(lines)):
        lines[x].insert(5 + leading_columns, accuracy[x-1])
        lines[x] = " ".join(lines[x])
    lines = "".join(lines)
   
    with open(result_file, "w") as results:
        results.truncate(0)
        results.write(lines)
    results.close()
    
    input('Press enter to exit')
